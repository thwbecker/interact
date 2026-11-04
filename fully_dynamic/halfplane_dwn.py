"""2-D in-plane (plane strain) stresses of a slip ramp on a straight segment in a full plane and
in a half plane with a traction-free surface, by the discrete wavenumber method (DWN): plane-wave
expansion in the horizontal wavenumber k and complex frequency omega, the free surface reflection
solved per k, the segment integral analytic per k. Prototype for exact half-plane kernels
(fd_plan.md M10). Frame: x horizontal, z depth (positive down), the half plane z > 0.
Source: segment of half-length a, centre (xs, zs), unit tangent s = (sx, sz), unit slip rate
along s since t = 0. Receiver at (xr, zr). Stresses returned in the (x, z) frame as
(sxx, sxz, szz)."""
import numpy as np


def _fields(k, w, alpha, beta, rho):
    lam = rho*(alpha**2 - 2*beta**2)
    mu = rho*beta**2
    ka = w/alpha
    kb = w/beta
    na = np.sqrt(k**2 - ka**2 + 0j)
    nb = np.sqrt(k**2 - kb**2 + 0j)
    na = np.where(na.real < 0, -na, na)
    nb = np.where(nb.real < 0, -nb, nb)
    return lam, mu, ka, kb, na, nb


def _incident(k, w, alpha, beta, rho, s, a, Zsign):
    """Plane-wave amplitudes of the stress (at depth offset Z = zr - zs with sign Zsign,
    excluding the factors exp(i k (xr - xs)) exp(-nu |Z|)), for the alpha and the beta waves,
    integrated over the segment: returns dict c -> (sxx, sxz, szz) amplitude arrays and nu_c."""
    lam, mu, ka, kb, na, nb = _fields(k, w, alpha, beta, rho)
    sx, sz = s
    nx, nz = sz, -sx                           # normal: the sign makes the full-plane part equal the closed form (inplane.segment_ramp_stress)
    m = np.array([[2*sx*nx, sx*nz + sz*nx], [sx*nz + sz*nx, 2*sz*nz]])*mu   # moment per unit slip and length
    out = {}
    for c, nu in (("a", na), ("b", nb)):
        Dr = [1j*k, -nu*Zsign]                 # receiver derivatives d/dx, d/dz on exp(i k X - nu |Z|)
        Ds = [-1j*k, nu*Zsign]                 # source derivatives d/dxi_x, d/dxi_z
        g = 1.0/(2*nu)
        # displacement amplitudes u_i = sum_pq m_pq Ds_q G_ip, G_ip = coef/(rho w^2)
        u = []
        for i in range(2):
            acc = 0
            for p in range(2):
                if c == "b":
                    Gip = (kb**2*(1.0 if i == p else 0.0) + Dr[i]*Dr[p])*g/(rho*w**2)
                else:
                    Gip = -(Dr[i]*Dr[p])*g/(rho*w**2)
                for q in range(2):
                    acc = acc + m[p, q]*Ds[q]*Gip
            u.append(acc)
        ux, uz = u
        exx = Dr[0]*ux
        ezz = Dr[1]*uz
        sxx = (lam + 2*mu)*exx + lam*ezz
        szz = lam*exx + (lam + 2*mu)*ezz
        sxz = mu*(Dr[1]*ux + Dr[0]*uz)
        # segment integral: source points xi = xs + tau s; exp(i k (xr - xi_x) - nu Zsign (zr - xi_z))
        cc = -1j*k*sx + nu*Zsign*sz
        seg = np.where(np.abs(cc*a) < 1e-8, 2*a + 0j, 2*np.sinh(cc*a)/np.where(cc == 0, 1, cc))
        out[c] = (sxx*seg, sxz*seg, szz*seg, nu)
    return out


def _incident_disp_traction_at_surface(k, w, alpha, beta, rho, s, a, zs):
    """Traction (sxz, szz) at z = 0 of the upgoing incident field, per k (factor exp(i k (x - xs)))."""
    inc = _incident(k, w, alpha, beta, rho, s, a, -1.0)
    txz = 0
    tzz = 0
    for c in ("a", "b"):
        sxx, sxz, szz, nu = inc[c]
        ph = np.exp(-nu*zs)                    # |Z| = zs at the surface
        txz = txz + sxz*ph
        tzz = tzz + szz*ph
    return txz, tzz


def spectrum(xr, zr, xs, zs, s, a, w, alpha, beta, rho, L, M, part="total"):
    """Frequency-domain stress (sxx, sxz, szz) at the receiver for unit slip (not rate), at one
    complex frequency w, wavenumbers k_m = 2 pi m / L, |m| <= M (sources repeated every L)."""
    lam, mu, ka, kb, na, nb = _fields(np.zeros(1), w, alpha, beta, rho)
    k = 2*np.pi*np.arange(-M, M + 1)/L
    lam, mu, ka, kb, na, nb = _fields(k, w, alpha, beta, rho)
    ph_x = np.exp(1j*k*(xr - xs))
    res = np.zeros(3, complex)
    if part in ("full", "total"):
        Zs = 1.0 if zr > zs else -1.0
        inc = _incident(k, w, alpha, beta, rho, s, a, Zs)
        for c in ("a", "b"):
            sxx, sxz, szz, nu = inc[c]
            e = np.exp(-nu*abs(zr - zs))*ph_x
            res += np.array([np.sum(sxx*e), np.sum(sxz*e), np.sum(szz*e)])/L
    if part in ("refl", "total"):
        txz, tzz = _incident_disp_traction_at_surface(k, w, alpha, beta, rho, s, a, zs)
        # reflected downgoing P (potential A) and S (potential B): displacements
        # P: (i k, -na) A exp(-na z); S: (nb, i k) B exp(-nb z)
        # tractions at z = 0: sxz = mu (d uz/dx + d ux/dz), szz = lam div + 2 mu d uz/dz
        PxzA = mu*(1j*k*(-na) + (-na)*1j*k)
        PzzA = lam*(1j*k*1j*k + na**2) + 2*mu*na**2
        PxzB = mu*(1j*k*1j*k - nb**2)
        PzzB = lam*(1j*k*nb - nb*1j*k) + 2*mu*(-nb*1j*k)
        det = PxzA*PzzB - PxzB*PzzA
        A = (-txz*PzzB + tzz*PxzB)/det
        B = (-tzz*PxzA + txz*PzzA)/det
        # reflected stresses at depth zr
        ea = np.exp(-na*zr)
        eb = np.exp(-nb*zr)
        uxA, uzA = 1j*k*A, -na*A
        uxB, uzB = nb*B, 1j*k*B
        sxx = ((lam + 2*mu)*(1j*k*uxA) + lam*(-na*uzA))*ea + ((lam + 2*mu)*(1j*k*uxB) + lam*(-nb*uzB))*eb
        szz = (lam*(1j*k*uxA) + (lam + 2*mu)*(-na*uzA))*ea + (lam*(1j*k*uxB) + (lam + 2*mu)*(-nb*uzB))*eb
        sxz = mu*((-na*uxA + 1j*k*uzA)*ea + (-nb*uxB + 1j*k*uzB)*eb)
        res += np.array([np.sum(sxx*ph_x), np.sum(sxz*ph_x), np.sum(szz*ph_x)])/L
    return res


def ramp_stress(xr, zr, xs, zs, s, a, times, alpha, beta, rho, part="total", T=None, N=1024, L=None, M=None, wI_T=6.0):
    """Time-domain stress for unit slip RATE since t = 0 (slip ramp), at the given times."""
    times = np.asarray(times, float)
    if T is None:
        T = 2.5*times.max()
    if L is None:
        L = 2*alpha*T + 4*(abs(xr - xs) + a)
    if M is None:
        M = int(L/(0.25*min(a, 1e9))) + 1
    wI = wI_T/T
    dw = 2*np.pi/T
    out = np.zeros((len(times), 3))
    acc = np.zeros((len(times), 3), complex)
    for j in range(N//2 + 1):
        w = j*dw + 1j*wI
        S = spectrum(xr, zr, xs, zs, s, a, w, alpha, beta, rho, L, M, part)*(-1.0/w**2)
        fac = 1.0 if j == 0 else 2.0
        acc += fac*np.outer(np.exp(-1j*(j*dw)*times), S)
    out = (np.exp(wI*times)[:, None]*acc.real)*dw/(2*np.pi)
    return out
