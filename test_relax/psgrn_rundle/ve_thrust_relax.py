#!/usr/bin/env python3
"""ve_thrust_relax.py: postseismic surface deformation of a 2-D (plane
strain) dipping fault in an elastic layer over a Maxwell half-space,
with or without gravity.

Problem.  Elastic layer 0 < z < H (z positive DOWN, free surface at
z = 0) over a half-space z > H that is elastic in bulk and Maxwell in
shear (relaxation time tM = eta/mu2).  A fault dipping at angle `dip`
towards +x from the surface trace at x = 0 (or from depth `top`) to
depth `extent` * H carries prescribed slip (uniform or tapered) applied
as a step at t = 0; thrust sense = hanging wall (the +x side at the
surface) moves up-dip.  Gravity enters as buoyancy at the two density
contrasts: sigma_zz(0) = rho1 g u_z(0) at the free surface and
[sigma_zz] = (rho2 - rho1) g u_z at the layer base (the standard
interface-buoyancy approximation; no self-gravitation, no advection of
pre-stress inside the layers).

Method.  Fourier transform in x; in z the exact exponential basis of
the plane-strain Navier system (eigenvalues +-|k|, Jordan chains, so
only decaying exponentials are ever formed and the system is well
conditioned at every k); the fault is a line of double couples whose
k-domain source is a jump of the displacement-stress vector; the
Maxwell half-space enters through mu2(s) = mu2 s tM / (1 + s tM) with
elastic bulk modulus (correspondence principle) and the time domain is
recovered with the fixed Talbot inversion.  The response is computed as
the homogeneous elastic half-space response (closed form in x per
source, exact; the step at the trace is resolved to the source
spacing) plus the layered/viscoelastic difference, which decays like
exp(-k H (2 - extent)) and is transformed on a graded k grid with a
Filon-trapezoid rule.

Usage: ve_thrust_relax.py [options]; -h lists them.  Lengths in km,
moduli in GPa, densities in kg/m3, slip in m, times in Maxwell times.
Output: <out>.npz with x, times and the displacement arrays, <out>.txt
with a far-field table, <out>_near.png, <out>_far.png.

Conventions in the output: x positive towards the hanging wall, u_x
positive towards +x, u_z positive UP (the code's internal z is down).
Displacements are given in units of the maximum slip.
"""
import sys, argparse
import numpy as np

# ----------------------------------------------------------------- basis
def basis(kappa, sgn, lam, mu):
    """Exponential basis of the plane-strain system for the state vector
    Y = (Ux, Uz, Sxz, Szz) with Fourier convention d/dx -> i k,
    k = sgn * kappa.  Returns (v1, v2, w1, w2):
      Y = e^{-kappa (z-za)} [a1 v1 + a2 (v2 + (z-za) v1)]       decaying down
        + e^{+kappa (z-zb)} [b1 w1 + b2 (w2 + (z-zb) w1)]       growing down
    i.e. A v1 = -kappa v1, (A + kappa) v2 = v1, A w1 = kappa w1,
    (A - kappa) w2 = w1 (verified against A in the test suite)."""
    l2m = lam + 2 * mu
    s = 1j * sgn
    v1 = np.array([-1 / (2 * kappa * mu), -s / (2 * kappa * mu), 1.0, s])
    v2 = np.array([1 / (2 * kappa**2 * (lam + mu)),
                   -s * l2m / (2 * kappa**2 * mu * (lam + mu)), 0.0, s / kappa])
    w1 = np.array([1 / (2 * kappa * mu), -s / (2 * kappa * mu), 1.0, -s])
    w2 = np.array([1 / (2 * kappa**2 * (lam + mu)),
                   s * l2m / (2 * kappa**2 * mu * (lam + mu)), 0.0, s / kappa])
    return v1, v2, w1, w2


def A_matrix(k, lam, mu):
    """dY/dz = A Y for the plane-strain system (used only for checks)."""
    ik = 1j * k
    l2m = lam + 2 * mu
    A = np.zeros((4, 4), complex)
    A[0, 1] = -ik; A[0, 2] = 1 / mu
    A[1, 0] = -ik * lam / l2m; A[1, 3] = 1 / l2m
    A[2, 0] = k * k * 4 * mu * (lam + mu) / l2m; A[2, 3] = -ik * lam / l2m
    A[3, 2] = -ik
    return A


def source_jump(k, xs, mxx, mxz, mzz, lam, mu):
    """Jump (below minus above) of the regular part of Y at the depth of a
    point double couple with moment density (mxx, mxz, mzz) per unit
    length at horizontal position xs: body force f_i = -d_j (m_ij delta),
    whose delta' part is absorbed as Y = Y_reg + a delta, giving
    [Y_reg] = b + A a with a = (0,0,mxz,mzz), b = (0,0,ik mxx, ik mxz)."""
    e = np.exp(-1j * k * xs)
    l2m = lam + 2 * mu
    return np.array([mxz / mu, mzz / l2m,
                     1j * k * (mxx - lam * mzz / l2m), 0.0]) * e


# ------------------------------------------------------------ one k, all s
class LayerSolver:
    """Surface response U(k; mu2) = (Ux, Uz) at z = 0 for one wavenumber.
    The layer part (independent of the half-space modulus) is reduced
    once to a particular solution plus a 2-d null space; each half-space
    modulus then costs a 4x4 solve."""

    def __init__(self, k, nodes, jumps, lam1, mu1, g1, dg, H):
        # nodes: sorted depths of the source planes (0 < z < H); jumps:
        # list of 4-vectors, one per node; g1 = rho1 g (dimensionless),
        # dg = (rho2 - rho1) g
        kappa, sgn = abs(k), np.sign(k)
        self.k, self.kappa, self.sgn = k, kappa, sgn
        v1, v2, w1, w2 = basis(kappa, sgn, lam1, mu1)
        zb = np.r_[0.0, nodes, H]
        M = len(zb) - 1                       # sub-layers
        n = 4 * M
        # unknown ordering per sub-layer j: a1, a2, b1, b2 at 4j..4j+3
        def top(j):                            # Y at z_j^+ as 4x4 block
            dz = zb[j + 1] - zb[j]
            e = np.exp(-kappa * dz)
            return np.column_stack([v1, v2, e * w1, e * (w2 - dz * w1)])
        def bot(j):                            # Y at z_{j+1}^- as 4x4 block
            dz = zb[j + 1] - zb[j]
            e = np.exp(-kappa * dz)
            return np.column_stack([e * v1, e * (v2 + dz * v1), w1, w2])
        E = np.zeros((n - 2, n), complex)
        rhs = np.zeros(n - 2, complex)
        T0 = top(0)
        E[0, 0:4] = T0[2]                      # Sxz(0) = 0
        E[1, 0:4] = T0[3] - g1 * T0[1]         # Szz(0) - rho1 g Uz(0) = 0
        r = 2
        for j in range(1, M):
            E[r:r + 4, 4 * j:4 * j + 4] = top(j)
            E[r:r + 4, 4 * (j - 1):4 * j] = -bot(j - 1)
            rhs[r:r + 4] = jumps[j - 1]
            r += 4
        # particular solution and null space (E has full row rank n-2):
        # E^H = Q R, null space = last two columns of Q, particular
        # solution from R^H y = rhs
        Q, R = np.linalg.qr(E.conj().T, mode='complete')
        self.N = Q[:, n - 2:]                              # (n, 2)
        y = np.linalg.solve(R[:n - 2, :].conj().T, rhs)
        xp = Q[:, :n - 2] @ y
        self.xp = xp
        B = np.zeros((4, n), complex)
        B[:, 4 * (M - 1):] = bot(M - 1)       # Y(H^-)
        IG = np.eye(4); IG[3, 1] = dg          # adds (rho2-rho1) g Uz to Szz
        self.BN = IG @ B @ self.N              # (4,2)
        self.Bp = IG @ B @ xp                  # (4,)
        self.T0N = T0[:2] @ self.N[0:4]        # surface displacement map
        self.T0p = T0[:2] @ xp[0:4]

    def surface(self, lam2, mu2):
        """(Ux, Uz) at z = 0 for half-space moduli lam2, mu2 (arrays of any
        shape, complex allowed)."""
        lam2 = np.asarray(lam2, complex); mu2 = np.asarray(mu2, complex)
        shp = np.broadcast(lam2, mu2).shape
        lam2 = np.broadcast_to(lam2, shp).ravel(); mu2 = np.broadcast_to(mu2, shp).ravel()
        kappa, s = self.kappa, 1j * self.sgn
        l2m = lam2 + 2 * mu2
        # decaying basis of the half-space, columns normalised
        V1 = np.stack([-1 / (2 * kappa * mu2), -s / (2 * kappa * mu2),
                       np.ones_like(mu2), s * np.ones_like(mu2)], -1)
        V2 = np.stack([1 / (2 * kappa**2 * (lam2 + mu2)),
                       -s * l2m / (2 * kappa**2 * mu2 * (lam2 + mu2)),
                       np.zeros_like(mu2), s / kappa * np.ones_like(mu2)], -1)
        V1 /= np.linalg.norm(V1, axis=-1, keepdims=True)
        V2 /= np.linalg.norm(V2, axis=-1, keepdims=True)
        # interface: (I+G) Y(H^-) = c1 V1 + c2 V2  ->  BN y - [V1 V2] c = -Bp
        Mx = np.zeros(shp + (4, 4), complex).reshape(-1, 4, 4)
        Mx[:, :, 0:2] = self.BN
        Mx[:, :, 2] = -V1
        Mx[:, :, 3] = -V2
        sol = np.linalg.solve(Mx, np.broadcast_to(-self.Bp, (len(mu2), 4))[..., None])[..., 0]
        y = sol[:, :2]
        u = self.T0p[None, :] + y @ self.T0N.T                # (n, 2)
        return u.reshape(shp + (2,))


# ------------------------------------------------------------------ Talbot
def talbot_weights(t, M):
    """Fixed Talbot nodes s_k and weights w_k (Abate & Valko 2004) such
    that f(t) = sum_k Re(w_k F_r(s_k)) for a REAL f with transform F_r."""
    r = 2.0 * M / (5.0 * t)
    th = np.arange(1, M) * np.pi / M
    cot = np.cos(th) / np.sin(th)
    s = r * th * (cot + 1j)
    sig = th + (th * cot - 1) * cot
    w = (r / M) * np.exp(t * s) * (1 + 1j * sig)
    return np.r_[r, s], np.r_[0.5 * (r / M) * np.exp(r * t), w]


def laplace_invert_complex(F, s_nodes, w, idx_conj):
    """Inverse transform of a COMPLEX-valued time function from samples
    F(s) at the Talbot nodes and at their conjugates: F_r(s) = (F(s) +
    conj F(conj s))/2, F_i(s) = (F(s) - conj F(conj s))/(2i), each a
    transform of a real function, then f = f_r + i f_i."""
    Fs, Fc = F[..., :len(s_nodes)], np.conj(F[..., idx_conj])
    Fr, Fi = 0.5 * (Fs + Fc), (Fs - Fc) / 2j
    fr = np.sum((w * Fr).real, axis=-1)
    fi = np.sum((w * Fi).real, axis=-1)
    return fr + 1j * fi


# ---------------------------------------------------------------- the model
def fault_sources(dip_deg, extent, top, taper, taper_width, nsrc, H=1.0):
    """Point double couples along the fault: positions, depths and moment
    densities per unit along-strike length, in units mu1 = 1, slip = 1.
    Returns xs, zs, (mxx, mxz, mzz) arrays and the slip profile."""
    d = np.deg2rad(dip_deg)
    z0, z1 = top * H, extent * H
    L0, L1 = z0 / np.sin(d), z1 / np.sin(d)            # down-dip distances
    xi = L0 + (np.arange(nsrc) + 0.5) * (L1 - L0) / nsrc
    dl = (L1 - L0) / nsrc
    if taper == 'uniform':
        s = np.ones(nsrc)
    elif taper == 'cos':                                # cosine taper over the bottom fraction
        s = np.ones(nsrc)
        wl = taper_width * (L1 - L0)
        m = xi > L1 - wl
        s[m] = 0.5 * (1 + np.cos(np.pi * (xi[m] - (L1 - wl)) / wl))
    elif taper == 'sin':                                # half sine over the whole fault
        s = np.sin(np.pi * (xi - L0) / (L1 - L0))
    else:
        raise ValueError(taper)
    xs, zs = xi * np.cos(d), xi * np.sin(d)
    # slip vector of the hanging wall relative to the footwall, thrust
    # sense: up-dip = (-cos d, -sin d) in (x, z-down); normal into the
    # hanging wall n = (sin d, -cos d); m = mu (b n^T + n b^T) dl
    mxx = -np.sin(2 * d) * s * dl
    mzz = np.sin(2 * d) * s * dl
    mxz = np.cos(2 * d) * s * dl
    return xs, zs, mxx, mxz, mzz, xi, s


def group_by_depth(zs, jumps):
    """Merge sources at the same depth (e.g. a horizontal fault)."""
    order = np.argsort(zs)
    nodes, J = [], []
    for i in order:
        if nodes and abs(zs[i] - nodes[-1]) < 1e-12:
            J[-1] = J[-1] + jumps[i]
        else:
            nodes.append(zs[i]); J.append(jumps[i].copy())
    return np.array(nodes), J


def surface_response(kgrid, src, lam1, mu1, lam2_of_s, mu2_of_s, s_nodes,
                     g1, dg, H=1.0):
    """U(k, s) for all k in kgrid and all s in s_nodes (s_nodes may be the
    single value None for a purely elastic evaluation with lam2_of_s,
    mu2_of_s constants).  Returns array (nk, ns, 2)."""
    xs, zs, mxx, mxz, mzz = src
    if s_nodes is None:
        lam2, mu2 = np.array([lam2_of_s]), np.array([mu2_of_s])
    else:
        lam2, mu2 = lam2_of_s(s_nodes), mu2_of_s(s_nodes)
    out = np.zeros((len(kgrid), len(mu2), 2), complex)
    for i, k in enumerate(kgrid):
        jumps = [source_jump(k, xs[j], mxx[j], mxz[j], mzz[j], lam1, mu1)
                 for j in range(len(xs))]
        nodes, J = group_by_depth(zs, jumps)
        ls = LayerSolver(k, nodes, J, lam1, mu1, g1, dg, H)
        out[i] = ls.surface(lam2, mu2)
    return out


def layered_response(kgrid, src, lam1, mu1, lam2_of_s, mu2_of_s, s_nodes, g1, dg,
                     H=1.0, chunk=256):
    """Vectorised (over k) version of surface_response: U(k, s) of shape
    (nk, ns, 2) for the layer over the half-space with moduli
    lam2_of_s(s), mu2_of_s(s) (s_nodes None: constants lam2_of_s,
    mu2_of_s).  Layer reduction by batched QR, half-space by batched 4x4
    solves."""
    if len(kgrid) > chunk:
        return np.concatenate([layered_response(kgrid[i:i + chunk], src, lam1, mu1, lam2_of_s,
                                                mu2_of_s, s_nodes, g1, dg, H, chunk)
                               for i in range(0, len(kgrid), chunk)])
    xs, zs, mxx, mxz, mzz = src
    if s_nodes is None:
        lam2, mu2 = np.array([lam2_of_s], complex), np.array([mu2_of_s], complex)
    else:
        lam2, mu2 = np.asarray(lam2_of_s(s_nodes), complex), np.asarray(mu2_of_s(s_nodes), complex)
    nk, ns = len(kgrid), len(mu2)
    k = kgrid
    kappa, sgn = np.abs(k), np.sign(k)
    s = 1j * sgn
    l2m = lam1 + 2 * mu1
    one = np.ones(nk)
    v1 = np.stack([-1 / (2 * kappa * mu1), -s / (2 * kappa * mu1), one, s], -1)
    v2 = np.stack([1 / (2 * kappa**2 * (lam1 + mu1)),
                   -s * l2m / (2 * kappa**2 * mu1 * (lam1 + mu1)), 0 * one, s / kappa], -1)
    w1 = np.stack([1 / (2 * kappa * mu1), -s / (2 * kappa * mu1), one, -s], -1)
    w2 = np.stack([1 / (2 * kappa**2 * (lam1 + mu1)),
                   s * l2m / (2 * kappa**2 * mu1 * (lam1 + mu1)), 0 * one, s / kappa], -1)
    # depth nodes and merged jumps
    nodes, inv = np.unique(np.round(zs, 12), return_inverse=True)
    ph = np.exp(-1j * np.outer(k, xs))                       # (nk, nsrc)
    Jsrc = np.stack([mxz / mu1 * ph, mzz / l2m * ph,
                     1j * k[:, None] * (mxx - lam1 * mzz / l2m) * ph, 0 * ph], -1)   # (nk, nsrc, 4)
    J = np.zeros((nk, len(nodes), 4), complex)
    for j in range(len(xs)):
        J[:, inv[j]] += Jsrc[:, j]
    zb = np.r_[0.0, nodes, H]
    M = len(zb) - 1
    n = 4 * M
    dz = np.diff(zb)                                          # (M,)
    e = np.exp(-np.outer(kappa, dz))                          # (nk, M)
    def top(j):
        return np.stack([v1, v2, e[:, j, None] * w1, e[:, j, None] * (w2 - dz[j] * w1)], -1)
    def bot(j):
        return np.stack([e[:, j, None] * v1, e[:, j, None] * (v2 + dz[j] * v1), w1, w2], -1)
    E = np.zeros((nk, n - 2, n), complex)
    rhs = np.zeros((nk, n - 2), complex)
    T0 = top(0)
    E[:, 0, 0:4] = T0[:, 2, :]
    E[:, 1, 0:4] = T0[:, 3, :] - g1 * T0[:, 1, :]
    r = 2
    for j in range(1, M):
        E[:, r:r + 4, 4 * j:4 * j + 4] = top(j)
        E[:, r:r + 4, 4 * (j - 1):4 * j] = -bot(j - 1)
        rhs[:, r:r + 4] = J[:, j - 1]
        r += 4
    Q, R = np.linalg.qr(np.conj(np.transpose(E, (0, 2, 1))), mode='complete')   # E^H = Q R
    N = Q[:, :, n - 2:]                                                           # (nk, n, 2)
    y = np.linalg.solve(np.conj(np.transpose(R[:, :n - 2, :], (0, 2, 1))), rhs[..., None])[..., 0]
    xp = np.einsum('kij,kj->ki', Q[:, :, :n - 2], y)                            # (nk, n)
    Bm = bot(M - 1)                                                               # (nk, 4, 4)
    IG = np.eye(4); IG[3, 1] = dg
    BN = np.einsum('ab,kbc,kcd->kad', IG, Bm, N[:, 4 * (M - 1):, :])            # (nk, 4, 2)
    Bp = np.einsum('ab,kbc,kc->ka', IG, Bm, xp[:, 4 * (M - 1):])                # (nk, 4)
    T0N = np.einsum('kab,kbc->kac', T0[:, :2, :], N[:, 0:4, :])                 # (nk, 2, 2)
    T0p = np.einsum('kab,kb->ka', T0[:, :2, :], xp[:, 0:4])                      # (nk, 2)
    # half-space decaying basis for every (k, s)
    ka = kappa[:, None]; sg = s[:, None]
    l2m2 = (lam2 + 2 * mu2)[None, :]
    mu2b, lam2b = mu2[None, :], lam2[None, :]
    on = np.ones((nk, ns))
    V1 = np.stack([-1 / (2 * ka * mu2b), -sg / (2 * ka * mu2b), on, sg * on], -1)
    V2 = np.stack([1 / (2 * ka**2 * (lam2b + mu2b)),
                   -sg * l2m2 / (2 * ka**2 * mu2b * (lam2b + mu2b)), 0 * on, sg / ka * on], -1)
    V1 /= np.linalg.norm(V1, axis=-1, keepdims=True)
    V2 /= np.linalg.norm(V2, axis=-1, keepdims=True)
    Mx = np.zeros((nk, ns, 4, 4), complex)
    Mx[:, :, :, 0:2] = BN[:, None, :, :]
    Mx[:, :, :, 2] = -V1
    Mx[:, :, :, 3] = -V2
    sol = np.linalg.solve(Mx, np.broadcast_to(-Bp[:, None, :], (nk, ns, 4))[..., None])[..., 0]
    yy = sol[..., :2]                                                             # (nk, ns, 2)
    return T0p[:, None, :] + np.einsum('kab,ksb->ksa', T0N, yy)


def halfspace_response(kgrid, src, lam, mu):
    """Surface response (nk, 2) of the homogeneous elastic half-space
    without gravity, closed form: for a point double couple (mxx, mxz,
    mzz) at depth d, with xi = |k| d, E = e^{-xi}, sg = sign k,
      mu Ux = E [ i sg mxx (xi/2 - (lam+2mu)/(2(lam+mu))) + mxz (xi - 1)
                  + i sg mzz (lam/(2(lam+mu)) - xi/2) ]
      mu Uz = E [ mxx (xi/2 - mu/(2(lam+mu))) - i sg mxz xi
                  - mzz (mu/(2(lam+mu)) + xi/2) ]
    times e^{-ik xs}.  Obtained from the 6x6 layer/half-space system
    (halfspace_response_solve) by noting that U e^{xi} is linear in xi
    and reading off the coefficients; both agree to roundoff."""
    chunk = 20000
    if len(kgrid) > chunk:
        return np.concatenate([halfspace_response(kgrid[i:i + chunk], src, lam, mu)
                               for i in range(0, len(kgrid), chunk)])
    xs, zs, mxx, mxz, mzz = src
    k = kgrid[:, None]
    kappa, sg = np.abs(k), np.sign(k)
    xi = kappa * zs[None, :]
    E = np.exp(-xi) * np.exp(-1j * k * xs[None, :])
    q = 1.0 / (2 * (lam + mu))
    Ux = E * (1j * sg * mxx * (xi / 2 - (lam + 2 * mu) * q) + mxz * (xi - 1)
              + 1j * sg * mzz * (lam * q - xi / 2))
    Uz = E * (mxx * (xi / 2 - mu * q) - 1j * sg * mxz * xi - mzz * (mu * q + xi / 2))
    return np.stack([Ux.sum(1), Uz.sum(1)], -1) / mu


def halfspace_surface_x(x, src, lam, mu):
    """Surface displacements (nx, 2) of the homogeneous elastic half-space
    directly in x, the inverse transform of halfspace_response done
    analytically: with X = x - xs, d the source depth and r2 = X^2 + d^2,
      e^{-|k| d}            <->  f0 = d / (pi r2)
      |k| d e^{-|k| d}      <->  f1 = d (d^2 - X^2) / (pi r2^2)
      i sgn(k) e^{-|k| d}   <->  g0 = -X / (pi r2)
      i sgn(k) |k| d e^{-|k| d} <-> g1 = -2 d^2 X / (pi r2^2)
    and the coefficients of halfspace_response.  Internal z-down
    convention (second component positive down)."""
    xs, zs, mxx, mxz, mzz = src
    X = x[:, None] - xs[None, :]
    d = zs[None, :]
    r2 = X**2 + d**2
    f0 = d / (np.pi * r2)
    f1 = d * (d**2 - X**2) / (np.pi * r2**2)
    g0 = -X / (np.pi * r2)
    g1 = -2 * d**2 * X / (np.pi * r2**2)
    q = 1.0 / (2 * (lam + mu))
    ux = mxx * (g1 / 2 - (lam + 2 * mu) * q * g0) + mxz * (f1 - f0) + mzz * (lam * q * g0 - g1 / 2)
    uz = mxx * (f1 / 2 - mu * q * f0) - mxz * g1 - mzz * (mu * q * f0 + f1 / 2)
    return np.stack([ux.sum(1), uz.sum(1)], -1) / mu


def halfspace_response_solve(kgrid, src, lam, mu, chunk=1500):
    """Fast path for the homogeneous elastic half-space without gravity:
    surface response (nk, 2) as the sum over sources of 6x6 systems
    (sub-layer above the source, decaying half-space below), vectorised
    over k and sources; chunked over k to bound memory."""
    if len(kgrid) > chunk:
        return np.concatenate([halfspace_response_solve(kgrid[i:i + chunk], src, lam, mu, chunk)
                               for i in range(0, len(kgrid), chunk)])
    xs, zs, mxx, mxz, mzz = src
    nk, ns = len(kgrid), len(xs)
    k = kgrid[:, None]
    kappa, sgn = np.abs(k), np.sign(k)
    s = 1j * sgn
    l2m = lam + 2 * mu
    one = np.ones((nk, 1))
    # basis vectors as (nk, 1, 4)
    v1 = np.stack([-1 / (2 * kappa * mu) * one, -s / (2 * kappa * mu) * one, one, s * one], -1)
    v2 = np.stack([1 / (2 * kappa**2 * (lam + mu)) * one,
                   -s * l2m / (2 * kappa**2 * mu * (lam + mu)) * one, 0 * one, s / kappa * one], -1)
    w1 = np.stack([1 / (2 * kappa * mu) * one, -s / (2 * kappa * mu) * one, one, -s * one], -1)
    w2 = np.stack([1 / (2 * kappa**2 * (lam + mu)) * one,
                   s * l2m / (2 * kappa**2 * mu * (lam + mu)) * one, 0 * one, s / kappa * one], -1)
    z = zs[None, :, None]                                   # (1, ns, 1)
    e = np.exp(-kappa[:, :, None] * z)                      # (nk, ns, 1)
    # unknowns (a1, a2, b1, b2, c1, c2); rows: Sxz(0), Szz(0), jump (4)
    M = np.zeros((nk, ns, 6, 6), complex)
    top = [v1, v2, e * w1, e * (w2 - z * w1)]                # Y(0) columns
    bot = [e * v1, e * (v2 + z * v1), w1, w2]                # Y(zs^-) columns
    for c in range(4):
        M[:, :, 0, c] = np.broadcast_to(top[c][..., 2], (nk, ns))
        M[:, :, 1, c] = np.broadcast_to(top[c][..., 3], (nk, ns))
        M[:, :, 2:, c] = -np.broadcast_to(bot[c], (nk, ns, 4))
    M[:, :, 2:, 4] = np.broadcast_to(v1, (nk, ns, 4))
    M[:, :, 2:, 5] = np.broadcast_to(v2, (nk, ns, 4))
    rhs = np.zeros((nk, ns, 6), complex)
    ph = np.exp(-1j * k * xs[None, :])                        # (nk, ns)
    rhs[:, :, 2] = mxz / mu * ph
    rhs[:, :, 3] = mzz / l2m * ph
    rhs[:, :, 4] = 1j * k * (mxx - lam * mzz / l2m) * ph
    sol = np.linalg.solve(M, rhs[..., None])[..., 0]          # (nk, ns, 6)
    U = np.zeros((nk, ns, 2), complex)
    for c in range(4):
        U += sol[:, :, c, None] * np.broadcast_to(top[c][..., :2], (nk, ns, 2))
    return U.sum(1)


def split_offset(kgrid, D):
    """The layered/viscoelastic horizontal field tends to a constant
    +-C far from the fault (the plate floats on the relaxed substrate),
    which the periodic Fourier sum would alias.  Estimate C from the two
    smallest wavenumbers (ik D_x e^{k} -> 2C) and remove 2C e^{-k}/(ik),
    the transform of C (2/pi) arctan(x/H), from D_x; the caller adds that
    function back after the inverse transform.  D: (nk, ..., 2)."""
    k0, k1 = kgrid[0], kgrid[1]
    y0 = 1j * k0 * D[0, ..., 0] * np.exp(k0) / 2
    y1 = 1j * k1 * D[1, ..., 0] * np.exp(k1) / 2
    C = (y0 - (y1 - y0) * k0 / (k1 - k0)).real
    Dr = D.copy()
    ek = np.exp(-kgrid)[(slice(None),) + (None,) * C.ndim]
    kk = kgrid[(slice(None),) + (None,) * C.ndim]
    Dr[..., 0] -= 2 * C[None, ...] * ek / (1j * kk)
    # the odd 1/x tail: D_x -> -i pi A sgn(k) as k -> 0, transform of
    # A x/(x^2 + H^2) is -i pi A sgn(k) e^{-|k|}; same extrapolation
    y0 = (Dr[0, ..., 0] * np.exp(k0) / (-1j * np.pi))
    y1 = (Dr[1, ..., 0] * np.exp(k1) / (-1j * np.pi))
    A = (y0 - (y1 - y0) * k0 / (k1 - k0)).real
    Dr[..., 0] -= -1j * np.pi * A[None, ...] * ek
    return Dr, C, A


def kgrid_graded(kmax, dk, kbreak=0.2, nsmall=400, kmin=1e-5):
    """Wavenumber grid for the inverse transform: dense and linear from
    kmin to kbreak (resolves the slow long-wavelength relaxation, whose
    structure sits at k H << 1 at late times), then uniform with spacing
    dk = pi / xmax (resolves e^{ikx} out to xmax).  Returns nodes and
    trapezoid weights; the integrand is finite at k -> 0 after the
    far-field terms have been removed (split_offset)."""
    ks = np.linspace(kmin, kbreak, nsmall, endpoint=False)
    ku = np.arange(kbreak, kmax, dk)
    k = np.r_[ks, ku]
    w = np.zeros_like(k)
    w[1:-1] = 0.5 * (k[2:] - k[:-2]); w[0] = 0.5 * (k[1] - k[0]); w[-1] = 0.5 * (k[-1] - k[-2])
    return k, w


def inverse_ft(kgrid, wk, U, x):
    """u(x) = (1/pi) Re int_0^inf U(k) e^{ikx} dk by a Filon-trapezoid
    rule: U is interpolated linearly between the nodes of kgrid and the
    oscillatory factor is integrated exactly on every interval, so the
    error depends on the smoothness of U at the grid scale and not on
    x (the plain trapezoid rule needs many nodes per period of e^{ikx},
    which at |x| near xmax produced 5e-4 wiggles).  wk is unused (kept
    for the call signature).  U: (nk, ...) complex; returns real
    array (len(x), ...)."""
    out = np.zeros((len(x),) + U.shape[1:])
    h = np.diff(kgrid)                                   # (nk-1,)
    step = 2000
    xx = x[:, None]
    for i0 in range(0, len(h), step):
        i1 = min(i0 + step, len(h))
        th = xx * h[None, i0:i1]                         # (nx, m)
        small = np.abs(th) < 1e-3
        ths = np.where(small, 1.0, th)
        eth = np.exp(1j * ths)
        I0 = (eth - 1) / (1j * ths)
        I1 = (eth - I0) / (1j * ths)
        W1 = np.where(small, 0.5 + 1j * th / 3, I1)
        W0 = np.where(small, 0.5 + 1j * th / 6, I0 - I1)
        ph = np.exp(1j * xx * kgrid[None, i0:i1]) * h[None, i0:i1]
        out += np.tensordot(ph * W0, U[i0:i1], axes=(1, 0)).real
        out += np.tensordot(ph * W1, U[i0 + 1:i1 + 1], axes=(1, 0)).real
    return out / np.pi


def run(p):
    H = p.H * 1e3
    mu1 = p.mu1 * 1e9; nu = p.nu
    lam1 = 2 * mu1 * nu / (1 - 2 * nu)
    mu2 = p.mu2_ratio * mu1; lam2 = lam1 * p.mu2_ratio     # same Poisson ratio
    K2 = lam2 + 2 * mu2 / 3
    # dimensionless: lengths / H, moduli / mu1, displacements / slip
    l1, m1 = lam1 / mu1, 1.0
    m2e, K2n = mu2 / mu1, K2 / mu1
    g1 = p.rho1 * p.g * H / mu1 if p.gravity else 0.0
    dg = (p.rho2 - p.rho1) * p.g * H / mu1 if p.gravity else 0.0
    xs, zs, mxx, mxz, mzz, xi, sprof = fault_sources(p.dip, p.extent, p.top,
                                                     p.taper, p.taper_width, p.nsrc)
    src = (xs, zs, mxx, mxz, mzz)                      # fine: homogeneous part
    srcc = fault_sources(p.dip, p.extent, p.top, p.taper, p.taper_width, p.nsrc_layer)[:5]
    # Maxwell half-space in the Laplace domain (time in Maxwell times)
    mu2_of_s = lambda s: m2e * s / (1.0 + s)
    if p.lam_const:                                  # Rundle (1982) convention
        lam2_of_s = lambda s: l1 * p.mu2_ratio + 0.0 * s
    else:                                            # elastic bulk modulus
        lam2_of_s = lambda s: K2n - 2.0 / 3.0 * mu2_of_s(s)
    # grids
    xmax = p.xmax                                   # in H
    dk = np.pi / xmax
    kmax_c = p.kmax_coarse                          # difference part
    # the layered/viscoelastic part has structure in k at the scale of
    # the relaxation front (0.1/H at 50 tM), finer than pi/xmax: oversample
    kc, wc = kgrid_graded(kmax_c, dk / p.oversample, nsmall=p.nsmall)
    x = np.linspace(-xmax, xmax, p.nx)
    times = [t for t in p.times if t > 0]
    print(f'# H {p.H:g} km, dip {p.dip:g}, extent {p.extent:g} H (bottom {p.extent*p.H:g} km), '
          f'top {p.top:g} H, slip {p.taper}, nsrc {p.nsrc}/{p.nsrc_layer}, gravity {"on" if p.gravity else "off"}, '
          f'mu2/mu1 {p.mu2_ratio:g}, nu {nu:g}')
    print(f'# layered k grid: {len(kc)} nodes to {kmax_c:g}/H; '
          f'x to +-{xmax:g} H, {p.nx} points; Talbot M {p.talbot}')
    # 2. layered elastic (t = 0+) and relaxed (mu2 -> 0) on the coarse grid
    Ue = layered_response(kc, srcc, l1, m1, l1 * p.mu2_ratio, m2e, None, g1, dg)[:, 0, :]
    lam_rel = l1 * p.mu2_ratio if p.lam_const else K2n - 2 / 3 * m2e * 1e-9
    Ur = layered_response(kc, srcc, l1, m1, lam_rel, m2e * 1e-9, None, g1, dg)[:, 0, :]
    Uh_c = halfspace_response(kc, srcc, l1, m1)
    # 3. viscoelastic times through Talbot (nodes and conjugates)
    s_all, w_all, idx_conj, slices = [], [], [], []
    for t in times:
        s, w = talbot_weights(t, p.talbot)
        n0 = len(s_all)
        s_all += list(s) + list(np.conj(s[1:]))
        idx_conj.append(np.r_[n0, n0 + len(s) + np.arange(len(s) - 1)])
        w_all.append(w); slices.append((n0, len(s)))
    s_all = np.array(s_all)
    Us = layered_response(kc, srcc, l1, m1, lam2_of_s, mu2_of_s, s_all, g1, dg)  # (nk, ns, 2)
    Ut = np.zeros((len(kc), len(times), 2), complex)
    for it, t in enumerate(times):
        n0, ns = slices[it]
        F = Us[:, n0:n0 + 2 * ns - 1, :] / s_all[None, n0:n0 + 2 * ns - 1, None]   # step source
        Fm = np.moveaxis(F, 1, -1)                                   # (nk, 2, nodes)
        Ut[:, it, :] = laplace_invert_complex(Fm, s_all[n0:n0 + ns], w_all[it],
                                              idx_conj[it] - n0)
    # 4. assemble in x: homogeneous (fine) + difference (coarse)
    # 1. homogeneous elastic half-space (no gravity, same material): closed form in x
    uh = halfspace_surface_x(x, src, l1, m1)             # (nx, 2)
    arct = (2 / np.pi) * np.arctan(x)
    tail = x / (x**2 + 1.0)
    De, Ce, Ae = split_offset(kc, Ue - Uh_c)
    Dr, Cr, Ar = split_offset(kc, Ur - Uh_c)
    Dt, Ct, At = split_offset(kc, Ut - Uh_c[:, None, :])
    d_el = inverse_ft(kc, wc, De, x); d_el[:, 0] += Ce * arct + Ae * tail
    d_rel = inverse_ft(kc, wc, Dr, x); d_rel[:, 0] += Cr * arct + Ar * tail
    d_t = inverse_ft(kc, wc, Dt, x)                      # (nx, nt, 2)
    d_t[:, :, 0] += Ct[None, :] * arct[:, None] + At[None, :] * tail[:, None]
    print('# far-field horizontal offset C/slip (u_x -> -+C): coseismic %.4f, ' % Ce
          + ' '.join(f't={tt:g}: {c:.4f}' for tt, c in zip(times, Ct)) + f', relaxed {Cr:.4f}')
    u0 = uh + d_el
    uinf = uh + d_rel
    ut = uh[:, None, :] + d_t
    # z down internally -> report uplift
    if not p.gravity:
        # a plate with an internal dislocation on an inviscid substrate
        # without gravity has no bounded static state (it kinks); the
        # mu2 -> 0 solution is then set by the k_min cutoff and is not
        # reported
        uinf[:] = np.nan
        print('# gravity off: no bounded relaxed limit (free plate kinks); relaxed curve not reported')
    res = dict(x=x, times=np.array([0.0] + times), H_km=p.H,
               ux=np.concatenate([u0[:, None, 0], ut[:, :, 0]], 1),
               uz=-np.concatenate([u0[:, None, 1], ut[:, :, 1]], 1),
               ux_relaxed=uinf[:, 0], uz_relaxed=-uinf[:, 1],
               ux_halfspace=uh[:, 0], uz_halfspace=-uh[:, 1],
               offset=np.r_[Ce, Ct], offset_relaxed=Cr,
               xi=xi, slip=sprof, params=str(vars(p)))
    # Talbot sanity: t -> 0 limit of the first (smallest) time vs elastic
    return res


def far_table(res, fh):
    x, t = res['x'], res['times']
    cols = [-20, -10, -5, -2, -1, -0.5, 0.5, 1, 2, 5, 10, 20]
    idx = [np.argmin(abs(x - c)) for c in cols]
    for comp in ('ux', 'uz'):
        fh.write(f'# {comp} / slip at x/H = ' + ' '.join(f'{c:7g}' for c in cols) + '\n')
        for it, tt in enumerate(t):
            fh.write(f'{("t=%g tM" % tt):>12s} ' + ' '.join(f'{res[comp][i, it]:7.4f}' for i in idx) + '\n')
        fh.write(f'{"relaxed":>12s} ' + ' '.join(f'{res[comp + "_relaxed"][i]:7.4f}' for i in idx) + '\n')
        a = np.abs(res[comp])
        fh.write(f'# max |{comp}|: ' + ' '.join(f't={tt:g}: {a[:, it].max():.4f} at x={x[a[:, it].argmax()]:.2f}'
                                               for it, tt in enumerate(t)) + '\n')


def plot(res, out, p):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x, t = res['x'], res['times']
    cmap = plt.get_cmap('viridis')
    title = (f'dip {p.dip:g}, bottom {p.extent:g} H, {p.taper} slip, gravity {"on" if p.gravity else "off"}, '
             f'H = {p.H:g} km, mu2/mu1 = {p.mu2_ratio:g}')
    for tag, xl in (('near', 3.0), ('far', p.xmax)):
        fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        for it, tt in enumerate(t):
            col = 'k' if it == 0 else cmap(it / (len(t) - 1))
            lab = 'coseismic' if it == 0 else f'{tt:g} tM'
            ax[0].plot(x, res['ux'][:, it], color=col, lw=1.2 if it else 1.8, label=lab)
            ax[1].plot(x, res['uz'][:, it], color=col, lw=1.2 if it else 1.8)
        if np.isfinite(res['ux_relaxed']).any():
            ax[0].plot(x, res['ux_relaxed'], 'r:', lw=1, label='relaxed (mu2 -> 0)')
            ax[1].plot(x, res['uz_relaxed'], 'r:', lw=1)
        ax[0].set_ylabel('horizontal u_x / slip  (+ towards hanging wall)')
        ax[1].set_ylabel('vertical u_z / slip  (+ up)')
        ax[1].set_xlabel('x / H  (fault trace at 0, dips towards +x)')
        for a in ax:
            a.axhline(0, color='0.6', lw=0.5); a.axvline(0, color='0.6', lw=0.5); a.grid(alpha=0.3)
            a.set_xlim(-xl, xl)
        if tag == 'far':
            for a in ax:
                a.set_yscale('symlog', linthresh=1e-3)
        ax[0].legend(fontsize=8, ncol=2)
        fig.suptitle(title, fontsize=10)
        fig.tight_layout()
        fig.savefig(f'{out}_{tag}.png', dpi=130)
        plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--H', type=float, default=40.0, help='elastic layer thickness [km]')
    ap.add_argument('--dip', type=float, default=30.0, help='fault dip [deg], dips towards +x')
    ap.add_argument('--extent', type=float, default=0.5, help='bottom depth of the fault / H (<= 1)')
    ap.add_argument('--top', type=float, default=0.0, help='top depth of the fault / H (0: surface breaking)')
    ap.add_argument('--taper', default='uniform', choices=['uniform', 'cos', 'sin'],
                    help='slip profile: uniform, cosine taper at the bottom, half sine')
    ap.add_argument('--taper-width', type=float, default=0.3, help='bottom taper as a fraction of fault length (cos)')
    ap.add_argument('--mu1', type=float, default=32.04, help='layer shear modulus [GPa]')
    ap.add_argument('--mu2-ratio', type=float, default=1.0, help='half-space / layer shear modulus (unrelaxed)')
    ap.add_argument('--nu', type=float, default=0.25, help="Poisson's ratio (both, unrelaxed)")
    ap.add_argument('--rho1', type=float, default=2800.0)
    ap.add_argument('--rho2', type=float, default=3300.0)
    ap.add_argument('--g', type=float, default=9.81)
    ap.add_argument('--gravity', type=int, default=1, help='1: buoyancy at surface and layer base; 0: off')
    ap.add_argument('--times', type=float, nargs='+', default=[0.5, 1, 2, 5, 10, 50], help='in Maxwell times')
    ap.add_argument('--nsrc', type=int, default=96, help='point sources along the fault (homogeneous part)')
    ap.add_argument('--nsrc-layer', type=int, default=24,
                    help='point sources for the layered difference (its discretisation error cancels)')
    ap.add_argument('--nx', type=int, default=4001)
    ap.add_argument('--xmax', type=float, default=50.0, help='half width of the x window / H')
    ap.add_argument('--kmax-coarse', type=float, default=40.0, help='k H cutoff of the layered difference')
    ap.add_argument('--talbot', type=int, default=32)
    ap.add_argument('--lam-const', action='store_true',
                    help='hold lambda of the half-space constant instead of the bulk modulus (Rundle 1982)')
    ap.add_argument('--nsmall', type=int, default=400, help='dense k points below k H = 0.2')
    ap.add_argument('--oversample', type=float, default=4.0, help='k spacing of the layered part = pi/(xmax*oversample)')
    ap.add_argument('--out', default='ve_relax')
    ap.add_argument('--no-plot', action='store_true')
    p = ap.parse_args(argv)
    if p.extent > 1.0 or p.extent <= p.top:
        sys.exit('extent must be in (top, 1]')
    if p.extent == 1.0 and p.gravity == 0:
        print('# note: a fault cutting the whole layer over an inviscid substrate without gravity '
              'has no bounded relaxed limit; results are finite at finite time only', file=sys.stderr)
    res = run(p)
    np.savez(p.out + '.npz', **res)
    with open(p.out + '.txt', 'w') as fh:
        fh.write('# ' + ' '.join(sys.argv) + '\n')
        far_table(res, fh)
    print(open(p.out + '.txt').read())
    if not p.no_plot:
        plot(res, p.out, p)
    return res


if __name__ == '__main__':
    main()
