#!/usr/bin/env python3
"""
inplane2d: standalone prototype for 2-D PLANE-STRAIN (in-plane, mode
I/II) quasi-static deformation of an elastic plate over a Maxwell
half-space, driven by dip-slip dislocations: the kernel generator
needed for thrust-fault-in-cross-section earthquake cycles
(implementation plan step 7b).  NO GRAVITY (see README caveats).

Method: wavenumber-domain (e^{ikx}) state-vector solution
Y = [ux, uz, sxz, szz] with the exact analytic basis of plane-strain
elasticity per uniform region, point-dislocation sources inserted as
state-vector jumps derived from the moment-tensor equivalent-force
representation, and (for the viscoelastic substrate) the
correspondence principle mu2 -> mu2(s) with numerical inverse
Laplace transform (Talbot contour).

Everything here is validated against an independent finite-element
reference (fem2d.py) before being trusted; see run_inplane_proto_test.
"""
import numpy as np
from scipy.linalg import expm

# ----------------------------------------------------------------------
# analytic basis for a uniform plane-strain region, fields ~ e^{ikx}
# Y(z) = [ux, uz, sxz, szz]; columns: [dec(A=1,B=0), dec(0,1),
# grow(1,0), grow(0,1)]; exponentials are referenced to zref (the
# caller picks zref per region so magnitudes stay O(1))
def basis(z, k, lam, mu, zref_dec, zref_grow):
    g = 2.0*mu/(lam + mu)
    lp2 = lam + 2.0*mu
    ed = np.exp(-k*(z - zref_dec))       # decaying downward
    eg = np.exp(+k*(z - zref_grow))      # growing downward
    kz = k*z
    B = np.zeros((4, 4), dtype=complex)
    # decaying, (A,B)=(1,0)
    B[:, 0] = [1j*ed, -ed, -2.0*mu*k*1j*ed, 2.0*mu*k*ed]
    # decaying, (A,B)=(0,1)
    B[:, 1] = [1j*kz*ed, (-(1.0 + g) - kz)*ed,
               mu*k*1j*(-g - 2.0*kz)*ed,
               k*((lam + 2.0*mu)*g + 2.0*mu*kz)*ed]
    # growing, (A,B)=(1,0)
    B[:, 2] = [1j*eg, eg, 2.0*mu*k*1j*eg, 2.0*mu*k*eg]
    # growing, (A,B)=(0,1)
    B[:, 3] = [1j*kz*eg, (-(1.0 + g) + kz)*eg,
               mu*k*1j*(-g + 2.0*kz)*eg,
               k*(-(lam + 2.0*mu)*g + 2.0*mu*kz)*eg]
    return B

# the first-order ODE matrix dY/dz = A Y (for basis verification)
def odemat(k, lam, mu):
    lp2 = lam + 2.0*mu
    A = np.zeros((4, 4), dtype=complex)
    A[0] = [0.0, -1j*k, 1.0/mu, 0.0]
    A[1] = [-1j*k*lam/lp2, 0.0, 0.0, 1.0/lp2]
    A[2] = [k*k*4.0*mu*(lam + mu)/lp2, 0.0, 0.0, -1j*k*lam/lp2]
    A[3] = [0.0, 0.0, -1j*k, 0.0]
    return A

# ----------------------------------------------------------------------
# state-vector jumps across z = d for a point dislocation at (xs, d):
# Burgers vector b*shat on a plane with unit normal nhat (2-D: line
# dislocation along strike).  Moment density m_ij = mu(s_i n_j +
# s_j n_i) + lam (s.n) delta_ij; jumps from the delta / delta'
# structure of the equivalent body force (see README derivation):
#   [ux]  = m_xz/(2 pi mu)
#   [uz]  = m_zz/(2 pi (lam+2mu))
#   [sxz] = (ik/2pi) (m_xx - lam m_zz/(lam+2mu))
#   [szz] = 0
def jumps(k, lam, mu, b, shat, nhat, xs):
    sdotn = shat[0]*nhat[0] + shat[1]*nhat[1]
    mxx = b*(2.0*mu*shat[0]*nhat[0] + lam*sdotn)
    mzz = b*(2.0*mu*shat[1]*nhat[1] + lam*sdotn)
    mxz = b*mu*(shat[0]*nhat[1] + shat[1]*nhat[0])
    ph = np.exp(-1j*k*xs)/(2.0*np.pi)
    # overall sign: [u] = b*shat on the +nhat side (checked against
    # the independent FEM split-node reference, gate T1)
    return -np.array([mxz/mu, mzz/(lam + 2.0*mu),
                      1j*k*(mxx - lam*mzz/(lam + 2.0*mu)), 0.0],
                     dtype=complex)*ph

# ----------------------------------------------------------------------
# conditioned region basis: the pure exponential columns degenerate
# as k*width -> 0 (e^{+-kz} -> 1), and the matcher then amplifies
# rounding into a long-wavelength garbage field (caught by the T1
# gate as a smooth ~10 percent stress error growing with depth).
# For small k*width use the matrix-exponential propagator columns
# expm(A (z-za)), which are perfectly conditioned there; switch to
# the exponential basis once the scales separate.
def basis_cs(zloc, k, lam, mu):
    """conditioned analytic basis for SMALL k*width: symmetric /
    antisymmetric (cosh/sinh) combinations of the growing and
    decaying families, in LOCAL coordinates zloc = z - za (the span
    is shift-invariant).  Exact solutions; columns stay independent
    as k -> 0 (pairwise conditioning ~ 1/(k w) instead of the
    catastrophic degeneracy of the pure exponentials)."""
    g = 2.0*mu/(lam + mu)
    lp2g = (lam + 2.0*mu)*g
    kz = k*zloc
    ch, sh = np.cosh(kz), np.sinh(kz)
    B = np.empty((4, 4), dtype=complex)
    B[:, 0] = [1j*ch, sh, 2j*mu*k*sh, 2.0*mu*k*ch]
    B[:, 1] = [1j*sh, ch, 2j*mu*k*ch, 2.0*mu*k*sh]
    B[:, 2] = [1j*kz*ch, -(1.0 + g)*ch + kz*sh,
               1j*mu*k*(-g*ch + 2.0*kz*sh),
               k*(-lp2g*sh + 2.0*mu*kz*ch)]
    B[:, 3] = [1j*kz*sh, -(1.0 + g)*sh + kz*ch,
               1j*mu*k*(-g*sh + 2.0*kz*ch),
               k*(-lp2g*ch + 2.0*mu*kz*sh)]
    return B

def region_basis(k, lam, mu, za, zb):
    w = zb - za
    if k*w < 2.0:
        return lambda z: basis_cs(z - za, k, lam, mu)
    return lambda z: basis(z, k, lam, mu, za, zb)

def basis_cs_vec(zloc, k, lam, mu):
    """vectorized basis_cs: zloc array (m,) -> (m, 4, 4)"""
    zloc = np.atleast_1d(zloc)
    g = 2.0*mu/(lam + mu)
    lp2g = (lam + 2.0*mu)*g
    kz = k*zloc
    ch, sh = np.cosh(kz), np.sinh(kz)
    m = zloc.size
    B = np.empty((m, 4, 4), dtype=complex)
    B[:, 0, 0] = 1j*ch;  B[:, 1, 0] = sh
    B[:, 2, 0] = 2j*mu*k*sh; B[:, 3, 0] = 2.0*mu*k*ch
    B[:, 0, 1] = 1j*sh;  B[:, 1, 1] = ch
    B[:, 2, 1] = 2j*mu*k*ch; B[:, 3, 1] = 2.0*mu*k*sh
    B[:, 0, 2] = 1j*kz*ch
    B[:, 1, 2] = -(1.0 + g)*ch + kz*sh
    B[:, 2, 2] = 1j*mu*k*(-g*ch + 2.0*kz*sh)
    B[:, 3, 2] = k*(-lp2g*sh + 2.0*mu*kz*ch)
    B[:, 0, 3] = 1j*kz*sh
    B[:, 1, 3] = -(1.0 + g)*sh + kz*ch
    B[:, 2, 3] = 1j*mu*k*(-g*sh + 2.0*kz*ch)
    B[:, 3, 3] = k*(-lp2g*ch + 2.0*mu*kz*sh)
    return B

def basis_vec(z, k, lam, mu, zref_dec, zref_grow):
    """vectorized exponential basis: z array (m,) -> (m, 4, 4)"""
    z = np.atleast_1d(z)
    g = 2.0*mu/(lam + mu)
    ed = np.exp(-k*(z - zref_dec))
    eg = np.exp(+k*(z - zref_grow))
    kz = k*z
    m = z.size
    B = np.empty((m, 4, 4), dtype=complex)
    B[:, 0, 0] = 1j*ed; B[:, 1, 0] = -ed
    B[:, 2, 0] = -2.0*mu*k*1j*ed; B[:, 3, 0] = 2.0*mu*k*ed
    B[:, 0, 1] = 1j*kz*ed
    B[:, 1, 1] = (-(1.0 + g) - kz)*ed
    B[:, 2, 1] = mu*k*1j*(-g - 2.0*kz)*ed
    B[:, 3, 1] = k*((lam + 2.0*mu)*g + 2.0*mu*kz)*ed
    B[:, 0, 2] = 1j*eg; B[:, 1, 2] = eg
    B[:, 2, 2] = 2.0*mu*k*1j*eg; B[:, 3, 2] = 2.0*mu*k*eg
    B[:, 0, 3] = 1j*kz*eg
    B[:, 1, 3] = (-(1.0 + g) + kz)*eg
    B[:, 2, 3] = mu*k*1j*(-g + 2.0*kz)*eg
    B[:, 3, 3] = k*(-(lam + 2.0*mu)*g + 2.0*mu*kz)*eg
    return B

def eval_k_vec(zarr, k, lam1, mu1, lam2, mu2, H, d, c):
    """vectorized eval_k: state vectors (m, 4) at depths zarr"""
    zarr = np.atleast_1d(zarr)
    Y = np.empty((zarr.size, 4), dtype=complex)
    m1 = zarr <= d
    m2 = (~m1) & (zarr <= H)
    m3 = zarr > H
    def rb(msk, za, zb, lm, mm, cc):
        if not np.any(msk):
            return
        w = zb - za
        if k*w < 2.0:
            B = basis_cs_vec(zarr[msk] - za, k, lm, mm)
        else:
            B = basis_vec(zarr[msk], k, lm, mm, za, zb)
        Y[msk] = np.einsum("mij,j->mi", B, cc)
    rb(m1, 0.0, d, lam1, mu1, c[0:4])
    rb(m2, d, H, lam1, mu1, c[4:8])
    if np.any(m3):
        B = basis_vec(zarr[m3], k, lam2, mu2, H, H)[:, :, :2]
        Y[m3] = np.einsum("mij,j->mi", B, c[8:10])
    return Y

# layered solve at one wavenumber k (>0) for one point source at
# depth d in the PLATE (0 < d < H): regions R1 [0,d], R2 [d,H],
# substrate [H, inf).  lam1,mu1 plate; lam2,mu2 substrate (may be
# complex for Laplace-domain viscoelasticity).  Returns the 10
# coefficients (c1(4), c2(4), c3(2)).
def solve_k(k, lam1, mu1, lam2, mu2, H, d, J,
            rho1=0.0, rho2=0.0, g=0.0):
    """gravity (g > 0) enters through the ADVECTED BOUNDARY
    conditions (local buoyancy; no self-gravitation, no internal
    buoyancy): with z positive DOWN and uz positive down,
      free surface:  szz(0) - rho1 g uz(0) = 0
      interface:     szz_below - szz_above - (rho2 - rho1) g uz = 0
    (sign convention validated by the T4 gates: bounded, buoyancy-
    suppressed relaxed limit; g -> 0 recovers the elastic gates)."""
    M = np.zeros((10, 10), dtype=complex)
    rhs = np.zeros(10, dtype=complex)
    B1 = region_basis(k, lam1, mu1, 0.0, d)
    B2 = region_basis(k, lam1, mu1, d, H)
    # substrate: decaying only, referenced to H (well conditioned)
    def B3(z):
        return basis(z, k, lam2, mu2, H, H)[:, :2]
    # free surface: sxz = 0 and szz - rho1 g uz = 0
    M[0, 0:4] = B1(0.0)[2, :]
    M[1, 0:4] = B1(0.0)[3, :] - rho1*g*B1(0.0)[1, :]
    # source jumps at z=d: Y2(d) - Y1(d) = J
    M[2:6, 4:8] = B2(d)
    M[2:6, 0:4] = -B1(d)
    rhs[2:6] = J
    # interface at z=H: ux, uz, sxz continuous; szz jump balances
    # the advected density contrast
    M[6:10, 8:10] = B3(H)
    M[6:10, 4:8] = -B2(H)
    M[9, 8:10] -= (rho2 - rho1)*g*B3(H)[1, :]
    # column equilibration for the solve
    cs = np.max(np.abs(M), axis=0)
    cs[cs == 0.0] = 1.0
    c = np.linalg.solve(M/cs, rhs)/cs
    return c

# evaluate the k-domain state vector at depth z given the coefficients
def eval_k(z, k, lam1, mu1, lam2, mu2, H, d, c):  # (c from solve_k)
    if z <= d:
        return region_basis(k, lam1, mu1, 0.0, d)(z) @ c[0:4]
    if z <= H:
        return region_basis(k, lam1, mu1, d, H)(z) @ c[4:8]
    return basis(z, k, lam2, mu2, H, H)[:, :2] @ c[8:10]

# sxx needs the constitutive relation (not in the state vector):
# sxx = (lam+2mu) ik ux + lam duz/dz with duz/dz = (szz - ik lam ux)/(lam+2mu)
def sxx_from_Y(Y, k, lam, mu):
    lp2 = lam + 2.0*mu
    return (lp2*1j*k*Y[0] + lam*(Y[3] - 1j*k*lam*Y[0])/lp2)

# ----------------------------------------------------------------------
# real-space fields: integrate over k with the reality symmetry
# u(x,z) = 2 Re int_0^inf Yhat(k) e^{ikx} dk  (Yhat(-k) = conj(Yhat(k))
# holds for these sources; verified in the T0 gate).  Gauss-Legendre
# panels on a log grid in k.
def kgrid(kmin, kmax, npanel=60, nquad=8):
    edges = np.geomspace(kmin, kmax, npanel + 1)
    xg, wg = np.polynomial.legendre.leggauss(nquad)
    ks, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        ks.append(0.5*(b - a)*xg + 0.5*(a + b))
        ws.append(0.5*(b - a)*wg)
    return np.concatenate(ks), np.concatenate(ws)

def fields_pairs_ref(xpts, zpts, sources, lam1, mu1, lam2, mu2, H,
                 kmin=None, kmax=None, npanel=80, nquad=8,
                 rho1=0.0, rho2=0.0, g=0.0):
    """REFERENCE implementation of fields at MATCHED points
    (xpts[i], zpts[i]), i.e. the diagonal
    of fields_xz, at O(n) per wavenumber instead of O(n^2): the
    k-domain state vector is evaluated once per receiver DEPTH and
    multiplied by that receiver's own phase e^{i k x_i}.  Same
    conventions and return keys as fields_xz; this is what kernel
    generation wants, since only matched source-receiver pairs enter
    an interaction matrix."""
    xpts = np.atleast_1d(np.asarray(xpts, float))
    zpts = np.atleast_1d(np.asarray(zpts, float))
    if xpts.size != zpts.size:
        raise ValueError("fields_pairs: xpts and zpts must match")
    cplx = any(np.iscomplexobj(np.asarray(m))
               for m in (lam1, mu1, lam2, mu2))
    dmin = min(s[4] for s in sources)
    if kmin is None: kmin = 1e-3/H
    if kmax is None: kmax = 60.0/max(np.min(np.abs(zpts - dmin)), 1e-2*H)
    ks, ws = kgrid(kmin, kmax, npanel, nquad)
    dt_ = complex if cplx else float
    out = {f: np.zeros(xpts.size, dtype=dt_) for f in
           ("ux", "uz", "sxx", "sxz", "szz")}
    Sp = np.array([1.0, -1.0, 1.0, -1.0])
    lamv = np.where(zpts <= H, lam1, lam2)
    muv = np.where(zpts <= H, mu1, mu2)
    lp2v = lamv + 2.0*muv
    for k, w in zip(ks, ws):
        Yk = np.zeros((zpts.size, 4), dtype=complex)
        sxxk = np.zeros(zpts.size, dtype=complex)
        Ym = np.zeros((zpts.size, 4), dtype=complex)
        sxxm = np.zeros(zpts.size, dtype=complex)
        for (b, sh, nh, x0, d) in sources:
            J = jumps(k, lam1, mu1, b, sh, nh, x0)
            c = solve_k(k, lam1, mu1, lam2, mu2, H, d, J, rho1, rho2, g)
            Yv = eval_k_vec(zpts, k, lam1, mu1, lam2, mu2, H, d, c)
            Yk += Yv
            sxxk += (lp2v*1j*k*Yv[:, 0] +
                     lamv*(Yv[:, 3] - 1j*k*lamv*Yv[:, 0])/lp2v)
            if cplx:
                Jm = jumps(-k, lam1, mu1, b, sh, nh, x0)
                cm = solve_k(k, lam1, mu1, lam2, mu2, H, d, Sp*Jm,
                             rho1, rho2, g)
                Y2 = eval_k_vec(zpts, k, lam1, mu1, lam2, mu2, H, d, cm)
                Ym += Sp[None, :]*Y2
                sxxm += -(lp2v*1j*k*Y2[:, 0] +
                          lamv*(Y2[:, 3] - 1j*k*lamv*Y2[:, 0])/lp2v)
        ph = np.exp(1j*k*xpts)
        if cplx:
            for i, f in enumerate(("ux", "uz", "sxz", "szz")):
                out[f] += w*(Yk[:, i]*ph + Ym[:, i]/ph)
            out["sxx"] += w*(sxxk*ph + sxxm/ph)
        else:
            for i, f in enumerate(("ux", "uz", "sxz", "szz")):
                out[f] += 2.0*w*np.real(Yk[:, i]*ph)
            out["sxx"] += 2.0*w*np.real(sxxk*ph)
    return out

def ve_fields_pairs(xpts, zpts, sources, lam1, mu1, lam2, mu2, H, tau2,
                    t, tau1=None, M=16, rho1=0.0, rho2=0.0, g=0.0,
                    lam_mode="bulk", ctx=None, **kw):
    """viscoelastic fields at matched points (Talbot inversion of
    fields_pairs); same signature idea as ve_fields_xz.

    ctx: a pair_ctx() for this geometry and k-grid.  The contour nodes
    differ only in the moduli, so the geometric part is built once
    here (or handed in, to share it across sample times as well)."""
    if ctx is None:
        ctx = pair_ctx(xpts, zpts, sources, H, **kw)
    elif kw:
        raise ValueError("ve_fields_pairs: pass either ctx or the "
                         "k-grid arguments, not both")
    def F(s):
        l2, m2 = maxwell_moduli(s, mu2, lam2, tau2, lam_mode)
        if tau1 is not None:
            l1, m1 = maxwell_moduli(s, mu1, lam1, tau1, lam_mode)
        else:
            l1, m1 = lam1, mu1
        out = fields_pairs_ctx(ctx, l1, m1, l2, m2,
                               rho1=rho1, rho2=rho2, g=g)
        return np.stack([out[c] for c in
                         ("ux", "uz", "sxx", "sxz", "szz")])/s
    R = talbot(F, t, M)
    return {c: R[i] for i, c in
            enumerate(("ux", "uz", "sxx", "sxz", "szz"))}

def fields_xz(xobs, zobs, sources, lam1, mu1, lam2, mu2, H,
              kmin=None, kmax=None, npanel=80, nquad=8,
              rho1=0.0, rho2=0.0, g=0.0):
    """sources: list of (b, shat, nhat, xs, d).  Returns dict with
    ux, uz, sxx, sxz, szz arrays over (len(zobs), len(xobs)).

    Real moduli: uses the conjugate symmetry Y(-k) = conj(Y(k)),
    real output.  COMPLEX moduli (Laplace-domain viscoelasticity):
    the symmetry fails; the negative-k half line is computed via the
    parity relation A(-k) = S A(k) S, S = diag(1,-1,1,-1), i.e. by a
    second solve with source S J(-k); complex output."""
    zs = np.atleast_1d(zobs); xs_ = np.atleast_1d(xobs)
    cplx = any(np.iscomplexobj(np.asarray(m))
               for m in (lam1, mu1, lam2, mu2))
    dmin = min(s[4] for s in sources)
    if kmin is None: kmin = 1e-3/H
    if kmax is None: kmax = 60.0/max(min(abs(z - dmin) for z in zs), 1e-2*H)
    ks, ws = kgrid(kmin, kmax, npanel, nquad)
    dt_ = complex if cplx else float
    out = {f: np.zeros((zs.size, xs_.size), dtype=dt_) for f in
           ("ux", "uz", "sxx", "sxz", "szz")}
    Sp = np.array([1.0, -1.0, 1.0, -1.0])
    for k, w in zip(ks, ws):
        Yk = np.zeros((4, zs.size), dtype=complex)
        sxxk = np.zeros(zs.size, dtype=complex)
        Ym = np.zeros((4, zs.size), dtype=complex)
        sxxm = np.zeros(zs.size, dtype=complex)
        for (b, sh, nh, x0, d) in sources:
            J = jumps(k, lam1, mu1, b, sh, nh, x0)
            c = solve_k(k, lam1, mu1, lam2, mu2, H, d, J,
                        rho1, rho2, g)
            if cplx:
                Jm = jumps(-k, lam1, mu1, b, sh, nh, x0)
                cm = solve_k(k, lam1, mu1, lam2, mu2, H, d, Sp*Jm,
                             rho1, rho2, g)
            Yv = eval_k_vec(zs, k, lam1, mu1, lam2, mu2, H, d, c)
            Yk += Yv.T
            lamv = np.where(zs <= H, lam1, lam2)
            muv = np.where(zs <= H, mu1, mu2)
            lp2v = lamv + 2.0*muv
            sxxk += (lp2v*1j*k*Yv[:, 0] +
                     lamv*(Yv[:, 3] - 1j*k*lamv*Yv[:, 0])/lp2v)
            if cplx:
                Y2 = eval_k_vec(zs, k, lam1, mu1, lam2, mu2, H, d, cm)
                Ym += (Sp[None, :]*Y2).T
                sxxm += -(lp2v*1j*k*Y2[:, 0] +
                          lamv*(Y2[:, 3] - 1j*k*lamv*Y2[:, 0])/lp2v)
        ph = np.exp(1j*np.outer(np.ones(zs.size), xs_)*k)
        if cplx:
            for i, f in enumerate(("ux", "uz", "sxz", "szz")):
                out[f] += w*(Yk[i][:, None]*ph + Ym[i][:, None]/ph)
            out["sxx"] += w*(sxxk[:, None]*ph + sxxm[:, None]/ph)
        else:
            for i, f in enumerate(("ux", "uz", "sxz", "szz")):
                out[f] += 2.0*w*np.real(Yk[i][:, None]*ph)
            out["sxx"] += 2.0*w*np.real(sxxk[:, None]*ph)
    return out

# ----------------------------------------------------------------------
# a uniform-slip fault segment = superposition of point dislocations
# along the fault trace (line integral discretized by ns points);
# fault from (x1,z1) to (x2,z2), slip b along the fault (shat), normal
# nhat perpendicular.  This IS the dislocation-density representation:
# uniform slip on a segment equals point dislocations at the two ENDS
# only in the antiplane case; in plane strain we keep the line
# distribution of the moment density (exact as ns -> inf).
def segment_sources(x1, z1, x2, z2, b, ns=40):
    L = np.hypot(x2 - x1, z2 - z1)
    shat = ((x2 - x1)/L, (z2 - z1)/L)
    nhat = (-shat[1], shat[0])
    t = (np.arange(ns) + 0.5)/ns
    return [(b*L/ns, shat, nhat, x1 + t_*(x2 - x1), z1 + t_*(z2 - z1))
            for t_ in t]

# ----------------------------------------------------------------------
# viscoelasticity: Maxwell in shear, ELASTIC IN BULK (constant bulk
# modulus K); plane strain uses the 3-D moduli, so
#   mu(s)  = mu s/(s + 1/tau),   tau = eta/mu
#   lam(s) = K - (2/3) mu(s)
def maxwell_moduli(s, mu, lam, tau, lam_mode="bulk"):
    """lam_mode 'bulk': constant bulk modulus (default);
    'lambda': constant first Lame modulus (Rundle 1982 figures)"""
    mus = mu*s/(s + 1.0/tau)
    if lam_mode == "lambda":
        return lam, mus
    K = lam + 2.0*mu/3.0
    return K - 2.0*mus/3.0, mus

# fixed-Talbot inverse Laplace transform (Abate & Valko), M nodes;
# F must accept a complex s and return a (possibly array) value.
def talbot(F, t, M=24):
    r = 2.0*M/(5.0*t)
    S = 0.5*np.exp(r*t)*F(r)
    for j in range(1, M):
        th = j*np.pi/M
        ct = 1.0/np.tan(th)
        s = r*th*(ct + 1j)
        sig = th + (th*ct - 1.0)*ct
        S = S + np.real(np.exp(t*s)*F(s)*(1.0 + 1j*sig))
    return S*r/M

# step-response stress fields for the layered problem with a MAXWELL
# SUBSTRATE (plate elastic), slip applied at t=0 and held: evaluates
# fields_xz with s-dependent substrate moduli on the Talbot contour.
# For a uniform Maxwell medium pass the same (lam, mu, tau) for the
# plate (tau1 not None).
def ve_fields_xz(xobs, zobs, sources, lam1, mu1, lam2, mu2, H, tau2,
                 t, tau1=None, M=16, rho1=0.0, rho2=0.0, g=0.0,
                 lam_mode="bulk", **kw):
    def F(s):
        l2, m2 = maxwell_moduli(s, mu2, lam2, tau2, lam_mode)
        if tau1 is not None:
            l1, m1 = maxwell_moduli(s, mu1, lam1, tau1, lam_mode)
        else:
            l1, m1 = lam1, mu1
        out = fields_xz(xobs, zobs, sources, l1, m1, l2, m2, H,
                        rho1=rho1, rho2=rho2, g=g, **kw)
        return np.stack([out[c] for c in
                         ("ux", "uz", "sxx", "sxz", "szz")])/s
    R = talbot(F, t, M)
    return {c: R[i] for i, c in
            enumerate(("ux", "uz", "sxx", "sxz", "szz"))}

# ----------------------------------------------------------------------
# BATCHED EVALUATION PATH
#
# The per-wavenumber loop above is correct but spends most of its time,
# at kernel-generation sizes, building a full 4x4 basis matrix per
# receiver and contracting it with the coefficient vector.  Both region
# bases are however SEPARABLE: with the four scalar functions of depth
# collected into f_m(z), each basis is
#
#     B(z) = sum_m f_m(z) C_m ,     C_m constant in z
#
#   cosh/sinh form:  f = (ch, sh, kz ch, kz sh),  kz = k (z - za)
#   exponential:     f = (ed, eg, kz ed, kz eg),  kz = k z,
#                    ed = e^{-k(z-zref_dec)},  eg = e^{+k(z-zref_grow)}
#
# so the state vector at a receiver is Y(z) = sum_m f_m(z) (C_m c) with
# four small matrix-vector products per (wavenumber, source) and no
# per-receiver 4x4 at all.  Summing the k-quadrature then turns into a
# matrix product: with
#
#     g_m[k, i] = w_k e^{i k x_i} f_m[k, z_i]
#
# the quadrature is  out[i, :] = sum_m (g_m^T V_m)[i, :],  V_m[k, :] =
# C_m(k) c(k).  g_m carries no material parameters, so it is built once
# per source geometry and reused for every Laplace-contour node and
# every sample time; the negative-k parity contribution reuses it as
# conj(g_m) because w and f are real.
#
# The C_m matrices below are read off the basis definitions; the
# t7_batched_fields.py gate checks them against basis()/basis_cs() and
# the batched fields against fields_pairs() to roundoff.

def _split_cs(k, lam, mu):
    """C_m for basis_cs: B(zloc) = ch C0 + sh C1 + (kz ch) C2 + (kz sh) C3.
    k may be an array (n,); returns four (n, 4, 4) arrays."""
    k = np.atleast_1d(np.asarray(k))
    z = np.zeros_like(k, dtype=complex)
    gm = 2.0*mu/(lam + mu)
    lp2g = (lam + 2.0*mu)*gm
    u = np.stack([1j + z, z, z, 2.0*mu*k + z])           # (4, n)
    v = np.stack([z, 1.0 + z, 2j*mu*k + z, z])
    p = np.stack([z, -(1.0 + gm) + z, -1j*mu*k*gm + z, z])
    q = np.stack([z, z, z, -k*lp2g + z])
    o = np.zeros_like(u)
    def M(c0, c1, c2, c3):                               # columns -> (n,4,4)
        return np.stack([c0, c1, c2, c3], axis=-1).transpose(1, 0, 2)
    return (M(u, v, p, q), M(v, u, q, p),
            M(o, o, u, v), M(o, o, v, u))

def _split_exp(k, lam, mu):
    """C_m for basis: B(z) = ed D0 + eg D1 + (kz ed) D2 + (kz eg) D3"""
    k = np.atleast_1d(np.asarray(k))
    z = np.zeros_like(k, dtype=complex)
    gm = 2.0*mu/(lam + mu)
    lp2g = (lam + 2.0*mu)*gm
    a = np.stack([1j + z, -1.0 + z, -2j*mu*k + z, 2.0*mu*k + z])
    b = np.stack([1j + z, 1.0 + z, 2j*mu*k + z, 2.0*mu*k + z])
    pd = np.stack([z, -(1.0 + gm) + z, -1j*mu*k*gm + z, k*lp2g + z])
    pg = np.stack([z, -(1.0 + gm) + z, -1j*mu*k*gm + z, -k*lp2g + z])
    o = np.zeros_like(a)
    def M(c0, c1, c2, c3):
        return np.stack([c0, c1, c2, c3], axis=-1).transpose(1, 0, 2)
    return (M(a, pd, o, o), M(o, o, b, pg),
            M(o, a, o, o), M(o, o, o, b))

def _basis_b(kb, z, za, zb, lam, mu, cs):
    """(nb, 4, 4) region basis at depth z for the region [za, zb], with
    the same cosh/sinh vs exponential choice as region_basis makes,
    per batch member (cs is the boolean mask k*(zb-za) < 2)"""
    kb = np.asarray(kb)
    nb = kb.size
    zv = np.broadcast_to(np.asarray(z, float), (nb,))
    B = np.empty((nb, 4, 4), dtype=complex)
    if cs.any():
        kk = kb[cs]
        C0, C1, C2, C3 = _split_cs(kk, lam, mu)
        kz = kk*(zv[cs] - za)
        ch, sh = np.cosh(kz), np.sinh(kz)
        f = (ch, sh, kz*ch, kz*sh)
        B[cs] = (f[0][:, None, None]*C0 + f[1][:, None, None]*C1 +
                 f[2][:, None, None]*C2 + f[3][:, None, None]*C3)
    if (~cs).any():
        kk = kb[~cs]
        D0, D1, D2, D3 = _split_exp(kk, lam, mu)
        zz = zv[~cs]
        ed = np.exp(-kk*(zz - za))
        eg = np.exp(kk*(zz - zb))
        kz = kk*zz
        B[~cs] = (ed[:, None, None]*D0 + eg[:, None, None]*D1 +
                  (kz*ed)[:, None, None]*D2 + (kz*eg)[:, None, None]*D3)
    return B

def jumps_b(kb, lam, mu, b, shat, nhat, xs):
    """vectorized jumps(): (nb, 4) for an array of wavenumbers"""
    kb = np.asarray(kb)
    sdotn = shat[0]*nhat[0] + shat[1]*nhat[1]
    mxx = b*(2.0*mu*shat[0]*nhat[0] + lam*sdotn)
    mzz = b*(2.0*mu*shat[1]*nhat[1] + lam*sdotn)
    mxz = b*mu*(shat[0]*nhat[1] + shat[1]*nhat[0])
    lp2 = lam + 2.0*mu
    J = np.zeros((kb.size, 4), dtype=complex)
    J[:, 0] = mxz/mu
    J[:, 1] = mzz/lp2
    J[:, 2] = 1j*kb*(mxx - lam*mzz/lp2)
    return -J*(np.exp(-1j*kb*xs)/(2.0*np.pi))[:, None]

def solve_k_b(kb, lam1, mu1, lam2, mu2, H, d, Jb, rho1=0.0, rho2=0.0,
              g=0.0):
    """batched solve_k over wavenumbers for ONE source depth d:
    kb (nb,), Jb (nb, 4) -> coefficients (nb, 10).  Same equations,
    same column equilibration, as solve_k."""
    kb = np.asarray(kb)
    nb = kb.size
    cs1 = kb*d < 2.0                     # region [0, d]
    cs2 = kb*(H - d) < 2.0               # region [d, H]
    B1_0 = _basis_b(kb, 0.0, 0.0, d, lam1, mu1, cs1)
    B1_d = _basis_b(kb, d, 0.0, d, lam1, mu1, cs1)
    B2_d = _basis_b(kb, d, d, H, lam1, mu1, cs2)
    B2_H = _basis_b(kb, H, d, H, lam1, mu1, cs2)
    # substrate: decaying pair only, referenced to H, so ed = eg = 1
    D0, _, D2, _ = _split_exp(kb, lam2, mu2)
    B3_H = (D0[:, :, :2] + (kb*H)[:, None, None]*D2[:, :, :2])
    M = np.zeros((nb, 10, 10), dtype=complex)
    rhs = np.zeros((nb, 10), dtype=complex)
    M[:, 0, 0:4] = B1_0[:, 2, :]
    M[:, 1, 0:4] = B1_0[:, 3, :] - rho1*g*B1_0[:, 1, :]
    M[:, 2:6, 4:8] = B2_d
    M[:, 2:6, 0:4] = -B1_d
    rhs[:, 2:6] = Jb
    M[:, 6:10, 8:10] = B3_H
    M[:, 6:10, 4:8] = -B2_H
    M[:, 9, 8:10] -= (rho2 - rho1)*g*B3_H[:, 1, :]
    csc = np.max(np.abs(M), axis=1)                  # per column
    csc[csc == 0.0] = 1.0
    c = np.linalg.solve(M/csc[:, None, :], rhs[..., None])[..., 0]/csc
    return c

class _PairCtx:
    """geometric part of a matched-pair evaluation: wavenumber grid,
    quadrature weights times receiver phases, and the depth functions
    f_m of the region bases.  None of it depends on the moduli, so one
    context serves every Laplace-contour node and every sample time of
    a given source and receiver set."""
    __slots__ = ("ks", "ws", "nz", "H", "zpts", "srcs")

def pair_ctx(xpts, zpts, sources, H, kmin=None, kmax=None, npanel=80,
             nquad=8):
    xpts = np.atleast_1d(np.asarray(xpts, float))
    zpts = np.atleast_1d(np.asarray(zpts, float))
    if xpts.size != zpts.size:
        raise ValueError("pair_ctx: xpts and zpts must match")
    dmin = min(s[4] for s in sources)
    if kmin is None:
        kmin = 1e-3/H
    if kmax is None:
        kmax = 60.0/max(np.min(np.abs(zpts - dmin)), 1e-2*H)
    ks, ws = kgrid(kmin, kmax, npanel, nquad)
    ctx = _PairCtx()
    ctx.ks, ctx.ws, ctx.nz, ctx.H, ctx.zpts = ks, ws, xpts.size, H, zpts
    wph = np.exp(1j*np.outer(xpts, ks))*ws[None, :]      # (nz, nk)
    ctx.srcs = []
    for (b, sh, nh, x0, d) in sources:
        # receivers grouped by region exactly as eval_k_vec groups them
        m1 = np.nonzero(zpts <= d)[0]
        m2 = np.nonzero((zpts > d) & (zpts <= H))[0]
        m3 = np.nonzero(zpts > H)[0]
        cs1 = ks*d < 2.0
        cs2 = ks*(H - d) < 2.0
        regions = []
        for rows, sl, za, zb, csm, sub in ((m1, (0, 4), 0.0, d, cs1, False),
                                           (m2, (4, 8), d, H, cs2, False),
                                           (m3, (8, 10), H, H, None, True)):
            if rows.size == 0:
                continue
            z = zpts[rows]
            groups = []
            if sub:
                # substrate: the decaying pair only, both references at
                # H, so the growing slots are absent: f = (ed, kz ed)
                idx = np.arange(ks.size)
                ed = np.exp(-np.outer(z - H, ks))
                f = (ed, np.outer(z, ks)*ed)
                groups.append(dict(kidx=idx, branch="sub",
                                   G=np.concatenate(
                                       [wph[np.ix_(rows, idx)]*fm
                                        for fm in f], axis=1)))
            else:
                for branch, sel in (("cs", csm), ("ex", ~csm)):
                    idx = np.nonzero(sel)[0]
                    if idx.size == 0:
                        continue
                    kk = ks[idx]
                    if branch == "cs":
                        kz = np.outer(z - za, kk)
                        ch, shh = np.cosh(kz), np.sinh(kz)
                        f = (ch, shh, kz*ch, kz*shh)
                    else:
                        ed = np.exp(-np.outer(z - za, kk))
                        eg = np.exp(np.outer(z - zb, kk))
                        kz = np.outer(z, kk)
                        f = (ed, eg, kz*ed, kz*eg)
                    groups.append(dict(kidx=idx, branch=branch,
                                       G=np.concatenate(
                                           [wph[np.ix_(rows, idx)]*fm
                                            for fm in f], axis=1)))
            regions.append(dict(rows=rows, sl=sl, sub=sub, groups=groups))
        ctx.srcs.append(dict(b=b, sh=sh, nh=nh, x0=x0, d=d,
                             regions=regions))
    return ctx

_SP = np.array([1.0, -1.0, 1.0, -1.0])

def fields_pairs_ctx(ctx, lam1, mu1, lam2, mu2, rho1=0.0, rho2=0.0,
                     g=0.0):
    """matched-pair fields for one set of moduli, given a pair_ctx.
    Same conventions and return keys as fields_pairs_ref."""
    ks, H = ctx.ks, ctx.H
    cplx = any(np.iscomplexobj(np.asarray(m))
               for m in (lam1, mu1, lam2, mu2))
    acc = np.zeros((ctx.nz, 5), dtype=complex)
    for ent in ctx.srcs:
        d = ent["d"]
        Jb = jumps_b(ks, lam1, mu1, ent["b"], ent["sh"], ent["nh"],
                     ent["x0"])
        runs = [(solve_k_b(ks, lam1, mu1, lam2, mu2, H, d, Jb,
                           rho1, rho2, g), 1.0)]
        if cplx:
            # negative-k half line by the parity relation, as in
            # fields_pairs_ref: a second solve with source S J(-k)
            Jm = jumps_b(-ks, lam1, mu1, ent["b"], ent["sh"], ent["nh"],
                         ent["x0"])
            runs.append((solve_k_b(ks, lam1, mu1, lam2, mu2, H, d,
                                   _SP[None, :]*Jm, rho1, rho2, g), -1.0))
        for cc, par in runs:
            for reg in ent["regions"]:
                a, b_ = reg["sl"]
                lm, mm = (lam2, mu2) if reg["sub"] else (lam1, mu1)
                for grp in reg["groups"]:
                    idx = grp["kidx"]
                    kk = ks[idx]
                    if grp["branch"] == "cs":
                        Cm = _split_cs(kk, lm, mm)
                    elif grp["branch"] == "ex":
                        Cm = _split_exp(kk, lm, mm)
                    else:
                        Ce = _split_exp(kk, lm, mm)
                        Cm = (Ce[0][:, :, :2], Ce[2][:, :, :2])
                    cs_ = cc[idx, a:b_]
                    nm, nki = len(Cm), idx.size
                    V = np.empty((nm*nki, 5), dtype=complex)
                    for m_ in range(nm):
                        Y = np.einsum("bij,bj->bi", Cm[m_], cs_)
                        r = slice(m_*nki, (m_ + 1)*nki)
                        V[r, 0] = Y[:, 0]
                        V[r, 1] = par*Y[:, 1]
                        V[r, 2] = Y[:, 2]
                        V[r, 3] = par*Y[:, 3]
                        V[r, 4] = par*1j*kk*Y[:, 0]
                    G = grp["G"] if par > 0.0 else np.conj(grp["G"])
                    acc[reg["rows"]] += G @ V
    lamv = np.where(ctx.zpts <= H, lam1, lam2)
    muv = np.where(ctx.zpts <= H, mu1, mu2)
    lp2v = lamv + 2.0*muv
    if cplx:
        col = [acc[:, i] for i in range(5)]
    else:
        col = [2.0*acc[:, i].real for i in range(5)]
    # sxx from the constitutive relation; the 1j k ux part is carried in
    # the fifth accumulator so that the k factor stays on the
    # wavenumber side of the quadrature
    sxx = (lp2v - lamv*lamv/lp2v)*col[4] + (lamv/lp2v)*col[3]
    return dict(ux=col[0], uz=col[1], sxx=sxx, sxz=col[2], szz=col[3])

def fields_pairs(xpts, zpts, sources, lam1, mu1, lam2, mu2, H,
                 kmin=None, kmax=None, npanel=80, nquad=8,
                 rho1=0.0, rho2=0.0, g=0.0):
    """matched-pair fields (batched path; see _PairCtx).  Identical
    interface and conventions to fields_pairs_ref, which it replaces in
    production use and is checked against by t7_batched_fields.py."""
    ctx = pair_ctx(xpts, zpts, sources, H, kmin, kmax, npanel, nquad)
    return fields_pairs_ctx(ctx, lam1, mu1, lam2, mu2, rho1, rho2, g)
