"""
2-D plane-strain (in-plane, mode II) kernels for straight constant-slip-rate
segments in a homogeneous full space.

Conventions (documented, to be matched to interact with a global sign flip)
- A segment has centre X, unit tangent s, unit normal n = rot90(s).
- Slip D = u(+n side) - u(-n side) along s. Positive slip on a segment
  lowers the shear traction on that segment (stress drop) in the static limit.
- Stresses are "stress from slip" (physical). interact's Is/In are in the
  stress-drop convention, i.e. opposite sign.
- sigma_point_ramp (generated, inplane_kernel.py) is the field of a point
  slip patch of unit length and unit slip rate created at t = 0, and carries
  the OPPOSITE sign to the physical field (see test_plane_wave_identity).
  The segment field is therefore  S_seg = - integral_{-a}^{a} S(x - x') dx'.
- Production kernel engine: closed-form primitives of Tada and Madariaga
  (2001) / Ando et al. (2007) (segment_ramp_stress_closed). The quadrature of
  the independently derived point field (segment_ramp_stress_quad) is kept
  as the cross-check; the two agree to about 1e-13 relative away from the
  fronts (test_inplane.py).
- Lag kernel K^(l), l >= 1: stress at time l*dt due to unit slip rate during
  [0, dt). K^(1) on the source element itself is approximately -mu/(2 beta)
  (exactly so when alpha*dt < a). History convolution:
      tau(t_{n+1}) = sum_{m <= n} K^(n+1-m) V^m.
"""
import numpy as np
from inplane_kernel import sigma_point_ramp

_GL_CACHE = {}


def _gauss(n):
    if n not in _GL_CACHE:
        v, w = np.polynomial.legendre.leggauss(n)
        _GL_CACHE[n] = (0.5 * (v + 1.0), 0.5 * w)  # on [0,1]
    return _GL_CACHE[n]


def _phi(v):
    """Cubic map [0,1] -> [0,1] with zero slope at both ends; removes
    inverse-square-root endpoint singularities of the integrand."""
    return 3 * v**2 - 2 * v**3, 6 * v * (1 - v)


def _integrate_S(x_lo, x_hi, xr, yr, t, alpha, beta, rho, ngauss, extra_bp=(), field=None):
    """integral_{x_lo}^{x_hi} S(xr - x', yr, t) dx' for arrays of equal shape.
    Breakpoints at the P and S front crossings and at any extra_bp arrays.
    Returns (sxx, sxy, syy) arrays. Intervals with x_hi <= x_lo give 0."""
    x_lo = np.asarray(x_lo, float); x_hi = np.asarray(x_hi, float)
    xr = np.asarray(xr, float); yr = np.asarray(yr, float); t = np.asarray(t, float)
    shp = np.broadcast(x_lo, x_hi, xr, yr, t).shape
    x_lo, x_hi, xr, yr, t = [np.broadcast_to(z, shp).ravel() for z in (x_lo, x_hi, xr, yr, t)]
    x_lo = np.minimum(x_lo, x_hi)                    # empty interval -> zero length
    bps = [x_lo, x_hi]
    for c in (alpha, beta):
        with np.errstate(invalid="ignore"):
            h = np.sqrt(np.maximum((c * t) ** 2 - yr**2, 0.0))
        ok = c * t > np.abs(yr)
        for sgn in (-1.0, 1.0):
            b = np.where(ok, xr + sgn * h, x_lo)
            bps.append(np.clip(b, x_lo, x_hi))
    for b in extra_bp:
        bps.append(np.clip(np.broadcast_to(np.asarray(b, float), shp).ravel(), x_lo, x_hi))
    B = np.sort(np.stack(bps, axis=1), axis=1)       # (M, nb)
    v, w = _gauss(ngauss)
    ph, dph = _phi(v)
    lo = B[:, :-1]; hi = B[:, 1:]                    # (M, nb-1)
    L = hi - lo
    xp = lo[:, :, None] + L[:, :, None] * ph[None, None, :]        # (M, ni, ng)
    wt = L[:, :, None] * (w * dph)[None, None, :]
    xi = xr[:, None, None] - xp
    if field is None:
        field = sigma_point_ramp
    comps = field(xi, yr[:, None, None], t[:, None, None], alpha, beta, rho)
    out = [np.sum(s * wt, axis=(1, 2)).reshape(shp) for s in comps]
    return out


def segment_ramp_stress_quad(xr, yr, t, a, alpha, beta, rho, ngauss=12, near=0.5, field=None, ncomp=3, pw_index=1):
    """Physical stress tensor (sxx, sxy, syy) at receiver point (xr, yr) in the
    source frame, time t, due to unit slip rate on the segment |x'| <= a, y' = 0
    since t = 0. Receivers near the segment line inside the slipping interval
    (|yr| < near*a and |xr| < a) are evaluated by the complement identity
        int_{-a}^{a} S = int_all S - int_{x' < -a} S - int_{x' > a} S,
    with int_all S = (0, mu/(2 beta), 0) (plane S wave; verified in tests),
    which avoids the hypersingular integrand at x' = xr.
    """
    xr = np.asarray(xr, float); yr = np.asarray(yr, float); t = np.asarray(t, float)
    shp = np.broadcast(xr, yr, t).shape
    xr, yr, t = [np.broadcast_to(z, shp).ravel().copy() for z in (xr, yr, t)]
    mu = rho * beta**2
    res = [np.zeros(xr.shape) for _ in range(ncomp)]
    # receivers close to the source line inside the slipping interval: use the
    # complement identity (the direct integrand has a near-hypersingular peak
    # of width |yr| at x' = xr that Gauss quadrature cannot resolve)
    onseg = (np.abs(yr) < near * a) & (np.abs(xr) < a)
    off = ~onseg
    if off.any():
        xo, yo, to = xr[off], yr[off], t[off]
        nearline = np.abs(yo) < 2 * a
        for sel, refine in ((~nearline, False), (nearline, True)):
            if not sel.any():
                continue
            xq, yq, tq = xo[sel], yo[sel], to[sel]
            bp = ()
            if refine:
                # geometric refinement around the nearest point: the integrand
                # has a peak of width |yr| at x' = xr
                d = np.maximum(np.abs(yq), 1e-3 * a)
                xc = np.clip(xq, -a, a)
                bp = [xc] + [xc + sgn * d * 2.0**k for k in range(0, 8) for sgn in (-1.0, 1.0)]
            o = _integrate_S(-a, a, xq, yq, tq, alpha, beta, rho, ngauss, extra_bp=bp, field=field)
            idx = np.flatnonzero(off)[sel]
            for k in range(ncomp):
                res[k][idx] = -o[k]
    if onseg.any():
        xs, ys, ts = xr[onseg], yr[onseg], t[onseg]
        Lc = alpha * ts * 1.000001
        # outer parts: x' < -a  (xi = xs - x' >= xs + a)  and  x' > a;
        # geometric breakpoints from the segment ends outward (integrand ~ 1/xi^2)
        dl = xs + a; dr = a - xs
        bpl = [-a - dl * 2.0**k for k in range(0, 24)]
        bpr = [a + dr * 2.0**k for k in range(0, 24)]
        left = _integrate_S(xs - Lc, np.full_like(xs, -a), xs, ys, ts, alpha, beta, rho, ngauss, extra_bp=bpl, field=field)
        right = _integrate_S(np.full_like(xs, a), xs + Lc, xs, ys, ts, alpha, beta, rho, ngauss, extra_bp=bpr, field=field)
        z = np.zeros_like(xs)
        # whole-cone integral of the (sign-flipped) point field = plane S wave,
        # +mu/(2 beta) in the fault-parallel shear component, zero elsewhere
        tot = [z] * ncomp
        tot[pw_index] = np.where(beta * ts > np.abs(ys), mu / (2 * beta), 0.0)
        for k in range(ncomp):
            res[k][onseg] = -(tot[k] - left[k] - right[k])
    return [r.reshape(shp) for r in res]


def segment_ramp_stress_closed(xr, yr, t, a, alpha, beta, rho):
    """Closed form: Tada and Madariaga (2001) eqs (A4)-(A6), which equal the
    Ando et al. (2007) eqs (4)-(6). Physical stress (sxx, sxy, syy) at (xr, yr)
    in the source frame, time t, unit slip rate on |x'| <= a since t = 0:
        -(mu/(2 beta)) [I(xr + a, yr, t) - I(xr - a, yr, t)].
    Valid for any receiver, including the source element itself (the
    endpoints are at finite distance from an interior collocation point).
    Loses relative accuracy for c t / r >> 1 through cancellation of the t^3
    terms between the P and S parts (see Tada and Madariaga, section 3.5);
    within the dynamic window used here this is not an issue."""
    from ref_ando2007 import _I_inplane
    mu = rho * beta**2
    xr = np.asarray(xr, float); yr = np.asarray(yr, float); t = np.asarray(t, float)
    Ip = _I_inplane(xr + a, yr, t, alpha, beta)
    Im = _I_inplane(xr - a, yr, t, alpha, beta)
    return [-(mu / (2 * beta)) * (u - v) for u, v in zip(Ip, Im)]


def segment_ramp_stress(xr, yr, t, a, alpha, beta, rho, method="closed", **kw):
    """Dispatcher. method="closed" (default, production): Tada-Madariaga
    primitives. method="quad": quadrature of the independently derived point
    field (inplane_kernel.py); used as the cross-check."""
    if method == "closed":
        return segment_ramp_stress_closed(xr, yr, t, a, alpha, beta, rho)
    return segment_ramp_stress_quad(xr, yr, t, a, alpha, beta, rho, **kw)


# ---------------- static references ----------------
def static_point_edge(x, y, mu, nu):
    """Hirth and Lothe edge dislocation, unit Burgers vector along x, cut on
    x < 0 of the plane y = 0, slip jump u(0+) - u(0-) = +1 on the cut."""
    D = mu / (2 * np.pi * (1 - nu))
    r2 = x * x + y * y
    return (-D * y * (3 * x * x + y * y) / r2**2,
            D * x * (x * x - y * y) / r2**2,
            D * y * (x * x - y * y) / r2**2)


def static_segment_stress(xr, yr, a, mu, nu):
    """Static stress (per unit slip) from uniform slip on |x'| <= a, y' = 0:
    P(x - a) - P(x + a)."""
    p1 = static_point_edge(xr - a, yr, mu, nu)
    p2 = static_point_edge(xr + a, yr, mu, nu)
    return [u - v for u, v in zip(p1, p2)]


# ---------------- geometry, projection, lag kernels ----------------
def rot90(s):
    return np.array([-s[1], s[0]])


def receiver_in_source_frame(Xs, ss, Xr):
    d = np.asarray(Xr, float) - np.asarray(Xs, float)
    return d @ ss, d @ rot90(ss)


def project(sxx, sxy, syy, ss, sr):
    """Rotate a source-frame stress tensor to global and return the shear and
    normal tractions on a receiver with tangent sr, normal rot90(sr)."""
    ns = rot90(ss); nr = rot90(sr)
    Q = np.stack([ss, ns], axis=1)        # columns: source frame axes in global
    S = np.empty(np.shape(sxx) + (2, 2))
    S[..., 0, 0] = sxx; S[..., 0, 1] = sxy; S[..., 1, 0] = sxy; S[..., 1, 1] = syy
    Sg = Q @ S @ Q.T
    tr = Sg @ nr
    return tr @ sr, tr @ nr


def lag_kernel(Xs, ss, a, Xr, sr, lags, dt, alpha, beta, rho, method="closed", ngauss=12):
    """K^(l) for the lags in `lags` (integers >= 1): shear and normal traction
    on the receiver at time l*dt due to unit slip rate on the source during
    [0, dt). Returns two arrays of shape lags.shape."""
    lags = np.asarray(lags)
    xr, yr = receiver_in_source_frame(Xs, ss, Xr)
    t1 = lags * dt
    t0 = (lags - 1) * dt
    s1 = segment_ramp_stress(xr, yr, t1, a, alpha, beta, rho, method=method, ngauss=ngauss)
    s0 = segment_ramp_stress(xr, yr, t0, a, alpha, beta, rho, method=method, ngauss=ngauss)
    d = [u - v for u, v in zip(s1, s0)]
    return project(d[0], d[1], d[2], ss, sr)


def static_traction(Xs, ss, a, Xr, sr, mu, nu):
    xr, yr = receiver_in_source_frame(Xs, ss, Xr)
    s = static_segment_stress(xr, yr, a, mu, nu)
    return project(s[0], s[1], s[2], ss, sr)


def lag_kernel_table(Xs, ss, a, Xr, sr, nlag, dt, alpha, beta, rho, method="closed", ngauss=12, chunk=200000):
    """K^(l) for l = 1..nlag for many receivers at once. Xr: (M,2) centres,
    sr: (M,2) tangents (or (2,) shared). Returns (ks, kn), each (M, nlag).
    Evaluates the ramp stress once per time level and differences consecutive
    levels, so the cost is one segment quadrature per (receiver, lag)."""
    Xr = np.atleast_2d(np.asarray(Xr, float))
    sr = np.asarray(sr, float)
    M = Xr.shape[0]
    xr, yr = receiver_in_source_frame(Xs, ss, Xr)       # (M,)
    tl = np.arange(0, nlag + 1) * dt                      # time levels
    ks = np.empty((M, nlag)); kn = np.empty((M, nlag))
    # process receivers in chunks to bound memory: each entry costs ~22*ngauss points
    per = max(1, chunk // ((nlag + 1) * (22 * ngauss if method == "quad" else 4)))
    for i0 in range(0, M, per):
        i1 = min(M, i0 + per)
        s = segment_ramp_stress(xr[i0:i1, None], yr[i0:i1, None], tl[None, :], a, alpha, beta, rho, method=method, ngauss=ngauss)
        d = [np.diff(c, axis=1) for c in s]              # (m, nlag)
        if sr.ndim == 1:
            t_s, t_n = project(d[0], d[1], d[2], ss, sr)
        else:
            t_s = np.empty_like(d[0]); t_n = np.empty_like(d[0])
            for j in range(i1 - i0):
                t_s[j], t_n[j] = project(d[0][j], d[1][j], d[2][j], ss, sr[i0 + j])
        ks[i0:i1] = t_s; kn[i0:i1] = t_n
    return ks, kn


def W_table(Xs, ss, a, Xr, sr, nlag, dt, alpha, beta, rho, gamma_switch=10.0, et=1.0):
    """Dynamic-correction lag table W[l] = K[l] - Ks dt, l = 1..nlag, for many
    receivers, shear and normal tractions, each (M, nlag). Column 0 is the
    current-step kernel S(et dt) - Ks dt (collocation at et dt after the start
    of the step; et = 1 end of step, et > 1 the Peirce-Siebrits epsilon
    scheme), column l >= 1 is S((et+l) dt) - S((et+l-1) dt) - Ks dt. Closed form while
    beta (l-1) dt < gamma_switch * r_far (r_far: farther segment endpoint),
    large-time series beyond (inplane_series.W_series), which removes the
    t^3 cancellation of the closed form. Also returns the static rows."""
    from inplane_series import W_series
    mu = rho * beta**2
    nu = (alpha**2 - 2 * beta**2) / (2 * (alpha**2 - beta**2))
    Xr = np.atleast_2d(np.asarray(Xr, float)); sr = np.asarray(sr, float)
    M = Xr.shape[0]
    xr, yr = receiver_in_source_frame(Xs, ss, Xr)
    rfar = np.maximum(np.hypot(xr + a, yr), np.hypot(xr - a, yr))
    lsw = np.minimum(np.ceil(gamma_switch * rfar / beta / dt - et).astype(int) + 2, nlag)   # first column using the series
    Ks_s, Ks_n = static_traction(Xs, ss, a, Xr, sr, mu, nu)
    Ws = np.empty((M, nlag)); Wn = np.empty((M, nlag))
    lags = np.arange(1, nlag + 1)
    # closed-form part (vectorised over all receivers for lags up to max(lsw))
    lc = int(lsw.max())
    tl = np.concatenate([[0.0], (et + np.arange(0, lc)) * dt])
    s = segment_ramp_stress_closed(xr[:, None], yr[:, None], tl[None, :], a, alpha, beta, rho)
    d = [np.diff(c, axis=1) for c in s]
    ts, tn = project(d[0], d[1], d[2], ss, sr) if sr.ndim == 1 else _project_rows(d, ss, sr)
    Ws[:, :lc] = ts - Ks_s[:, None] * dt
    Wn[:, :lc] = tn - Ks_n[:, None] * dt
    # series part, per receiver from its own switch lag
    if lc < nlag or (lsw < lc).any():
        for i in range(M):
            l0 = int(lsw[i])
            if l0 >= nlag:
                continue
            tt = (lags[l0 - 1:] - 2 + et) * dt          # start time of the interval for column l-1: (et + l - 2) dt
            c = W_series(xr[i], yr[i], tt, dt, a, alpha, beta, rho)
            cs, cn = project(c[0], c[1], c[2], ss, sr if sr.ndim == 1 else sr[i])
            Ws[i, l0 - 1:] = cs; Wn[i, l0 - 1:] = cn
    return Ws, Wn, Ks_s, Ks_n


def _project_rows(d, ss, sr):
    t_s = np.empty_like(d[0]); t_n = np.empty_like(d[0])
    for j in range(d[0].shape[0]):
        t_s[j], t_n[j] = project(d[0][j], d[1][j], d[2][j], ss, sr[j])
    return t_s, t_n
