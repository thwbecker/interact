"""
Antiplane (mode III) kernels built with the same pipeline as the in-plane
ones: 2-D SH Green's function, double time integral, equivalent double
couple, stresses, then segment quadrature from inplane.segment_ramp_stress.

Point field: slip patch of unit length and unit slip rate created at t = 0
on the plane y = 0 at the origin, antiplane slip u3. Components returned:
(s31, s32). As for the in-plane field, the sign is opposite to the
physical one (the segment routine applies the minus sign), and the
whole-cone integral equals the plane SH wave, (0, +mu/(2 beta)).

Derivation (sympy, run once at import; expressions are small):
    G33(r, t)  = H(t - r/beta) / (2 pi mu sqrt(t^2 - r^2/beta^2))
    Ghat       = double time integral = (u acosh(u/r) - sqrt(u^2 - r^2)) / (2 pi mu beta), u = beta t
    u3         = mu d/dy Ghat   (double couple m_32 = m_23 = mu * slip)
    s31 = mu d/dx u3,  s32 = mu d/dy u3
"""
import numpy as np
import sympy as sp

_x, _y, _t = sp.symbols("x y t", real=True)
_be, _rho = sp.symbols("beta rho", positive=True)
_mu = _rho * _be**2
_R = sp.sqrt(_x**2 + _y**2)
_u = _be * _t
_Z = sp.sqrt(_u**2 - _R**2)
_Ghat = (_u * sp.acosh(_u / _R) - _Z) / (2 * sp.pi * _mu * _be)
_u3 = _mu * sp.diff(_Ghat, _y)
_s31 = _mu * sp.diff(_u3, _x)
_s32 = _mu * sp.diff(_u3, _y)
_f = sp.lambdify((_x, _y, _t, _be, _rho), (_s31, _s32), modules="numpy")


def sigma_point_ramp_ap(x, y, t, alpha, beta, rho):
    """(s31, s32) at (x, y), time t; alpha is accepted for interface
    compatibility and unused (no P wave in antiplane)."""
    x = np.asarray(x, float); y = np.asarray(y, float); t = np.asarray(t, float)
    x, y, t = np.broadcast_arrays(x, y, t)
    R = np.hypot(x, y)
    m = beta * t > R
    s31 = np.zeros(x.shape); s32 = np.zeros(x.shape)
    if m.any():
        with np.errstate(all="ignore"):
            a, b = _f(x[m], y[m], t[m], beta, rho)
        s31[m] = a; s32[m] = b
    return s31, s32


def segment_ramp_stress_ap(xr, yr, t, a, alpha, beta, rho, method="closed", ngauss=12):
    """(s31, s32). method="closed": Ando et al. (2007) eqs (30)-(31)
    (production). method="quad": quadrature of the point field derived here
    (cross-check)."""
    if method == "closed":
        from ref_ando2007 import antiplane_segment_ramp_stress_ref
        return antiplane_segment_ramp_stress_ref(xr, yr, t, a, beta, rho)
    from inplane import segment_ramp_stress_quad
    return segment_ramp_stress_quad(xr, yr, t, a, alpha, beta, rho, ngauss=ngauss,
                                    field=sigma_point_ramp_ap, ncomp=2, pw_index=1)


def static_segment_stress_ap(xr, yr, a, mu):
    """Static (s31, s32) per unit slip from uniform antiplane slip on |x'|<=a:
    screw dislocation difference P(x-a) - P(x+a), P = (-mu y/(2 pi r^2), mu x/(2 pi r^2))."""
    def P(x, y):
        r2 = x * x + y * y
        return (-mu * y / (2 * np.pi * r2), mu * x / (2 * np.pi * r2))
    p1 = P(xr - a, yr); p2 = P(xr + a, yr)
    return [u - v for u, v in zip(p1, p2)]


def W_table_ap(Xs, ss, a, Xr, sr, nlag, dt, alpha, beta, rho, et=1.0):
    """Antiplane analogue of inplane.W_table: shear traction W[l] (M, nlag)
    and static row; the normal entries are zero (returned for interface
    compatibility)."""
    import inplane as ip
    from ref_ando2007 import antiplane_segment_ramp_stress_ref
    mu = rho * beta**2
    Xr = np.atleast_2d(np.asarray(Xr, float)); sr = np.asarray(sr, float)
    xr, yr = ip.receiver_in_source_frame(Xs, ss, Xr)
    nr = ip.rot90(sr)                               # receiver normal (global)
    # components of the receiver normal in the source frame
    ns_src = np.array([np.dot(nr, ss), np.dot(nr, ip.rot90(ss))])
    tl = np.concatenate([[0.0], (et + np.arange(0, nlag)) * dt])
    s31, s32 = antiplane_segment_ramp_stress_ref(xr[:, None], yr[:, None], tl[None, :], a, beta, rho)
    K = np.diff(ns_src[0] * s31 + ns_src[1] * s32, axis=1)
    p1, p2 = static_segment_stress_ap(xr, yr, a, mu)
    Ks = ns_src[0] * p1 + ns_src[1] * p2
    Ws = K - Ks[:, None] * dt
    z = np.zeros_like(Ws)
    return Ws, z, Ks, np.zeros_like(Ks)
