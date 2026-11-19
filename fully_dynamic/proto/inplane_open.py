"""2-D plane-strain kernels for an OPENING (tensile) constant-rate segment,
companion of inplane.py: physical stress (sxx, sxy, syy) in the source frame
at (xr, yr), time t, from unit opening rate on |x'| <= a, y' = 0 since t = 0
(opening = u_y(+) - u_y(-)). Quadrature of the generated point field
(inplane_kernel_open.py, same sign convention as the shear one), with the
complement identity on the segment: the whole-line integral of the physical
field is the plane P wave (sxx, sxy, syy) = (-lambda/(2 alpha), 0, -rho alpha/2)."""
import numpy as np
import inplane as ip
from inplane_kernel_open import sigma_point_ramp_open

def segment_ramp_stress_open(xr, yr, t, a, alpha, beta, rho, ngauss=12, near=0.5, tshift=0.0):
    """tshift > 0: return the average of the values at t - tshift and t + tshift (tables use
    a small shift, 1e-4 dt: when a wave front crosses a segment end exactly at a sample time,
    two quadrature breakpoints coincide and the point field cancels catastrophically; the
    stress itself is continuous in time, so the symmetric average changes it by O(tshift^2)
    away from fronts and O(tshift) at a kink)"""
    if tshift > 0.0:
        lo = segment_ramp_stress_open(xr, yr, np.asarray(t, float) - tshift, a, alpha, beta, rho, ngauss, near)
        hi = segment_ramp_stress_open(xr, yr, np.asarray(t, float) + tshift, a, alpha, beta, rho, ngauss, near)
        return [0.5*(u + v) for u, v in zip(lo, hi)]
    xr = np.asarray(xr, float); yr = np.asarray(yr, float); t = np.asarray(t, float)
    shp = np.broadcast(xr, yr, t).shape
    xr, yr, t = [np.broadcast_to(z, shp).ravel().copy() for z in (xr, yr, t)]
    lam = rho*(alpha**2 - 2*beta**2)
    res = [np.zeros(xr.shape) for _ in range(3)]
    onseg = (np.abs(yr) < near*a) & (np.abs(xr) < a)
    off = ~onseg
    if off.any():
        xo, yo, to = xr[off], yr[off], t[off]
        nearline = np.abs(yo) < 2*a
        for sel, refine in ((~nearline, False), (nearline, True)):
            if not sel.any(): continue
            xq, yq, tq = xo[sel], yo[sel], to[sel]
            bp = ()
            if refine:
                d = np.maximum(np.abs(yq), 1e-3*a); xc = np.clip(xq, -a, a)
                bp = [xc] + [xc + sgn*d*2.0**k for k in range(0, 8) for sgn in (-1.0, 1.0)]
            o = ip._integrate_S(-a, a, xq, yq, tq, alpha, beta, rho, ngauss, extra_bp=bp, field=sigma_point_ramp_open)
            idx = np.flatnonzero(off)[sel]
            for k in range(3): res[k][idx] = -o[k]
    if onseg.any():
        xs, ys, ts = xr[onseg], yr[onseg], t[onseg]
        Lc = alpha*ts*1.000001
        dl = xs + a; dr = a - xs
        bpl = [-a - dl*2.0**k for k in range(0, 24)]; bpr = [a + dr*2.0**k for k in range(0, 24)]
        left = ip._integrate_S(xs - Lc, np.full_like(xs, -a), xs, ys, ts, alpha, beta, rho, ngauss, extra_bp=bpl, field=sigma_point_ramp_open)
        right = ip._integrate_S(np.full_like(xs, a), xs + Lc, xs, ys, ts, alpha, beta, rho, ngauss, extra_bp=bpr, field=sigma_point_ramp_open)
        pw = alpha*ts > np.abs(ys)
        tot = [np.where(pw, lam/(2*alpha), 0.0), np.zeros_like(xs), np.where(pw, rho*alpha/2, 0.0)]   # point-field sign
        for k in range(3): res[k][onseg] = -(tot[k] - left[k] - right[k])
    return [r.reshape(shp) for r in res]
