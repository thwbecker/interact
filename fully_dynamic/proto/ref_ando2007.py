"""
Reference kernels transcribed from Ando, Kame and Yamashita (2007, EPS 59,
363-373; CC BY), equations (3)-(8) (in-plane, after Tada and Madariaga 2001)
and (30)-(31) (antiplane, after Cochard and Madariaga 1994).

Their convention (their eq. 1):
    sigma_pq = -(mu/(2 beta)) * D  -(mu/(2 beta)) * sum K_pq D
with the primitive I_pq(x1, x2, t) giving, through
    -(mu/(2 beta)) * [I(x1 + ds/2, x2, t) - I(x1 - ds/2, x2, t)],
the stress at (x1, x2) and time t due to unit slip velocity on the element
|x1| <= ds/2 of the x1 axis acting from time 0 to t (ramp response).
Source element on the x1 axis, normal x2 to the left of increasing x1.

This module is used only as an independent check.

STATUS
- Antiplane (eqs 30-31): transcription verified (static limit reproduces the
  screw-dislocation field to 1e-7; agrees with antiplane.py to machine
  precision away from the fronts). Used by the tests.
- In-plane (Tada and Madariaga 2001 eqs A4-A6, identical to Ando et al.
  2007 eqs 4-6): transcribed from the typeset PDF pages (all factors of
  p = beta/alpha in the numerator). Verified against the static edge
  dislocation limit; used by the tests.
"""
import numpy as np


def _I_inplane(x1, x2, t, alpha, beta):
    """Tada and Madariaga (2001) eqs (A4)-(A6), read from the typeset PDF
    (identical to Ando et al. 2007 eqs 4-6). All factors of p = beta/alpha
    are in the numerator."""
    x1 = np.asarray(x1, float); x2 = np.asarray(x2, float); t = np.asarray(t, float)
    x1, x2, t = np.broadcast_arrays(x1, x2, t)
    r = np.hypot(x1, x2)
    p = beta / alpha
    HL = (t > r / alpha).astype(float)
    HT = (t > r / beta).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        gL = np.sqrt(np.maximum((alpha * t / r) ** 2 - 1.0, 0.0))
        gT = np.sqrt(np.maximum((beta * t / r) ** 2 - 1.0, 0.0))
        gL3 = gL**3; gT3 = gT**3
        q = 2 * (3 * x1**2 - x2**2) / (3 * r**2)
        q2 = 2 * (3 * x2**2 - x1**2) / (3 * r**2)
        a2 = x2**2 / r**2; a1 = x1**2 / r**2
        # (A4)
        I11 = (-(1 / np.pi) * HL * (2 * x2 / r) * p * (q * p**2 * gL3 + (1 - 2 * a2 * p**2) * gL)
               + (1 / np.pi) * HT * (2 * x2 / r) * (q * gT3 + (1 - 2 * a2) * gT))
        # (A5)
        I22 = ((1 / np.pi) * HL * (2 * x2 / r) * p * (q * p**2 * gL3 + (2 * a1 * p**2 - 1) * gL)
               - (1 / np.pi) * HT * (2 * x2 / r) * (q * gT3 + (2 * a1 - 1) * gT))
        # (A6)
        Hpl = (x1 > 0).astype(float) * (t > np.abs(x2) / beta).astype(float)
        arg = np.abs(x1) / np.sqrt(np.maximum((beta * t) ** 2 - x2**2, 1e-300))
        arc = np.arccos(np.clip(arg, -1.0, 1.0))
        I12 = (Hpl
               - (1 / np.pi) * np.sign(x1) * HL * (2 * np.abs(x1) / r) * p * (q2 * p**2 * gL3 + 2 * a2 * p**2 * gL)
               + (1 / np.pi) * np.sign(x1) * HT * ((2 * np.abs(x1) / r) * (q2 * gT3 + 2 * a2 * gT) - arc))
    return I11, I12, I22


def segment_ramp_stress_ref(xr, yr, t, a, alpha, beta, rho):
    """Physical stress (sxx, sxy, syy) at (xr, yr), time t, unit slip rate on
    |x1| <= a since t = 0, from the Ando et al. (2007) primitives."""
    mu = rho * beta**2
    Ip = _I_inplane(xr + a, yr, t, alpha, beta)
    Im = _I_inplane(xr - a, yr, t, alpha, beta)
    return [-(mu / (2 * beta)) * (u - v) for u, v in zip(Ip, Im)]


def _I_antiplane(x1, x2, t, beta):
    x1 = np.asarray(x1, float); x2 = np.asarray(x2, float); t = np.asarray(t, float)
    x1, x2, t = np.broadcast_arrays(x1, x2, t)
    r = np.hypot(x1, x2)
    Hs = (t > r / beta).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        gb = np.sqrt(np.maximum((beta * t / r) ** 2 - 1.0, 0.0))
        I31 = -(x2 / (np.pi * r)) * Hs * gb                                    # eq. (30)
        Hpl = (x1 > 0).astype(float) * (t > np.abs(x2) / beta).astype(float)
        arg = np.abs(x1) / np.sqrt(np.maximum((beta * t) ** 2 - x2**2, 1e-300))
        arc = np.arccos(np.clip(arg, -1.0, 1.0))
        I32 = Hpl + (1 / np.pi) * np.sign(x1) * Hs * ((np.abs(x1) / r) * gb - arc)   # eq. (31)
    return I31, I32


def antiplane_segment_ramp_stress_ref(xr, yr, t, a, beta, rho):
    """Physical (s31, s32) at (xr, yr), time t, unit antiplane slip rate on
    |x1| <= a since t = 0."""
    mu = rho * beta**2
    Ip = _I_antiplane(xr + a, yr, t, beta)
    Im = _I_antiplane(xr - a, yr, t, beta)
    return [-(mu / (2 * beta)) * (u - v) for u, v in zip(Ip, Im)]
