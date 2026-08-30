#!/usr/bin/env python3
"""T3: elastic plate over Maxwell substrate (the thrust-cycle target,
no gravity).  Gates:
 (a) t -> 0+ reproduces the layered ELASTIC solution (unrelaxed);
 (b) t -> infty approaches the elastic solution with the substrate
     shear modulus relaxed away (mu2 -> 0, bulk kept);
 (c) interface tractions (sxz, szz) stay continuous across z = H at
     all sampled times, while ux, sxx may jump;
 (d) the plate-point stress evolves monotonically between the two
     limits (single held dislocation, no reloading)."""
import numpy as np
from inplane2d import fields_xz, ve_fields_xz, segment_sources

yr = 3.15576e7
mu1 = lam1 = 30e9
mu2, lam2 = 30e9, 30e9
K2 = lam2 + 2*mu2/3
tau = 10.0*yr
d1, d2, H = 4e3, 12e3, 20e3
src = segment_sources(0.0, d1, 0.0, d2, 1.0, ns=12)
kw = dict(kmin=2e-7, kmax=4e-3, npanel=32, nquad=5)
xs = np.array([6e3, 14e3])
zp, zsb = 15e3, 30e3          # plate / substrate points

el0 = fields_xz(xs, np.array([zp, zsb]), src, lam1, mu1, lam2, mu2, H, **kw)
mur = 1e-8*mu2
elinf = fields_xz(xs, np.array([zp, zsb]), src, lam1, mu1,
                  K2 - 2*mur/3, mur, H, **kw)
ok = True
# (a) unrelaxed limit
ve = ve_fields_xz(xs, np.array([zp, zsb]), src, lam1, mu1, lam2, mu2,
                  H, tau, 1e-3*tau, M=12, **kw)
for c in ("sxx", "sxz", "szz"):
    e = np.max(np.abs(np.real(ve[c]) - el0[c]))/np.max(np.abs(el0[c]))
    print(f"T3a t=0+   {c}: rel dev from unrelaxed elastic {e:8.2e}")
    ok = ok and e < 5e-3
# (b) relaxed limit (slow small-k modes: generous time and tolerance)
ve9 = ve_fields_xz(xs, np.array([zp, zsb]), src, lam1, mu1, lam2, mu2,
                   H, tau, 300.0*tau, M=14, **kw)
for c in ("sxx", "sxz", "szz"):
    sc = np.max(np.abs(elinf[c])) + 0.02*np.max(np.abs(el0[c]))
    e = np.max(np.abs(np.real(ve9[c]) - elinf[c]))/sc
    print(f"T3b t=300tau {c}: rel dev from relaxed elastic  {e:8.2e}")
    ok = ok and e < 0.05
# (c) interface traction continuity through time
# exactly at the interface from both sides (the evaluator
# branches at z <= H); sxz, szz, ux, uz must be continuous,
# sxx legitimately jumps once the substrate has relaxed
zint = np.array([H, H*(1.0 + 1e-12)])
for tfac in (0.3, 1.0, 10.0):
    vi = ve_fields_xz(xs, zint, src, lam1, mu1, lam2, mu2, H, tau,
                      tfac*tau, M=12, **kw)
    for c in ("sxz", "szz"):
        num = np.max(np.abs(vi[c][0] - vi[c][1]))
        den = np.max(np.abs(vi[c])) + 1e-30
        print(f"T3c t={tfac:4.1f}tau {c}: interface jump/scale "
              f"{num/den:8.2e}")
        ok = ok and num/den < 1e-6
# (d) monotonic evolution at the plate point
ts = np.array([0.1, 1.0, 10.0, 30.0])
v = []
for tf in ts:
    vv = ve_fields_xz(xs[:1], np.array([zp]), src, lam1, mu1, lam2,
                      mu2, H, tau, tf*tau, M=12, **kw)
    v.append(float(np.real(vv["sxz"][0, 0])))
v = np.array(v)
mono = np.all(np.diff(v) > 0) or np.all(np.diff(v) < 0)
print("T3d plate sxz(t):", " ".join(f"{q:9.1f}" for q in v),
      "monotonic" if mono else "NOT monotonic")
ok = ok and mono
print("T3 " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
