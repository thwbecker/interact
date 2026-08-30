#!/usr/bin/env python3
"""T4: gravity (advected-boundary buoyancy) gates.

 (a) elastic fields with g on change only at the rho g/(mu k) level
     (sub-percent at seismogenic wavelengths), so the Okada anchor
     still holds;
 (b) plate over Maxwell substrate WITH gravity: the t -> infty limit
     matches the direct elastic solve with the substrate shear
     modulus relaxed away AND the same gravity rows (the buoyant
     flexure state);
 (c) the through-plate degeneracy is CURED: for a fault cutting the
     whole plate, the relaxed surface uplift is kmin-sensitive
     without gravity (neutral mode leaking through the long-
     wavelength cutoff) and kmin-INsensitive with gravity."""
import numpy as np
from inplane2d import fields_xz, ve_fields_xz, segment_sources

yr = 3.15576e7
mu1 = lam1 = mu2 = lam2 = 30e9
K2 = lam2 + 2*mu2/3
rho1, rho2, g = 2800.0, 3300.0, 9.81
tau = 10.0*yr
H = 20e3
ok = True

# (a) contained fault, elastic, g on vs off
d1, d2 = 4e3, 12e3
src = segment_sources(0.0, d1, 0.0, d2, 1.0, ns=12)
kw = dict(kmin=2e-7, kmax=4e-3, npanel=32, nquad=5)
xs = np.array([6e3, 14e3]); zz = np.array([15e3, 30e3])
e0 = fields_xz(xs, zz, src, lam1, mu1, lam2, mu2, H, **kw)
eg = fields_xz(xs, zz, src, lam1, mu1, lam2, mu2, H,
               rho1=rho1, rho2=rho2, g=g, **kw)
for c in ("sxx", "sxz", "szz"):
    dev = np.max(np.abs(eg[c] - e0[c]))/np.max(np.abs(e0[c]))
    print(f"T4a elastic {c}: |g-on - g-off|/scale = {dev:8.2e}")
    ok = ok and dev < 0.02
# (b) relaxed limit with gravity
mur = 1e-8*mu2
elinf = fields_xz(xs, zz, src, lam1, mu1, K2 - 2*mur/3, mur, H,
                  rho1=rho1, rho2=rho2, g=g, **kw)
ve9 = ve_fields_xz(xs, zz, src, lam1, mu1, lam2, mu2, H, tau,
                   300.0*tau, M=14, rho1=rho1, rho2=rho2, g=g, **kw)
for c in ("sxx", "sxz", "szz"):
    sc = np.max(np.abs(elinf[c])) + 0.02*np.max(np.abs(e0[c]))
    e = np.max(np.abs(np.real(ve9[c]) - elinf[c]))/sc
    print(f"T4b relaxed-with-g {c}: dev from direct solve {e:8.2e}")
    ok = ok and e < 0.05
# (c) through-plate degeneracy cured: k-domain RELAXED surface
# response of a fault cutting the whole plate.  Without gravity the
# long-wavelength response is limited only by the (vanishing)
# substrate rigidity and blows up as k -> 0; buoyancy caps it.
from inplane2d import jumps, solve_k, eval_k
mur = 1e-10*mu2
lr = K2 - 2*mur/3
dsrc = 10e3                      # a mid-plate point source
J1 = jumps(1e-6, lam1, mu1, 1.0, (0, 1), (1, 0), 0.0)
rat = {}
for kk in (1e-6, 1e-5):
    J = jumps(kk, lam1, mu1, 1.0, (0, 1), (1, 0), 0.0)
    u = {}
    for gg in (0.0, g):
        c = solve_k(kk, lam1, mu1, lr, mur, H, dsrc, J,
                    rho1, rho2, gg)
        u[gg] = abs(eval_k(0.0, kk, lam1, mu1, lr, mur, H, dsrc, c)[1])
    rat[kk] = u[0.0]/u[g]
    print(f"T4c k={kk:.0e}: relaxed |uz(k)| surface, g=0 / g-on = "
          f"{rat[kk]:10.1f}")
ok = ok and rat[1e-6] > 10.0 and rat[1e-6] > rat[1e-5]
print("T4 " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
