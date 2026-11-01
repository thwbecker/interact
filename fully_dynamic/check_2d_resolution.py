#!/usr/bin/env python3
"""Resolution of the rate-and-state length scales on a 2-D interact section (in-plane, plane
strain), for deciding on an element size before QD and FD comparisons (fd_plan.md M9).

Per element (geometry geom_*.in, rsf_ab.in with a b, rsf_dc.in with Dc, rsf_sigma.in with the
effective normal stress, all in the section directory):
  h       element length (2 x half-length)
  Lb      mu' Dc / (b sigma), mu' = mu / (1 - nu) for plane strain
  Lambda0 (9 pi / 32) Lb, the quasi-static estimate of the process zone at rest (rupture shrinks
          it further, so the dynamic zone is smaller)
  h*      (2 / pi) mu' Dc / ((b - a) sigma) for velocity-weakening elements (a < b)
Prints per group the smallest Lambda0 / h and h* / h, with the counts of elements below the
thresholds (defaults 3 and 20, of the order recommended in the SEAS benchmark descriptions;
cautious lower bounds rather than convergence guarantees), and the largest element size that
would meet them everywhere.
Usage: check_2d_resolution.py section_dir [mu] [nu] [min Lambda0/h] [min h*/h]"""
import sys
import os
import glob
import numpy as np

d = sys.argv[1]
mu = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0e10
nu = float(sys.argv[3]) if len(sys.argv) > 3 else 0.25
lam_min = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
hst_min = float(sys.argv[5]) if len(sys.argv) > 5 else 20.0
g = np.loadtxt(sorted(glob.glob(os.path.join(d, "geom_*.in")))[0], ndmin=2)
ab = np.loadtxt(os.path.join(d, "rsf_ab.in"), ndmin=2)
dc = np.loadtxt(os.path.join(d, "rsf_dc.in"), ndmin=1).ravel()
sig = np.abs(np.loadtxt(os.path.join(d, "rsf_sigma.in"), ndmin=1).ravel())
n = len(g)
if not (len(ab) == len(dc) == len(sig) == n):
    sys.exit(f"check_2d_resolution: {n} elements but {len(ab)}, {len(dc)}, {len(sig)} lines in rsf_ab, rsf_dc, rsf_sigma")
h = 2.0*g[:, 5]
grp = g[:, 7].astype(int)
a = ab[:, 0]
b = ab[:, 1]
mup = mu/(1.0 - nu)
Lb = mup*dc/(b*sig)
lam = 9.0*np.pi/32.0*Lb
vw = a < b
hst = np.full(n, np.inf)
hst[vw] = 2.0/np.pi*mup*dc[vw]/((b[vw] - a[vw])*sig[vw])
print(f"check_2d_resolution: {d}: {n} elements of {h.min():.1f} to {h.max():.1f} m, mu' {mup:.3e} Pa, "
      f"{vw.sum()} velocity weakening; thresholds Lambda0/h >= {lam_min:g}, h*/h >= {hst_min:g}")
worst = 0.0
for k in np.unique(grp):
    m = grp == k
    mv = m & vw
    r_lam = lam[m]/h[m]
    line = (f"  group {k}: {m.sum()} elements, h {h[m].min():.1f} to {h[m].max():.1f} m, Lambda0 {lam[m].min():.1f} to "
            f"{lam[m].max():.1f} m, min Lambda0/h {r_lam.min():.2f} ({(r_lam < lam_min).sum()} below)")
    if mv.any():
        r_h = hst[mv]/h[mv]
        line += f", VW {mv.sum()}: h* {hst[mv].min():.0f} to {hst[mv].max():.0f} m, min h*/h {r_h.min():.1f} ({(r_h < hst_min).sum()} below)"
    print(line)
hmax_lam = np.min(lam/lam_min)
hmax_hst = np.min(hst[vw]/hst_min) if vw.any() else np.inf
print(f"  element size meeting both thresholds everywhere: {min(hmax_lam, hmax_hst):.1f} m "
      f"(Lambda0: {hmax_lam:.1f} m, h*: {hmax_hst:.1f} m)")
i = int(np.argmin(lam/h))
print(f"  least resolved element: {i} (group {grp[i]}, x {g[i, 0]/1e3:.2f} km, depth {-g[i, 1]/1e3:.2f} km): "
      f"h {h[i]:.1f} m, Dc {dc[i]:.3g} m, sigma {sig[i]/1e6:.1f} MPa, a {a[i]:.4f}, b {b[i]:.4f}, Lambda0 {lam[i]:.1f} m")
