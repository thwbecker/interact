#!/usr/bin/env python3
"""Coarsen a 2-D interact geometry with its rate-and-state input files by joining m consecutive
elements of each group into one (for fully dynamic tests at a coarser element size, fd_plan.md
M7): new element = chord from the start of the first to the end of the last element, strike from
the chord, half-length half the chord; per-element parameters averaged (rsf_ab.in, rsf_dc.in,
rsf_sigma.in, rsf_vpl.in, rsf_ic.in: shear stress and slip rate averaged).
Usage: coarsen_2d.py input_dir output_dir m [geometry file name]"""
import sys
import os
import glob
import numpy as np

src, dst, m = sys.argv[1], sys.argv[2], int(sys.argv[3])
gname = sys.argv[4] if len(sys.argv) > 4 else os.path.basename(glob.glob(os.path.join(src, "geom_*.in"))[0])
os.makedirs(dst, exist_ok=True)
g = np.loadtxt(os.path.join(src, gname))
grp = g[:, 7].astype(int)
strike = np.radians(g[:, 3])
t = np.stack([np.sin(strike), np.cos(strike)], axis=1)       # unit tangent in (x, y)
hl = g[:, 5]
blocks = []
for k in np.unique(grp):
    idx = np.where(grp == k)[0]
    nb = max(1, int(round(len(idx)/m)))         # nearly equal blocks (element sizes within a few percent)
    for bl in np.array_split(idx, nb):
        blocks.append(bl)
rows = []
for bl in blocks:
    i0, i1 = bl[0], bl[-1]
    p0 = g[i0, :2] - hl[i0]*t[i0]
    p1 = g[i1, :2] + hl[i1]*t[i1]
    c = 0.5*(p0 + p1)
    d = p1 - p0
    L = np.hypot(d[0], d[1])
    st = np.degrees(np.arctan2(d[0], d[1])) % 360.0
    rows.append(f"{c[0]:.10e} {c[1]:.10e} 0.0 {st:.6f} 90.0 {0.5*L:.10e} 0.0 {grp[i0]:d}")
open(os.path.join(dst, gname), "w").write("\n".join(rows) + "\n")
for fn in ("rsf_ab.in", "rsf_dc.in", "rsf_sigma.in", "rsf_vpl.in", "rsf_ic.in"):
    p = os.path.join(src, fn)
    if not os.path.exists(p):
        continue
    a = np.loadtxt(p, ndmin=2)
    out = np.array([a[bl].mean(axis=0) for bl in blocks])
    np.savetxt(os.path.join(dst, fn), out, fmt="%.6e")
print(f"coarsen_2d: {len(g)} elements in {len(np.unique(grp))} groups -> {len(blocks)} elements (m = {m})")
