#!/usr/bin/env python3
"""Resample a 2-D interact geometry with its rate-and-state input files to coarser, nearly equal
elements (for fully dynamic tests at a coarser element size, fd_plan.md M7). Each group is taken
as a polyline through its elements in file order; its arc length S is cut into nb = round(S/h)
equal pieces, h the target element length; a new element is the chord between the ends of a
piece. Per-element inputs (rsf_ab.in, rsf_dc.in, rsf_sigma.in, rsf_vpl.in, rsf_ic.in) are
averaged over the old elements, weighted by their overlap with the piece along the arc.
Usage: coarsen_2d.py input_dir output_dir h [geometry file name]
  h: target element length in m, or, if below 50, a factor times the median element length"""
import sys
import os
import glob
import numpy as np

src, dst, h = sys.argv[1], sys.argv[2], float(sys.argv[3])
gname = sys.argv[4] if len(sys.argv) > 4 else os.path.basename(sorted(glob.glob(os.path.join(src, "geom_*.in")))[0])
os.makedirs(dst, exist_ok=True)
g = np.loadtxt(os.path.join(src, gname))
grp = g[:, 7].astype(int)
strike = np.radians(g[:, 3])
t = np.stack([np.sin(strike), np.cos(strike)], axis=1)       # unit tangent in (x, y)
hl = g[:, 5]
if h < 50.0:
    h = h*np.median(2.0*hl)
files = {}
for fn in ("rsf_ab.in", "rsf_dc.in", "rsf_sigma.in", "rsf_vpl.in", "rsf_ic.in"):
    p = os.path.join(src, fn)
    if os.path.exists(p):
        files[fn] = np.loadtxt(p, ndmin=2)
rows = []
out = {fn: [] for fn in files}
for k in np.unique(grp):
    idx = np.where(grp == k)[0]
    p0 = g[idx, :2] - hl[idx, None]*t[idx]
    p1 = g[idx, :2] + hl[idx, None]*t[idx]
    vert = np.vstack([p0[:1], p1])                           # polyline through the element ends
    seg = 2.0*hl[idx]
    s = np.concatenate([[0.0], np.cumsum(seg)])
    S = s[-1]
    nb = max(1, int(round(S/h)))
    edges = np.linspace(0.0, S, nb + 1)

    def at(x):
        i = min(np.searchsorted(s, x, side="right") - 1, len(seg) - 1)
        f = (x - s[i])/seg[i]
        return vert[i] + f*(vert[i + 1] - vert[i])
    for b in range(nb):
        a, c = at(edges[b]), at(edges[b + 1])
        mid = 0.5*(a + c)
        d = c - a
        L = np.hypot(d[0], d[1])
        st = np.degrees(np.arctan2(d[0], d[1])) % 360.0
        rows.append(f"{mid[0]:.10e} {mid[1]:.10e} 0.0 {st:.6f} 90.0 {0.5*L:.10e} 0.0 {k:d}")
        w = np.clip(np.minimum(s[1:], edges[b + 1]) - np.maximum(s[:-1], edges[b]), 0.0, None)
        for fn, arr in files.items():
            out[fn].append((w[:, None]*arr[idx]).sum(axis=0)/w.sum())
open(os.path.join(dst, gname), "w").write("\n".join(rows) + "\n")
for fn in files:
    np.savetxt(os.path.join(dst, fn), np.array(out[fn]), fmt="%.6e")
gn = np.loadtxt(os.path.join(dst, gname), ndmin=2)
print(f"coarsen_2d: {len(g)} elements in {len(np.unique(grp))} groups -> {len(gn)} elements of "
      f"{2*gn[:, 5].min():.1f} to {2*gn[:, 5].max():.1f} m (target {h:.1f} m)")
