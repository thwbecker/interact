#!/usr/bin/env python3
"""
box-counting estimate of the basin boundary dimension from basin_rs2
output grids. counts boxes of side delta (in pixels) that contain
both outcomes; N(delta) ~ delta^(-D0) for a boundary of box-counting
dimension D0 within the plane. prints N(delta) and local slopes.

for the off-center patch grids (b74e_patch etc.), the pixel scale
lies in the asymptotic regime and a smooth boundary gives D0 = 1;
on origin-centered grids, expect a crossover from D0 ~ 2 at the
pixel scale (sub-resolution winding) toward 1 at coarse delta.

usage: basin_boxcount.py <prefix> [<prefix> ...]
"""
import sys
import glob
import numpy as np


def load(prefix):
    files = sorted(glob.glob(prefix + ".*.dat"))
    if not files:
        sys.exit("no files for " + prefix)
    dat = np.vstack([np.loadtxt(f) for f in files])
    x = np.unique(dat[:, 0]); y = np.unique(dat[:, 1])
    ii = np.searchsorted(x, dat[:, 0]); jj = np.searchsorted(y, dat[:, 1])
    P = np.full((len(y), len(x)), -9, int)
    P[jj, ii] = np.where(dat[:, 6] == 1, -1,
                         np.where(dat[:, 6] == 2, -9, dat[:, 5])).astype(int)
    return P


for prefix in sys.argv[1:]:
    P = load(prefix)
    B = P != 4                  # minority indicator (unclassified rare)
    n0, n1 = P.shape
    print("== %s (%d x %d, minority fraction %.4f)" %
          (prefix, n0, n1, B.mean()))
    deltas = [d for d in (1, 2, 4, 8, 16, 32, 64) if d <= min(n0, n1)//4]
    Ns = []
    for d in deltas:
        m0, m1 = n0//d, n1//d
        blocks = B[:m0*d, :m1*d].reshape(m0, d, m1, d)
        smin = blocks.min(axis=(1, 3))
        smax = blocks.max(axis=(1, 3))
        N = int(np.sum(smin != smax))   # boxes containing both outcomes
        Ns.append(N)
    print("   delta(px)   N(delta)   local D0")
    for i, (d, N) in enumerate(zip(deltas, Ns)):
        if i > 0 and Ns[i-1] > 0 and N > 0:
            D = -np.log(N/Ns[i-1])/np.log(d/deltas[i-1])
            print("   %8d %10d      %.2f" % (d, N, D))
        else:
            print("   %8d %10d" % (d, N))
    good = [(d, N) for d, N in zip(deltas, Ns) if N > 0]
    if len(good) > 2:
        ld = np.log([g[0] for g in good]); lN = np.log([g[1] for g in good])
        print("   overall D0 (fit over all delta): %.2f" %
              (-np.polyfit(ld, lN, 1)[0]))
    print()
