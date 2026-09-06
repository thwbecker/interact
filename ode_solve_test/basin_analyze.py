#!/usr/bin/env python3
"""
analysis of basin_rs2 output grids: class/period tallies, minority
basin fractions, Grebogi-style uncertainty exponent (globally and per
radial annulus), and, for eigenplane grids, the log-polar radial
autocorrelation test of spiral self-similarity against the linear
prediction Delta ln r = 2 pi Re(lambda) / Im(lambda).

usage: basin_analyze.py <prefix> [<prefix> ...]

reads <prefix>.*.dat rank files; grid geometry and knd are taken from
the data and header. requires numpy and scipy.
"""
import sys
import glob
import numpy as np
from scipy.ndimage import map_coordinates

B1, B2, RHO = 1.0, 0.84, 0.048


def eig_growth(knd):
    """ln of the once-around radial growth factor of the fixed point"""
    kcr = ((B1-1) + RHO*(2*B1+(B2-1)*(2+RHO)) +
           np.sqrt(4*RHO*RHO*((B1-1)+B2) + ((B1-1)+RHO*RHO*(B2-1))**2))/(2+2*RHO)
    k = knd*kcr
    J = np.array([[B1-1-k+RHO*B2, 1, RHO-1], [-k, 0, 0], [-RHO*B2, 0, -RHO]])
    w = np.linalg.eig(J)[0]
    wl = w[np.argmax(w.real)]
    if abs(wl.imag) < 1e-12:
        return np.nan
    return 2*np.pi*wl.real/abs(wl.imag)


def load(prefix):
    files = sorted(glob.glob(prefix + ".*.dat"))
    if not files:
        sys.exit("no files for " + prefix)
    knd, eig = np.nan, 0
    with open([f for f in files if f.endswith(".0.dat")][0]) as f:
        for line in f:
            if not line.startswith("#"):
                break
            t = line.split()
            if "knd" in t:
                knd = float(t[t.index("knd")+1])
            if "eig_plane" in t:
                eig = int(t[t.index("eig_plane")+1])
    dat = np.vstack([np.loadtxt(f) for f in files])
    x = np.unique(dat[:, 0]); y = np.unique(dat[:, 1]); n = len(x)
    ii = np.searchsorted(x, dat[:, 0]); jj = np.searchsorted(y, dat[:, 1])
    P = np.full((len(y), n), -9, int)
    P[jj, ii] = np.where(dat[:, 6] == 1, -1,
                         np.where(dat[:, 6] == 2, -9, dat[:, 5])).astype(int)
    return x, y, P, knd, eig


def uncertainty(P, ok, mask, ds):
    f = []
    for d in ds:
        diff = tot = 0
        for A, O, M in ((P, ok, mask), (P.T, ok.T, mask.T)):
            m = O[:, :-d] & O[:, d:] & M[:, :-d]
            diff += np.sum((A[:, :-d] != A[:, d:]) & m)
            tot += np.sum(m)
        f.append(diff/tot if tot else np.nan)
    return np.array(f)


for prefix in sys.argv[1:]:
    x, y, P, knd, eig = load(prefix)
    ok = P != -9
    print("== %s (knd %.6f, eig_plane %d, %d x %d)" %
          (prefix, knd, eig, len(x), len(y)))
    vals, cnts = np.unique(P[ok], return_counts=True)
    print("   classes: " + ", ".join(
        "%s %d" % ("chaotic" if v == -1 else "period-%d" % v, c)
        for v, c in zip(vals, cnts)))
    X, Y = np.meshgrid(x, y)
    R = np.hypot(X, Y)
    rmax = min(x[-1], y[-1])
    sel = ok & (R > 0.5*rmax)
    print("   minority basin fraction (r > %.1e): %.3f%%" %
          (0.5*rmax, 100*np.mean(P[sel] != 4)))
    ds = np.array([1, 2, 4, 8])
    dgrid = x[1] - x[0]
    print("   uncertainty fraction f(d), d in pixels of %.3e:" % dgrid)
    fg = uncertainty(P, ok, np.ones_like(ok, bool), ds)
    a = np.polyfit(np.log(ds), np.log(fg), 1)[0]
    print("     global:  f = %s  alpha = %.2f" %
          (" ".join("%.4f" % v for v in fg), a))
    for r0, r1 in ((0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0)):
        fa = uncertainty(P, ok, (R >= r0*rmax) & (R < r1*rmax), ds)
        aa = np.polyfit(np.log(ds), np.log(np.maximum(fa, 1e-9)), 1)[0]
        print("     r/rmax [%.2f, %.2f]: f = %s  alpha = %.2f" %
              (r0, r1, " ".join("%.4f" % v for v in fa), aa))
    # spiral self-similarity (eigenplane grids)
    lng = eig_growth(knd)
    if eig and np.isfinite(lng):
        Bm = (P != 4).astype(float)   # minority basin indicator
        nth, nlr = 360, 240
        lr = np.linspace(np.log(0.05*rmax), np.log(rmax), nlr)
        th = np.linspace(0, 2*np.pi, nth, endpoint=False)
        Rg, Tg = np.meshgrid(np.exp(lr), th)
        ci = (Rg*np.cos(Tg) - x[0])/(x[1]-x[0])
        cj = (Rg*np.sin(Tg) - y[0])/(y[1]-y[0])
        LP = map_coordinates(Bm, [cj, ci], order=0)
        LPm = LP - LP.mean(axis=1, keepdims=True)
        dlr = lr[1]-lr[0]
        maxlag = int(min(2.2*lng, lr[-1]-lr[0]-dlr)/dlr)
        ac = np.array([np.mean(np.sum(LPm[:, :-l]*LPm[:, l:], axis=1)/(nlr-l))
                       for l in range(1, maxlag)])
        lags = np.arange(1, maxlag)*dlr
        w = (lags > 0.5*lng) & (lags < 1.5*lng)
        ipk = np.where(w)[0][np.argmax(ac[w])]
        print("   spiral test: predicted Delta ln r = %.3f, "
              "autocorr peak at %.3f (norm. value %.3f)" %
              (lng, lags[ipk], ac[ipk]/ac[0]))
    print()
