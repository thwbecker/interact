#!/usr/bin/env python3
"""
complexity.py: measures of earthquake-cycle complexity, following
Sun & Becker (2026, GJI 247) and Sun, Gilpin & Becker (submitted):

  recurrence     intervals, coefficient of variation and its class,
                 cluster-N count of repeated intervals, return map
  extrema        distinct extrema of a stress time series (prominence and
                 small-swing filtering, tolerance for equal values)
  perm_entropy   topological permutation entropy (log2 of the number of
                 distinct ordinal patterns) of an event sequence
  d2             Grassberger-Procaccia correlation dimension with a
                 Theiler window
  twonn          TwoNN intrinsic dimension (Facco et al. 2017)
  resample       states at random times by monotone cubic interpolation
                 along the trajectory, against the over-weighting of
                 densely sampled coseismic phases
  pca_pr         PCA of a double-centred field (time x cells) and the
                 participation ratio

usage:
  complexity.py --test            validation on the Lorenz system and a
                                  closed orbit
  complexity.py rundir [...]      measures for rsf_solve run directories
                                  written by run_model (series.bin and
                                  catalogs), first fault group only
"""
import argparse
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks
from scipy.spatial.distance import pdist

SEC_PER_YEAR = 365.25 * 86400.0


def cv_class(cv):
    """classes of Sun & Becker (2026)"""
    if cv < 0.1:
        return "periodic"
    if cv < 0.3:
        return "semiperiodic"
    if cv < 0.6:
        return "irregular"
    return "highly irregular"


def recurrence(onsets):
    """intervals, CV, CV class and cluster-N (intervals within 5 percent of
    the median counted as the same value)"""
    dt = np.diff(np.sort(onsets))
    if len(dt) < 2:
        return dict(n=len(dt), mean=np.nan, cv=np.nan, cls="n/a", cluster_n=0, intervals=dt)
    cv = dt.std(ddof=1) / dt.mean()
    s = np.sort(dt)
    cluster_n = 1 + int(np.sum(np.diff(s) > 0.05 * np.median(dt)))
    return dict(n=len(dt), mean=dt.mean(), cv=cv, cls=cv_class(cv), cluster_n=cluster_n, intervals=dt)


def extrema(t, s, prom_frac=0.05, swing_frac=0.2):
    """local maxima and minima of s with a prominence of prom_frac of the
    range; alternating extremum pairs whose swing is below swing_frac of
    the range are removed.  Returns (times, values, kind) with kind +1
    for maxima, -1 for minima, in time order"""
    rng = s.max() - s.min()
    ip, _ = find_peaks(s, prominence=prom_frac * rng)
    im, _ = find_peaks(-s, prominence=prom_frac * rng)
    idx = np.r_[ip, im]
    kind = np.r_[np.ones(len(ip)), -np.ones(len(im))]
    o = np.argsort(idx)
    idx, kind = idx[o], kind[o]
    keep = np.ones(len(idx), dtype=bool)
    for i in range(len(idx) - 1):
        if keep[i] and keep[i + 1] and kind[i] != kind[i + 1] and \
           abs(s[idx[i]] - s[idx[i + 1]]) < swing_frac * rng:
            keep[i] = keep[i + 1] = False
    idx, kind = idx[keep], kind[keep]
    return t[idx], s[idx], kind


def n_distinct(v, tol):
    """number of distinct values, values closer than tol counted as one"""
    if len(v) == 0:
        return 0
    s = np.sort(v)
    return 1 + int(np.sum(np.diff(s) > tol))


def perm_entropy(seq, order=6, tol=0.0):
    """topological permutation entropy: log2 of the number of distinct
    ordinal patterns of length `order`; values closer than tol rank equal"""
    seq = np.asarray(seq, dtype=float)
    if len(seq) < order:
        return np.nan
    pats = set()
    for i in range(len(seq) - order + 1):
        w = seq[i:i + order]
        # rank with ties: values within tol share the rank of the smallest
        o = np.argsort(w, kind="stable")
        ranks = np.empty(order, dtype=int)
        r = 0
        for j in range(order):
            if j > 0 and w[o[j]] - w[o[j - 1]] > tol:
                r = j
            ranks[o[j]] = r
        pats.add(tuple(ranks))
    return np.log2(len(pats))


def resample(t, fields, n, rng, t0=None, t1=None):
    """fields (list of 1-D arrays on t) at n random times in [t0, t1],
    monotone cubic interpolation along the trajectory"""
    t0 = t[0] if t0 is None else t0
    t1 = t[-1] if t1 is None else t1
    tr = np.sort(rng.uniform(t0, t1, n))
    # drop repeated time stamps (restarts, event-located steps)
    keep = np.r_[True, np.diff(t) > 0]
    out = [PchipInterpolator(t[keep], f[keep])(tr) for f in fields]
    return tr, np.column_stack(out)


def d2(points, times=None, theiler=0.0, nfit=(0.03, 0.3), nr=40):
    """Grassberger-Procaccia correlation dimension: slope of log C(r)
    against log r, C the fraction of pairs (coordinates standardized)
    closer than r, excluding pairs with |t_i - t_j| < theiler.  The fit
    range runs from the nfit[0]**2 to the nfit[1]**2 quantile of the pair
    distances; with the defaults the Lorenz system gives 2.04 (published
    2.05) and a closed orbit 1.00 for 6000 points.  Returns (D2, r, C)"""
    x = (points - points.mean(axis=0)) / points.std(axis=0)
    d = pdist(x)
    if times is not None and theiler > 0:
        dtp = pdist(times[:, None])
        d = d[dtp >= theiler]
    d = d[d > 0]
    lo, hi = np.quantile(d, nfit[0] ** 2), np.quantile(d, nfit[1] ** 2)
    r = np.logspace(np.log10(lo), np.log10(hi), nr)
    ds = np.sort(d)
    c = np.searchsorted(ds, r) / len(ds)
    ok = c > 0
    slope = np.polyfit(np.log(r[ok]), np.log(c[ok]), 1)[0]
    return slope, r, c


def twonn(points, discard=0.1):
    """TwoNN intrinsic dimension (Facco et al. 2017) from the ratio of
    second to first nearest-neighbour distances"""
    from scipy.spatial import cKDTree
    x = (points - points.mean(axis=0)) / points.std(axis=0)
    dd, _ = cKDTree(x).query(x, k=3)
    r1, r2 = dd[:, 1], dd[:, 2]
    ok = r1 > 0
    mu = np.sort(r2[ok] / r1[ok])
    n = len(mu)
    m = int(n * (1 - discard))
    f = np.arange(1, n + 1) / n
    xx, yy = np.log(mu[:m]), -np.log(1 - f[:m])
    return float(np.sum(xx * yy) / np.sum(xx * xx))


def pca_pr(field):
    """PCA of a (time x cells) field after removing the temporal mean of
    each cell and the spatial mean of each snapshot; returns (variance
    fractions, spatial modes, participation ratio)"""
    f = field - field.mean(axis=0)
    f = f - f.mean(axis=1, keepdims=True)
    u, s, vt = np.linalg.svd(f, full_matrices=False)
    lam = s ** 2
    frac = lam / lam.sum()
    pr = lam.sum() ** 2 / np.sum(lam ** 2)
    return frac, vt, pr


def lorenz_test(n=6000, seed=1):
    from scipy.integrate import solve_ivp

    def rhs(t, u):
        x, y, z = u
        return [10 * (y - x), x * (28 - z) - y, x * y - 8.0 / 3.0 * z]
    s = solve_ivp(rhs, (0, 200), [1, 1, 1], max_step=0.01, rtol=1e-9, atol=1e-9)
    sel = s.t > 20
    rng = np.random.default_rng(seed)
    tr, p = resample(s.t[sel], [s.y[0][sel], s.y[1][sel], s.y[2][sel]], n, rng)
    dl, _, _ = d2(p, tr, theiler=1.0)
    th = np.linspace(0, 200 * np.pi, 20000)
    tr, p = resample(th, [np.cos(th), np.sin(th), 0.3 * np.cos(2 * th)], n, rng)
    dc, _, _ = d2(p, tr, theiler=4 * np.pi)
    return dl, dc, twonn(p)


def run_measures(rdir, tmin_frac=0.25, n=6000, seed=1):
    """measures of one run directory (first fault group): CV and cluster-N
    of the per-group catalog, distinct extrema and PE of mean tau/sigma,
    D2 and TwoNN of (max log10 V, mean tau/sigma, mean log10 theta) on
    random times after tmin_frac of the run, PCA participation ratio of
    log10 V"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from rsf_series import Run, cat_columns
    r = Run(rdir)
    sel = r.group == 0
    t = r.t
    tmin = t[0] + tmin_frac * (t[-1] - t[0])
    cat = r.gcat[0] if r.gcat[0] is not None else r.cat
    gc = cat_columns(None, group=r.gcat[0] is not None)
    on = cat[:, gc["onset"]] if cat is not None else np.zeros(0)
    on = on[on >= tmin]
    rec = recurrence(on)
    mu = (r.tau[:, sel] / r.sigma[:, sel]).mean(axis=1)
    w = t >= tmin
    te, ve, ke = extrema(t[w], mu[w])
    peaks = ve[ke > 0]
    rngm = mu[w].max() - mu[w].min()
    ndist = n_distinct(ve, 3e-5)
    pe = perm_entropy(peaks, 6, 3e-5)
    rng = np.random.default_rng(seed)
    tr, p = resample(t[w], [r.logv[w][:, sel].max(axis=1), mu[w], r.logtheta[w][:, sel].mean(axis=1)], n, rng)
    dim2, _, _ = d2(p, tr, theiler=2 * rec["mean"] if np.isfinite(rec["mean"]) else 0.0)
    dtn = twonn(p)
    # PCA on uniformly resampled log10 V (per cell), 4000 times
    tu = np.linspace(tmin, t[-1], 4000)
    keep = np.r_[True, np.diff(t) > 0]
    lv = np.column_stack([np.interp(tu, t[keep], r.logv[keep][:, i]) for i in np.nonzero(sel)[0]])
    frac, _, pr = pca_pr(lv)
    return dict(run=os.path.basename(os.path.abspath(rdir)), events=rec["n"] + 1, mean_T=rec["mean"],
                cv=rec["cv"], cls=rec["cls"], cluster_n=rec["cluster_n"], n_extrema=ndist,
                pe=pe, d2=dim2, twonn=dtn, pr=pr, stress_range=rngm)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rundirs", nargs="*")
    p.add_argument("--test", action="store_true")
    p.add_argument("--tmin_frac", type=float, default=0.25, help="fraction of the run discarded as spin-up")
    a = p.parse_args()
    if a.test:
        dl, dc, tn = lorenz_test()
        print("Lorenz D2 %.3f (published 2.05); closed orbit D2 %.3f, TwoNN %.3f (both should be 1)" % (dl, dc, tn))
    for d in a.rundirs:
        m = run_measures(d, a.tmin_frac)
        print("%-18s events %4i T %7.1f yr CV %.3f %-16s cluster-N %3i extrema %4i PE %.2f D2 %.2f TwoNN %.2f PR %.1f" %
              (m["run"], m["events"], m["mean_T"], m["cv"], m["cls"], m["cluster_n"], m["n_extrema"],
               m["pe"], m["d2"], m["twonn"], m["pr"]))


if __name__ == "__main__":
    main()
