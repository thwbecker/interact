#!/usr/bin/env python3
"""
compare the BS98 rate-and-state runs (rundir of run_bs98) with the run_checks
resolution, tolerance and initial-condition variants: event statistics
per fault, synchrony in model C, small main-fault events in model B,
and event rasters per variant

usage: compare_checks.py bs98_rundir checks_rundir
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from rsf_series import _load_cat, cat_columns  # noqa: E402

GC = cat_columns(None, group=True)
LAB = {3: ("2", "1", "3"), 5: ("5", "2", "1", "3", "4")}


def load(d):
    """per-fault (onset, M0) arrays from the per-group catalogs, and the cell count"""
    cats = []
    k = 0
    while os.path.exists(os.path.join(d, "rsf_catalog.g%03i.dat" % k)):
        c = _load_cat(os.path.join(d, "rsf_catalog.g%03i.dat" % k))
        if c is None:
            cats.append((np.zeros(0), np.zeros(0)))
        else:
            cats.append((c[:, GC["onset"]], c[:, GC["M0"]] / 1e20))
        k += 1
    with open(os.path.join(d, "geom.in")) as f:
        ncell = sum(1 for _ in f)
    return cats, ncell


def stats(t, m):
    """interval mean and CV and moment range, skipping the start-up event"""
    t = t[t > 10.0]
    m = m[-len(t):] if len(t) else m[:0]
    if len(t) < 3:
        return "n %2i" % len(t)
    dt = np.diff(t)
    return "n %2i T %5.1f cv %.2f M0 %.2f-%.2f" % (len(t), dt.mean(), dt.std() / dt.mean(), m.min(), m.max())


def sync_c(cats):
    """model C: fraction of fault-1 events with events on faults 5 and 4
    (our 0 and 4) within 1 yr; skips the start-up event"""
    t1 = cats[2][0]
    t1 = t1[t1 > 10.0]
    n = 0
    for t in t1:
        ok = True
        for k in (0, 4):
            tk = cats[k][0]
            if len(tk) == 0 or np.min(np.abs(tk - t)) > 1.0:
                ok = False
        n += ok
    return n, len(t1)


def raster(ax, cats, labels, tmax, title):
    for k, (t, m) in enumerate(cats):
        ax.scatter(t, np.full(len(t), k), s=30 * m / 0.5 + 2, c="k", marker="|")
    ax.set_yticks(range(len(cats)))
    ax.set_yticklabels(["fault " + s for s in labels])
    ax.set_ylim(-0.7, len(cats) - 0.3)
    ax.set_xlim(0, tmax)
    ax.set_title(title, fontsize=9, loc="left")


def main():
    base = sys.argv[1]
    var = sys.argv[2]
    singles = (("beta60.5", ("b60.5_res", "b60.5_tol")), ("beta54.6", ("b54.6_res", "b54.6_tol")),
               ("beta29.5", ("b29.5_res", "b29.5_tol")))
    for b, vs in singles:
        for name, d in [(b + " (base)", os.path.join(base, b))] + [(v, os.path.join(var, v)) for v in vs]:
            if not os.path.exists(os.path.join(d, "rsf_catalog.g000.dat")):
                continue
            cats, n = load(d)
            t = cats[0][0]
            dt = np.diff(t[t > 10.0])
            print("%-18s cells %4i events %2i intervals %s last event %.0f yr" %
                  (name, n, len(t), " ".join("%.0f" % x for x in dt), t[-1]))
    models = (("modelA", "A", 500), ("modelB", "B", 1500), ("modelC", "C", 300))
    for mb, mv, tmax in models:
        runs = [("base", os.path.join(base, mb))]
        cand = ["%s_%s" % (mv, v) for v in ("res", "tol", "seed2", "vpert0")]
        if mv == "C":
            cand += ["Cs_vpert0", "Cs_seed1", "Cs_seed2"]
        for v in cand:
            d = os.path.join(var, v)
            if os.path.exists(os.path.join(d, "rsf_catalog.g000.dat")) and \
               os.path.exists(os.path.join(d, "rsf_catalog.dat")):
                runs.append((v, d))
        fig, axs = plt.subplots(len(runs), 1, figsize=(10, 1.3 * len(runs) * (1.6 if mv == "C" else 1.0) + 0.6),
                                sharex=True, squeeze=False)
        print("model %s" % mv)
        for i, (v, d) in enumerate(runs):
            cats, n = load(d)
            labels = LAB[len(cats)]
            line = "  %-7s cells %4i | " % (v, n) + " | ".join(
                "%s: %s" % (labels[k], stats(*cats[k])) for k in range(len(cats)))
            extra = ""
            if mv == "C":
                ns, nt = sync_c(cats)
                extra = " | fault 1 events with 5 and 4 within 1 yr: %i of %i" % (ns, nt)
                pairs = []
                for ka, kb in ((0, 4), (1, 3)):
                    ta, tb = cats[ka][0], cats[kb][0]
                    if len(ta) == len(tb) and len(ta):
                        pairs.append("%.3g" % np.abs(ta - tb).max())
                    else:
                        pairs.append("n %i/%i" % (len(ta), len(tb)))
                extra += " | max |t5-t4|, |t2-t3| [yr]: " + ", ".join(pairs)
            if mv == "B":
                t, m = cats[1]
                extra = " | main-fault events with M0 < 0.25e20: %i" % np.sum((m < 0.25) & (t > 10.0))
            print(line + extra)
            raster(axs[i, 0], cats, labels, tmax, "model %s, %s (%i cells)" % (mv, v, n))
        axs[-1, 0].set_xlabel("time [yr]")
        fig.tight_layout()
        fig.savefig(os.path.join(var, "checks_model%s.png" % mv), dpi=120)


if __name__ == "__main__":
    main()
