#!/usr/bin/env python3
"""
plot_antiplane.py [rundir]: analysis of the 2-D antiplane runs of
run_antiplane_sweep (directories <model>_<law>_X<X>)

Events from rsf_catalog.dat (threshold 1e-3 m/s); events starting within
60 s of the end of the previous one are merged; the first quarter of each
run is discarded. Rupture length = number of cells above the threshold
times the cell size; full rupture: at least 0.9 W. The nucleation point
of an event is the cell of maximum slip rate in the first series record
after onset with max |v| above 1e-3 m/s, given as position along the
velocity-weakening part (0 = top or one end, 1 = bottom or other end).

Writes antiplane_bifurcation.png (interevent time / mean and rupture
length / W against X), antiplane_nucleation.png (nucleation position
against X), antiplane_summary.png (events per full-rupture cycle, CV,
fraction of central nucleation), antiplane_spacetime.png and
antiplane_summary.txt.

X is W/hc for Wu & Chen (2014), with their boundaries at W/h* = 1.2, 5.5,
33, 105, i.e. W/hc = 2.4, 11, 66, 210 for b/a = 2, and W/L_inf for
Cattania (2019), with her estimate alpha = 0.45 (W/L_inf)^1/2 and onset
of partial ruptures at W/L_inf = 5 +- 2.
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

YR = 365.25 * 86400.0
SETS = (("wuchen", "aging", "Wu & Chen model A, aging"), ("wuchen", "slip", "Wu & Chen model A, slip law"),
        ("cattania", "aging", "Cattania (2019), aging"), ("cattania", "slip", "Cattania (2019), slip law"))
COL = {"aging": "C0", "slip": "C3"}


def params(d):
    w = open(os.path.join(d, "params.txt")).read().split()
    return {w[i]: (w[i + 1] if w[i] in ("model", "law") else float(w[i + 1])) for i in range(0, len(w) - 1, 2)}


def read_series(d, nc):
    # memory-mapped: the series of the larger runs do not fit in memory
    fn = os.path.join(d, "series.bin")
    h = np.fromfile(fn, "<i4", 3)
    rec = np.dtype([("t", "<f8"), ("f", "<f4", (h[2], h[1]))])
    nrec = (os.path.getsize(fn) - 12) // rec.itemsize      # whole records only
    s = np.memmap(fn, dtype=rec, mode="r", offset=12, shape=(nrec,))
    return np.array(s["t"]) / YR, s["f"][:, 0, :]


def analyze(d):
    p = params(d)
    c = np.loadtxt(os.path.join(d, "rsf_catalog.dat"), ndmin=2)
    tmax = float(open(os.path.join(d, "rsf_solve.log")).read().split("stop ")[1].split()[0]) \
        if "stop " in open(os.path.join(d, "rsf_solve.log")).read() else None
    t, lv = read_series(d, int(p["nvw"] + p["nvs"]))
    if tmax is None:
        tmax = t[-1]
    if c.size == 0:
        return p, dict(on=np.zeros(0), ln=np.zeros(0), nuc=np.zeros(0))
    on, ar, nc = c[:, 1], c[:, 2], c[:, 4]
    mon, mln, last = [], [], -1e30
    for i in range(len(on)):
        if mon and (on[i] - last) * YR < 60.0:
            mln[-1] = max(mln[-1], nc[i] * p["ds"])
            last = max(last, ar[i])
            continue
        mon.append(on[i])
        mln.append(nc[i] * p["ds"])
        last = ar[i]
    mon, mln = np.array(mon), np.array(mln)
    keep = mon >= 0.25 * t[-1]
    mon, mln = mon[keep], mln[keep]
    nvw = int(p["nvw"])
    nuc = []
    for t0 in mon:
        i = np.searchsorted(t, t0 - 1e-12)
        while i < len(t) and np.max(lv[i]) < -3:
            i += 1
        nuc.append((np.argmax(lv[i][:nvw]) + 0.5) / nvw if i < len(t) else np.nan)
    return p, dict(on=mon, ln=mln, nuc=np.array(nuc), t=t, lv=lv)


def stats(p, r):
    on, ln = r["on"], r["ln"]
    full = ln >= 0.9 * p["W"]
    s = dict(n=len(on), nfull=int(full.sum()), alpha=np.nan, cv=np.nan, cvfull=np.nan,
             ffull=full.mean() if len(on) else np.nan,
             fcenter=np.mean(np.abs(r["nuc"] - 0.5) < 0.25) if len(on) else np.nan)
    if len(on) > 2:
        dt = np.diff(on)
        s["cv"] = dt.std(ddof=1) / dt.mean()
    if full.sum() > 2:
        s["cvfull"] = np.diff(on[full]).std(ddof=1) / np.diff(on[full]).mean()
        s["alpha"] = len(on) / full.sum()
    return s


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
    res = {}
    for d in sorted(glob.glob(os.path.join(rd, "*_X*"))):
        if not os.path.exists(os.path.join(d, "rsf_catalog.dat")):
            continue
        lines_log = open(os.path.join(d, "rsf_solve.log")).read().splitlines()
        if not lines_log or not lines_log[-1].startswith("wall"):
            continue                    # run not finished
        p, r = analyze(d)
        res[(p["model"], p["law"], p["X"])] = (p, r, stats(p, r))
    lines = ["%-9s %-6s %6s %6s %5s %6s %6s %6s %6s %6s" % ("model", "law", "X", "events", "full", "alpha", "CV",
                                                          "CVfull", "f_full", "f_cent")]
    for k in sorted(res):
        s = res[k][2]
        lines.append("%-9s %-6s %6.1f %6i %5i %6.2f %6.3f %6.3f %6.2f %6.2f" % (k[0], k[1], k[2], s["n"], s["nfull"],
                     s["alpha"], s["cv"], s["cvfull"], s["ffull"], s["fcenter"]))
    open(os.path.join(rd, "antiplane_summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    rng = np.random.default_rng(1)

    def bounds(a, model):
        xb = (2.4, 11, 66, 210) if model == "wuchen" else (3, 5, 7)
        for x in xb:
            a.axvline(x, color="0.8", ls=":" if model == "wuchen" else "--", lw=1)

    # bifurcation diagrams
    fig, ax = plt.subplots(2, 4, figsize=(18, 7), sharey="row")
    for j, (model, law, lab) in enumerate(SETS):
        for k in [k for k in res if k[0] == model and k[1] == law]:
            p, r, s = res[k]
            if len(r["on"]) < 3:
                continue
            dt = np.diff(r["on"])
            ax[0, j].plot(k[2] * np.exp(rng.uniform(-0.06, 0.06, len(dt))), dt / dt.mean(), ".", ms=3,
                          color=COL[law], alpha=0.6)
            ax[1, j].plot(k[2] * np.exp(rng.uniform(-0.06, 0.06, len(r["ln"]))), r["ln"] / p["W"], ".", ms=3,
                          color=COL[law], alpha=0.6)
        for i in range(2):
            ax[i, j].set_xscale("log")
            ax[i, j].set_yscale("log")
            bounds(ax[i, j], model)
        ax[0, j].set_title(lab, fontsize=9)
        ax[1, j].set_xlabel("W/hc" if model == "wuchen" else "W/L_inf")
    ax[0, 0].set_ylabel("interevent time / mean")
    ax[1, 0].set_ylabel("rupture length / W")
    fig.text(0.01, 0.005, "dotted: Wu & Chen boundaries (W/hc 2.4, 11, 66, 210); dashed: Cattania onset of partial "
             "ruptures (W/L_inf 5 +- 2)", fontsize=8, color="0.4")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(os.path.join(rd, "antiplane_bifurcation.png"), dpi=130)

    # nucleation position
    fig, ax = plt.subplots(1, 4, figsize=(18, 3.6), sharey=True)
    for j, (model, law, lab) in enumerate(SETS):
        for k in [k for k in res if k[0] == model and k[1] == law]:
            p, r, s = res[k]
            full = r["ln"] >= 0.9 * p["W"]
            xj = k[2] * np.exp(rng.uniform(-0.06, 0.06, len(r["nuc"])))
            ax[j].plot(xj[full], r["nuc"][full], "o", ms=3, mfc="none", color=COL[law])
            ax[j].plot(xj[~full], r["nuc"][~full], ".", ms=3, color="k")
        ax[j].set_xscale("log")
        bounds(ax[j], model)
        ax[j].set_title(lab, fontsize=9)
        ax[j].set_xlabel("W/hc" if model == "wuchen" else "W/L_inf")
    ax[0].set_ylabel("nucleation position / W\n(open: full ruptures, dots: partial)")
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "antiplane_nucleation.png"), dpi=130)

    # summary
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for j, model in enumerate(("wuchen", "cattania")):
        for law in ("aging", "slip"):
            ks = sorted(k for k in res if k[0] == model and k[1] == law)
            ax[j].plot([k[2] for k in ks], [res[k][2]["alpha"] for k in ks], "o-", color=COL[law],
                       label="%s law, events per cycle" % law)
            ax[j].plot([k[2] for k in ks], [res[k][2]["cv"] + 1 for k in ks], "s--", color=COL[law], mfc="none",
                       label="%s law, 1 + CV" % law)
        if model == "cattania":
            xx = np.logspace(0, 2, 50)
            ax[j].plot(xx, np.maximum(1, 0.45 * np.sqrt(xx)), "k-", lw=1, label="Cattania: 0.45 (W/L_inf)^1/2")
        bounds(ax[j], model)
        ax[j].set_xscale("log")
        ax[j].set_yscale("log")
        ax[j].set_xlabel("W/hc" if model == "wuchen" else "W/L_inf")
        ax[j].set_title("Wu & Chen (2014) model A" if model == "wuchen" else "Cattania (2019)", fontsize=10)
        ax[j].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "antiplane_summary.png"), dpi=130)

    # space-time for the largest X of each set
    fig, ax = plt.subplots(1, 4, figsize=(18, 3.8))
    for j, (model, law, lab) in enumerate(SETS):
        ks = sorted(k for k in res if k[0] == model and k[1] == law)
        if not ks:
            ax[j].axis("off")
            continue
        p, r, s = res[ks[-1]]
        lv = r["lv"]
        w = slice(max(0, len(lv) - 4000), len(lv), max(1, len(lv[-4000:]) // 2000))
        im = ax[j].pcolormesh(np.arange(lv[w].shape[0]), (np.arange(lv.shape[1]) + 0.5) * p["ds"], lv[w].T,
                              cmap="magma", vmin=-12, vmax=0, shading="nearest", rasterized=True)
        ax[j].invert_yaxis()
        ax[j].axhline(p["W"], color="w", lw=0.5, ls="--")
        ax[j].set_title("%s, X %g (last 4000 records)" % (lab, ks[-1][2]), fontsize=8)
        ax[j].set_ylabel("position along fault [m]")
    fig.colorbar(im, ax=ax, fraction=0.02, label="log10 V [m/s]")
    fig.savefig(os.path.join(rd, "antiplane_spacetime.png"), dpi=120, bbox_inches="tight")


if __name__ == "__main__":
    main()
