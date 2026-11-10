#!/usr/bin/env python3
"""
compare_papers.py [rundir]: quantitative comparison of the
run_parameter_maps runs with Wu & Chen (2014) and Cattania (2019), per
state evolution law, from the classification of plot_parameter_maps.py.

Wu & Chen: their regime boundaries W/h* = 0.6, 1.2, 5.5, 33, 105 with
W/h* = (W/hc)(1 - a/b) predict a regime for every run (stable, slow,
centre, side, narrow, wide); reported are the fraction of runs in the
predicted regime, the mismatches, and per b/a the W/h* interval (between
neighbouring grid points) of each transition seen here, against their
value.

Cattania: her estimate alpha = 0.45 (W/L_inf)^1/2 events per full-rupture
cycle (at least 1) and her onset of partial ruptures at W/L_inf = 5 +- 2;
reported are alpha(model)/alpha(estimate) per run and per a/b the
W/L_inf interval in which partial ruptures (alpha > 1.2) set in.

Writes comparison.txt, comparison_wuchen.png, comparison_cattania.png.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_parameter_maps as pm  # noqa: E402

ORDER = ("stable", "slow", "centre", "side", "narrow", "wide")
# a periodic run without series (nucleation point unknown) matches either
# periodic regime
EQUIV = {"periodic": {"centre", "side"}}
BOUNDS = (0.6, 1.2, 5.5, 33.0, 105.0)
LAWCOL = {"aging": "C0", "slip": "C3", "prz": "C2"}


def wc_predict(x, aob):
    """predicted regime(s): both neighbours if W/h* lies within 10 percent
    of one of their boundaries"""
    wh = x * (1.0 - aob)
    i = int(np.searchsorted(BOUNDS, wh))
    pred = {ORDER[i]}
    for j, b in enumerate(BOUNDS):
        if abs(np.log(wh / b)) < np.log(1.1):
            pred |= {ORDER[j], ORDER[j + 1]}
    return pred, wh


def collect(rd):
    import glob
    res = {}
    for d in sorted(glob.glob(os.path.join(rd, "*_ab*_X*"))):
        log = os.path.join(d, "rsf_solve.log")
        if not os.path.exists(log):
            continue
        ll = open(log).read().splitlines()
        if not ll or not ll[-1].startswith("wall"):
            continue
        name = os.path.basename(d)
        model, law = name.split("_")[:2]
        aob = float(name.split("_ab")[1].split("_X")[0])
        try:
            p, s, cat, maxlv = pm.classify(d)
        except Exception:               # noqa: BLE001
            continue
        res[(model, law, aob, p["X"])] = (s, cat)
    return res


def transitions(xs, cats):
    """W/hc intervals of the first change into each regime group"""
    rank = {"stable": 0, "slow": 0, "centre": 1, "periodic": 1, "side": 2, "narrow": 3, "wide": 3}
    out = {}
    for lev, name in ((1, "seismic"), (2, "side"), (3, "partial")):
        lo, hi = None, None
        for x, c in zip(xs, cats):
            if rank[c] >= lev:
                hi = x
                break
            lo = x
        out[name] = (lo, hi)
    return out


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
    res = collect(rd)
    laws = [l for l in ("aging", "slip", "prz") if any(k[1] == l for k in res)]
    L = []
    # Wu & Chen
    L.append("Wu & Chen (2014) setup: regime predicted by their boundaries against the regime found")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    ref = {"seismic": 1.2, "side": 5.5, "partial": 33.0}
    for law in laws:
        ks = sorted(k for k in res if k[0] == "wuchen" and k[1] == law)
        if not ks:
            continue
        match = [bool(({res[k][1]} | EQUIV.get(res[k][1], set())) & wc_predict(k[3], k[2])[0]) for k in ks]
        L.append("  %s law: %i of %i runs in the predicted regime (%.0f %%; within 10 %% of a boundary "
                 "either neighbour counts)" % (law, sum(match), len(ks),
                 100.0 * np.mean(match)))
        for k, m in zip(ks, match):
            if not m:
                L.append("    b/a %.2f W/hc %5.1f (W/h* %6.2f): predicted %-14s found %s" %
                         (1 / k[2], k[3], wc_predict(k[3], k[2])[1], "/".join(sorted(wc_predict(k[3], k[2])[0])),
                          res[k][1]))
        for aob in sorted(set(k[2] for k in ks)):
            kk = sorted(k for k in ks if k[2] == aob)
            tr = transitions([k[3] for k in kk], [res[k][1] for k in kk])
            for j, name in enumerate(("seismic", "side", "partial")):
                lo, hi = tr[name]
                f = 1.0 - aob
                if hi is None:
                    txt = "not reached up to W/h* %.2g" % (kk[-1][3] * f)
                elif lo is None:
                    txt = "already at the smallest W/h*, %.2g" % (hi * f)
                else:
                    txt = "between W/h* %.2g and %.2g" % (lo * f, hi * f)
                L.append("    b/a %.2f, onset of %-7s: %-32s (Wu & Chen %g)" % (1 / aob, name, txt, ref[name]))
                xpos = 1 / aob + {"aging": -0.04, "slip": 0, "prz": 0.04}[law]
                if hi is None:          # not reached: lower bound, arrow up
                    ax[j].annotate("", xy=(xpos, kk[-1][3] * f * 2.0), xytext=(xpos, kk[-1][3] * f),
                                   arrowprops=dict(arrowstyle="->", color=LAWCOL[law], lw=2))
                if hi is not None:
                    y = (lo * f if lo is not None else hi * f / 2.0, hi * f)
                    ax[j].plot([1 / aob + {"aging": -0.04, "slip": 0, "prz": 0.04}[law]] * 2, y, "-",
                               color=LAWCOL[law], lw=4, alpha=0.7, label=law if aob == sorted(set(
                                   k[2] for k in ks))[0] else None)
    for j, name in enumerate(("seismic", "side", "partial")):
        ax[j].axhline(ref[name], color="k", ls="--", lw=1, label="Wu & Chen")
        ax[j].set_yscale("log")
        ax[j].set_xlabel("b/a")
        ax[j].set_title("onset of %s behaviour" % name if name != "side" else "nucleation moves to the ends",
                        fontsize=10)
        ax[j].legend(fontsize=7)
    ax[0].set_ylabel("W/h* interval of the transition\n(arrows: not reached up to that value)")
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "comparison_wuchen.png"), dpi=130)
    # Cattania
    L.append("")
    L.append("Cattania (2019) setup: events per full-rupture cycle against 0.45 (W/L_inf)^1/2")
    fig, ax = plt.subplots(1, len(laws), figsize=(5 * len(laws), 4), sharey=True, squeeze=False)
    xx = np.logspace(0, 1.7, 50)
    for j, law in enumerate(laws):
        ks = sorted(k for k in res if k[0] == "cattania" and k[1] == law and res[k][0]["n"] > 0)
        a = ax[0, j]
        a.plot(xx, np.maximum(1, 0.45 * np.sqrt(xx)), "k-", lw=1, label="Cattania estimate")
        a.axvspan(3, 7, color="0.9")
        for aob, col in zip((0.5, 0.625, 0.75, 0.875), ("C0", "C1", "C2", "C3")):
            kk = [k for k in ks if k[2] == aob]
            a.plot([k[3] for k in kk], [res[k][0]["alpha"] for k in kk], "o-", color=col, label="a/b %g" % aob)
            if kk:
                onset = next((k[3] for k in kk if res[k][0]["alpha"] > 1.2), None)
                before = [k[3] for k in kk if onset is None or k[3] < onset]
                ratios = ["%.2f" % (res[k][0]["alpha"] / max(1.0, 0.45 * np.sqrt(k[3]))) for k in kk]
                otxt = ("not reached up to W/L_inf %g" % kk[-1][3] if onset is None else
                        "between W/L_inf %s and %g" % (max(before) if before else "<", onset))
                L.append("  %s law, a/b %.3f: W/L_inf %s, alpha/estimate %s; onset of partial ruptures %s "
                         "(Cattania 5 +- 2)" % (law, aob, [k[3] for k in kk], ratios, otxt))
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlabel("W/L_inf")
        a.set_title("%s law" % law, fontsize=10)
        a.legend(fontsize=7)
    ax[0, 0].set_ylabel("events per full-rupture cycle")
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "comparison_cattania.png"), dpi=130)
    open(os.path.join(rd, "comparison.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
