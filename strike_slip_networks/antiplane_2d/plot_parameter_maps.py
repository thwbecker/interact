#!/usr/bin/env python3
"""
plot_parameter_maps.py [rundir]: regime maps from the runs of
run_parameter_maps (directories <model>_<law>_ab<a/b>_X<X>), drawn like
Wu & Chen (2014, their Figure 11: W/hc against b/a, regime by symbol)
and Cattania (2019, her Figure 2a: W/L_inf against a/b, number of
events per full-rupture cycle by colour), one panel per state
evolution law.

Regimes (Wu & Chen style), from the events after the first quarter of
each run (merged within one hour, as in plot_antiplane.py):
  stable     no catalog event (1e-3 m/s) and max slip rate below 10 vpl
  slow       no catalog event, but max slip rate above 10 vpl
  periodic   all events full ruptures (>= 0.9 W), interval CV < 0.05,
             nucleation point unknown (run without series)
  centre     all events full ruptures (>= 0.9 W), interval CV < 0.05,
             nucleation within 0.25 W of the centre
  side       as centre, nucleation farther out
  narrow     partial ruptures with at most three distinct rupture
             lengths (to 0.05 W)
  wide       partial ruptures with more distinct lengths
The curves are Wu & Chen's boundaries W/h* = 0.6, 1.2, 5.5, 33, 105
with h* = hc/(1 - a/b) (this reproduces the curve shapes of their
figure; at b/a = 2, h* = 2 hc).

Writes parameter_map_wuchen.png, parameter_map_cattania.png,
parameter_maps.txt.
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_antiplane as pa  # noqa: E402

CATS = (("stable", "C0", "o"), ("slow", "orange", "o"), ("periodic", "0.5", "o"), ("centre", "red", "o"),
        ("side", "limegreen", "o"),
        ("narrow", "cyan", "o"), ("wide", "magenta", "o"))
LAWS = ("aging", "slip", "prz")


def classify(d):
    p, r = pa.analyze(d)
    s = pa.stats(p, r)
    mon = np.loadtxt(os.path.join(d, "rsf_monitor.dat"), ndmin=2)
    tmon = mon[:, 2]
    maxlv = mon[tmon >= 0.25 * tmon[-1], 4].max()
    vpl = p["vpl"]
    if s["n"] == 0:
        cat = "slow" if maxlv > np.log10(10 * vpl) else "stable"
    elif s["ffull"] == 1.0 and s["cv"] < 0.05:
        if np.all(np.isnan(r["nuc"])):
            cat = "periodic"            # run without series: nucleation point unknown
        else:
            cat = "centre" if np.nanmedian(np.abs(r["nuc"] - 0.5)) < 0.25 else "side"
    else:
        nd = len(np.unique(np.round(r["ln"] / p["W"] / 0.05)))
        cat = "narrow" if nd <= 3 else "wide"
    return p, s, cat, maxlv


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
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
            p, s, cat, maxlv = classify(d)
        except Exception as e:          # noqa: BLE001
            print("skipping", d, e)
            continue
        res[(model, law, aob, p["X"])] = (s, cat, maxlv)
    lines = ["%-9s %-6s %6s %6s %5s %5s %6s %7s %6s %-7s" % ("model", "law", "a/b", "X", "n", "full", "alpha", "CV",
                                                           "maxlv", "regime")]
    for k in sorted(res):
        s, cat, maxlv = res[k]
        lines.append("%-9s %-6s %6.3f %6.1f %5i %5i %6.2f %7.3f %6.2f %-7s" % (k[0], k[1], k[2], k[3], s["n"],
                     s["nfull"], s["alpha"], s["cv"], maxlv, cat))
    open(os.path.join(rd, "parameter_maps.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    laws = [l for l in LAWS if any(k[1] == l for k in res)]

    # Wu & Chen style
    fig, ax = plt.subplots(1, len(laws), figsize=(6 * len(laws), 4.6), squeeze=False)
    ba = np.linspace(1.02, 4.2, 200)
    for j, law in enumerate(laws):
        a = ax[0, j]
        for c, col in zip((0.6, 1.2, 5.5, 33, 105), ("C0", "orange", "red", "limegreen", "cyan")):
            a.plot(c * ba / (ba - 1.0), ba, color=col, lw=1.5)
        for cat, col, mk in CATS:
            ks = [k for k in res if k[0] == "wuchen" and k[1] == law and res[k][1] == cat]
            a.scatter([k[3] for k in ks], [1.0 / k[2] for k in ks], c=col, marker=mk, s=70, edgecolor="k",
                      label=cat, zorder=3)
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlim(0.8, 600)
        a.set_ylim(1.0, 4.3)
        a.set_yticks([1, 1.33, 2, 3, 4])
        a.set_yticklabels(["1", "1.33", "2", "3", "4"])
        a.set_xlabel("W/hc")
        a.set_ylabel("b/a")
        a.set_title("Wu & Chen (2014) setup, %s law" % law, fontsize=10)
    ax[0, 0].legend(fontsize=7, loc="upper right")
    fig.text(0.01, 0.01, "curves: Wu & Chen boundaries W/h* = 0.6, 1.2, 5.5, 33, 105 with h* = hc/(1 - a/b)",
             fontsize=8, color="0.4")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(rd, "parameter_map_wuchen.png"), dpi=130)

    # Cattania style
    fig, ax = plt.subplots(1, len(laws), figsize=(6 * len(laws), 4.6), squeeze=False)
    norm = matplotlib.colors.LogNorm(1, 10)
    for j, law in enumerate(laws):
        a = ax[0, j]
        a.axvspan(3, 7, color="0.9")
        a.axvline(5, color="0.6", lw=1)
        a.axhline(0.378, color="k", lw=0.8)
        ks = [k for k in res if k[0] == "cattania" and k[1] == law]
        seis = [k for k in ks if res[k][0]["n"] > 0 and np.isfinite(res[k][0]["alpha"])]
        asei = [k for k in ks if res[k][0]["n"] == 0]
        sc = a.scatter([k[3] for k in seis], [k[2] for k in seis], c=[res[k][0]["alpha"] for k in seis],
                       cmap="viridis_r", norm=norm, s=80, edgecolor="k", zorder=3)
        a.scatter([k[3] for k in asei], [k[2] for k in asei], c="w", s=80, edgecolor="k", zorder=3)
        a.set_xscale("log")
        a.set_xlim(1, 50)
        a.set_ylim(0.1, 1.0)
        a.set_xlabel("W/L_inf")
        a.set_ylabel("a/b")
        a.set_title("Cattania (2019) setup, %s law" % law, fontsize=10)
        plt.colorbar(sc, ax=a, label="events per full-rupture cycle")
    fig.text(0.01, 0.01, "grey band: Cattania's onset of partial ruptures, W/L_inf = 5 +- 2; line: a/b = 0.378 "
             "(above it L_inf sets the nucleation size); white: no catalog events", fontsize=8, color="0.4")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(rd, "parameter_map_cattania.png"), dpi=130)


if __name__ == "__main__":
    main()
