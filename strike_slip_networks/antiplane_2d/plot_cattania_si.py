#!/usr/bin/env python3
"""
plot_cattania_si.py [rundir]: the runs of run_cattania_si.py (Cattania
setup with her S1 settings, b 0.02 and 0.01, aging, slip and two-state
slip laws), one column per law:

  row 1  W/L_inf against a/b, events per full-rupture cycle by colour,
         circles b = 0.02, diamonds b = 0.01, white: no events (as in her
         Figure 2a)
  row 2  events per full-rupture cycle against W/L_inf with her estimate
         0.45 (W/L_inf)^1/2 (her Figure 2b)
  row 3  CV of the interevent times against W/L_inf (her Figure 4b)

Events as in plot_antiplane.py (merged within one hour, first quarter
discarded). Writes cattania_si.png and cattania_si.txt.
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

LAWS = (("aging", "aging law"), ("slip", "slip law"), ("two", "two state variables, slip law"))
MK = {0.02: "o", 0.01: "D"}
ABCOL = {0.3: "C7", 0.5: "C0", 0.625: "C1", 0.75: "C2"}


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
    res = {}
    for d in sorted(glob.glob(os.path.join(rd, "si_*_b*_ab*_X*"))):
        try:
            ll = open(os.path.join(d, "rsf_solve.log")).read().splitlines()
        except OSError:
            continue
        if not ll or not ll[-1].startswith("wall"):
            continue
        name = os.path.basename(d)
        law = name.split("_")[1]
        b = float(name.split("_b")[1].split("_ab")[0])
        ab = float(name.split("_ab")[1].split("_X")[0])
        x = float(name.split("_X")[1])
        try:
            p, r = pa.analyze(d)
        except Exception:               # noqa: BLE001
            continue
        res[(law, b, ab, x)] = pa.stats(p, r)
    lines = ["%-6s %5s %6s %6s %6s %5s %6s %6s" % ("law", "b", "a/b", "W/Linf", "events", "full", "alpha", "CV")]
    for k in sorted(res):
        s = res[k]
        lines.append("%-6s %5.2f %6.3f %6.0f %6i %5i %6.2f %6.3f" % (k + (s["n"], s["nfull"], s["alpha"], s["cv"])))
    open(os.path.join(rd, "cattania_si.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    fig, ax = plt.subplots(3, 3, figsize=(16, 13))
    norm = matplotlib.colors.LogNorm(1, 10)
    xx = np.logspace(0, 1.6, 50)
    for j, (law, lab) in enumerate(LAWS):
        ks = [k for k in res if k[0] == law]
        a = ax[0, j]
        a.axvspan(3, 7, color="0.9")
        a.axhline(0.378, color="k", lw=0.8)
        for b, mk in MK.items():
            kk = [k for k in ks if k[1] == b]
            seis = [k for k in kk if res[k]["n"] > 0 and np.isfinite(res[k]["alpha"])]
            asei = [k for k in kk if res[k]["n"] == 0]
            yoff = 0.012 if b == 0.01 else -0.012
            sc = a.scatter([k[3] for k in seis], [k[2] + yoff for k in seis],
                           c=[res[k]["alpha"] for k in seis], cmap="viridis_r", norm=norm, marker=mk, s=70,
                           edgecolor="k", zorder=3)
            a.scatter([k[3] for k in asei], [k[2] + yoff for k in asei], c="w", marker=mk, s=70, edgecolor="k",
                      zorder=3)
        plt.colorbar(sc, ax=a, label="events per full-rupture cycle")
        a.set_ylim(0.2, 0.85)
        a.set_ylabel("a/b")
        a.set_title(lab, fontsize=11)
        ax[1, j].plot(xx, np.maximum(1, 0.45 * np.sqrt(xx)), "k-", lw=1, label="0.45 (W/L_inf)^1/2")
        for i in (1, 2):
            ax[i, j].axvspan(3, 7, color="0.9")
        for b, mk in MK.items():
            for ab, col in ABCOL.items():
                kk = sorted(k for k in ks if k[1] == b and k[2] == ab and res[k]["n"] > 2)
                if not kk:
                    continue
                kw = dict(marker=mk, color=col, ls="-" if b == 0.02 else "--", lw=0.8, mfc=col if b == 0.02 else "w",
                          label="b %g, a/b %g" % (b, ab))
                ax[1, j].plot([k[3] for k in kk], [res[k]["alpha"] for k in kk], **kw)
                ax[2, j].plot([k[3] for k in kk], [res[k]["cv"] for k in kk], **kw)
        ax[1, j].set_yscale("log")
        ax[1, j].set_ylim(0.8, 6)
        ax[1, j].set_ylabel("events per full-rupture cycle")
        ax[2, j].set_ylabel("CV of interevent times")
        ax[2, j].set_ylim(-0.05, 1.2)
        ax[1, j].legend(fontsize=6, ncol=2)
        for i in range(3):
            ax[i, j].set_xscale("log")
            ax[i, j].set_xlim(0.7, 45)
            ax[i, j].set_xticks([1, 2, 4, 8, 16, 32])
            ax[i, j].set_xticklabels(["1", "2", "4", "8", "16", "32"])
            ax[i, j].minorticks_off()
            ax[i, j].set_xlabel("W/L_inf")
    fig.text(0.01, 0.003, "Settings of Cattania's S1: events above 0.1 m/s, min(L_inf, L_b)/dx = 7.7 (two state "
             "variables: G dc2/(b2 sigma)/32), creeping region of depth W. Grey band: her onset W/L_inf 5 +- 2; "
             "line in row 1: a/b 0.378.", fontsize=8, color="0.4")
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    fig.savefig(os.path.join(rd, "cattania_si.png"), dpi=110)


if __name__ == "__main__":
    main()
