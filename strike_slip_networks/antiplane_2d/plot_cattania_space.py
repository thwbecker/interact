#!/usr/bin/env python3
"""
plot_cattania_space.py [rundir]: the Cattania (2019) setup for the
aging law, the slip law and two state variables (slip law), from the
runs of run_cattania_space.py, laid out after her figures:

  row 1  W/L_inf against a/b, events per full-rupture cycle by colour
         (her Figure 2a; white: no catalog events)
  row 2  events per full-rupture cycle against W/L_inf with her estimate
         0.45 (W/L_inf)^1/2 (her Figure 2b)
  row 3  coefficient of variation of all interevent times against
         W/L_inf (her Figure 4b)
  row 4  survival function of rupture length / W for the a/b = 0.75 runs
         (her Figure 3a; a/b 0.625 where 0.75 is missing)

Events as in plot_antiplane.py (merged within one hour, first quarter
discarded). Writes cattania_space.png and cattania_space.txt.
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
ABCOL = {0.5: "C0", 0.625: "C1", 0.75: "C2", 0.875: "C3"}


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
    res = {}
    for d in sorted(glob.glob(os.path.join(rd, "cattania_*_ab*_X*"))):
        try:
            ll = open(os.path.join(d, "rsf_solve.log")).read().splitlines()
        except OSError:
            continue
        if not ll or not ll[-1].startswith("wall"):
            continue
        name = os.path.basename(d)
        law = name.split("_")[1]
        ab = float(name.split("_ab")[1].split("_X")[0])
        try:
            p, r = pa.analyze(d)
        except Exception:               # noqa: BLE001
            continue
        res[(law, ab, p["X"])] = (p, r, pa.stats(p, r))
    lines = ["%-6s %6s %6s %6s %5s %6s %6s %8s" % ("law", "a/b", "W/Linf", "events", "full", "alpha", "CV",
                                                  "alpha/est")]
    for k in sorted(res):
        s = res[k][2]
        est = max(1.0, 0.45 * np.sqrt(k[2]))
        lines.append("%-6s %6.3f %6.1f %6i %5i %6.2f %6.3f %8.2f" % (k[0], k[1], k[2], s["n"], s["nfull"],
                     s["alpha"], s["cv"], s["alpha"] / est))
    open(os.path.join(rd, "cattania_space.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    fig, ax = plt.subplots(4, 3, figsize=(16, 16))
    norm = matplotlib.colors.LogNorm(1, 10)
    xx = np.logspace(0.2, 1.9, 50)
    for j, (law, lab) in enumerate(LAWS):
        ks = sorted(k for k in res if k[0] == law)
        a = ax[0, j]
        a.axvspan(3, 7, color="0.9")
        a.axhline(0.378, color="k", lw=0.8)
        seis = [k for k in ks if res[k][2]["n"] > 0 and np.isfinite(res[k][2]["alpha"])]
        asei = [k for k in ks if res[k][2]["n"] == 0]
        sc = a.scatter([k[2] for k in seis], [k[1] for k in seis], c=[res[k][2]["alpha"] for k in seis],
                       cmap="viridis_r", norm=norm, s=90, edgecolor="k", zorder=3)
        a.scatter([k[2] for k in asei], [k[1] for k in asei], c="w", s=90, edgecolor="k", zorder=3)
        a.set_xscale("log")
        a.set_xlim(1.5, 90)
        a.set_ylim(0.3, 0.9)
        a.set_ylabel("a/b")
        a.set_title(lab, fontsize=11)
        plt.colorbar(sc, ax=a, label="events per full-rupture cycle")
        for i in (1, 2):
            ax[i, j].axvspan(3, 7, color="0.9")
        ax[1, j].plot(xx, np.maximum(1, 0.45 * np.sqrt(xx)), "k-", lw=1, label="0.45 (W/L_inf)^1/2")
        for ab in sorted(set(k[1] for k in ks)):
            kk = [k for k in ks if k[1] == ab and res[k][2]["n"] > 2]
            ax[1, j].plot([k[2] for k in kk], [res[k][2]["alpha"] for k in kk], "o-", color=ABCOL.get(ab, "k"),
                          label="a/b %g" % ab)
            ax[2, j].plot([k[2] for k in kk], [res[k][2]["cv"] for k in kk], "o-", color=ABCOL.get(ab, "k"),
                          label="a/b %g" % ab)
        ax[1, j].set_yscale("log")
        ax[1, j].set_ylim(0.8, 12)
        ax[1, j].set_ylabel("events per full-rupture cycle")
        ax[2, j].set_ylabel("CV of interevent times")
        ax[2, j].set_ylim(-0.05, 1.6)
        ab_sf = 0.75 if any(k[1] == 0.75 for k in ks) else 0.625
        cm = plt.cm.viridis(np.linspace(0, 1, 11))
        for i, k in enumerate(sorted(k for k in ks if k[1] == ab_sf)):
            ln = np.sort(res[k][1]["ln"] / res[k][0]["W"])[::-1]
            if len(ln) > 1:
                ax[3, j].step(ln, np.arange(1, len(ln) + 1) / len(ln), where="post", color=cm[i % 11],
                              label="W/L_inf %g" % k[2])
        r_ = np.logspace(-2, 0, 10)
        ax[3, j].plot(r_, 0.1 * r_ ** -0.5, "k:", lw=1)
        ax[3, j].set_xscale("log")
        ax[3, j].set_yscale("log")
        ax[3, j].set_xlabel("rupture length / W (a/b %g)" % ab_sf)
        ax[3, j].set_ylabel("N/N_tot")
        for i in range(3):
            ax[i, j].set_xscale("log")
            ax[i, j].set_xlim(1.5, 90)
            ax[i, j].set_xticks([2, 4, 8, 16, 32, 64])
            ax[i, j].set_xticklabels(["2", "4", "8", "16", "32", "64"])
            ax[i, j].minorticks_off()
            ax[i, j].set_xlabel("W/L_inf")
        ax[1, j].legend(fontsize=7)
        ax[3, j].legend(fontsize=6)
    fig.text(0.01, 0.003, "grey band: Cattania's onset of partial ruptures, W/L_inf 5 +- 2; line in row 1: a/b = "
             "0.378; dotted in row 4: slope -1/2. For two state variables W/L_inf refers to the one-state "
             "L_inf with the same spring-slider critical stiffness.", fontsize=8, color="0.4")
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    fig.savefig(os.path.join(rd, "cattania_space.png"), dpi=110)


if __name__ == "__main__":
    main()
