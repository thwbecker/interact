#!/usr/bin/env python3
"""
plot_slip_space.py [rundir]: complexity measures for all si_<law>_b<b>_
ab<a/b>_X<X> runs in rundir (run_cattania_si.py, run_slip_space.py), one
column per a/b and one line per law (including two-state variants
two-r<rho>-q<b2/b1>) and b:

  row 1  events per full-rupture cycle
  row 2  CV of the interevent times
  row 3  number of distinct rupture lengths (rounded to 0.05 W)

Events as in plot_antiplane.py (merged within one hour, first quarter
discarded, full rupture at least 0.9 W or 0.9 times the longest event).
Writes slip_space.png and slip_space.txt.
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
        s = pa.stats(p, r)
        lmax = r["ln"].max() if len(r["ln"]) else 1.0
        s["nlen"] = len(np.unique(np.round(r["ln"] / min(p["W"], lmax) / 0.05))) if len(r["ln"]) else 0
        res[(law, b, ab, x)] = s
    lines = ["%-16s %5s %6s %6s %6s %5s %6s %6s %5s" % ("law", "b", "a/b", "W/Linf", "events", "full", "alpha", "CV",
                                                      "nlen")]
    for k in sorted(res):
        s = res[k]
        lines.append("%-16s %5.2f %6.3f %6.0f %6i %5i %6.2f %6.3f %5i" % (k + (s["n"], s["nfull"], s["alpha"],
                                                                           s["cv"], s["nlen"])))
    open(os.path.join(rd, "slip_space.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    abs_ = sorted(set(k[2] for k in res))
    series = sorted(set((k[0], k[1]) for k in res))
    cols = {s: plt.cm.tab10(i % 10) for i, s in enumerate(series)}
    fig, ax = plt.subplots(3, len(abs_), figsize=(4.2 * len(abs_), 10), squeeze=False, sharex=True)
    for j, ab in enumerate(abs_):
        for law, b in series:
            ks = sorted(k for k in res if k[0] == law and k[1] == b and k[2] == ab and res[k]["n"] > 2)
            if not ks:
                continue
            x = [k[3] for k in ks]
            kw = dict(marker="o" if b == 0.02 else "D", color=cols[(law, b)], lw=1, label="%s, b %g" % (law, b))
            ax[0, j].plot(x, [res[k]["alpha"] for k in ks], **kw)
            ax[1, j].plot(x, [res[k]["cv"] for k in ks], **kw)
            ax[2, j].plot(x, [res[k]["nlen"] for k in ks], **kw)
        ax[0, j].set_title("a/b %g" % ab, fontsize=11)
        ax[2, j].set_xlabel("W/L_inf")
        for i in range(3):
            ax[i, j].set_xscale("log")
    ax[0, 0].set_ylabel("events per full-rupture cycle")
    ax[1, 0].set_ylabel("CV of interevent times")
    ax[2, 0].set_ylabel("distinct rupture lengths")
    ax[0, -1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "slip_space.png"), dpi=110)


if __name__ == "__main__":
    main()
