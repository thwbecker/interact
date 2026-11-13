#!/usr/bin/env python3
"""
slip_budget.py [rundir ...]: Cattania's (2019) slip budget at the loading
point for runs of the Cattania setup that have the station file
fltst_load.dat (slip of the first creeping cell below the velocity-
weakening part; written by run_antiplane since stations were added).

Cattania's argument: an event nucleates each time the creeping region has
slipped S_n since the previous event, and the whole fault breaks once it
has slipped S_full since the last full rupture, so that the number of
events per full-rupture cycle is alpha = S_full/S_n. Measured here from
the slip S(t) at the loading point at the onsets of the (merged) events,
after the first quarter of the run:

  S_1     median slip from a full rupture (>= 0.9 W) to the next event,
          the first nucleation after a full rupture (her S_n)
  S_full  median slip between consecutive full ruptures
  S_ev    median slip between consecutive events of any size (mixes
          intervals of different kinds once partial ruptures occur)

Both are given in units of s* = pi (b-a) sigma L_inf / G, the form of her
S_n = pi dtau L_n / mu' with dtau and L_n replaced by (b-a) sigma and L_inf
(a scale, not her estimate, which needs the stress drop and L_n of
each event).

Writes slip_budget.txt and slip_budget.png in the first directory given
(or "."), from all cattania_*_ab*_X* run directories below it.
"""
import glob
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_antiplane as pa  # noqa: E402

YR = 365.25 * 86400.0
G, SIG = 30e9, 50e6
LAWS = (("aging", "C0", "o"), ("slip", "C3", "s"), ("two", "C2", "^"))


def budget(d):
    p, r = pa.analyze(d)
    st = np.loadtxt(os.path.join(d, "fltst_load.dat"), ndmin=2)
    t, slip = st[:, 0] / YR, st[:, 1]
    on = r["on"]
    full = r["ln"] >= 0.9 * p["W"]
    s_on = np.interp(on, t, slip)
    linf = p["lscale"]
    sstar = math.pi * (p["b"] - p["a"]) * SIG * linf / G
    dn = np.diff(s_on)
    df = np.diff(s_on[full])
    ifull = np.nonzero(full[:-1])[0]
    d1 = s_on[ifull + 1] - s_on[ifull]
    return dict(p=p, n=len(on), nfull=int(full.sum()),
                sev=np.median(dn) / sstar if len(dn) else np.nan,
                sn=np.median(d1) / sstar if len(d1) else np.nan,
                sfull=np.median(df) / sstar if len(df) else np.nan,
                sstar=sstar)


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
    res = {}
    for d in sorted(glob.glob(os.path.join(rd, "**", "cattania_*_ab*_X*"), recursive=True)):
        if not os.path.exists(os.path.join(d, "fltst_load.dat")):
            continue
        try:
            ll = open(os.path.join(d, "rsf_solve.log")).read().splitlines()
            if not ll or not ll[-1].startswith("wall"):
                continue
            b = budget(d)
        except Exception as e:          # noqa: BLE001
            print("skipping", d, e)
            continue
        name = os.path.basename(d)
        law = name.split("_")[1]
        ab = float(name.split("_ab")[1].split("_X")[0])
        res[(law, ab, b["p"]["X"], name)] = b
    lines = ["%-6s %6s %6s %6s %5s %8s %9s %8s %8s %10s" % ("law", "a/b", "W/Linf", "events", "full", "S_1/s*",
                                                           "S_full/s*", "S_ev/s*", "ratio", "s* [m]")]
    for k in sorted(res):
        b = res[k]
        lines.append("%-6s %6.3f %6.1f %6i %5i %8.1f %9.1f %8.1f %8.2f %10.3e" % (k[0], k[1], k[2], b["n"],
                     b["nfull"], b["sn"], b["sfull"], b["sev"], b["sfull"] / b["sn"], b["sstar"]))
    open(os.path.join(rd, "slip_budget.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    xx = np.logspace(0.2, 2.5, 50)
    for law, col, mk in LAWS:
        for ab, fill in ((0.5, col), (0.625, "w"), (0.75, "0.8")):
            ks = sorted(k for k in res if k[0] == law and k[1] == ab)
            if not ks:
                continue
            x = [k[2] for k in ks]
            kw = dict(marker=mk, color=col, mfc=fill, ls="-", lw=0.8, label="%s, a/b %g" % (law, ab))
            ax[0].plot(x, [res[k]["sn"] for k in ks], **kw)
            ax[1].plot(x, [res[k]["sfull"] for k in ks], **kw)
            ax[2].plot(x, [res[k]["sfull"] / res[k]["sn"] for k in ks], **kw)
    ax[1].plot(xx, 0.5 * np.sqrt(xx), "k:", lw=1, label="slope 1/2")
    ax[2].plot(xx, np.maximum(1, 0.45 * np.sqrt(xx)), "k-", lw=1, label="0.45 (W/L_inf)^1/2")
    for a, t in zip(ax, ("S_1: slip from a full rupture to the next event", "S_full: slip between full ruptures",
                         "S_full / S_1")):
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlabel("W/L_inf")
        a.set_title(t, fontsize=10)
    ax[0].set_ylabel("slip at the loading point / (pi (b-a) sigma L_inf / G)")
    ax[2].legend(fontsize=6)
    ax[1].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "slip_budget.png"), dpi=130)


if __name__ == "__main__":
    main()
