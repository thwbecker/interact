#!/usr/bin/env python3
"""
cattania_checks.py rundir [rundir ...]: further comparisons of Cattania
(2019) setup runs (aging law, any directory with params.txt and
rsf_catalog.dat; events as in plot_antiplane.py, merged within one hour,
first quarter of each run discarded):

  distributions  survival function of rupture length / W (her Figure 3a,
                 with the slope -1/2 that she derives for large W/L_inf)
                 and probability density of the normalized interevent
                 time lambda*dt (her Figure 4a, lambda the mean event
                 rate), with an exponential (Poisson) and a 1/x decay
                 for reference
  full threshold events per full-rupture cycle when a full rupture is
                 defined as at least 0.9, 0.8 or 0.7 W

Writes cattania_distributions.png and cattania_checks.txt in the current
directory.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_antiplane as pa  # noqa: E402


def main():
    dirs = sys.argv[1:]
    runs = []
    for d in dirs:
        try:
            p, r = pa.analyze(d)
        except Exception as e:          # noqa: BLE001
            print("skipping", d, e)
            continue
        runs.append((p, r, os.path.basename(os.path.normpath(d))))
    runs.sort(key=lambda z: (z[0]["a"] / z[0]["b"], z[0]["X"]))
    lines = ["%-34s %5s %6s %6s %16s %16s %16s" % ("run", "a/b", "W/Linf", "events", "alpha(0.9W)", "alpha(0.8W)",
                                                  "alpha(0.7W)")]
    for p, r, name in runs:
        row = "%-34s %5.3f %6.0f %6i" % (name, p["a"] / p["b"], p["X"], len(r["on"]))
        for f in (0.9, 0.8, 0.7):
            nf = int(np.sum(r["ln"] >= f * p["W"]))
            row += " %16s" % ("%.2f (%i full)" % (len(r["on"]) / nf, nf) if nf else "-")
        lines.append(row)
    lines.append("")
    lines.append("rupture lengths / W (sorted) and interevent-time CV:")
    for p, r, name in runs:
        ln = np.sort(r["ln"] / p["W"])
        dt = np.diff(r["on"])
        lines.append("%-34s CV %.3f  partial lengths / W: %s" % (
            name, dt.std(ddof=1) / dt.mean() if len(dt) > 2 else np.nan,
            " ".join("%.2f" % v for v in np.unique(np.round(ln[ln < 0.9], 2)))))
    open("cattania_checks.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    cm = plt.cm.viridis(np.linspace(0, 0.9, max(1, len(runs))))
    bins = np.logspace(-4, 1, 31)
    for i, (p, r, name) in enumerate(runs):
        lab = "W/L_inf %g, a/b %g" % (p["X"], p["a"] / p["b"])
        ln = np.sort(r["ln"] / p["W"])[::-1]
        ax[0].step(ln, np.arange(1, len(ln) + 1) / len(ln), where="post", color=cm[i], label=lab)
        dt = np.diff(r["on"])
        if len(dt) > 2:
            h, e = np.histogram(dt / dt.mean(), bins=bins, density=True)
            c = np.sqrt(e[1:] * e[:-1])
            ok = h > 0
            ax[1].plot(c[ok], h[ok] * 10.0 ** i, "o-", color=cm[i], ms=3, label=lab + " (x 10^%i)" % i)
    x = np.logspace(-2, 0, 10)
    ax[0].plot(x, 0.1 * x ** -0.5, "k:", label="slope -1/2")
    xx = np.logspace(-4, 1, 50)
    ax[1].plot(xx, np.exp(-xx), "k--", lw=0.8, label="exp(-x)")
    ax[1].plot(xx, 0.1 / xx, "k:", lw=0.8, label="1/x")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("rupture length / W")
    ax[0].set_ylabel("N / N_tot")
    ax[0].set_title("rupture lengths (cf. Cattania Fig. 3a)", fontsize=10)
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("lambda dt")
    ax[1].set_ylabel("pdf (offset by powers of 10)")
    ax[1].set_title("normalized interevent times (cf. Cattania Fig. 4a)", fontsize=10)
    for a in ax:
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig("cattania_distributions.png", dpi=130)


if __name__ == "__main__":
    main()
