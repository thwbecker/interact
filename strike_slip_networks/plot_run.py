#!/usr/bin/env python3
"""
plot_run.py: overview figure for one run directory.

Left column: mean (line) and +- one standard deviation (band) over the
cells of each fault of log10 slip rate, shear stress, normal stress and
log10 state variable against time.  Right column: cumulative mean slip
per fault, event moment magnitude against time (per-fault catalog if
present, else the global one), a space-time image of log10 slip rate
along strike for the first fault, and a map of the fault system
colored by a - b.

usage: plot_run.py rundir [out.png] [--tmin yr] [--tmax yr]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rsf_series import Run, cat_columns


def band(ax, t, m, s, color, label):
    ax.plot(t, m, color=color, lw=0.8, label=label)
    ax.fill_between(t, m - s, m + s, color=color, alpha=0.25, lw=0)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rundir")
    p.add_argument("out", nargs="?", default=None)
    p.add_argument("--tmin", type=float, default=None)
    p.add_argument("--tmax", type=float, default=None)
    args = p.parse_args()

    r = Run(args.rundir)
    par = r.par
    out = args.out or os.path.join(args.rundir, "overview.png")
    tmin = args.tmin if args.tmin is not None else r.t[0]
    tmax = args.tmax if args.tmax is not None else r.t[-1]
    sel = (r.t >= tmin) & (r.t <= tmax)
    t = r.t[sel]
    colors = plt.cm.tab10(np.arange(max(r.ngroup, 1)) % 10)

    fig, ax = plt.subplots(4, 2, figsize=(14, 12))
    for a in list(ax[:, 0]) + list(ax[:3, 1]):
        a.sharex(ax[0, 0])
    for fld, lab, a in (("logv", "log10 V [m/s]", ax[0, 0]),
                        ("tau", "tau [MPa]", ax[1, 0]),
                        ("sigma", "sigma [MPa]", ax[2, 0]),
                        ("logtheta", "log10 theta [s]", ax[3, 0])):
        m, s = r.group_stats(fld)
        for k in range(r.ngroup):
            band(a, t, m[k][sel], s[k][sel], colors[k], "fault %i" % k)
        if fld == "logv":
            # maximum over all cells, for the coseismic peaks the mean hides
            a.plot(t, r.logv[sel].max(axis=1), color="k", lw=0.4, alpha=0.6, label="max")
        a.set_ylabel(lab)
    ax[0, 0].legend(loc="upper right", fontsize=8, ncol=min(r.ngroup, 3))
    ax[0, 0].set_title("mean and std over cells per fault")

    # cumulative mean slip per fault
    m, s = r.group_stats("slip")
    for k in range(r.ngroup):
        ax[0, 1].plot(t, m[k][sel], color=colors[k], lw=1.0)
    ax[0, 1].set_ylabel("mean slip [m]")
    ax[0, 1].set_title("cumulative slip and events")

    # event magnitudes
    have_g = any(c is not None for c in r.gcat)
    if have_g:
        cc = cat_columns(None, group=True)
        for k in range(r.ngroup):
            c = r.gcat[k]
            if c is None:
                continue
            w = (c[:, cc["onset"]] >= tmin) & (c[:, cc["onset"]] <= tmax)
            ax[1, 1].vlines(c[w, cc["onset"]], 0, c[w, cc["Mw"]], color=colors[k], lw=1.2)
            ax[1, 1].plot(c[w, cc["onset"]], c[w, cc["Mw"]], "o", color=colors[k], ms=3)
    elif r.cat is not None:
        cc = cat_columns(None)
        c = r.cat
        w = (c[:, cc["onset"]] >= tmin) & (c[:, cc["onset"]] <= tmax)
        ax[1, 1].vlines(c[w, cc["onset"]], 0, c[w, cc["Mw"]], color="k", lw=1.2)
        ax[1, 1].plot(c[w, cc["onset"]], c[w, cc["Mw"]], "ko", ms=3)
    if r.cat is not None:
        cc = cat_columns(None)
        mw = r.cat[:, cc["Mw"]]
        ax[1, 1].set_ylim(max(0.0, mw.min() - 1.0), mw.max() + 0.5)
    ax[1, 1].set_ylabel("Mw")

    # space-time plot of log10 V for fault 0 against record index mapped to time
    g0 = r.group == 0
    yy = r.y[g0] / 1e3
    im = ax[2, 1].pcolormesh(t, yy, r.logv[sel][:, g0].T, cmap="magma", vmin=-12, vmax=0,
                             shading="nearest", rasterized=True)
    ax[2, 1].set_ylabel("fault 0 along strike [km]")
    cb = fig.colorbar(im, ax=ax[2, 1], pad=0.01)
    cb.set_label("log10 V")
    ax[3, 0].set_xlabel("time [yr]")
    ax[2, 1].set_xlabel("time [yr]")
    for a in list(ax[:, 0]) + list(ax[:3, 1]):
        a.set_xlim(tmin, tmax)
    # map view, cells colored by a - b, S1 axis drawn at the origin
    amb = r.a - r.b
    sc = ax[3, 1].scatter(r.x / 1e3, r.y / 1e3, c=amb, s=4, cmap="coolwarm",
                          vmin=-np.abs(amb).max(), vmax=np.abs(amb).max())
    fig.colorbar(sc, ax=ax[3, 1], pad=0.01).set_label("a - b")
    az = np.radians(par.get("sr_az", 135.0))
    span = 0.15 * max(np.ptp(r.x), np.ptp(r.y), 1.0) / 1e3
    ax[3, 1].plot([-span * np.sin(az), span * np.sin(az)], [-span * np.cos(az), span * np.cos(az)],
                  "k-", lw=2, label="S1 axis")
    ax[3, 1].set_aspect("equal", adjustable="datalim")
    ax[3, 1].set_xlabel("x [km]")
    ax[3, 1].set_ylabel("y [km]")
    ax[3, 1].legend(loc="best", fontsize=8)

    fig.suptitle("%s: %i fault(s), L %.0f km, W %.0f km, a %.4f b %.4f dc %.3f m, L/h* %.1f, "
                 "tau_dot %.1f kPa/yr, %i cells, %i events" %
                 (os.path.basename(os.path.abspath(args.rundir)), par["nfault"], par["L"] / 1e3,
                  par["W"] / 1e3, par["a"], par["b"], par["dc"], par["L_over_hstar"],
                  par["taudot_kPa_per_yr"], par["N"], 0 if r.cat is None else len(r.cat)),
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
