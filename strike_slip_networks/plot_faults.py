#!/usr/bin/env python3
"""
plot_faults.py: per-fault diagnostic figure for one run directory.

Top: map view of the fault system (cells colored by a - b), with the
sense of shear on each fault shown by half arrows (from the sign of
the resolved far-field shear rate in params.json) and the S1 axis.
Below: one column per fault, rows
  1  log10 slip rate: mean (line) and +- std (band) over the fault's
     cells, maximum over its cells (thin black)
  2  shear stress tau (mean, band) and normal stress sigma (dashed,
     right axis)
  3  log10 state variable theta (mean, band)
  4  cumulative mean slip (line) and event Mw on the fault (stems,
     right axis, per-fault catalog)
Optionally (--spacetime) a fifth row with log10 V along strike
against time.

usage: plot_faults.py rundir [rundir2 ...] [-o out.png] [--tmin yr] [--tmax yr] [--spacetime]\nSeveral runs are placed side by side (one map and one column per fault each).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rsf_series import Run, cat_columns


def band(ax, t, m, s, color, lw=0.8):
    ax.plot(t, m, color=color, lw=lw)
    ax.fill_between(t, m - s, m + s, color=color, alpha=0.25, lw=0)


def fault_frame(r, k):
    """centre, along-strike unit vector and normal of fault k, and the
    along-strike coordinate of its cells [m]"""
    sel = r.group == k
    xc, yc = r.x[sel].mean(), r.y[sel].mean()
    strike = r.par.get("strikes", [0.0] * r.ngroup)[k]
    phi = np.radians(strike)
    sx, sy = np.sin(phi), np.cos(phi)
    along = (r.x[sel] - xc) * sx + (r.y[sel] - yc) * sy
    return (xc, yc), (sx, sy), (np.cos(phi), -np.sin(phi)), along, sel


def draw_map(ax, r):
    amb = r.a - r.b
    vmax = max(np.abs(amb).max(), 1e-4)
    sc = ax.scatter(r.x / 1e3, r.y / 1e3, c=amb, s=6, cmap="coolwarm", vmin=-vmax, vmax=vmax)
    plt.colorbar(sc, ax=ax, pad=0.01, fraction=0.05).set_label("a - b")
    rates = r.par.get("fault_rates", [])
    Lk = r.par["L"]
    span = 0.12 * Lk / 1e3
    for k in range(r.ngroup):
        (xc, yc), (sx, sy), (nx, ny), along, sel = fault_frame(r, k)
        taud = rates[k]["tau_dot"] if k < len(rates) else 1.0
        sense = 1.0 if taud >= 0 else -1.0      # +: left-lateral
        # left-lateral: standing on either side, the other side moves to
        # the left, i.e. the block on the +n side (right of strike) moves
        # along +s and the block on the -n side along -s
        off = 0.06 * Lk / 1e3
        for side in (+1.0, -1.0):
            px = xc / 1e3 + side * off * nx
            py = yc / 1e3 + side * off * ny
            d = side * sense * span
            ax.annotate("", xy=(px + d * sx, py + d * sy), xytext=(px - d * sx, py - d * sy),
                        arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2))
        lab = "%s, %.0f kPa/yr" % ("LL" if sense > 0 else "RL", abs(taud) * 365.25 * 86400 / 1e3)
        ax.text(xc / 1e3, yc / 1e3 - 2.4 * off, "fault %i (%s)" % (k, lab), ha="center", va="top",
                fontsize=8)
    # S1 axis as a thin double arrow in the lower left, clear of the faults
    az = np.radians(r.par.get("sr_az", 135.0))
    ax.set_xlim(r.x.min() / 1e3 - 0.1 * Lk / 1e3, r.x.max() / 1e3 + 0.1 * Lk / 1e3)
    x0 = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    y0 = r.y.min() / 1e3 - 0.15 * Lk / 1e3
    s1 = 0.5 * span
    ax.annotate("", xy=(x0 + s1 * np.sin(az), y0 + s1 * np.cos(az)),
                xytext=(x0 - s1 * np.sin(az), y0 - s1 * np.cos(az)),
                arrowprops=dict(arrowstyle="<->", color="0.4", lw=1.0))
    ax.text(x0 + 1.2 * s1, y0, "S1 axis", ha="left", va="center", fontsize=8, color="0.4")
    # fill the allotted box at equal scale: x range from the data, y
    # range from the box shape
    xmin, xmax = r.x.min() / 1e3 - 0.1 * Lk / 1e3, r.x.max() / 1e3 + 0.1 * Lk / 1e3
    fig = ax.figure
    bb = ax.get_position()
    w_in = bb.width * fig.get_figwidth()
    h_in = bb.height * fig.get_figheight()
    yhalf = 0.5 * (xmax - xmin) * h_in / w_in
    ymid = 0.5 * (r.y.min() + r.y.max()) / 1e3
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(min(ymid - 1.4 * yhalf, y0 - s1), ymid + 0.8 * yhalf)
    ax.set_aspect("auto")
    ax.set_ylabel("y [km]")
    ax.set_title("map view [km], cells colored by a - b; arrows: sense of shear from the resolved far-field rate",
                 fontsize=9)
    ax.set_xlabel("x [km]", fontsize=8)


def plot_column(ax, r, k, tmin, tmax, col, first_col, last_col, spacetime, title):
    """one column of panels for fault k of run r"""
    cc = cat_columns(None, group=True)
    sel_t = (r.t >= tmin) & (r.t <= tmax)
    t = r.t[sel_t]
    (xc, yc), svec, nvec, along, sel = fault_frame(r, k)
    # 1 slip rate
    f = r.logv[sel_t][:, sel]
    band(ax[0], t, f.mean(axis=1), f.std(axis=1), col)
    ax[0].plot(t, f.max(axis=1), color="k", lw=0.4, alpha=0.7)
    ax[0].set_title(title, color=col, fontsize=10)
    # 2 stresses
    f = r.tau[sel_t][:, sel]
    band(ax[1], t, f.mean(axis=1), f.std(axis=1), col)
    ax2 = ax[1].twinx()
    g = r.sigma[sel_t][:, sel]
    ax2.plot(t, g.mean(axis=1), color="0.3", lw=0.8, ls="--")
    ax2.fill_between(t, g.mean(axis=1) - g.std(axis=1), g.mean(axis=1) + g.std(axis=1),
                     color="0.5", alpha=0.15, lw=0)
    ax2.tick_params(axis="y", labelsize=7, colors="0.3")
    if last_col:
        ax2.set_ylabel("sigma [MPa] (dashed)", color="0.3", fontsize=8)
    # 3 state
    f = r.logtheta[sel_t][:, sel]
    band(ax[2], t, f.mean(axis=1), f.std(axis=1), col)
    # 4 slip and events
    f = r.slip[sel_t][:, sel]
    ax[3].plot(t, np.abs(f.mean(axis=1)), color=col, lw=1.0)
    ax[3].fill_between(t, np.abs(f).min(axis=1), np.abs(f).max(axis=1), color=col, alpha=0.2, lw=0)
    ax4 = ax[3].twinx()
    c = r.gcat[k] if k < len(r.gcat) else None
    if c is not None:
        w = (c[:, cc["onset"]] >= tmin) & (c[:, cc["onset"]] <= tmax)
        ax4.vlines(c[w, cc["onset"]], 0, c[w, cc["Mw"]], color="k", lw=0.8, alpha=0.6)
        ax4.plot(c[w, cc["onset"]], c[w, cc["Mw"]], "ko", ms=2.5)
        mw_all = r.cat[:, cat_columns(None)["Mw"]] if r.cat is not None else c[:, cc["Mw"]]
        ax4.set_ylim(max(0.0, mw_all.min() - 1.0), mw_all.max() + 0.4)
    ax4.tick_params(axis="y", labelsize=7)
    if last_col:
        ax4.set_ylabel("Mw (stems)", fontsize=8)
    # 5 space-time
    if spacetime:
        order = np.argsort(along)
        im = ax[4].pcolormesh(t, along[order] / 1e3, r.logv[sel_t][:, sel][:, order].T, cmap="magma",
                              vmin=-12, vmax=0, shading="nearest", rasterized=True)
        if last_col:
            plt.colorbar(im, ax=ax[4], pad=0.01).set_label("log10 V")
        ax[4].set_ylabel("along strike [km]")
    ax[-1].set_xlabel("time [yr]")
    if first_col:
        ax[0].set_ylabel("log10 V [m/s]")
        ax[1].set_ylabel("tau [MPa]")
        ax[2].set_ylabel("log10 theta [s]")
        ax[3].set_ylabel("|slip| [m] (mean, min-max)")
    return ax2, ax4


def run_title(r, name):
    par = r.par
    rates = par.get("fault_rates", [])
    extra = ""
    if len(rates) == 1:
        extra = ", strike %.0f, sigma_dot %.1f kPa/yr" % (rates[0]["strike"], rates[0]["sigma_dot_kPa_per_yr"])
    return ("%s: %i fault(s), L %.0f km, W %.0f km, a %.4f b %.4f dc %.3f m, L/h* %.1f, "
            "tau_ref %.1f kPa/yr, %i cells, %i events%s" %
            (name, par["nfault"], par["L"] / 1e3, par["W"] / 1e3, par["a"], par["b"], par["dc"],
             par["L_over_hstar"], par["taudot_kPa_per_yr"], par["N"],
             0 if r.cat is None else len(r.cat), extra))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rundirs", nargs="+", help="one or more run directories; runs are placed side by side")
    p.add_argument("-o", "--out", default=None)
    p.add_argument("--tmin", type=float, default=None)
    p.add_argument("--tmax", type=float, default=None)
    p.add_argument("--spacetime", action="store_true")
    args = p.parse_args()

    runs = [Run(d) for d in args.rundirs]
    out = args.out or os.path.join(args.rundirs[0], "faults.png")
    tmin = args.tmin if args.tmin is not None else min(r.t[0] for r in runs)
    tmax = args.tmax if args.tmax is not None else max(r.t[-1] for r in runs)
    nrow = 5 if args.spacetime else 4
    ncol = sum(r.ngroup for r in runs)
    colors = plt.cm.tab10(np.arange(max(ncol, 1)) % 10)

    fig = plt.figure(figsize=(4.6 * ncol + 1.5, 2.3 * nrow + 3.0))
    gs = fig.add_gridspec(nrow + 1, ncol, height_ratios=[1.3] + [1.0] * nrow, hspace=0.45, wspace=0.3,
                          top=0.95, bottom=0.05)
    axes = np.empty((nrow, ncol), dtype=object)
    sig_axes = []
    j = 0
    for r, d in zip(runs, args.rundirs):
        name = os.path.basename(os.path.abspath(d))
        axm = fig.add_subplot(gs[0, j:j + r.ngroup])
        draw_map(axm, r)
        if len(runs) > 1:
            rates = r.par.get("fault_rates", [])
            sub = ", ".join("strike %.0f: tau_dot %.1f, sigma_dot %.1f kPa/yr" %
                            (q["strike"], q["tau_dot_kPa_per_yr"], q["sigma_dot_kPa_per_yr"]) for q in rates)
            axm.set_title("%s (%i events)\n%s" % (name, 0 if r.cat is None else len(r.cat), sub), fontsize=8)
        for k in range(r.ngroup):
            for i in range(nrow):
                axes[i, j] = fig.add_subplot(gs[i + 1, j], sharex=axes[0, 0] if (i or j) else None)
            title = "fault %i" % k if len(runs) == 1 else "%s, fault %i" % (name, k)
            ax2, ax4 = plot_column(axes[:, j], r, k, tmin, tmax, colors[j], j == 0, j == ncol - 1,
                                   args.spacetime, title)
            sig_axes.append(ax2)
            j += 1
    # common y ranges per row across columns, for comparability
    for i in range(min(nrow, 4)):
        lo = min(a.get_ylim()[0] for a in axes[i, :])
        hi = max(a.get_ylim()[1] for a in axes[i, :])
        for a in axes[i, :]:
            a.set_ylim(lo, hi)
    # and a common normal stress axis
    lo = min(a.get_ylim()[0] for a in sig_axes)
    hi = max(a.get_ylim()[1] for a in sig_axes)
    for a in sig_axes:
        a.set_ylim(lo, hi)
    axes[0, 0].set_xlim(tmin, tmax)
    if len(runs) == 1:
        fig.suptitle(run_title(runs[0], os.path.basename(os.path.abspath(args.rundirs[0]))),
                     fontsize=10, y=0.995)
    else:
        par = runs[0].par
        fig.suptitle("L %.0f km, W %.0f km, a %.4f b %.4f dc %.3f m, L/h* %.1f, tau_ref %.1f kPa/yr "
                     "(first run; others may differ)" %
                     (par["L"] / 1e3, par["W"] / 1e3, par["a"], par["b"], par["dc"],
                      par["L_over_hstar"], par["taudot_kPa_per_yr"]), fontsize=10, y=0.995)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
