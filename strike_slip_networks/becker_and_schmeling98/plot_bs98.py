#!/usr/bin/env python3
"""
plot_bs98.py: figures for the rate-and-state analogues of Becker and
Schmeling (1998, GJI 135, 165-176), from the run directories written
by run_bs98.

usage: plot_bs98.py [rundir]

writes
  bs98_fig4_5.png   mean Coulomb stress (tau - f0 sigma) and moment of
                    the single faults with beta = 45 and 29.5 degrees
                    over the first 550 yr (their Figs 4 and 5)
  bs98_fig6.png     event time against event number for the four
                    single faults, with the prediction of their eq. (11)
                    in the form t_(n+1) = t_n + dmu sigma(t_n)/(tau_dot C(beta)),
                    dmu from the recurrence time of the beta = 45 run
                    (their Fig. 6)
  bs98_fig7_8.png   maps and per-fault moment against time for models A
                    and B (their Figs 7 and 8)
  bs98_fig9_10.png  maps and space-time pattern of slipping cells
                    (|v| > 1 cm/s) for models C and D (their Figs 9, 10)

Fault numbers in the figures follow the paper: for A and B the middle
fault is 1 and the side faults are 2 (west) and 3 (east); for C and D
the order from west to east is 5, 2, 1, 3, 4.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from rsf_series import Run, cat_columns  # noqa: E402

SEC_PER_YEAR = 365.25 * 86400.0
GC = cat_columns(None, group=True)
BETAS = (("beta29.5", 29.5, "s"), ("beta45", 45.0, "^"), ("beta54.6", 54.6, "x"), ("beta60.5", 60.5, "+"))


def c_beta(beta, mu):
    """C(beta) = |sin 2 beta| + mu cos 2 beta, their eq. (5)"""
    b2 = np.radians(2.0 * beta)
    return abs(np.sin(b2)) + mu * np.cos(b2)


def onsets(run, k=0):
    c = run.gcat[k] if run.gcat[k] is not None else run.cat
    if c is None:
        return np.zeros(0), np.zeros(0)
    return c[:, GC["onset"]], c[:, GC["M0"]]


def fig4_5(rd):
    fig, ax = plt.subplots(2, 2, figsize=(11, 5), sharex=True)
    for j, (name, beta) in enumerate((("beta45", 45.0), ("beta29.5", 29.5))):
        r = Run(os.path.join(rd, name))
        f0 = r.par["f0"]
        sel = r.t <= 550.0
        cs = r.tau[sel].mean(axis=1) - f0 * r.sigma[sel].mean(axis=1)
        ax[0, j].plot(r.t[sel], cs, "k-", lw=0.8)
        ax[0, j].set_title("beta = %g deg" % beta)
        t, m0 = onsets(r)
        w = t <= 550.0
        ax[1, j].vlines(t[w], 0, m0[w] / 1e20, color="k", lw=1.2)
        ax[1, j].set_xlabel("time [yr]")
    ax[0, 0].set_ylabel("mean tau - f0 sigma [MPa]")
    ax[1, 0].set_ylabel("M0 [1e20 N m]")
    ymax = max(a.get_ylim()[1] for a in ax[1])
    for a in ax[1]:
        a.set_ylim(0, ymax)
    fig.suptitle("rate-and-state analogue of Becker and Schmeling (1998) Figs 4 and 5", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "bs98_fig4_5.png"), dpi=130)


def fig6(rd):
    r45 = Run(os.path.join(rd, "beta45"))
    tau_dot = r45.par["taudot"] * SEC_PER_YEAR / 1e6       # MPa/yr
    sig0 = r45.par["sigma"] / 1e6
    mu = r45.par["f0"]
    t45, _ = onsets(r45)
    dmu = np.mean(np.diff(t45[1:])) * tau_dot * c_beta(45.0, mu) / sig0
    fig, ax = plt.subplots(figsize=(7, 5))
    lines = ["dmu_eff %.3f from the beta = 45 run" % dmu]
    for name, beta, mk in BETAS:
        r = Run(os.path.join(rd, name))
        t, _ = onsets(r)
        if len(t) == 0:
            continue
        n = np.arange(len(t))
        ax.plot(n, t, mk, color="k", ms=6, mfc="none", label="beta = %g" % beta)
        # prediction from the first event on
        cb = c_beta(beta, mu)
        tp = [t[0]]
        while len(tp) < len(t) + 3:
            sig = sig0 - tau_dot * np.cos(np.radians(2.0 * beta)) * tp[-1]
            if sig <= 0.0:
                break
            tp.append(tp[-1] + dmu * sig / (tau_dot * cb))
        ax.plot(np.arange(len(tp)), tp, "k-", lw=0.7)
        lines.append("beta %g: %i events, last interval %.0f yr (first %.0f)" %
                     (beta, len(t), t[-1] - t[-2] if len(t) > 1 else 0, t[1] - t[0] if len(t) > 1 else 0))
    ax.set_xlabel("event no.")
    ax.set_ylabel("time [yr]")
    ax.set_ylim(0, 2000)
    ax.legend(loc="center right")
    ax.set_title("rate-and-state analogue of Becker and Schmeling (1998) Fig. 6\n"
                 "symbols: model, lines: their eq. (11) with dmu from beta = 45", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "bs98_fig6.png"), dpi=130)
    return lines


def draw_map(ax, r, labels):
    for k in range(r.ngroup):
        s = r.group == k
        o = np.argsort(r.x[s])
        ax.plot(r.x[s][o] / 1e3, r.y[s][o] / 1e3, "k-", lw=2)
        ax.text(r.x[s].mean() / 1e3, r.y[s].mean() / 1e3 + 8, labels[k], ha="center", fontsize=9)
    ax.set_aspect("equal")
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")


def fig7_8(rd):
    fig = plt.figure(figsize=(11, 10))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.55, 1, 1, 1], hspace=0.45)
    labels = ("2", "1", "3")
    order = (1, 0, 2)           # paper faults 1, 2, 3 = our 1, 0, 2
    for j, name in enumerate(("modelA", "modelB")):
        r = Run(os.path.join(rd, name))
        axm = fig.add_subplot(gs[0, j])
        draw_map(axm, r, labels)
        axm.set_title("%s (%s)" % (name, "leftstep" if j == 0 else "rightstep"))
        ymax = 0.0
        axs = []
        for i, k in enumerate(order):
            a = fig.add_subplot(gs[i + 1, j])
            t, m0 = onsets(r, k)
            a.vlines(t, 0, m0 / 1e20, color="k", lw=1.0)
            a.text(0.02, 0.85, "fault %s" % labels[k], transform=a.transAxes)
            a.set_ylabel("M0 [1e20 N m]")
            ymax = max(ymax, m0.max() / 1e20 if len(m0) else 0.0)
            a.set_xlim(0, r.t[-1])
            axs.append(a)
        for a in axs:
            a.set_ylim(0, 1.1 * ymax)
        axs[-1].set_xlabel("time [yr]")
    fig.suptitle("rate-and-state analogue of Becker and Schmeling (1998) Figs 7 and 8", fontsize=10)
    fig.savefig(os.path.join(rd, "bs98_fig7_8.png"), dpi=130, bbox_inches="tight")


def fig9_10(rd):
    fig = plt.figure(figsize=(11, 7.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.28, 1], hspace=0.35)
    labels = ("5", "2", "1", "3", "4")
    for j, name in enumerate(("modelC", "modelD")):
        r = Run(os.path.join(rd, name))
        axm = fig.add_subplot(gs[0, j])
        draw_map(axm, r, labels)
        axm.set_title("%s (%s)" % (name, "extended leftstep" if j == 0 else "extended rightstep"))
        a = fig.add_subplot(gs[1, j])
        hot = r.logv > -2.0
        ii, cc = np.nonzero(hot)
        a.plot(r.x[cc] / 1e3, r.t[ii], "k.", ms=1.0)
        a.set_xlabel("x [km]")
        a.set_ylabel("t [yr]")
        a.set_ylim(0, r.t[-1])
        for k in range(r.ngroup):
            a.text(r.x[r.group == k].mean() / 1e3, 0.01 * r.t[-1], labels[k], ha="center", fontsize=9)
    fig.suptitle("rate-and-state analogue of Becker and Schmeling (1998) Figs 9 and 10\n"
                 "dots: cells with |v| > 1 cm/s at series records", fontsize=10)
    fig.savefig(os.path.join(rd, "bs98_fig9_10.png"), dpi=130, bbox_inches="tight")


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
    fig4_5(rd)
    for line in fig6(rd):
        print(line)
    fig7_8(rd)
    fig9_10(rd)
    print("wrote bs98_fig4_5.png bs98_fig6.png bs98_fig7_8.png bs98_fig9_10.png in", rd)


if __name__ == "__main__":
    main()
