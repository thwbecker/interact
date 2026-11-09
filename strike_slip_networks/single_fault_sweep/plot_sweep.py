#!/usr/bin/env python3
"""
plot_sweep.py [rundir]: figures and a table for the runs of
run_single_sweep (one state variable, aging law, against two state
variables, slip law; L/h* 10, 20, 30; strike -20 ... 20 degrees relative
to the optimal orientation)

  sweep_timeseries.png  mean tau/sigma of the velocity-weakening cells
                        against time, one panel per run
  sweep_recurrence.png  interval to the next event against event time,
                        symbol size by Mw, colour by orientation
  sweep_summary.png     per run: number of events with Mw >= 6, CV of
                        their intervals after removing a linear trend in
                        time, order-6 permutation entropy of the peak
                        tau/sigma sequence, fraction of those events that
                        are partial ruptures (area below half of the
                        velocity-weakening area)

Recurrence statistics use events with Mw >= 6 only, as in Sun & Becker
(2026): the catalog also lists small re-ruptures of one or two cells
seconds after a large event (Mw about 4.5), which would otherwise
dominate the interval statistics.
  sweep_summary.txt     the same numbers as a table
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from complexity import extrema, perm_entropy  # noqa: E402
from rsf_series import Run, cat_columns  # noqa: E402

LAWS = (("one", "one state, aging"), ("two", "two states, slip"))
XS = (10, 20, 30)
STRIKES = (-20, -10, 0, 10, 20)
CC = cat_columns(None)


def load(rd, law, x, s):
    d = os.path.join(rd, "%s_X%i_s%i" % (law, x, s))
    if not os.path.exists(os.path.join(d, "series.bin")):
        return None
    try:
        return Run(d)
    except (IOError, ValueError):
        return None


def measures(r):
    vw = r.a < r.b + (r.par.get("b2", 0.0) or 0.0)
    mu = (r.tau[:, vw] / r.sigma[:, vw]).mean(axis=1)
    out = dict(mu=mu)
    c = r.cat
    if c is None or len(c) < 3:
        out.update(n=0 if c is None else len(c), nsmall=0, cv=np.nan, pe=np.nan, partial=np.nan, on=np.zeros(0),
                   dt=np.zeros(0), mw=np.zeros(0))
        return out
    big = c[:, CC["Mw"]] >= 6.0
    nsmall = int(np.sum(~big))
    c = c[big]
    if len(c) < 3:
        out.update(n=len(c), nsmall=nsmall, cv=np.nan, pe=np.nan, partial=np.nan, on=np.zeros(0),
                   dt=np.zeros(0), mw=np.zeros(0))
        return out
    on = c[:, CC["onset"]]
    mw = c[:, CC["Mw"]]
    area = c[:, CC["area"]]
    a_vw = np.sum(4.0 * r.par["dx"] / 2.0 * r.par["W"] / 2.0 * np.ones(np.sum(vw)))
    dt = np.diff(on)
    tt = on[:-1]
    if len(dt) >= 3:
        p = np.polyfit(tt, dt, 1)
        res = dt - np.polyval(p, tt)
        cv = res.std(ddof=1) / dt.mean()
    else:
        cv = np.nan
    te, ve, ke = extrema(r.t, mu)
    pe = perm_entropy(ve[ke > 0], 6, 3e-5)
    out.update(n=len(on), nsmall=nsmall, cv=cv, pe=pe, partial=np.mean(area < 0.5 * a_vw), on=on, dt=dt, mw=mw)
    return out


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
    res = {}
    for law, _ in LAWS:
        for x in XS:
            for s in STRIKES:
                r = load(rd, law, x, s)
                if r is not None:
                    res[(law, x, s)] = (r, measures(r))
    cols = [(law, x) for law, _ in LAWS for x in XS]
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(STRIKES)))
    # time series grid
    fig, ax = plt.subplots(len(STRIKES), len(cols), figsize=(3.0 * len(cols), 1.5 * len(STRIKES)),
                           sharex=True, squeeze=False)
    for j, (law, x) in enumerate(cols):
        for i, s in enumerate(STRIKES):
            a = ax[i, j]
            if (law, x, s) in res:
                r, m = res[(law, x, s)]
                a.plot(r.t, m["mu"], "k-", lw=0.5)
                a.text(0.02, 0.05, "%i ev >= M6, PE %.1f" % (m["n"], m["pe"]) if m["n"] >= 3 else "%i ev" % m["n"],
                       transform=a.transAxes, fontsize=7)
            if i == 0:
                a.set_title("%s, L/h* %i" % (dict(LAWS)[law], x), fontsize=8)
            if j == 0:
                a.set_ylabel("strike %+i\ntau/sigma" % s, fontsize=8)
            a.tick_params(labelsize=7)
    for a in ax[-1]:
        a.set_xlabel("time [yr]", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "sweep_timeseries.png"), dpi=120)
    # recurrence against time
    fig, ax = plt.subplots(2, 3, figsize=(13, 7), squeeze=False)
    for j, (law, x) in enumerate(cols):
        a = ax[j // 3, j % 3]
        for i, s in enumerate(STRIKES):
            if (law, x, s) in res:
                m = res[(law, x, s)][1]
                if len(m["dt"]):
                    a.scatter(m["on"][:-1], m["dt"], s=4 + 8 * (m["mw"][:-1] - 5.5).clip(0), color=colors[i],
                              label="%+i" % s, edgecolor="k", lw=0.2)
        a.set_yscale("log")
        a.set_title("%s, L/h* %i" % (dict(LAWS)[law], x), fontsize=9)
        a.set_xlabel("event time [yr]", fontsize=8)
        a.set_ylabel("interval to next event [yr]", fontsize=8)
    ax[0, 0].legend(title="strike rel. to optimal", fontsize=7, title_fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "sweep_recurrence.png"), dpi=120)
    # summary heat maps
    names = (("n", "events with Mw >= 6"), ("cv", "CV of detrended intervals"), ("pe", "PE of peak tau/sigma [bits]"),
             ("partial", "fraction of partial ruptures"))
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.8))
    lines = ["%-4s %3s %4s %5s %6s %6s %5s %7s" % ("law", "L/h*", "strk", "n>=6", "n<6", "CV_dt", "PE", "partial")]
    for k, (key, title) in enumerate(names):
        z = np.full((len(cols), len(STRIKES)), np.nan)
        for j, (law, x) in enumerate(cols):
            for i, s in enumerate(STRIKES):
                if (law, x, s) in res:
                    z[j, i] = res[(law, x, s)][1][key]
        im = ax[k].imshow(z, cmap="viridis", aspect="auto")
        for j in range(len(cols)):
            for i in range(len(STRIKES)):
                if np.isfinite(z[j, i]):
                    ax[k].text(i, j, "%.2g" % z[j, i], ha="center", va="center", fontsize=7, color="w")
        ax[k].set_xticks(range(len(STRIKES)))
        ax[k].set_xticklabels(["%+i" % s for s in STRIKES])
        ax[k].set_yticks(range(len(cols)))
        ax[k].set_yticklabels(["%s L/h* %i" % (law, x) for law, x in cols], fontsize=8)
        ax[k].set_xlabel("strike rel. to optimal [deg] (+: clamping)", fontsize=8)
        ax[k].set_title(title, fontsize=9)
        plt.colorbar(im, ax=ax[k], fraction=0.05)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "sweep_summary.png"), dpi=120)
    for law, x in cols:
        for s in STRIKES:
            if (law, x, s) in res:
                m = res[(law, x, s)][1]
                lines.append("%-4s %3i %+4i %5i %6i %6.3f %5.2f %7.2f" % (law, x, s, m["n"], m["nsmall"], m["cv"],
                                                                        m["pe"], m["partial"]))
    with open(os.path.join(rd, "sweep_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
