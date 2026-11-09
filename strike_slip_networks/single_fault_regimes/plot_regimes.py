#!/usr/bin/env python3
"""
plot_regimes.py [rundir]: analysis of the runs of run_regimes (single
strip, backslip loading, W/h* = 1 ... 64, laws one, slip, two)

Events are taken from rsf_catalog.dat; an event starting within one day
of the end of the previous one is merged into it (this removes the
short re-ruptures of one or two cells at the strip tips seconds to
minutes after large slip-law events); the first 20 percent of the run
is discarded. Rupture length = rupture area / strip depth; a full
rupture has a length of at least 0.9 W (W = 48 km, the velocity-
weakening length).

Writes
  regimes_summary.png   events per full-rupture cycle, CV of the
                        intervals (all events and full ruptures only) and
                        fraction of full ruptures against W/h*, with the
                        estimate of Cattania (2019), alpha = 0.45
                        (W/L_inf)^1/2 (L_inf = h*/2), and the regime
                        boundaries of Wu & Chen (2014) at W/h* = 1.2, 5.5,
                        33, 105 (both for 1-D antiplane faults)
  regimes_bifurcation.png
                        bifurcation-style diagrams: every normalized
                        interevent time lambda*dt (lambda the mean event
                        rate) and every rupture length / W against W/h*,
                        one column per law
  regimes_intervals.png histograms of lambda*dt, one panel per law and
                        W/h*, common axes
  regimes_lengths.png   survival functions of rupture length / W per law
  regimes_spacetime.png log10 slip rate along strike against time for
                        selected runs
  regimes_summary.txt   the numbers
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from rsf_series import Run, cat_columns, _load_cat  # noqa: E402

LAWS = (("one", "one state, aging", "C0", "o"), ("slip", "one state, slip", "C3", "s"),
        ("two", "two states, slip", "C2", "^"))
XS = (1, 2, 4, 8, 16, 32, 64)
WVW = 48e3
DEPTH = 15e3
CC = cat_columns(None)
YR = 365.25 * 86400.0


def events(d):
    """merged events after spin-up: arrays onset [yr], length [m], Mw; and tmax"""
    c = _load_cat(os.path.join(d, "rsf_catalog.dat"))
    with open(os.path.join(d, "run.cmd")) as f:
        tmax = float(f.read().split("-stop_time_yr")[1].split()[0])
    if c is None:
        return np.zeros(0), np.zeros(0), np.zeros(0), tmax
    on, ar = c[:, CC["onset"]], c[:, CC["arrest"]]
    ln = c[:, CC["area"]] / DEPTH
    mw = c[:, CC["Mw"]]
    mon, mln, mmw, last_ar = [], [], [], -1e30
    for i in range(len(on)):
        if mon and (on[i] - last_ar) * YR < 86400.0:
            mln[-1] = max(mln[-1], ln[i])
            mmw[-1] = max(mmw[-1], mw[i])
            last_ar = max(last_ar, ar[i])
            continue
        mon.append(on[i])
        mln.append(ln[i])
        mmw.append(mw[i])
        last_ar = ar[i]
    mon, mln, mmw = np.array(mon), np.array(mln), np.array(mmw)
    keep = mon >= 0.2 * tmax
    return mon[keep], mln[keep], mmw[keep], tmax


def stats(on, ln):
    full = ln >= 0.9 * WVW
    out = dict(n=len(on), nfull=int(full.sum()), alpha=np.nan, cv=np.nan, cvfull=np.nan,
               ffull=full.mean() if len(on) else np.nan)
    if len(on) > 2:
        dt = np.diff(on)
        out["cv"] = dt.std(ddof=1) / dt.mean()
    if full.sum() > 2:
        dtf = np.diff(on[full])
        out["cvfull"] = dtf.std(ddof=1) / dtf.mean()
        out["alpha"] = len(on) / full.sum()
    return out


def main():
    rd = sys.argv[1] if len(sys.argv) > 1 else "."
    res = {}
    for law, _, _, _ in LAWS:
        for x in XS:
            d = os.path.join(rd, "%s_W%i" % (law, x))
            if os.path.exists(os.path.join(d, "rsf_catalog.dat")):
                on, ln, mw, tmax = events(d)
                res[(law, x)] = (on, ln, mw, stats(on, ln))
    lines = ["%-5s %4s %6s %6s %7s %7s %7s %7s" % ("law", "W/h*", "events", "full", "alpha", "CV", "CVfull", "f_full")]
    for (law, x), (on, ln, mw, s) in sorted(res.items(), key=lambda k: ([l[0] for l in LAWS].index(k[0][0]), k[0][1])):
        lines.append("%-5s %4i %6i %6i %7.2f %7.3f %7.3f %7.2f" % (law, x, s["n"], s["nfull"], s["alpha"], s["cv"],
                                                                s["cvfull"], s["ffull"]))
    with open(os.path.join(rd, "regimes_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # summary against W/h*
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    xx = np.logspace(-0.2, 2.3, 100)
    ax[0].plot(xx, np.maximum(1.0, 0.45 * np.sqrt(2 * xx)), "k-", lw=1, label="Cattania (2019), 0.45 (W/L_inf)^1/2")
    for a in ax:
        for xb in (1.2, 5.5, 33, 105):
            a.axvline(xb, color="0.7", ls=":", lw=1)
        a.set_xscale("log")
        a.set_xlabel("W/h*  (h* = 2 L_inf)")
    for law, lab, col, mk in LAWS:
        xs = [x for x in XS if (law, x) in res]
        for k, key in enumerate(("alpha", "cv", "ffull")):
            ax[k].plot(xs, [res[(law, x)][3][key] for x in xs], marker=mk, color=col, label=lab, lw=1)
        ax[1].plot(xs, [res[(law, x)][3]["cvfull"] for x in xs], marker=mk, color=col, ls="--", lw=1, mfc="none")
    for law, lab, col, mk in LAWS:
        st = [x for x in XS if (law, x) in res and res[(law, x)][3]["n"] == 0]
        ax[0].plot(st, [0.9] * len(st), "x", color=col, ms=8)
    ax[0].text(0.75, 0.82, "x: no events (stable)", fontsize=7)
    ax[0].set_yscale("log")
    ax[0].set_ylabel("events per full-rupture cycle")
    ax[1].set_ylabel("CV of intervals (solid: all events, dashed: full ruptures)")
    ax[2].set_ylabel("fraction of full ruptures")
    ax[0].legend(fontsize=8)
    ax[2].text(0.02, 0.02, "dotted: Wu & Chen (2014) boundaries\nW/h* = 1.2, 5.5, 33, 105", fontsize=7,
               color="0.4", transform=ax[2].transAxes)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "regimes_summary.png"), dpi=130)

    # bifurcation-style diagrams
    fig, ax = plt.subplots(2, len(LAWS), figsize=(5 * len(LAWS), 7), sharex=True, sharey="row")
    rng = np.random.default_rng(1)
    for i, (law, lab, col, mk) in enumerate(LAWS):
        for x in XS:
            if (law, x) not in res or len(res[(law, x)][0]) < 3:
                continue
            on, ln = res[(law, x)][0], res[(law, x)][1]
            dt = np.diff(on)
            jit = x * np.exp(rng.uniform(-0.08, 0.08, len(dt)))
            ax[0, i].plot(jit, dt / dt.mean(), ".", color=col, ms=3, alpha=0.6)
            jit = x * np.exp(rng.uniform(-0.08, 0.08, len(ln)))
            ax[1, i].plot(jit, ln / WVW, ".", color=col, ms=3, alpha=0.6)
        ax[0, i].set_title(lab, fontsize=10)
        for k in range(2):
            ax[k, i].set_xscale("log")
            ax[k, i].set_yscale("log")
            for xb in (1.2, 5.5, 33, 105):
                ax[k, i].axvline(xb, color="0.8", ls=":", lw=1)
        ax[1, i].set_xlabel("W/h*")
    ax[0, 0].set_ylabel("interevent time / mean")
    ax[1, 0].set_ylabel("rupture length / W")
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "regimes_bifurcation.png"), dpi=130)

    # interevent-time histograms
    xs_all = [x for x in XS if any((l[0], x) in res for l in LAWS)]
    fig, ax = plt.subplots(len(LAWS), len(xs_all), figsize=(2.3 * len(xs_all), 2.2 * len(LAWS)),
                           sharex=True, squeeze=False)
    bins = np.logspace(-6, 1, 36)
    for i, (law, lab, col, mk) in enumerate(LAWS):
        for j, x in enumerate(xs_all):
            a = ax[i, j]
            if (law, x) in res and len(res[(law, x)][0]) > 2:
                on = res[(law, x)][0]
                dt = np.diff(on)
                lam = 1.0 / dt.mean()
                h, e = np.histogram(lam * dt, bins=bins, density=True)
                a.bar(e[:-1], h, width=np.diff(e), align="edge", color=col, alpha=0.7)
                a.text(0.03, 0.85, "%i ev" % len(on), transform=a.transAxes, fontsize=7)
            a.set_xscale("log")
            a.set_yscale("log")
            a.set_ylim(1e-2, 1e3)
            a.axvline(1.0, color="k", lw=0.5, ls=":")
            if i == 0:
                a.set_title("W/h* = %i" % x, fontsize=9)
            if j == 0:
                a.set_ylabel("%s\npdf" % lab, fontsize=8)
            if i == len(LAWS) - 1:
                a.set_xlabel("lambda dt", fontsize=8)
            a.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "regimes_intervals.png"), dpi=120)

    # rupture length survival functions
    fig, ax = plt.subplots(1, len(LAWS), figsize=(5 * len(LAWS), 4), sharey=True)
    cm = plt.cm.viridis(np.linspace(0, 1, len(XS)))
    for i, (law, lab, col, mk) in enumerate(LAWS):
        for j, x in enumerate(XS):
            if (law, x) in res and len(res[(law, x)][1]) > 1:
                ln = np.sort(res[(law, x)][1] / WVW)[::-1]
                ax[i].step(ln, np.arange(1, len(ln) + 1) / len(ln), color=cm[j], label="W/h* %i" % x, where="post")
        r = np.logspace(-2, 0, 10)
        ax[i].plot(r, 0.1 * r ** -0.5, "k:", lw=1)
        ax[i].text(0.012, 0.9, "slope -1/2", fontsize=7)
        ax[i].set_xscale("log")
        ax[i].set_yscale("log")
        ax[i].set_xlabel("rupture length / W")
        ax[i].set_title(lab, fontsize=9)
    ax[0].set_ylabel("N / N_tot (survival)")
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(rd, "regimes_lengths.png"), dpi=120)

    # space-time plots for selected runs
    sel = [x for x in (4, 16, 64) if any((l[0], x) in res for l in LAWS)]
    fig, ax = plt.subplots(len(LAWS), len(sel), figsize=(5 * len(sel), 2.6 * len(LAWS)), squeeze=False)
    for i, (law, lab, col, mk) in enumerate(LAWS):
        for j, x in enumerate(sel):
            a = ax[i, j]
            d = os.path.join(rd, "%s_W%i" % (law, x))
            if (law, x) not in res:
                a.axis("off")
                continue
            r = Run(d)
            w = r.t >= r.t[-1] - 1500.0
            xs_ = (r.x - r.x.min()) / 1e3
            o = np.argsort(xs_)
            im = a.pcolormesh(np.arange(w.sum()), xs_[o], r.logv[w][:, o].T, cmap="magma", vmin=-12, vmax=0,
                              shading="nearest", rasterized=True)
            a.set_title("%s, W/h* %i (last 1500 yr, by record)" % (lab, x), fontsize=8)
            a.set_ylabel("along strike [km]", fontsize=8)
            a.tick_params(labelsize=7)
    fig.colorbar(im, ax=ax, fraction=0.02, label="log10 V [m/s]")
    fig.savefig(os.path.join(rd, "regimes_spacetime.png"), dpi=120, bbox_inches="tight")


if __name__ == "__main__":
    main()
