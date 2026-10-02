#!/usr/bin/env python3
"""
Compare BP1 pseudo-transient runs: event catalogues, interevent times, iteration
statistics, and a summary figure.

usage: python3 compare_runs.py run_a.csv run_b.csv ...
Labels are the file stems. Logs written by bp1_pt.py; a log appended to on
restart may lack a header, which is handled.

Event definition: a contiguous set of logged samples with Vmax > v_seis
(default 1e-3 m/s); samples separated by less than gap_s seconds (default 1 day)
are merged into one event. Onset is the first sample above v_seis, duration the
time until the last. With dense logging (bp1_pt.py --log_dense_Vmax) these are
resolved to the time step; with sparse logging they are only indicative, and a
note is printed.
"""
import sys
import os
import numpy as np

COLS = ["step", "t_yr", "dt_s", "it", "resid_pa", "Vmax", "V_0km", "V_7p5km",
        "V_15km", "tau_7p5km_MPa", "slip_7p5km"]
V_SEIS = 1e-3
YR = 365.25 * 86400.0
FULL_SLIP = 0.5          # m at 7.5 km to call an event "full"
GAP_S = 86400.0          # merge fast-slip samples closer than this


def load(path):
    with open(path) as f:
        first = f.readline()
    if first.startswith("step"):
        d = np.genfromtxt(path, delimiter=",", names=True)
    else:
        d = np.genfromtxt(path, delimiter=",", names=COLS)
    return np.atleast_1d(d)


def events(d, v_seis=V_SEIS, gap_s=GAP_S):
    t, V = d["t_yr"], d["Vmax"]
    fast = np.where(V > v_seis)[0]
    if len(fast) == 0:
        return []
    groups = [[fast[0]]]
    for i in fast[1:]:
        if (t[i] - t[groups[-1][-1]]) * YR > gap_s:
            groups.append([i])
        else:
            groups[-1].append(i)
    out = []
    for g in groups:
        s, e = g[0], g[-1]
        pre = max(s - 1, 0)       # last sample before the fast phase
        slip = d["slip_7p5km"][e] - d["slip_7p5km"][pre]
        # dense sampling if the median step between samples within the event is one step
        dense = len(g) > 2 and np.median(np.diff(d["step"][g])) <= 1.5
        out.append(dict(t_on=t[s], dur_s=(t[e] - t[s]) * YR, peakV=V[s:e + 1].max(),
                        slip=slip, tau_before=d["tau_7p5km_MPa"][pre],
                        tau_after=d["tau_7p5km_MPa"][e], full=slip > FULL_SLIP,
                        nsamp=len(g), dense=dense))
    return out


def summarize(label, d):
    ev = events(d)
    print(f"\n{label}: {len(d)} samples, {d['t_yr'][-1]:.1f} yr, {int(d['step'][-1])} steps")
    inter = d["Vmax"] < 1e-8
    if np.any(inter):
        print(f"  interseismic PT iterations per step: mean {d['it'][inter].mean():.0f}")
    cos = d["Vmax"] > V_SEIS
    if np.any(cos):
        print(f"  coseismic PT iterations per step: mean {d['it'][cos].mean():.0f}, max {d['it'][cos].max():.0f}")
    if ev and not all(e["dense"] for e in ev):
        print("  note: some events are sparsely sampled; durations, peaks and 'before' stresses are indicative only")
    print("  events (onset yr, duration s, peak V m/s, slip at 7.5 km m, stress at 7.5 km MPa, samples):")
    for e in ev:
        kind = "full   " if e["full"] else "partial"
        flag = "" if e["dense"] else " (sparse)"
        print(f"   {kind} {e['t_on']:9.3f}  {e['dur_s']:8.1f}  {e['peakV']:5.2f}  "
              f"{e['slip']:5.2f}  {e['tau_before']:6.2f} -> {e['tau_after']:6.2f}  {e['nsamp']:5d}{flag}")
    fulls = [e["t_on"] for e in ev if e["full"]]
    if len(fulls) > 1:
        iv = np.diff(fulls)
        print("  full-to-full intervals (yr):", np.round(iv, 2))
        if len(iv) > 1:
            print(f"  mean interval excluding the first: {iv[1:].mean():.1f} yr "
                  f"(first excluded: unstable initial state)")
    return ev


def main(paths):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = [(os.path.splitext(os.path.basename(p))[0], load(p)) for p in paths if os.path.exists(p)]
    for p in paths:
        if not os.path.exists(p):
            print("missing:", p)
    if not runs:
        return
    for label, d in runs:
        summarize(label, d)

    fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for label, d in runs:
        ax[0].semilogy(d["t_yr"], d["Vmax"], lw=0.8, label=label)
        ax[1].plot(d["t_yr"], d["slip_7p5km"], lw=0.8, label=label)
        ax[2].plot(d["t_yr"], d["tau_7p5km_MPa"], lw=0.8, label=label)
    ax[0].axhline(V_SEIS, color="k", ls=":", lw=0.6)
    ax[0].set_ylabel("max slip rate [m/s]")
    ax[1].set_ylabel("slip at 7.5 km [m]")
    ax[2].set_ylabel("shear stress at 7.5 km [MPa]")
    ax[2].set_xlabel("time [yr]")
    ax[0].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("compare_bp1.png", dpi=130)
    print("\nfigure written to compare_bp1.png")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
