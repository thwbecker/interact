#!/usr/bin/env python3
"""
plot_bp3ve_demo.py: large-event recurrence vs Maxwell time from a
run_bp3ve_demo directory (dipping thrust in an elastic plate over a
Maxwell half-space, gravity on, sigma_n fixed).  Contained-fault
backslip loading: recurrence lengthens with relaxation strength.
usage: plot_bp3ve_demo.py rundir [slip_min_m] [nskip] [gap_yr]
       nskip: intervals to drop from the START of each sequence
       (default 1).  The first cycles of these runs are a start-up
       transient and, at 1500 yr, the sequence can still be drifting
       toward its limit; the last interval is printed alongside the
       mean so that a drift is visible rather than averaged away.
       gap_yr: catalog rows closer together than this are ONE rupture
       episode (default 1 yr).  A rupture that arrests and immediately
       re-nucleates is logged as two rows a fraction of a day apart,
       and counting rows rather than episodes halves the recurrence.
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rundir = sys.argv[1]
smin = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
nskip = int(sys.argv[3]) if len(sys.argv) > 3 else 1
gap = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
nrow = {}
os.chdir(rundir)

def big_events(f):
    """onset times of the large rupture EPISODES: rows within gap
    years of one another belong to one episode (see gap_yr above)"""
    tt = []
    for l in open(f):
        if l.startswith("#"):
            continue
        c = l.split()
        if float(c[7]) > smin:
            tt.append(float(c[1]))
    tt = sorted(tt)
    out = []
    for t in tt:
        if not out or t - out[-1] > gap:
            out.append(t)
    nrow[f] = (len(tt), len(out))
    return out

def trec(f):
    """mean interval after dropping the first nskip intervals, with
    the last interval alongside it as a drift indicator"""
    t = big_events(f)
    if len(t) < nskip + 3:
        return (None, len(t), None)
    d = np.diff(t)[nskip:]
    return (np.mean(d), len(t), d[-1])

te, ne, tel = trec("cat_el.dat")
if te is None:
    raise SystemExit(f"plot_bp3ve_demo: only {ne} large events in "
                     "cat_el.dat: run longer (-stop_time_yr) or lower "
                     "the slip threshold (2nd argument)")
def rowinfo(f):
    r, e = nrow[f]
    return f" [{r} rows -> {e} episodes]" if r != e else ""
print(f"elastic: t_rec {te:.1f} yr, last interval {tel:.1f} yr "
      f"({ne} large events, first {nskip} interval(s) dropped)"
      + rowinfo("cat_el.dat"))
tms, trs = [], []
for f in sorted(glob.glob("cat_tm*.dat")):
    tm = float(f.replace("cat_tm", "").replace(".dat", ""))
    r, n, rl = trec(f)
    if r:
        tms.append(tm); trs.append(r)
        print(f"tM = {tm:6.1f} yr: t_rec {r:6.1f} yr, last {rl:6.1f} yr "
              f"({n} events, {100*(r/te - 1):+5.1f} percent of elastic)"
              + rowinfo(f))
    else:
        print(f"tM = {tm:6.1f} yr: only {n} large events, skipped")
tms, trs = np.array(tms), np.array(trs)
if tms.size == 0:
    raise SystemExit("plot_bp3ve_demo: no viscoelastic run had enough "
                     "large events to measure a recurrence")

fig, ax = plt.subplots(figsize=(7, 5))
ax.axhline(1.0, color="0.6", lw=0.8, ls=":", label="elastic")
o = np.argsort(tms)
ax.plot(te/tms[o], trs[o]/te, "ks-", ms=7, mfc="none",
        label="rsf_solve -ve_mode 3 (plate over Maxwell, gravity)")
ax.set_xlabel("T_rec(elastic) / t_M   (relaxation strength)")
ax.set_ylabel("T_rec / T_rec(elastic)")
import glob as _g
_lg = sorted(_g.glob("run_tm*.log"))
_sig = "shear + normal tractions"
if _lg and not any("SHEAR + NORMAL" in open(l, errors="ignore").read()
                   for l in _lg):
    _sig = "shear only, sigma_n frozen"
ax.set_title("Contained dipping thrust: recurrence change by\n"
             f"substrate relaxation ({_sig})", fontsize=11)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("bp3ve_demo.png", dpi=150, bbox_inches="tight")
print("wrote bp3ve_demo.png")
