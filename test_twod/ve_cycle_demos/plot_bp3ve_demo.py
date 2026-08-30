#!/usr/bin/env python3
"""
plot_bp3ve_demo.py: large-event recurrence vs Maxwell time from a
run_bp3ve_demo directory (dipping thrust in an elastic plate over a
Maxwell half-space, gravity on, sigma_n fixed).  Contained-fault
backslip loading: recurrence lengthens with relaxation strength.
usage: plot_bp3ve_demo.py rundir [slip_min_m]
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rundir = sys.argv[1]
smin = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
os.chdir(rundir)

def big_events(f):
    tt = []
    for l in open(f):
        if l.startswith("#"):
            continue
        c = l.split()
        if float(c[7]) > smin:
            tt.append(float(c[1]))
    return sorted(set(tt))

def trec(f):
    t = big_events(f)
    return (np.mean(np.diff(t[1:])), len(t)) if len(t) > 3 else (None, len(t))

te, ne = trec("cat_el.dat")
print(f"elastic: t_rec {te:.1f} yr ({ne} large events)")
tms, trs = [], []
for f in sorted(glob.glob("cat_tm*.dat")):
    tm = float(f.replace("cat_tm", "").replace(".dat", ""))
    r, n = trec(f)
    if r:
        tms.append(tm); trs.append(r)
        print(f"tM = {tm:6.1f} yr: t_rec {r:.1f} yr ({n} events)")
tms, trs = np.array(tms), np.array(trs)

fig, ax = plt.subplots(figsize=(7, 5))
ax.axhline(1.0, color="0.6", lw=0.8, ls=":", label="elastic")
o = np.argsort(tms)
ax.plot(te/tms[o], trs[o]/te, "ks-", ms=7, mfc="none",
        label="rsf_solve -ve_mode 3 (plate over Maxwell, gravity)")
ax.set_xlabel("T_rec(elastic) / t_M   (relaxation strength)")
ax.set_ylabel("T_rec / T_rec(elastic)")
ax.set_title("Contained dipping thrust: recurrence lengthening by\n"
             "substrate relaxation (loading-pathway mechanism)",
             fontsize=11)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("bp3ve_demo.png", dpi=150, bbox_inches="tight")
print("wrote bp3ve_demo.png")
