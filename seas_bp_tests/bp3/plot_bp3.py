#!/usr/bin/env python3
"""
plot_bp3.py: wrap-up figure for a SEAS BP3-QD run directory
(run_bp3 output): interevent times, station slip-rate histories,
cumulative station slip, and fault stresses; prints the
characteristic-interval summary for comparison with Erickson et al.
(BSSA 2023) Figure 8 (dip 90: one event every ~90 yr; dip 60 thrust:
~60/87/90/95 yr; dip 30 thrust: ~65/80 yr).

usage: plot_bp3.py [rundir] [tmin_yr]
       rundir   default .
       tmin_yr  exclude events before this time from the interval
                summary (spin-up; default 150)
writes <rundir>/bp3_summary.png
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rundir = sys.argv[1] if len(sys.argv) > 1 else "."
tmin = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0
yr = 3.15576e7
os.chdir(rundir)

# ---- catalog: dedupe restart-boundary rows, sort by onset ------------------
cat = []
with open("rsf_catalog.dat") as f:
    for l in f:
        if l.startswith("#"):
            continue
        c = l.split()
        cat.append((float(c[1]), float(c[7])))    # onset [yr], mean slip [m]
cat = sorted(set(cat))
tev = np.array([c[0] for c in cat])
sev = np.array([c[1] for c in cat])
iv = np.diff(tev)
msel = tev[1:] >= tmin
print(f"plot_bp3: {tev.size} events; intervals after t = {tmin:.0f} yr:")
print("  " + " ".join(f"{d:7.2f}" for d in iv[msel]))
uq = sorted(set(np.round(iv[msel], 1)))
print("  distinct (0.1-yr rounding): " +
      " ".join(f"{u:.1f}" for u in uq) + "  [yr]")
print("  mean slip per event: " +
      " ".join(f"{s:5.2f}" for s in sev[1:][msel]) + "  [m]")

fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.5))
ax = axs[0, 0]
ax.plot(np.arange(2, tev.size + 1), iv, "ks-", ms=5, mfc="none")
ax.axvline(1 + np.searchsorted(tev[1:], tmin) + 0.5, color="0.8", lw=0.8)
ax.set_xlabel("event number")
ax.set_ylabel("time between events [yr]")
ax.set_title("interevent times (cf. BSSA 2023 Fig. 8d-f)", fontsize=10)

# ---- station files ----------------------------------------------------------
stf = sorted(glob.glob("fltst_dp*.dat"))
show = [s for s in stf if any(k in s for k in
                              ("dp000", "dp075", "dp125", "dp250"))]
ax = axs[0, 1]
axs2 = axs[1, 0]
ax3 = axs[1, 1]
if stf:
    for s in show:
        d = np.loadtxt(s, comments="#")
        lab = s.replace("fltst_", "").replace(".dat", "")
        ax.plot(d[:, 0]/yr, d[:, 2], lw=0.7, label=lab)
        axs2.plot(d[:, 0]/yr, d[:, 1], lw=1.0, label=lab)
        ax3.plot(d[:, 0]/yr, d[:, 3], lw=0.8, label=lab + " tau")
    # normal stress: only where it moves (sigma-coupled dips)
    d0 = np.loadtxt(show[0], comments="#")
    if np.ptp(d0[:, 4]) > 1e-3:
        for s in show:
            d = np.loadtxt(s, comments="#")
            lab = s.replace("fltst_", "").replace(".dat", "")
            ax3.plot(d[:, 0]/yr, d[:, 4], "--", lw=0.8,
                     label=lab + " sigma")
    ax.set_ylabel("log10 slip rate [m/s]")
    ax.axhline(-9, color="0.8", lw=0.5)
    axs2.set_ylabel("slip [m]")
    ax3.set_ylabel("stress [MPa]")
    for a in (ax, axs2, ax3):
        a.set_xlabel("t [yr]")
        a.legend(fontsize=7, ncol=2)
    ax.set_title("station slip rates", fontsize=10)
    axs2.set_title("station cumulative slip", fontsize=10)
    ax3.set_title("station shear (and normal, dashed) stress",
                  fontsize=10)
else:
    for a in (ax, axs2, ax3):
        a.text(0.5, 0.5, "no fltst_dp*.dat station files",
               ha="center", transform=a.transAxes)
fig.suptitle(f"BP3-QD summary: {os.path.basename(os.getcwd())}",
             fontsize=12)
fig.tight_layout()
fig.savefig("bp3_summary.png", dpi=150, bbox_inches="tight")
print("plot_bp3: wrote bp3_summary.png")
