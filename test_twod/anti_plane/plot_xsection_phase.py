#!/usr/bin/env python3
"""
plot_xsection_phase.py: stress at fixed off-fault probe points versus
cycle phase, one line per model, from the probe time series that
plot_xsection.py stores in each run's xsect_fields.npz (re-run
plot_xsection.py per directory if the npz predates the probes).

Each panel is one probe (plate points at mid-fault depth, substrate
points below the plate base; x = 0.25 and 0.75 of the plate
thickness).  Curves show sigma minus its own cycle mean, so
amplitude AND phase lag of the loading at that point read off
directly: elastic is the zero-lag reference, small tM/T_rec kills
the far-field amplitude, tM/T_rec ~ 1 shifts the substrate response
by about a Maxwell time.

usage: plot_xsection_phase.py [dir1 dir2 ...]
       default: xsect_r0.1 xsect_r0.5 xsect_r1 xsect_r2 xsect_r10 xsect_el
writes xsect_phase_sxz.png and xsect_phase_syz.png
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dirs = sys.argv[1:] if len(sys.argv) > 1 else \
    ["xsect_r0.1", "xsect_r0.5", "xsect_r1", "xsect_r2", "xsect_r10",
     "xsect_el"]
dirs = sorted(dirs, key=lambda d: d == "xsect_el")
runs = []
for d in dirs:
    f = os.path.join(d, "xsect_fields.npz")
    if not os.path.exists(f):
        print(f"plot_xsection_phase: skipping {d} (no xsect_fields.npz)")
        continue
    z = np.load(f, allow_pickle=True)
    if "probe_t" not in z.files:
        print(f"plot_xsection_phase: skipping {d} (npz has no probe "
              "series; re-run plot_xsection.py there)")
        continue
    runs.append((d, z))
if not runs:
    sys.exit("plot_xsection_phase: nothing to plot")

def phaseplot(comp, cbl, outname):
    zx = runs[0][1]["probe_xz"]
    npr = zx.shape[0]
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    axs = axs.ravel()
    nve = sum(1 for d, _ in runs if d != "xsect_el")
    cmap = plt.get_cmap("viridis")
    ive = 0
    for d, z in runs:
        t = z["probe_t"]; ta, tb = z["t_cyc"][0], z["t_cyc"][1]
        ph = (t - ta)/(tb - ta)
        m = (ph >= 0.0) & (ph <= 1.0)
        F = z[comp]/1e6                    # (nt, npr) [MPa]
        lab = str(z["label"])
        if d == "xsect_el":
            kw = dict(color="k", ls="--", lw=2.0, zorder=10)
        else:
            kw = dict(color=cmap(ive/max(nve - 1, 1)), lw=1.5)
            ive += 1
        for k in range(npr):
            y = F[m, k] - np.mean(F[m, k])
            axs[k].plot(ph[m], y, label=lab, **kw)
    for k in range(npr):
        x_km, z_km = zx[k, 0]/1e3, zx[k, 1]/1e3
        where = "plate" if z_km < float(runs[0][1]["H"])/1e3 else "substrate"
        axs[k].set_title(f"x = {x_km:.0f} km, z = {z_km:.0f} km ({where})",
                         fontsize=10)
        axs[k].axhline(0.0, color="0.7", lw=0.5)
        axs[k].set_ylabel(cbl + " - cycle mean [MPa]", fontsize=9)
    for k in (2, 3):
        axs[k].set_xlabel("cycle phase")
    axs[0].legend(fontsize=8, ncol=2)
    fig.suptitle(f"{cbl} at fixed points through the cycle: amplitude "
                 "and phase lag of the loading vs tM/T_rec (elastic "
                 "dashed)", fontsize=11)
    fig.tight_layout()
    fig.savefig(outname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"plot_xsection_phase: wrote {outname} ({len(runs)} models)")

phaseplot("probe_sxz", "sigma_xz", "xsect_phase_sxz.png")
phaseplot("probe_syz", "sigma_yz", "xsect_phase_syz.png")
