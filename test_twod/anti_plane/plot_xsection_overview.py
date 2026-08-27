#!/usr/bin/env python3
"""
plot_xsection_overview.py: assemble the run_xsection series into one
overview figure, model rows (elastic on top, then increasing
tM/T_rec), snapshot columns (the ACTUAL post-event field, then the
change through the cycle at phases 0.25, 0.50, 0.75, pre), so the
effect of viscous relaxation reads directly down the columns.  Uses
the per-run xsect_fields.npz files written by plot_xsection.py (run
that first), a shared DISCRETE color scale (21 subdivisions) for all
change panels and a second shared discrete scale for the post-event
column.

usage: plot_xsection_overview.py [dir1 dir2 ...] (SAT env overrides the
       0.99 color-scale saturation quantile)
       default: xsect_el xsect_r0.1 xsect_r0.5 xsect_r1 xsect_r2 xsect_r10
writes xsect_overview_sxy.png (and _syz)
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

NLEV = 21          # target discrete color subdivisions (nice-number levels)
from matplotlib.ticker import MaxNLocator
def discrete(vmax):
    """symmetric discrete levels with round values, about NLEV bins"""
    lev = MaxNLocator(nbins=NLEV, symmetric=True).tick_values(-vmax, vmax)
    return lev, BoundaryNorm(lev, len(lev)-1), plt.get_cmap("RdBu_r", len(lev)-1)

sat = float(os.environ.get("SAT", 0.99))
dirs = sys.argv[1:] if len(sys.argv) > 1 else \
    ["xsect_el", "xsect_r0.1", "xsect_r0.5", "xsect_r1", "xsect_r2",
     "xsect_r10"]
yr = 3.15576e7
runs = []
for d in dirs:
    f = os.path.join(d, "xsect_fields.npz")
    if not os.path.exists(f):
        print(f"plot_xsection_overview: skipping {d} (no xsect_fields.npz; "
              "run plot_xsection.py there first)")
        continue
    runs.append((d, np.load(f, allow_pickle=True)))
if not runs:
    sys.exit("plot_xsection_overview: nothing to plot")

def overview(comp, cbl, outname):
    nr = len(runs); ncol = 5
    fig, axs = plt.subplots(nr, ncol, figsize=(3.1*ncol + 1.6, 2.6*nr + 1.6),
                            sharex=True, sharey=True, squeeze=False,
                            constrained_layout=True)
    # shared scales over ALL runs
    va = vd = 1e-12
    for _, z in runs:
        F = z[comp]                       # (5, nz, nx): post, 0.25, ..., pre
        va = max(va, np.quantile(np.abs(F[0]), sat)/1e6)
        vd = max(vd, np.quantile(np.abs(F[1:] - F[0]), sat)/1e6)
    _, nrm0, cmp0 = discrete(va)
    _, nrmd, cmpd = discrete(vd)
    im0 = imd = None
    for r, (d, z) in enumerate(runs):
        xg, zg = z["xg"]/1e3, z["zg"]/1e3
        H, D = float(z["H"])/1e3, float(z["D"])
        F = z[comp]; ts = z["t_snap"]/yr; ph = z["phases"]
        lab = str(z["label"])
        im0 = axs[r, 0].pcolormesh(xg, zg, F[0]/1e6, cmap=cmp0, norm=nrm0,
                                   shading="auto")
        for c in range(1, ncol):
            imd = axs[r, c].pcolormesh(xg, zg, (F[c] - F[0])/1e6, cmap=cmpd,
                                       norm=nrmd, shading="auto")
        for c in range(ncol):
            ax = axs[r, c]
            ax.axhline(H, color="k", lw=0.7, ls="--")
            ax.plot([0, 0], [0, D], color="k", lw=2.0)
            if r == 0:
                ax.set_title("post-event state" if c == 0 else
                             f"change, phase {ph[c]}", fontsize=9)
            if r == nr - 1:
                ax.set_xlabel("x [km]")
        axs[r, 0].set_ylabel(f"{lab}\ndepth [km]", fontsize=9)
    axs[0, 0].set_ylim(zg[-1], 0)
    fig.colorbar(im0, ax=axs[-1, 0], location="bottom", fraction=0.15,
                 label=cbl + " (state) [MPa]", extend="both")
    fig.colorbar(imd, ax=axs[-1, 1:], location="bottom", fraction=0.15,
                 shrink=0.5, label="d " + cbl + " [MPa]", extend="both")
    fig.suptitle(f"{cbl} through the last cycle: post-event state and the "
                 "interseismic change at cycle phases 0.25/0.50/0.75/pre; "
                 "dashed: plate base; thick: fault", fontsize=11)
    fig.savefig(outname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"plot_xsection_overview: wrote {outname} ({nr} models)")

overview("sxy", "sigma_xy", "xsect_overview_sxy.png")
overview("syz", "sigma_yz", "xsect_overview_syz.png")
