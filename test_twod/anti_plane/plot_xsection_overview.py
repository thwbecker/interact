#!/usr/bin/env python3
"""
plot_xsection_overview.py: assemble the run_xsection series into one
overview figure, model rows ordered by INCREASING Maxwell time top to
bottom, with the elastic model as the BOTTOM row: elastic is the
tM -> infinity limit, so the r = tM/T_rec = 10 row should visibly
converge to it, and the effect of viscous relaxation grows reading
UP the figure.  All panels show the ACTUAL reconstructed stress field
(the loading-relevant sigma_xz, and sigma_yz) at the snapshot times
through the last cycle (right after the event, cycle phases
0.25/0.50/0.75, and right before the next event), on ONE shared
discrete color scale (~21 nice-number subdivisions, quantile
saturated) so values read off the plot and rows compare directly.

Uses the per-run xsect_fields.npz files written by plot_xsection.py
(run that first).

usage: plot_xsection_overview.py [mode] [sat] [dir1 dir2 ...]
       mode: cyc (default) subtracts each run's own cycle-mean field;
             abs plots the actual reconstructed fields (see caveat);
             dif plots each VE model MINUS THE ELASTIC model, both
             cycle-mean-referenced (the pure viscous perturbation of
             the cycle pattern; elastic must be among the dirs and
             gets no row of its own)
       sat:  saturation quantile of the shared color scale (0.99)
       dirs: xsect_r0.1 xsect_r0.5 xsect_r1 xsect_r2 xsect_r10 xsect_el
       (all optional; a leading cyc/abs is the mode, a leading number
       the quantile, everything else directories)
writes xsect_overview_sxz.png (and _syz); mode abs appends _abs, and
its row labels carry each run's standing deficit (lag) at the plotted
snapshots

CAVEAT on mode abs: the absolute field includes each run's STANDING
slip-deficit level, which is genuine solver output but only partly
comparable across models: the elastic run's level is pinned by the
seed/initial stress forever, and the viscoelastic level equilibrates
on the slow tail of the relaxation ladder (tens of Maxwell times),
far beyond these run durations for tM >~ T_rec.  So in MODE=abs the
large-tM rows need not visually converge to the elastic row even
though the physics does in the tM -> inf limit.  MODE=cyc removes
that run-specific DC and shows the actual stress pattern of the cycle
itself; there the large-tM rows do converge to elastic.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

NLEV = 21          # target discrete color subdivisions (nice-number levels)
args = sys.argv[1:]
MODE = "cyc"
if args and args[0] in ("cyc", "abs", "dif"):
    MODE = args.pop(0)
SAT = 0.99
if args:
    try:
        SAT = float(args[0]); args.pop(0)
    except ValueError:
        pass
def discrete(vmax):
    """symmetric discrete levels with round values, about NLEV bins"""
    lev = MaxNLocator(nbins=NLEV, symmetric=True).tick_values(-vmax, vmax)
    return BoundaryNorm(lev, len(lev)-1), plt.get_cmap("RdBu_r", len(lev)-1)

def getcomp(z, comp):
    """npz key with fallback to the pre-rename key (sxz was sxy)"""
    if comp in z.files:
        return z[comp]
    return z[{"sxz": "sxy"}.get(comp, comp)]

dirs = args if args else \
    ["xsect_r0.1", "xsect_r0.5", "xsect_r1", "xsect_r2", "xsect_r10",
     "xsect_el"]
# elastic (tM -> inf) belongs at the bottom regardless of argument order
dirs = sorted(dirs, key=lambda d: d == "xsect_el")
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
    Fel = None
    if MODE == "dif":
        el = [z for d, z in runs if str(z["label"]) == "elastic"]
        if not el:
            sys.exit("plot_xsection_overview: mode dif needs the elastic "
                     "run among the directories")
        Fel = np.array(getcomp(el[0], comp), dtype=float)
        Fel = Fel - np.mean(Fel, axis=0)
    def prep(z):
        F = np.array(getcomp(z, comp), dtype=float)
        if MODE in ("cyc", "dif"):
            F = F - np.mean(F, axis=0)
        if MODE == "dif":
            F = F - Fel
        return F
    rows = runs if MODE != "dif" else \
        [(d, z) for d, z in runs if str(z["label"]) != "elastic"]
    nr = len(rows); ncol = 5
    fig, axs = plt.subplots(nr, ncol, figsize=(3.1*ncol + 1.4, 2.6*nr + 1.0),
                            sharex=True, sharey=True, squeeze=False)
    # one shared, saturated scale over ALL runs and ALL snapshots
    va = 1e-12
    for _, z in rows:
        va = max(va, np.quantile(np.abs(prep(z)), SAT)/1e6)
    nrm, cmp_ = discrete(va)
    im = None
    for r, (d, z) in enumerate(rows):
        xg, zg = z["xg"]/1e3, z["zg"]/1e3
        H, D = float(z["H"])/1e3, float(z["D"])
        F = prep(z); ts = z["t_snap"]/yr; ph = z["phases"]
        lab = str(z["label"])
        for c in range(ncol):
            ax = axs[r, c]
            im = ax.pcolormesh(xg, zg, F[c]/1e6, cmap=cmp_, norm=nrm,
                               shading="auto")
            ax.axhline(H, color="k", lw=0.7, ls="--")
            ax.plot([0, 0], [0, D], color="k", lw=2.0)
            if r == 0:
                ax.set_title("post-event" if c == 0 else
                             ("pre-event" if c == ncol - 1 else
                              f"phase {ph[c]}"), fontsize=9)
            if r == nr - 1:
                ax.set_xlabel("x [km]")
        if MODE == "abs" and "dmean" in z.files:
            lab += f"\nlag {-float(z['dmean']):.0f} m"
        axs[r, 0].set_ylabel(f"{lab}\ndepth [km]", fontsize=9)
    axs[0, 0].set_ylim(zg[-1], 0)
    fig.colorbar(im, ax=axs, label=cbl + " [MPa]", shrink=0.55, pad=0.015,
                 aspect=40, extend="both")
    what = {"abs": "actual reconstructed field",
            "cyc": "actual field around each run's own cycle mean",
            "dif": "cycle pattern minus the elastic cycle pattern"}[MODE]
    rowtxt = ("rows: increasing tM/T_rec (elastic reference subtracted)"
              if MODE == "dif" else
              "rows: increasing tM/T_rec, elastic = tM -> inf limit at bottom")
    fig.suptitle(f"{cbl} through the cycle ({what}; {rowtxt}; dashed: "
                 "plate base; thick: fault)", fontsize=11)
    fig.savefig(outname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"plot_xsection_overview: wrote {outname} ({nr} models)")

tag = {"cyc": "", "abs": "_abs", "dif": "_dif"}[MODE]
overview("sxz", "sigma_xz", f"xsect_overview_sxz{tag}.png")
overview("syz", "sigma_yz", f"xsect_overview_syz{tag}.png")
