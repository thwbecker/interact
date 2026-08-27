#!/usr/bin/env python3
"""
plot_cycles2d.py: render figures from a run_cycles2d output directory.

Reads the field frames tmp_rsf/rsf_vel.gNNN.FFFFFF.bin and (if present)
tmp_rsf/rsf_tau.gNNN.FFFFFF.bin (float32 (x1, x2, value) triples, one
frame per accepted-step interval; value = log10|v| [m/s] or tau [MPa]),
frame times from rsf_vel.times, per-fault event catalogs from
rsf_catalog.gNNN.dat, and geometry from rsf_geom.gNNN.dat.

Figures written into the run directory:
  spacetime_v.png     log10 slip rate, depth vs time (per fault, up to 8)
  spacetime_tau.png   shear stress, depth vs time (per fault, up to 8)
  stress_profiles.png tau(z) cross-section evolution through the last
                      full cycle of fault 0 (interseismic sequence,
                      just before, and just after the event)
  event_raster.png    (nfault > 1) event onsets per fault vs time,
                      marker size by mean coseismic slip

usage: plot_cycles2d.py rundir [tmin_yr]
       tmin_yr crops early spin-up from the spacetime plots (default 0)
"""
import sys, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rundir = sys.argv[1] if len(sys.argv) > 1 else "."
tmin   = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
os.chdir(rundir)

# frame times: frame step time[yr] time[s] log10vmax ...
T = np.loadtxt("rsf_vel.times", comments="#")
if T.ndim == 1: T = T[None, :]
frames = T[:, 0].astype(int); t_yr = T[:, 2]; l10vmax = T[:, 4]

groups = sorted(int(f.split(".g")[1].split(".")[0])
                for f in glob.glob("rsf_geom.g*.dat"))
ngr = len(groups)

def read_frames(pref, g):
    """(nframe, ncell) value array for group g, plus the depth axis [km].
    Frame rows follow the geometry-file row order, so the depth of cell k
    is taken from rsf_geom.gNNN.dat column yc (antiplane: y = -depth)."""
    fl = sorted(glob.glob(f"tmp_rsf/{pref}.g{g:03d}.*.bin"))
    if not fl: return None, None
    geo = np.loadtxt(f"rsf_geom.g{g:03d}.dat", comments="#")
    if geo.ndim == 1: geo = geo[None, :]
    dep = np.abs(geo[:, 4])/1e3         # -yc -> depth [km]
    V = np.empty((len(fl), dep.size), dtype=np.float32)
    for i, f in enumerate(fl):
        V[i] = np.fromfile(f, dtype="<f4").reshape(-1, 3)[:, 2]
    # sort by depth so pcolormesh gets a monotonic axis
    order = np.argsort(dep)
    return V[:, order], dep[order]

def spacetime(pref, cmap, label, outname, vlim=None, xaxis="time"):
    """xaxis 'time': calendar time (coseismic phases collapse to lines);
    xaxis 'frame': uniform frame index (the solver's own step density,
    so nucleation and rupture are resolved; time ticks annotated)"""
    nshow = min(ngr, 8)
    fig, axs = plt.subplots(nshow, 1, figsize=(11, 2.2*nshow + 1.2),
                            sharex=True, squeeze=False)
    im = None
    for k in range(nshow):
        g = groups[k]
        V, dep = read_frames(pref, g)
        if V is None: plt.close(fig); return False
        n = min(V.shape[0], len(t_yr))
        sel = t_yr[:n] >= tmin
        tt = t_yr[:n][sel]; VV = V[:n][sel]
        if xaxis == "time":
            xe = np.concatenate([[tt[0]], 0.5*(tt[1:]+tt[:-1]), [tt[-1]]])
        else:
            xe = np.arange(len(tt)+1) - 0.5
        de = np.concatenate([[dep[0]-(dep[1]-dep[0])/2],
                             0.5*(dep[1:]+dep[:-1]),
                             [dep[-1]+(dep[-1]-dep[-2])/2]])
        ax = axs[k, 0]
        kw = dict(cmap=cmap)
        if vlim: kw.update(vmin=vlim[0], vmax=vlim[1])
        im = ax.pcolormesh(xe, de, VV.T, shading="flat", **kw)
        ax.invert_yaxis()
        ax.set_ylabel(f"fault {g}\ndepth [km]", fontsize=8)
        if xaxis == "frame" and k == nshow-1:
            # annotate calendar time at a few frame positions
            ii = np.linspace(0, len(tt)-1, 8).astype(int)
            ax.set_xticks(ii)
            ax.set_xticklabels([f"{tt[i]:.1f}" for i in ii], fontsize=7)
    axs[-1, 0].set_xlabel("time [yr]" if xaxis == "time"
                          else "frame index (time [yr] annotated)")
    fig.colorbar(im, ax=axs[:, 0], label=label, shrink=0.8, pad=0.01)
    fig.suptitle(f"{os.path.basename(os.getcwd())}: {label}", fontsize=11)
    fig.savefig(outname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True

ok_v = spacetime("rsf_vel", "magma", "log10 |v| [m/s]", "spacetime_v.png",
                 vlim=(-12, 0))
spacetime("rsf_vel", "magma", "log10 |v| [m/s]", "spacetime_v_frames.png",
          vlim=(-12, 0), xaxis="frame")
ok_t = spacetime("rsf_tau", "viridis", "shear stress [MPa]", "spacetime_tau.png")
spacetime("rsf_tau", "viridis", "shear stress [MPa]", "spacetime_tau_frames.png",
          xaxis="frame")

# stress cross-section evolution through the last full cycle of fault 0
if ok_t:
    V, dep = read_frames("rsf_tau", groups[0])
    n = min(V.shape[0], len(t_yr))
    # event frames on fault 0: use global log10 vmax crossings of -3
    seismic = l10vmax[:n] > -3.0
    onset_idx = np.where(seismic & ~np.roll(seismic, 1))[0]
    if len(onset_idx) >= 2:
        i0, i1 = onset_idx[-2], onset_idx[-1]
        # profiles: just after event at i0, evenly through the cycle,
        # just before and just after the event at i1
        after0 = i0
        while after0 < n-1 and seismic[after0]: after0 += 1
        pick = list(np.linspace(after0, i1-1, 6).astype(int)) + [i1-1]
        after1 = i1
        while after1 < n-1 and seismic[after1]: after1 += 1
        fig, ax = plt.subplots(figsize=(6.4, 6.2))
        cm = plt.cm.plasma(np.linspace(0.05, 0.85, len(pick)))
        for c, i in zip(cm, pick):
            ax.plot(V[i], dep, color=c, lw=1.2,
                    label=f"t = {t_yr[i]:.1f} yr")
        ax.plot(V[after1], dep, "k--", lw=1.4,
                label=f"post-event, t = {t_yr[after1]:.1f} yr")
        ax.invert_yaxis()
        ax.set_xlabel("shear stress [MPa]"); ax.set_ylabel("depth [km]")
        ax.set_title("fault 0: stress cross-section through the last cycle")
        ax.legend(fontsize=7.5)
        fig.savefig("stress_profiles.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

# multi-fault event raster
cats = sorted(glob.glob("rsf_catalog.g*.dat"))
if ngr > 1 and cats:
    fig, ax = plt.subplots(figsize=(11, 0.28*ngr + 2.2))
    for f in cats:
        g = int(f.split(".g")[1].split(".")[0])
        rows = [l.split() for l in open(f) if not l.startswith("#")]
        tt = np.array([float(r[1]) for r in rows if int(r[5]) > 3])
        us = np.array([float(r[7]) for r in rows if int(r[5]) > 3])
        if len(tt):
            ax.scatter(tt, np.full(len(tt), g), s=8 + 25*us/max(us.max(), 1e-9),
                       c="crimson", alpha=0.8, edgecolors="none")
    ax.set_xlabel("time [yr]"); ax.set_ylabel("fault index")
    ax.set_title(f"{os.path.basename(os.getcwd())}: event onsets "
                 "(marker size by mean coseismic slip)")
    fig.savefig("event_raster.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

print("plot_cycles2d: done")
