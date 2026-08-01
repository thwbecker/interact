#!/usr/bin/env python3
"""
Generate consistent per-patch rsf_solve input files from a
faults2geom-style geometry (group ID in the last geom.in column, with
the companion .idx file) and per-group, depth-dependent friction
profiles:

  rsf_ab.in     a b per patch            (-rsf_file)
  rsf_dc.in     D_c [m] per patch        (-rsf_dc_file)
  rsf_sigma.in  sigma_eff [Pa] per patch (-rsf_sigma_file)
  rsf_vpl.in    plate rate [m/s] per patch (-rsf_vpl_file); by default
                only the megathrust (group 0) is loaded by backslip and
                the splays are driven purely by stress transfer
  rsf_ic.in     initial "tau vel" per patch (-rsf_ic_file): steady
                state shear stress at v_ic under each patch's own a, b,
                and sigma, tau = sigma (f0 + (a-b) ln(v_ic/v0)); without
                this, a global tau0 combined with per-cell sigma leaves
                shallow patches far above failure and the integrator
                storms at t = 0. f0 and v0 below MUST match the -f0 and
                -v0 run flags

Friction profiles are per-group files of control points "depth_km a b"
(linear interpolation, end values extended); assign groups to profile
files in the group_profiles table below. Effective normal stress is a
capped gradient min(dsigma_dz * depth, sigma_cap). The script prints,
per group, the a-b range, the stability classification census in the
sense of critical-stiffness arguments (velocity strengthening b-a < 0
stable; weakening patches classified by the ratio of the process zone
scale to the element size), and the discretization check h* / dx with
h* = G D_c / ((b-a) sigma), warning where elements are coarser than
h*/3. All parameters are plain assignments below; no environment
variables.
"""
import os
import sys
import numpy as np

profile = sys.argv[1] if len(sys.argv) > 1 else "PD06"
geom = "geom_%s.in" % profile
idxfile = geom + ".idx"
variants = ["nosigma", "sigma"]
                          # one run directory per section AND per normal
                          # stress treatment: <profile>_nosigma holds sigma
                          # fixed at the values in rsf_sigma.in, <profile>_sigma
                          # evolves it through the In interaction matrix.  The
                          # per-patch input files are IDENTICAL in the two, so
                          # they are written twice from the same arrays and the
                          # only difference between the runs is the flag set in
                          # run_cycle.  Keeping the inputs byte identical is the
                          # point: it makes the two catalogs directly comparable
outdirs = [profile + "_" + v + "/" for v in variants]
outdir = outdirs[0]
outplot = outdir + "rsf_input_%s.png" % profile
make_plot = 1     # write outplot (needs matplotlib); 0: off
group_profiles = {}                              # megathrust (group 0) and the
                                                 # splays both take the default
                                                 # below: one geotherm, one law.
                                                 # To give the splays their own,
                                                 # put the megathrust file here
                                                 # as {0: "..."} and change the
                                                 # default
default_profile = "profile_cascadia_TzT.dat"     # a-b(z) from the subduction
                                                 # geotherm, via ab2profile.py
#default_profile = "profile_cascadia.dat"        # the shallower k_c/k estimate
                                                 # read off the figure; same
                                                 # transition depths to within
                                                 # 0.3 and 1.6 km, but a-b
                                                 # reaches -0.0060 rather than
                                                 # -0.0035, so events come out
                                                 # larger and less frequent
#default_profile = "profile_cascadia_deep.dat"   # the deeper figure estimate
bandfile = "stability_bands.dat"   # digitized k_c/k bands, drawn on the plot
                                   # and used for the coverage check below;
                                   # None: skip
band_estimate = "shallow"          # which estimate bandfile row set to use
dc_by_group = {}              # D_c [m] per group; others use dc_default
dc_default = 0.05             # sets h* with a-b and sigma; see the census
dsigma_dz = 18.0e6            # effective normal stress gradient [Pa/km],
                              # lithostatic minus hydrostatic
sigma_cap = 50.0e6            # cap [Pa]
sigma_min = 2.0e6             # floor [Pa] for the shallowest patches

slip_sense = 1.0              # thrust, the sense of a subduction megathrust
#slip_sense = -1.0            # sign of slip in interact's 2D convention that
                              # corresponds to normal (extensional) faulting on
                              # this geometry.  vpl, v_ic and tau_ic below all
                              # carry it, so the initial state sits on the same
                              # branch of the (odd) friction law that the
                              # backslip loading drives towards.  Confirm the
                              # sense with a static forward run before trusting
                              # a normal-stress-evolving cycle: with
                              # -calc_sigma_dot off the two signs are related by
                              # an exact symmetry, with it on they are not

vpl_by_group = {0: slip_sense*1.2e-9}   # ~38 mm/yr, Juan de Fuca convergence,
                                        # on the interface only; the splays are
                                        # driven purely by stress transfer.
                                        # Same order as the 1e-9 of the Uphoff
                                        # runs, so recurrence times stay
                                        # comparable
shear_modulus = 3.0e10        # [Pa], only for the h* diagnostics
f0 = 0.6                      # reference friction, must match -f0
v0 = 1.0e-6                   # reference velocity, must match -v0
v_ic = slip_sense*1.0e-9      # initial sliding velocity for the ICs
ic_uniform_tau = None         # per-patch steady state; a uniform tau0 over
                              # a depth-varying sigma freezes the fault
coord_unit_km = 0.001         # geometry length unit in km (0.001: meters)

rows = np.loadtxt(geom, usecols=(0, 1, 5, 7))
gids = rows[:, 3].astype(int)
depth_km = -rows[:, 1] * coord_unit_km
hl_km = rows[:, 2] * coord_unit_km
n = len(rows)

profiles = {}
def get_profile(gid):
    fname = group_profiles.get(gid, default_profile)
    if fname not in profiles:
        p = np.loadtxt(fname)
        profiles[fname] = p[np.argsort(p[:, 0])]
    return profiles[fname]

a = np.zeros(n); b = np.zeros(n)
dc = np.zeros(n); sig = np.zeros(n); vpl = np.zeros(n)
for i in range(n):
    p = get_profile(gids[i])
    a[i] = np.interp(depth_km[i], p[:, 0], p[:, 1])
    b[i] = np.interp(depth_km[i], p[:, 0], p[:, 2])
    dc[i] = dc_by_group.get(gids[i], dc_default)
    sig[i] = min(max(dsigma_dz*depth_km[i], sigma_min), sigma_cap)
    vpl[i] = vpl_by_group.get(gids[i], 0.0)

for od in outdirs:
    os.makedirs(od, exist_ok=True)
    np.savetxt(od + "rsf_ab.in", np.column_stack([a, b]), fmt="%.6f %.6f")
    np.savetxt(od + "rsf_dc.in", dc, fmt="%.6e")
    np.savetxt(od + "rsf_sigma.in", sig, fmt="%.6e")
    np.savetxt(od + "rsf_vpl.in", vpl, fmt="%.6e")
if ic_uniform_tau is None:
    tau_ic = slip_sense*sig*(f0 + (a - b)*np.log(abs(v_ic)/v0))
else:
    tau_ic = np.full(n, ic_uniform_tau)
for od in outdirs:
    np.savetxt(od + "rsf_ic.in", np.column_stack([tau_ic, np.full(n, v_ic)]),
               fmt="%.6e %.6e")

print("wrote rsf_ab.in rsf_dc.in rsf_sigma.in rsf_vpl.in rsf_ic.in for "
      f"{n} patches in " + " ".join(outdirs))
print(f"  sigma {sig.min()/1e6:.2f} to {sig.max()/1e6:.2f} MPa; in the sigma "
      "variant the limiter bounds in run_cycle have to bracket this")
print(f"{'group':>5s} {'npatch':>6s} {'depth km':>12s} {'a-b range':>18s} "
      f"{'VW':>4s} {'VS':>4s} {'min h*/dx':>9s}")
for g in sorted(set(gids)):
    sel = gids == g
    amb = a[sel] - b[sel]
    vw = amb < 0
    hstar = np.full(sel.sum(), np.inf)
    if vw.any():
        hstar[vw] = shear_modulus*dc[sel][vw] / (-amb[vw]*sig[sel][vw]) / 1e3
    dx = 2.0*hl_km[sel]
    ratio = hstar/dx
    flag = "  <-- coarser than h*/3" if vw.any() and ratio[vw].min() < 3.0 else ""
    print(f"{g:5d} {sel.sum():6d} {depth_km[sel].min():5.1f}-{depth_km[sel].max():5.1f} "
          f"{amb.min():+8.4f}..{amb.max():+8.4f} {vw.sum():4d} {(~vw).sum():4d} "
          f"{(ratio[vw].min() if vw.any() else float('inf')):9.1f}{flag}")
print("h* = G D_c / ((b-a) sigma) on velocity-weakening patches; keep")
print("h*/dx well above ~3 (raise D_c or refine dx where flagged)")
lb = shear_modulus*dc/(b*sig)/1e3
print(f"L_b = G D_c / (b sigma) is {lb.min():.2f} to {lb.max():.2f} km; "
      f"element length {2*hl_km.min():.2f} to {2*hl_km.max():.2f} km")
for g in sorted(set(gids)):
    print(f"  group {g:2d}: {group_profiles.get(g, default_profile)}")

# the figure's velocity-weakening window against what the mesh covers.  A
# fault whose downdip tip is still inside the weakening window has no
# stabilizing region below the seismogenic zone, so ruptures nucleate on
# and run into the tip, which is a numerical edge and not a physical one
bands = None
if bandfile:
    try:
        rec = [ln.split() for ln in open(bandfile)
               if ln.strip() and not ln.startswith("#")]
        bands = [(r[1], float(r[2]), float(r[3])) for r in rec if r[0] == band_estimate]
    except (OSError, IndexError, ValueError):
        bands = None
if bands:
    red = [q for q in bands if q[0] == "red"]
    if red:
        rtop, rbot = red[0][1], red[-1][2]
        zmax = depth_km.max()
        print(f"figure ({band_estimate} estimate): weakening band {rtop:.2f} to {rbot:.2f} km depth")
        if zmax < rbot:
            print(f"  WARNING: the mesh bottoms out at {zmax:.1f} km, inside that band. "
                  f"There is no downdip velocity-strengthening region in this model, so "
                  f"ruptures will run into the fault tip. Extend the interface past "
                  f"{rbot:.0f} km (ext_depth_km in make_faults_cascadia.py) or accept "
                  f"that the downdip limit here is the mesh, not the friction")
        else:
            print(f"  mesh reaches {zmax:.1f} km, so the downdip transition is inside the model")

# ---------------------------------------------------------------- plot
# quick look at the per-patch fields actually written.  The signs matter
# here: tau_ic, v_ic and vpl all have to sit on the same branch of the
# (odd) friction law, so a sign slip shows up as tau_ic and vpl on
# opposite sides of zero in the right-hand panel
if make_plot:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping " + outplot)
    else:
        cols = plt.cm.tab10(np.arange(10))
        fig, ax = plt.subplots(1, 3, figsize=(13, 4.6), sharey=True)
        for k, g in enumerate(sorted(set(gids))):
            s = gids == g
            o = np.argsort(depth_km[s])
            z = depth_km[s][o]
            ax[0].plot(a[s][o], z, "-", color=cols[k % 10], lw=1.4, label=f"group {g}")
            ax[0].plot(b[s][o], z, "--", color=cols[k % 10], lw=1.0)
            ax[1].plot((a - b)[s][o], z, "-", color=cols[k % 10], lw=1.4)
            ax[2].plot(tau_ic[s][o]/1e6, z, "-", color=cols[k % 10], lw=1.4)
            ax[2].plot(sig[s][o]/1e6, z, ":", color=cols[k % 10], lw=1.0)
        vw = (a - b) < 0
        if vw.any():
            ax[1].axhspan(depth_km[vw].min(), depth_km[vw].max(),
                          color="0.85", zorder=0)
            ax[1].set_title("a - b (solid), velocity weakening shaded\n"
                            f"{depth_km[vw].min():.2f} to {depth_km[vw].max():.2f} km")
        # the digitized k_c/k bands, so the profile can be checked against
        # the figure it came from
        if bands:
            face = {"darkblue": "#202657", "blue": "#25518f",
                    "cream": "#ece7e2", "red": "#a83b1f", "darkred": "#4e1318"}
            for nm, z0, z1 in bands:
                ax[0].axhspan(z0, z1, color=face.get(nm, "0.8"), alpha=0.35,
                              lw=0, zorder=0)
        ax[1].axvline(0.0, color="k", lw=0.8)
        ax[2].axvline(0.0, color="k", lw=0.8)
        ax[0].set_xlabel("a (solid), b (dashed)")
        ax[0].set_ylabel("depth [km]")
        ax[0].legend(fontsize=6, ncol=2)
        ax[0].set_title(f"{n} patches, D_c = {dc.min():g} to {dc.max():g} m")
        ax[1].set_xlabel("a - b")
        ax[2].set_xlabel("tau_ic (solid), sigma (dotted) [MPa]")
        ax[2].set_title(f"v_ic = {v_ic:+.1e} m/s, vpl = "
                        f"{min(vpl):+.1e} to {max(vpl):+.1e} m/s")
        ax[0].invert_yaxis()
        for q in ax:
            q.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(outplot, dpi=140)
        print(f"wrote {outplot}")
