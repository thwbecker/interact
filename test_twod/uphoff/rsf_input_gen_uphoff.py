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
import sys
import numpy as np

geom = "geom_uphoff.in"
idxfile = "geom_uphoff.in.idx"
group_profiles = {}                              # every other group:
default_profile = "profile_uphoff.dat"           # this profile (same law on all faults)
dc_by_group = {}              # D_c [m] per group; others use dc_default
dc_default = 0.05             # their L
dsigma_dz = 0.0               # their sigma_n is CONSTANT
sigma_cap = 50.0e6            # cap [Pa]
sigma_min = 50.0e6            # = cap: constant 50 MPa everywhere
vpl_by_group = {0: -1.0e-9}   # their Vp (extension, normal sense); main fault only
shear_modulus = 3.204e10      # [Pa], only for the h* diagnostics
f0 = 0.6                      # reference friction, must match -f0
v0 = 1.0e-6                   # reference velocity, must match -v0
v_ic = 1.0e-9                 # initial sliding velocity for the ICs
ic_uniform_tau = 26.5461e6    # their uniform T0 [Pa]; None: per-patch steady state
coord_unit_km = 1.0           # geometry length unit in km (1.0: km)

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

np.savetxt("rsf_ab.in", np.column_stack([a, b]), fmt="%.6f %.6f")
np.savetxt("rsf_dc.in", dc, fmt="%.6e")
np.savetxt("rsf_sigma.in", sig, fmt="%.6e")
np.savetxt("rsf_vpl.in", vpl, fmt="%.6e")
if ic_uniform_tau is None:
    tau_ic = sig*(f0 + (a - b)*np.log(v_ic/v0))
else:
    tau_ic = np.full(n, ic_uniform_tau)
np.savetxt("rsf_ic.in", np.column_stack([tau_ic, np.full(n, v_ic)]),
           fmt="%.6e %.6e")

print(f"wrote rsf_ab.in rsf_dc.in rsf_sigma.in rsf_vpl.in rsf_ic.in for {n} patches")
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
