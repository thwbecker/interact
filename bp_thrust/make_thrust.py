#!/usr/bin/env python3
"""
Generate a surface-breaking, listric (curved) thrust for rsf_solve, to exercise
and validate the normal-stress path (-calc_sigma_dot) with dip slip
(-rsf_slip_mode 1).  The geometry follows Ozawa et al. (2023) in spirit: a thrust
that is about 50 km along strike and 20 km down dip, breaking the surface, with
the dip decreasing with depth (listric), so that dip slip produces genuine
normal-stress changes concentrated near the bend and the free surface.

interact geometry convention (verified against calc_quad_base_vecs), for
strike = 0:
  strike direction   = (0, 1, 0)          (+y, along strike)
  down-dip direction = (cos t, 0, -sin t) (+x and downward, t = dip)
  fault normal       = (sin t, 0, cos t)
So the thrust breaks the surface along y at x = 0, z = 0, and dips toward +x.

Outputs (interact native formats, geometry order):
  geom_thrust.in    x y z strike dip L W group   (meters, degrees; L,W HALF-lengths)
  rsf_thrust.dat    a  b                          (per patch)
  ic_thrust.in      tau[Pa]  v[m/s]               (per patch, for -rsf_ic_file)
  sigma_thrust.in   sigma0[Pa]                    (per patch, for -rsf_sigma_file;
                                                   depth-dependent, OPTIONAL demo)

All parameters are hardcoded below.  Note the resolution: at ds = 1 km with the
Ozawa-like D_c the cohesive zone is only marginally resolved, so this is a
machinery / path validation, not a converged benchmark; use a finer ds for
quantitative work.
"""
import math

# ---- parameters (edit here) ------------------------------------------------
ds       = 1.0        # cell size along strike and along the dip path [km]
Lstrike  = 50.0       # along-strike length [km]
Ldip     = 20.0       # down-dip length measured along the (curved) dip path [km]
dip_top  = 30.0       # dip at the surface [deg]
dip_bot  = 10.0       # dip at the deepest row [deg]

# rate-and-state (Ozawa-like); VW interior with a VS frame to keep nucleation in
a_vw     = 0.015      # a in the velocity-weakening interior
a_vs     = 0.025      # a in the velocity-strengthening frame (a > b)
b0       = 0.020      # b everywhere
frame_km = 3.0        # width of the VS frame on each edge [km]

sigma0   = 58.0e6     # uniform initial normal stress [Pa]
f0       = 0.6        # reference friction
v0       = 1.0e-6     # reference velocity [m/s]
vpl      = 1.0e-9     # loading (plate) rate [m/s]

# seeded nucleation patch (interior), to trigger slip so sigma actually evolves
nuc_halfstrike = 4.0  # half-size along strike [km]
nuc_halfdip    = 4.0  # half-size along dip   [km]
nuc_center_y   = 0.0  # along-strike center [km]
nuc_center_s   = 10.0 # down-dip path distance of the patch center [km]
Vnuc           = 1.0e-2   # nucleation-patch initial velocity [m/s]

# optional depth-dependent sigma0 for the -rsf_sigma_file demo (not used by the
# default validation run, which keeps sigma uniform so any spread is purely from
# calc_sigma_dot).  Linear increase with depth from sigma_surf to sigma_deep.
sigma_surf = 40.0e6   # [Pa] at the surface
sigma_deep = 70.0e6   # [Pa] at the deepest row
# ---------------------------------------------------------------------------

KM = 1.0e3
Ns = int(round(Lstrike/ds))     # columns along strike
Nd = int(round(Ldip/ds))        # rows down dip
ds_strike = Lstrike/Ns          # actual along-strike cell size [km]
w = Ldip/Nd                     # down-dip width per row [km]
Lhalf = ds_strike/2.0*KM        # half strike-length [m]
Whalf = w/2.0*KM                # half down-dip width [m]

# steady-state prestress at v = vpl in the VW interior: mu_ss = f0 + (a-b)ln(vpl/v0)
def tau_ss(a):
    return (f0 + (a-b0)*math.log(vpl/v0))*sigma0

rows = []          # (theta_deg, x_center_km, z_center_km, s_center_km)
x_top, z_top = 0.0, 0.0     # top edge of current row in the (x,z) dip plane [km]
s_top = 0.0                 # down-dip path distance to the top edge [km]
for r in range(Nd):
    frac = 0.0 if Nd == 1 else r/(Nd-1)
    theta = dip_top + (dip_bot-dip_top)*frac      # dip decreases with depth
    tr = math.radians(theta)
    cx, sx = math.cos(tr), math.sin(tr)
    x_c = x_top + 0.5*w*cx
    z_c = z_top - 0.5*w*sx
    s_c = s_top + 0.5*w
    rows.append((theta, x_c, z_c, s_c))
    x_top += w*cx
    z_top -= w*sx
    s_top += w

fg = open("geom_thrust.in","w")
fr = open("rsf_thrust.dat","w")
fi = open("ic_thrust.in","w")
fsig = open("sigma_thrust.in","w")

nvw = nnuc = 0
# geometry order: row-major (down-dip outer, along-strike inner), matching a
# natural grid; any consistent order is fine as long as all files agree
for r in range(Nd):
    theta, x_c, z_c, s_c = rows[r]
    depth = -z_c                                   # positive down [km]
    for c in range(Ns):
        y_c = (c + 0.5)*ds_strike - Lstrike/2.0    # centered along strike [km]

        # VS frame near the outer edges (in strike and in down-dip path)
        edge_strike = min(y_c + Lstrike/2.0, Lstrike/2.0 - y_c)
        edge_dip    = min(s_c, Ldip - s_c)
        in_frame = (edge_strike < frame_km) or (edge_dip < frame_km)
        a = a_vs if in_frame else a_vw
        if not in_frame:
            nvw += 1

        # seeded nucleation patch (interior)
        is_nuc = (abs(y_c - nuc_center_y) < nuc_halfstrike) and \
                 (abs(s_c - nuc_center_s) < nuc_halfdip) and (not in_frame)
        if is_nuc:
            tau = f0*sigma0          # at reference friction, above VW steady state
            v   = Vnuc
            nnuc += 1
        else:
            tau = tau_ss(a)
            v   = vpl

        fg.write(f"{x_c*KM:.6e} {y_c*KM:.6e} {z_c*KM:.6e} 0.0 {theta:.6f} {Lhalf:.1f} {Whalf:.1f} 0\n")
        fr.write(f"{a:.6e} {b0:.6e}\n")
        fi.write(f"{tau:.8e} {v:.6e}\n")
        # depth-dependent sigma (optional demo); depth runs 0 at surface to max
        sig = sigma_surf + (sigma_deep-sigma_surf)*(depth/max(1e-9,(-rows[-1][2])))
        fsig.write(f"{sig:.8e}\n")

for f in (fg,fr,fi,fsig):
    f.close()

zmax = -rows[-1][2]
xmax = rows[-1][1] + 0.5*w*math.cos(math.radians(rows[-1][0]))
print(f"thrust: {Ns} x {Nd} = {Ns*Nd} patches ; ds_strike={ds_strike:g} km w={w:g} km")
print(f"  dip {dip_top:g}->{dip_bot:g} deg ; surface-breaking ; max depth ~{zmax:.2f} km ; +x extent ~{xmax:.2f} km")
print(f"  VW cells={nvw} ; nucleation cells={nnuc} ; uniform sigma0={sigma0:g} Pa")
print(f"  files: geom_thrust.in rsf_thrust.dat ic_thrust.in sigma_thrust.in")
