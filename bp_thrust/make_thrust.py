#!/usr/bin/env python3
"""
Generate a surface-breaking, listric (curved) thrust for rsf_solve.

Two cases are provided, selected by the `case` parameter below (or by a second
command-line argument):

  "ozawa"      Faithful to Ozawa et al. (2023, GJI 232, 1471; arXiv 2110.12165),
               their Fig 4 thrust: 50 km along strike, 20 km along the dip path,
               dip tapering from 30 deg at the surface to 10 deg at depth,
               surface-breaking, half-space. b = 0.020 fixed, with a-b = -0.01 in
               the velocity-weakening interior and a-b = +0.01 on the
               velocity-strengthening edges (Fig 4), d_c = 0.02 m uniform,
               initial normal and
               shear tractions uniform at 58 MPa and 100 MPa, backslip loading at
               V_pl = 1e-9 m/s for both shear and normal (the normal loading and
               evolution come from running with -calc_sigma_dot). Run with dip
               slip (-rsf_slip_mode 1). Two interpretation choices are noted
               below where the paper does not pin them down: the initial slip
               velocity, and the exact a-b map of their Fig 4.

  "validation" The original normal-stress-path check: a near-steady-state
               background with a seeded nucleation patch, used to show that
               sigma evolves with -calc_sigma_dot and stays pinned without it.

interact geometry convention (verified against calc_quad_base_vecs), strike 0:
  strike direction   = (0, 1, 0)
  down-dip direction = (cos t, 0, -sin t)   (t = dip)
  fault normal       = (sin t, 0, cos t)

Outputs (interact native formats, geometry order):
  geom_thrust.in    x y z strike dip L W group   (meters, degrees; L,W HALF-lengths)
  rsf_thrust.dat    a  b
  ic_thrust.in      tau[Pa]  v[m/s]              (for -rsf_ic_file)
  sigma_thrust.in   sigma0[Pa]                   (uniform 58 MPa in the ozawa case,
                                                  depth-dependent in the validation
                                                  case; only used with -rsf_sigma_file)

Optional command line: make_thrust.py [ds_km] [case] overrides the two most
useful knobs so a run script can sweep resolution for the mesh-convergence study
without editing this file.

Interpretation notes for the ozawa case:
  1. The paper gives a uniform initial shear traction of 100 MPa, which is well
     above the velocity-weakening steady-state level (about 37 MPa here), so the
     fault starts strongly overstressed and the first event is an initialization
     transient. Analyze the stabilized later cycles, not the first event. The
     initial slip velocity is set to V_pl (a locked, high-stress start); the
     paper does not state it explicitly.
  2. The a-b map follows Fig 4: a velocity-weakening interior at a-b = -0.01 and
     velocity-strengthening edges at a-b = +0.01. The one remaining choice is the
     width of the strengthening edge band (frame_km); set it to match the figure.
"""
import sys, math

# ---- parameters (edit here) ------------------------------------------------
case     = "ozawa"    # "ozawa" or "validation"
ds       = 1.0        # cell size along strike and along the dip path [km]

# geometry (both cases)
Lstrike  = 50.0       # along-strike length [km]
Ldip     = 20.0       # down-dip length along the curved dip path [km]
dip_top  = 30.0       # dip at the surface [deg]
dip_bot  = 10.0       # dip at the deepest row [deg]

# rate-and-state (both cases): b fixed, a-b from Ozawa Fig 4
#   interior (velocity-weakening):     a-b = -0.01  -> a = 0.010
#   edges    (velocity-strengthening): a-b = +0.01  -> a = 0.030
# (the paper text quotes a/b = 0.75, but Fig 4 shows a-b = -0.01 in the interior;
#  the figure values are used here)
b0       = 0.020      # b everywhere
a_vw     = 0.010      # a in the velocity-weakening interior (a-b = -0.01)
a_vs     = 0.030      # a in the velocity-strengthening edges (a-b = +0.01)
frame_km = 3.0        # width of the VS edge band on each side [km]; set to match Fig 4
dc       = 0.02       # characteristic slip distance [m] (used only in the printout)

f0       = 0.6        # reference friction
v0       = 1.0e-6     # reference velocity [m/s]
vpl      = 1.0e-9     # loading (plate) rate [m/s]

# ozawa case: uniform initial tractions, no seed
sigma0_ozawa = 58.0e6  # uniform initial normal stress [Pa]
tau0_ozawa   = 100.0e6 # uniform initial shear stress  [Pa]

# validation case: near-steady-state background with a seeded nucleation patch
sigma0_val   = 58.0e6
nuc_halfstrike = 4.0   # [km]
nuc_halfdip    = 4.0   # [km]
nuc_center_y   = 0.0   # [km]
nuc_center_s   = 10.0  # down-dip path distance of the patch center [km]
Vnuc           = 1.0e-2

# depth-dependent sigma for the validation -rsf_sigma_file demo
sigma_surf = 40.0e6
sigma_deep = 70.0e6
# ---------------------------------------------------------------------------

# optional command-line overrides: ds [km], case
if len(sys.argv) > 1:
    ds = float(sys.argv[1])
if len(sys.argv) > 2:
    case = sys.argv[2]
if case not in ("ozawa", "validation"):
    sys.exit(f"make_thrust.py: unknown case '{case}' (use ozawa or validation)")

KM = 1.0e3
Ns = int(round(Lstrike/ds))
Nd = int(round(Ldip/ds))
ds_strike = Lstrike/Ns
w = Ldip/Nd
Lhalf = ds_strike/2.0*KM
Whalf = w/2.0*KM

def tau_ss(a, sigma):
    return (f0 + (a-b0)*math.log(vpl/v0))*sigma

# build the curved rows: (dip_deg, x_center_km, z_center_km, s_center_km)
rows = []
x_top, z_top, s_top = 0.0, 0.0, 0.0
for r in range(Nd):
    frac = 0.0 if Nd == 1 else r/(Nd-1)
    theta = dip_top + (dip_bot-dip_top)*frac
    tr = math.radians(theta)
    cx, sx = math.cos(tr), math.sin(tr)
    rows.append((theta, x_top + 0.5*w*cx, z_top - 0.5*w*sx, s_top + 0.5*w))
    x_top += w*cx
    z_top -= w*sx
    s_top += w

fg = open("geom_thrust.in","w")
fr = open("rsf_thrust.dat","w")
fi = open("ic_thrust.in","w")
fsig = open("sigma_thrust.in","w")

nvw = nnuc = 0
for r in range(Nd):
    theta, x_c, z_c, s_c = rows[r]
    depth = -z_c
    for c in range(Ns):
        y_c = (c + 0.5)*ds_strike - Lstrike/2.0

        edge_strike = min(y_c + Lstrike/2.0, Lstrike/2.0 - y_c)
        edge_dip    = min(s_c, Ldip - s_c)
        in_frame = (edge_strike < frame_km) or (edge_dip < frame_km)
        a = a_vs if in_frame else a_vw
        if not in_frame:
            nvw += 1

        if case == "ozawa":
            tau = tau0_ozawa           # uniform overstressed start, no seed
            v   = vpl
            sig = sigma0_ozawa
        else:
            is_nuc = (abs(y_c - nuc_center_y) < nuc_halfstrike) and \
                     (abs(s_c - nuc_center_s) < nuc_halfdip) and (not in_frame)
            if is_nuc:
                tau = f0*sigma0_val
                v   = Vnuc
                nnuc += 1
            else:
                tau = tau_ss(a, sigma0_val)
                v   = vpl
            sig = sigma_surf + (sigma_deep-sigma_surf)*(depth/max(1e-9,(-rows[-1][2])))

        fg.write(f"{x_c*KM:.6e} {y_c*KM:.6e} {z_c*KM:.6e} 0.0 {theta:.6f} {Lhalf:.1f} {Whalf:.1f} 0\n")
        fr.write(f"{a:.6e} {b0:.6e}\n")
        fi.write(f"{tau:.8e} {v:.6e}\n")
        fsig.write(f"{sig:.8e}\n")

for f in (fg,fr,fi,fsig):
    f.close()

zmax = -rows[-1][2]
xmax = rows[-1][1] + 0.5*w*math.cos(math.radians(rows[-1][0]))
print(f"case={case} : {Ns} x {Nd} = {Ns*Nd} patches ; ds_strike={ds_strike:g} km w={w:g} km")
print(f"  dip {dip_top:g}->{dip_bot:g} deg ; surface-breaking ; max depth ~{zmax:.2f} km ; +x extent ~{xmax:.2f} km")
if case == "ozawa":
    Lb = 3.204e10*dc/(b0*sigma0_ozawa)
    print(f"  uniform tau0={tau0_ozawa:g} Pa sigma0={sigma0_ozawa:g} Pa ; VW cells={nvw}")
    print(f"  d_c={dc:g} m ; cohesive zone Lb = G dc/(b sigma0) ~ {Lb:.0f} m (keep ds below ~Lb/3 for convergence)")
else:
    print(f"  VW cells={nvw} ; nucleation cells={nnuc} ; validation start (seeded)")
