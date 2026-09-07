#!/usr/bin/env python3
"""
gen_bp3.py: generate the SEAS BP3-QD configuration (Erickson et al.,
BSSA 2023) for rsf_solve: 2-D PLANE-STRAIN dipping fault in a
half-plane (interact's TWO_DIM_HALFPLANE_PLANE_STRAIN Crouch &
Starfield segments with free surface), quasi-dynamic rate-and-state
cycles.

BP3: fault frictional from the surface trace to 40 km DOWN-DIP;
a = 0.010 above 15 km (down-dip), linear to amax = 0.025 over
15-18 km, amax to 40 km; b = 0.015, sigma_n = 50 MPa, Dc = 8 mm,
V0 = 1e-6, f0 = 0.6, Vp = 1e-9 m/s, G = 32.04 GPa (rho 2670,
cs 3464 m/s), nu = 0.25; steady creep at Vp below 40 km (implicit
in the backslip formulation); BP1's uniform initial condition
(prestress at the amax steady state of Vinit = Vp).

Geometry convention: segments in the (x, y<=0) half-plane; a 2-D
in-plane element has w = 0, z = 0, dip = +90 (mode flag); the
segment's in-plane direction is set by strike: alpha = 90 - strike
(CCW from +x), so a fault dipping at theta toward +x uses
strike = 90 + theta (tangent (cos theta, -sin theta), pointing
down-dip).  Positive along-segment (STRIKE-component) slip with this
tangent is hanging-wall-UP, i.e. THRUST (reverse) sense, which is
also the spec's sign convention (eq 6: thrust = positive slip; Table
1: Vp, VL, Vinit positive for thrust).  Checked with interact itself
on the dip-60 geometry at 2 km cells: a unit of positive slip on all
segments uplifts the hanging wall (the side the fault dips toward)
by 0.51 m at 5 km from the trace and subsides the footwall by 0.29
m; the shear and normal stress changes on the fault agree with an
independent Okada (1992) evaluation to 0.1 percent at every cell,
and positive slip UNclamps the shallow fault (Oglesby et al., 1998),
as thrust slip must.  An earlier version of this docstring had the
sense inverted (a prototype with the hanging wall on the wrong side);
runs made under that label are physically the opposite sense, see
README_bp3.md.

BRANCH.  The sense of faulting is set by the SIGN OF THE SLIP against
this fixed tangent, i.e. by the sign of vpl, because that is what
reverses the slip-induced normal-stress change.  branch +1 (vpl > 0)
is THRUST, branch -1 (vpl < 0) is NORMAL.  It matters: at 100 m the
elastic sigma-coupled cycle is period 2 with 66.8/87.8 yr and
sigma_n rising interseismically (the thrust deficit clamps the
shallow fault) on the thrust branch, and near a single 94 yr interval
with sigma_n falling to 35 MPa on the normal branch.

Two traps this argument exists to close:

  * flipping the segment tangent does NOT switch the branch.  The
    normal is derived from the tangent, so both flip together, and a
    dislocation with Burgers vector +b along t across a plane with
    normal n is the same physical dislocation as -b along -t across
    normal -n.  Both interaction matrices are invariant and the run
    is unchanged (checked: identical event times to four decimals
    over 950 yr).  An earlier version of this script flipped the
    tangent for sense = -1, which is why that had no effect.
  * the initial condition has to be mirrored with the loading.  The
    prestress and Vinit below are the steady state at +Vp; used with
    -vpl -1e-9 they start the fault on the opposite branch, and at
    100 m the result is NO events in 1500 yr, with sigma_n creeping
    down to 27 MPa until the velocity-weakening patch is subcritical.

usage: gen_bp3.py ds_km dip_deg [prefix] [branch]
       branch +1 (default) THRUST, or -1 NORMAL: -1 negates the
       prestress and Vinit.  The caller MUST pass the matching -vpl,
       which the script prints; run_bp3 and run_bp3ve_demo derive
       both from one flag so they cannot disagree.
writes <prefix>_{geom.in,rsf.dat,ic.in}
"""
import sys
import numpy as np

ds = float(sys.argv[1])*1e3 if len(sys.argv) > 1 else 200.0
dip = np.deg2rad(float(sys.argv[2])) if len(sys.argv) > 2 else np.deg2rad(60.0)
pref = sys.argv[3] if len(sys.argv) > 3 else "bp3"
branch = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
if branch not in (1.0, -1.0):
    sys.exit("gen_bp3: branch must be +1 (thrust) or -1 (normal)")

Wf = 40e3                 # frictional down-dip extent
H, h = 15e3, 3e3          # VW extent, transition (down-dip distances)
a0, amax, b0 = 0.010, 0.025, 0.015
sig0, dc, V0, f0 = 50e6, 0.008, 1e-6, 0.6
Vp = 1e-9
G, cs = 32.04e9, 3464.0
eta = G/(2.0*cs)
n = int(round(Wf/ds))
strike = 90.0 + np.rad2deg(dip)      # tangent down-dip, both branches

fg = open(pref + "_geom.in", "w")
fr = open(pref + "_rsf.dat", "w")
fi = open(pref + "_ic.in", "w")
# prestress and initial velocity carry the sign of the branch: they
# are the steady state of the loading the run will actually apply
tau0 = branch*(sig0*amax*np.arcsinh(Vp/(2*V0)*np.exp((f0 + b0*np.log(V0/Vp))/amax))
               + eta*Vp)
for i in range(n):
    d = (i + 0.5)*ds                       # down-dip center distance
    x, y = d*np.cos(dip), -d*np.sin(dip)
    fg.write(f"{x:.6e} {y:.6e} 0 {strike:.8f} 90 {ds/2:.6e} 0 0\n")
    if d < H:
        a = a0
    elif d < H + h:
        a = a0 + (amax - a0)*(d - H)/h
    else:
        a = amax
    fr.write(f"{a:.6f} {b0:.6f}\n")
    fi.write(f"{tau0:.8e} {branch*Vp:.6e}\n")
# SEAS station file: 12 stations by down-dip distance (spec sec. 6)
with open(pref + "_stations.dat", "w") as fs:
    for xd_km in (0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35):
        idx = min(max(int(round(xd_km*1e3/ds - 0.5)), 0), n - 1)
        fs.write(f"dp{int(round(xd_km*10)):03d} {idx}\n")

print(f"gen_bp3: {n} segments, ds {ds:.0f} m, dip {np.rad2deg(dip):.0f}, "
      f"tau0 {tau0/1e6:.4f} MPa, "
      f"{'THRUST' if branch > 0 else 'NORMAL'} branch")
print(f"opts: -dc {dc} -sigma_init {sig0:.3e} -f0 {f0} -v0 {V0} "
      f"-shear_modulus {G:.4e} -s_wave_speed {cs} -vpl {branch*Vp:.3e}")
