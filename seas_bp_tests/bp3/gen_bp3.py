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
tangent is hanging-wall-DOWN, i.e. NORMAL-sense; run with
-vpl -1e-9 for the thrust branch (or +1e-9 for normal) and check
the surface-displacement sense.

usage: gen_bp3.py ds_km dip_deg [prefix] [sense]
       sense +1 (default) or -1: flips the segment tangent
       (positive slip direction), switching thrust vs normal
writes <prefix>_{geom.in,rsf.dat,ic.in}
"""
import sys
import numpy as np

ds = float(sys.argv[1])*1e3 if len(sys.argv) > 1 else 200.0
dip = np.deg2rad(float(sys.argv[2])) if len(sys.argv) > 2 else np.deg2rad(60.0)
pref = sys.argv[3] if len(sys.argv) > 3 else "bp3"
sense = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0

Wf = 40e3                 # frictional down-dip extent
H, h = 15e3, 3e3          # VW extent, transition (down-dip distances)
a0, amax, b0 = 0.010, 0.025, 0.015
sig0, dc, V0, f0 = 50e6, 0.008, 1e-6, 0.6
Vp = 1e-9
G, cs = 32.04e9, 3464.0
eta = G/(2.0*cs)
n = int(round(Wf/ds))
strike = 90.0 + np.rad2deg(dip) + (180.0 if sense < 0 else 0.0)

fg = open(pref + "_geom.in", "w")
fr = open(pref + "_rsf.dat", "w")
fi = open(pref + "_ic.in", "w")
tau0 = sig0*amax*np.arcsinh(Vp/(2*V0)*np.exp((f0 + b0*np.log(V0/Vp))/amax)) + eta*Vp
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
    fi.write(f"{tau0:.8e} {Vp:.6e}\n")
print(f"gen_bp3: {n} segments, ds {ds:.0f} m, dip {np.rad2deg(dip):.0f}, "
      f"tau0 {tau0/1e6:.4f} MPa")
print(f"opts: -dc {dc} -sigma_init {sig0:.3e} -f0 {f0} -v0 {V0} "
      f"-shear_modulus {G:.4e} -s_wave_speed {cs}")
