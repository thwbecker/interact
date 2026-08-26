#!/usr/bin/env python3
"""
gen_bp2d.py: generate a 2-D antiplane BP5-style depth-dependent
rate-and-state setup for rsf_solve (out-of-plane slip on one or more
vertical shear zones, properties depending on depth only).

Writes geom_bp2d.in (x y z strike dip L 0 group; antiplane elements,
y = -depth [m], L = half-extent [m]), rsf_bp2d.dat (a b per patch),
and ic_bp2d.in (tau vel per patch, HBI/BP5-style steady-state tau at
v_init including radiation damping).

Depth profile (BP5 flavor, SEAS BP5 / HBI conventions): VS from the
surface to h_s, linear VW transition over h_t, VW core, transition
back, VS to the bottom; a varies between a0 (VW) and amax (VS), b
uniform.

usage: gen_bp2d.py [ds_km] [nfault] [fault_dx_km] [depth_km] [H_vw_km]
       defaults    0.5      1        20            40         12
(remember to pass a matching -dc to rsf_solve; the discretization
warnings there check ds against Lb = G dc/(b sigma))
"""
import sys
import numpy as np

ds       = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5   # cell [km]
nfault   = int(sys.argv[2])   if len(sys.argv) > 2 else 1
fault_dx = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0  # spacing [km]
depth    = float(sys.argv[4]) if len(sys.argv) > 4 else 40.0  # fault depth extent [km]
H_vw_in  = float(sys.argv[5]) if len(sys.argv) > 5 else 12.0  # VW core extent [km]

# BP5-QD friction/medium parameters (SEAS BP5 / HBI)
a0, amax = 0.004, 0.04
b0   = 0.03
f0   = 0.6
V0   = 1e-6
Vpl  = 1e-9
sig  = 25.0e6          # [Pa]
G    = 32.04e9         # [Pa]
cs   = 3464.0          # [m/s]
eta  = G/(2.0*cs)
# depth structure [km]: VS cap, transition, VW core, transition, VS below
h_s, h_t = 4.0, 2.0
H_vw = H_vw_in                      # VW core from h_s+h_t to h_s+h_t+H_vw
# nucleation patch (BP5 style): elevated initial velocity in a small
# depth window of the VW core on fault 0, so the first event starts
# immediately instead of waiting decades for the linear instability
# of the (otherwise exactly steady) initial state to express itself
zn_c = h_s + h_t + 0.5*H_vw
zn1, zn2, Vnuc = zn_c - 1.0, zn_c + 1.0, 3e-2

nz = int(round(depth/ds))
half = ds*1e3/2.0
fg = open("geom_bp2d.in","w")
fr = open("rsf_bp2d.dat","w")
fi = open("ic_bp2d.in","w")
for kf in range(nfault):
    x = kf*fault_dx*1e3
    for j in range(nz):
        zmid = (j+0.5)*ds        # depth [km]
        # BP5-style a(z): a0 in the VW core, ramping to amax in the VS regions
        if zmid < h_s:
            a = amax
        elif zmid < h_s+h_t:
            a = amax + (a0-amax)*(zmid-h_s)/h_t
        elif zmid < h_s+h_t+H_vw:
            a = a0
        elif zmid < h_s+2*h_t+H_vw:
            a = a0 + (amax-a0)*(zmid-(h_s+h_t+H_vw))/h_t
        else:
            a = amax
        fg.write(f"{x:.6e} {-zmid*1e3:.6e} 0 0 -90 {half:.6e} 0 {kf}\n")
        fr.write(f"{a:.6f} {b0:.6f}\n")
        if (kf == 0) and (zn1 <= zmid <= zn2):
            v = Vnuc            # nucleation patch, BP5-style IC
        else:
            v = Vpl
        tau = sig*(f0 + (a-b0)*np.log(v/V0)) + eta*v
        fi.write(f"{tau:.8e} {v:.6e}\n")
fg.close(); fr.close(); fi.close()
print(f"gen_bp2d: {nfault} fault(s), {nz} cells each, ds {ds} km, "
      f"VW {h_s+h_t}-{h_s+h_t+H_vw} km, Lb {G*0.14/(b0*sig)/1e3:.2f} km")
