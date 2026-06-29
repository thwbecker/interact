#!/usr/bin/env python3
"""
Generate a TWO-fault test case for interact's H-matrix / b = A x benchmarking
(compress_interaction_matrix, sweep_hmat_bAx_2fault.sh) and, secondarily, for
rsf_solve. Two identical vertical strike-slip faults are laid down as parallel
planes with a controlled geometry, so the resulting interaction operator has two
spatially separated clusters and a coherent cross-fault far-field.

Controlled geometry (all inter-fault lengths in units of the fault width W):
  - each fault has a 2:1 aspect ratio, L = 2 W (enforced by imax = 2*jmax square
    cells)
  - the two faults are parallel (same strike) and separated by 0.5 W in the
    fault-normal direction
  - the second fault is offset along strike by 0.05 W relative to the first
The absolute scale (W in km) is arbitrary for the compression test: the H-matrix
structure depends on the geometric proportions, not the absolute size, so W only
sets the kernel magnitude, not the relative accuracy/compression.

Resolution is set by jmax (down-dip cells); total patches = 2*(2*jmax)*jmax =
4*jmax^2, chosen to match the single-fault sweeps:
  jmax = 63  -> 126x63  per fault, 15876 patches total (~ single-fault 0.5km, N=16000)
  jmax = 126 -> 252x126 per fault, 63504 patches total (~ single-fault 0.25km, N=64000)

The per-fault friction/medium/initial-condition recipe is the BP5-QD one from
bp5/make_bp5.py, scaled to this fault, and the nucleation patches sit on opposite
along-strike ends so an rsf_solve cycle would couple asymmetrically. Those
fields are unused by the b = A x test (compress reads only the geometry).

Outputs (interact native formats), tagged by jmax:
  geom_2f_j<jmax>.in   x y z strike dip L W group  (L,W HALF-lengths [m]; group 0/1)
  rsf_2f_j<jmax>.dat    a  b
  ic_2f_j<jmax>.in      tau[Pa]  v[m/s]
  dc_2f_j<jmax>.in      D_c[m]

Usage:
  python3 make_two_fault.py [jmax] [W_km] [sep_frac] [offset_frac]
Defaults: jmax=63  W_km=40  sep_frac=0.5  offset_frac=0.05
"""
import sys, numpy as np

jmax    = int(sys.argv[1])   if len(sys.argv) > 1 else 63     # down-dip cells (resolution)
Wkm     = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0   # fault width W [km] (absolute scale)
sep_fr  = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5    # fault-normal separation / W
off_fr  = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05   # along-strike offset / W

imax = 2*jmax                 # 2:1 aspect ratio (L = 2 W), square cells
ds   = Wkm/jmax               # cell size [km]
Lkm  = imax*ds                # along-strike length [km] = 2 W
sep  = sep_fr*Wkm             # fault-normal separation [km]
off  = off_fr*Wkm             # along-strike offset [km]

# --- BP5-QD friction / medium parameters (unchanged from bp5/make_bp5.py) ---
a0, a_max = 0.004, 0.04
b0  = 0.03
dc0 = 0.14
dc_nuc = 0.13
f0  = 0.6
V0  = 1e-6
Vpl = 1e-9
sig = 25.0
G   = 32.04
cs  = 3.464
eta = G/(2.0*cs)
Vnuc = 3e-2

# VW core / nucleation geometry, scaled to this fault (caps from the 100x40 BP5)
vw_hx  = min(30.0, 0.30*Lkm)
vw_dc  = min(10.0, 0.25*Wkm)
vw_hd  = min(6.0,  0.30*Wkm)
taper  = 2.0
nuc_h  = min(6.0, 0.5*vw_hx)
nuc_off= min(24.0, 0.8*vw_hx)

half = ds*1e3/2.0             # interact half-length [m] (square cell)

tag = f"j{jmax}"
fg = open(f"geom_2f_{tag}.in", "w")
fr = open(f"rsf_2f_{tag}.dat", "w")
fi = open(f"ic_2f_{tag}.in",  "w")
fd = open(f"dc_2f_{tag}.in",  "w")

nvw = nnuc = npatch = 0
for fault in (0, 1):
    xn    = (fault - 0.5)*sep                  # fault-normal offset [km]: -sep/2, +sep/2
    xoff  = fault*off                          # along-strike offset [km]: 0, off
    nuc_x = -nuc_off if fault == 0 else +nuc_off
    for i in range(imax):
        for j in range(jmax):
            xc  = (i + 0.5 - imax/2)*ds         # along-strike coord, centered on this fault
            dep = (j + 0.5)*ds                  # depth [km], 0 at free surface
            x   = xc + xoff                     # absolute along-strike position [km]
            r = max(abs(dep-vw_dc)-vw_hd, abs(xc)-vw_hx)/taper
            a = min(a0 + r*(a_max-a0), a_max); a = max(a, a0)
            v = Vpl; dc = dc0
            if abs(xc-nuc_x) < nuc_h and abs(dep-vw_dc) < nuc_h:
                v = Vnuc; dc = dc_nuc; nnuc += 1
            if a < 0.01: nvw += 1
            tau = sig*(a*np.arcsinh(0.5*v/V0*np.exp((f0 + b0*np.log(V0/Vpl))/a)) + eta*v)
            # strike=0 => strike dir +y, dip=90 => normal +x; plane at x=xn, along-strike->y, depth->z
            # NOTE: write the half-lengths at full precision, not %.1f: for non-integer
            # ds the shallowest patch center and half-width must cancel to z_top<=0, and
            # rounding the half up (e.g. 666.6667->666.7 vs z=-666.667) breaches the free
            # surface and read_geometry rejects the whole geometry.
            fg.write(f"{xn*1e3:.6e} {x*1e3:.6e} {-dep*1e3:.6e} 0.0 90.0 "
                     f"{half:.6e} {half:.6e} {fault}\n")
            fr.write(f"{a:.6e} {b0:.6e}\n")
            fi.write(f"{tau*1e6:.8e} {v:.6e}\n")
            fd.write(f"{dc:.6e}\n")
            npatch += 1

fg.close(); fr.close(); fi.close(); fd.close()
print(f"two-fault: jmax={jmax} W={Wkm:g}km L={Lkm:g}km (2:1)  sep={sep:g}km (={sep_fr:g}W) "
      f"offset={off:g}km (={off_fr:g}W)  {imax}x{jmax} per fault, {npatch} patches total ; "
      f"VW={nvw} nucleation={nnuc}")
print(f"wrote geom_2f_{tag}.in rsf_2f_{tag}.dat ic_2f_{tag}.in dc_2f_{tag}.in")
