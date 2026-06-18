#!/usr/bin/env python3
"""
Generate a TWO-fault test case for interact's rsf_solve, built from two
BP5-QD-style vertical strike-slip faults placed next to each other.

The single-fault BP5-QD setup (mirrored from bp5/make_bp5.py, see
https://strike.scec.org/cvws/seas/download/SEAS_BP5.pdf) is a 100 km (strike)
x 40 km (depth) vertical strike-slip fault with a velocity-weakening (VW)
rectangle embedded in a velocity-strengthening (VS) region and a high-velocity
nucleation patch.  Here we lay down two such faults as parallel planes
separated by `sep` km in the fault-normal direction, and we put the nucleation
patch on opposite along-strike ends of the two faults so that the two systems
couple through elastic stress transfer rather than firing in lockstep.  This
gives an interacting two-fault sequence whose event timing is sensitive to the
elastic operator, which is what we want when probing H-matrix accuracy.

The friction, medium, and initial-condition recipe is copied unchanged from
bp5/make_bp5.py so the per-fault physics is the familiar BP5 one; only the
geometry (two offset planes) and the nucleation side differ.

Outputs (interact native formats), tagged by resolution and separation:
  geom_2f_<ds>km_<sep>km.in   x y z strike dip L W group  (L,W HALF-lengths [m];
                                                           group 0 / 1 per fault)
  rsf_2f_<ds>km_<sep>km.dat    a  b                        (per patch)
  ic_2f_<ds>km_<sep>km.in      tau[Pa]  v[m/s]             (per patch)
  dc_2f_<ds>km_<sep>km.in      D_c[m]                      (per patch)

Usage:
  python3 make_two_fault.py [ds_km] [sep_km] [Lx_km] [Ld_km]
Defaults: ds=2.0  sep=4.0  Lx=100  Ld=40   (two 50x20 = 1000-cell faults,
2000 patches total: a moderate sampling that resolves the BP5 process zone at
roughly the 1 to 2 km level and keeps a dense reference run tractable on a
multi-core node).  Use ds=1.0 for a finer production sampling, or a smaller
Lx/Ld and larger ds for a quick plumbing check.

Note: the cross-fault coupling strength depends on `sep` relative to the VW
patch size; sep of a few km gives strong interaction for these 12 km-tall VW
cores.  The value is a modeling choice, not a fixed BP5 quantity, so treat the
resulting sequence as illustrative of an interacting pair rather than as a
standardized benchmark.
"""
import sys, numpy as np

ds  = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0     # cell size [km]
sep = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0     # fault-normal separation [km]
Lx  = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0   # along-strike length [km]
Ld  = float(sys.argv[4]) if len(sys.argv) > 4 else 40.0    # depth extent [km]
imax, jmax = int(round(Lx/ds)), int(round(Ld/ds))

# --- BP5-QD friction / medium parameters (unchanged from bp5/make_bp5.py) ---
a0, a_max = 0.004, 0.04   # direct-effect a in VW core / VS exterior
b0  = 0.03                # evolution parameter b (uniform)
dc0 = 0.14                # characteristic slip distance D_c [m] (VW/VS bulk)
dc_nuc = 0.13             # reduced D_RS inside the nucleation patch
f0  = 0.6                 # reference friction
V0  = 1e-6                # reference velocity [m/s]
Vpl = 1e-9                # plate rate = initial loading velocity [m/s]
sig = 25.0                # effective normal stress [MPa]
G   = 32.04               # shear modulus [GPa]
cs  = 3.464               # shear-wave speed [km/s]
eta = G/(2.0*cs)          # radiation-damping factor
Vnuc = 3e-2               # nucleation-patch initial velocity [m/s]

# VW core and nucleation geometry, scaled so they still fit when Lx/Ld shrink.
# BP5 uses a 60x12 km VW core (|x|<=30, |dep-10|<=6) on a 100x40 fault and a
# 12x12 km nucleation patch.  We keep those absolute sizes when the fault is
# large enough, and shrink them proportionally for small test faults.
vw_hx  = min(30.0, 0.30*Lx)        # VW half-extent along strike [km]
vw_dc  = min(10.0, 0.25*Ld)        # VW center depth [km]
vw_hd  = min(6.0,  0.30*Ld)        # VW half-extent in depth [km]
taper  = 2.0                       # linear a-taper width [km]
nuc_h  = min(6.0, 0.5*vw_hx)       # nucleation half-size [km]
nuc_off= min(24.0, 0.8*vw_hx)      # nucleation offset from center along strike [km]

half = ds*1e3/2.0          # interact half-length [m] (L=W -> ds-km square cell)

tag = f"{ds:g}km_{sep:g}km"
fg = open(f"geom_2f_{tag}.in", "w")
fr = open(f"rsf_2f_{tag}.dat", "w")
fi = open(f"ic_2f_{tag}.in",  "w")
fd = open(f"dc_2f_{tag}.in",  "w")

nvw = nnuc = npatch = 0
for fault in (0, 1):
    xn = (fault - 0.5) * sep                 # fault-normal offset [km]: -sep/2, +sep/2
    nuc_x = -nuc_off if fault == 0 else +nuc_off   # opposite ends -> asymmetric coupling
    for i in range(imax):
        for j in range(jmax):
            x   = (i + 0.5 - imax/2)*ds      # along-strike coordinate [km], centered
            dep = (j + 0.5)*ds               # depth [km], 0 at the free surface
            # a-distribution: VW core with a linear taper to the VS exterior
            r = max(abs(dep-vw_dc)-vw_hd, abs(x)-vw_hx)/taper
            a = min(a0 + r*(a_max-a0), a_max); a = max(a, a0)
            v = Vpl
            dc = dc0
            if abs(x-nuc_x) < nuc_h and abs(dep-vw_dc) < nuc_h:
                v = Vnuc; dc = dc_nuc; nnuc += 1
            if a < 0.01: nvw += 1
            # steady-state-at-Vpl initial shear stress evaluated at local v
            tau = sig*(a*np.arcsinh(0.5*v/V0*np.exp((f0 + b0*np.log(V0/Vpl))/a)) + eta*v)
            # interact geometry: strike=0 => strike dir +y, dip=90 => normal +x.
            # Fault plane at x = xn (km), along-strike -> y, depth -> z=-dep.
            fg.write(f"{xn*1e3:.6e} {x*1e3:.6e} {-dep*1e3:.6e} 0.0 90.0 "
                     f"{half:.1f} {half:.1f} {fault}\n")
            fr.write(f"{a:.6e} {b0:.6e}\n")
            fi.write(f"{tau*1e6:.8e} {v:.6e}\n")          # MPa -> Pa
            fd.write(f"{dc:.6e}\n")                        # per-cell D_c [m]
            npatch += 1

fg.close(); fr.close(); fi.close(); fd.close()
print(f"two-fault: ds={ds:g} km sep={sep:g} km  {imax}x{jmax} per fault, "
      f"{npatch} patches total ; VW={nvw} nucleation={nnuc}")
print(f"wrote geom_2f_{tag}.in rsf_2f_{tag}.dat ic_2f_{tag}.in dc_2f_{tag}.in")
