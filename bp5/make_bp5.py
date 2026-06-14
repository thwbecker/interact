#!/usr/bin/env python3
"""
Generate SEAS BP5-QD (rectangular) input files for interact's rsf_solve.

Mirrors the HBI `bp5r` setup (https://strike.scec.org/cvws/seas/download/SEAS_BP5.pdf).
Vertical strike-slip fault, 100 km (along strike) x 40 km (depth), square cells of
side `ds` km.  A velocity-weakening (VW) rectangle is embedded in a velocity-
strengthening (VS) region, and a high-velocity nucleation patch seeds the first event.

Outputs three files in interact's native formats:
  geom_bp5_<ds>km.in   x y z strike dip L W group   (L,W are HALF-lengths [m])
  rsf_bp5_<ds>km.dat    a  b                         (per patch, geometry order)
  ic_bp5_<ds>km.in      tau[Pa]  v[m/s]              (per patch, for -rsf_ic_file)

Usage:  python3 make_bp5.py [ds_km]      (default ds = 1.0 km  -> 100x40 = 4000 cells)
"""
import sys, numpy as np

ds = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0    # cell size [km]
Lx, Ld = 100.0, 40.0                                     # strike length, depth [km]
imax, jmax = int(round(Lx/ds)), int(round(Ld/ds))

# --- BP5-QD friction / medium parameters (SEAS / HBI conventions) ---
a0, a_max = 0.004, 0.04   # direct-effect a in VW core / VS exterior
b0  = 0.03                # evolution parameter b (uniform)
dc0 = 0.14                # characteristic slip distance D_c [m] (uniform here)
f0  = 0.6                 # reference friction
V0  = 1e-6                # reference velocity [m/s]   (HBI: vref)
Vpl = 1e-9                # plate rate = initial loading velocity [m/s]
sig = 25.0                # effective normal stress [MPa]
G   = 32.04               # shear modulus [GPa]
cs  = 3.464               # shear-wave speed [km/s]
eta = G/(2.0*cs)          # radiation-damping factor (HBI numeric convention)
Vnuc = 3e-2               # nucleation-patch initial velocity [m/s]

half = ds*1e3/2.0         # interact half-length [m]  (L=W -> ds-km square cell)

fg = open(f"geom_bp5_{ds:g}km.in", "w")
fr = open(f"rsf_bp5_{ds:g}km.dat", "w")
fi = open(f"ic_bp5_{ds:g}km.in",  "w")
nvw = nnuc = 0
for i in range(imax):
    for j in range(jmax):
        x   = (i + 0.5 - imax/2)*ds      # along-strike coordinate [km], centered
        dep = (j + 0.5)*ds               # depth [km], 0 at the free surface
        # a-distribution: VW core |x|<=30, |dep-10|<=6 ; linear 2 km taper to VS
        r = max(abs(dep-10)-6, abs(x)-30)/2.0
        a = min(a0 + r*(a_max-a0), a_max); a = max(a, a0)
        v = Vpl
        if abs(x+24) < 6 and abs(dep-10) < 6:     # 12x12 km nucleation patch
            v = Vnuc; nnuc += 1
        if a < 0.01: nvw += 1
        # steady-state-at-Vpl initial stress evaluated at the local velocity v
        tau = sig*(a*np.arcsinh(0.5*v/V0*np.exp((f0 + b0*np.log(V0/Vpl))/a)) + eta*v)
        # interact geometry: strike=0 => strike dir +y, dip=90 => normal +x, dip +z(up)
        #   so BP5 along-strike -> interact y, depth -> z=-dep, fault plane at x=0
        fg.write(f"0.0 {x*1e3:.6e} {-dep*1e3:.6e} 0.0 90.0 {half:.1f} {half:.1f} 0\n")
        fr.write(f"{a:.6e} {b0:.6e}\n")
        fi.write(f"{tau*1e6:.8e} {v:.6e}\n")          # MPa -> Pa
fg.close(); fr.close(); fi.close()
print(f"ds={ds:g} km : {imax}x{jmax} = {imax*jmax} cells ; VW={nvw} nucleation={nnuc}")
