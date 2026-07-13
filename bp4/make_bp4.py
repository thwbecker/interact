#!/usr/bin/env python3
"""
Generate SEAS BP4-QD (rectangular, whole-space) input files for interact's
rsf_solve, modeled on make_bp5.py.

BP4 is the whole-space sibling of BP5 (Erickson and Jiang, SEAS BP4, 2019).
A vertical strike-slip fault is embedded in a homogeneous whole space (run with
-full_space), symmetric about x3 = 0.  The rate-and-state region is |x3| <= Wf
and is meshed here; outside it creeps at the plate rate.  A velocity-weakening
rectangle (|x3| <= H, |x2| <= l/2) sits at the center, surrounded by a transition
strip of width h to velocity strengthening.  The first event is seeded in a square
favorable zone at one corner of the VW patch.  Two seed options are provided (see
nuc_seed below): a higher initial slip rate Vi with the background state unchanged,
which matches Jiang et al. 2022 and gives immediate initiation (default), or the
BP4 description PDF variant of a 1.1 tau0 prestress bump at Vinit, which nucleates
slowly near 200 yr.  The recurrence interval is the same either way.

BP4 differs from BP5 in nearly every number (see Table 1 of each benchmark), so
the parameters below are BP4's, not BP5's.

Outputs (interact native formats, geometry order):
  geom_bp4_<ds>km_L<Lstrike>.in   x y z strike dip L W group   (L,W HALF-lengths [m])
  rsf_bp4_<ds>km_L<Lstrike>.dat    a  b                          (per patch)
  ic_bp4_<ds>km_L<Lstrike>.in      tau[Pa]  v[m/s]               (per patch, for -rsf_ic_file)
D_c is uniform (L = 0.04 m), so no per-cell dc file is written; pass -dc 0.04.

Usage:  python3 make_bp4.py [ds_km] [Lstrike_km]
  ds_km      cell size (default 0.5 km, the suggested BP4 resolution)
  Lstrike_km along-strike extent of the meshed domain (default 100 km).  BP4's
             frictional domain is infinite along strike, so this is a numerical
             truncation; the benchmark asks for two domain sizes to check that
             results are domain-size independent.
"""
import sys, numpy as np

ds      = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5     # cell size [km]
Lstrike = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0   # along-strike domain [km]

# With the full-space fixes to eval_okada and dc3d.F, the Okada kernel accepts
# patches at any z, so the fault is placed at its true position, symmetric about
# x3 = 0 (z_off = 0), and the BP4 down-dip coordinate x3 is recovered directly as
# x3 = -z with no offset. z_off is kept only as a fallback: on a kernel that still
# rejects z > 0, set z_off > Wf to push the whole fault below z = 0 (a whole-space
# shift changes nothing physically), and then x3 = depth - z_off.
z_off = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0     # fault-center depth [km]

# --- BP4-QD Table 1 parameters (whole space) ---
a0, a_max = 0.0065, 0.025   # a in VW core / VS exterior
b0   = 0.013                # b (uniform)
dc   = 0.04                 # L = D_c [m] (uniform, no reduced nucleation value)
f0   = 0.6                  # reference friction
V0   = 1e-6                 # reference velocity [m/s]
Vpl  = 1e-9                 # plate rate
Vinit= 1e-9                 # initial slip rate (uniform, including nucleation)
sig  = 50.0                 # effective normal stress [MPa]
G    = 32.04                # shear modulus [GPa]
cs   = 3.464                # shear-wave speed [km/s]
eta  = G/(2.0*cs)           # radiation-damping factor; eta*V is in MPa here

# geometry [km]
H   = 15.0                  # VW half-height (VW spans |x3| <= H = 30 km tall)
lVW = 60.0                  # VW length (|x2| <= lVW/2)
h   = 3.0                   # VW-VS transition width
Wf  = 40.0                  # fault half-height (RSF region |x3| <= Wf = 80 km tall)
wi  = 12.0                  # favorable nucleation square width
hi  = 1.5                   # offset of nucleation square from the VW corner
# nucleation seed. "velocity" matches Jiang et al. 2022 and the other benchmark
# participants: a higher initial slip rate Vi is imposed in the square with the
# background state left unchanged, giving immediate initiation of the first event
# so simulated event times line up with the published runs. "prestress" is the
# BP4 description PDF variant (a 1.1 tau0 bump at Vinit), which nucleates slowly,
# near 200 yr. The recurrence interval is the same either way; only the timing of
# the first event differs.
nuc_seed = "velocity"
Vi       = 0.03            # nucleation slip rate [m/s] (velocity seed)
tau_bump = 1.1             # nucleation prestress factor (prestress seed)

# nucleation square, at one VW corner offset hi inward; centered here at
# (x2, x3) = (-lVW/2 + hi + wi/2, -H + hi + wi/2) = (-22.5, -7.5), matching the
# BP4 fault station strk-225dp-750.  By whole-space symmetry the mirror corner
# (+x3) is equivalent.
nuc_x2 = -lVW/2 + hi + wi/2.0     # -22.5
nuc_x3 = -H     + hi + wi/2.0     # -7.5

imax = int(round(Lstrike/ds))
jmax = int(round(2.0*Wf/ds))      # x3 from -Wf to +Wf
half = ds*1e3/2.0

# uniform scalar prestress tau0 (BP4 eq 18, using amax), in MPa
tau0 = sig*a_max*np.arcsinh(0.5*Vinit/V0*np.exp((f0 + b0*np.log(V0/Vinit))/a_max)) + eta*Vinit

fg = open(f"geom_bp4_{ds:g}km_L{Lstrike:g}.in", "w")
fr = open(f"rsf_bp4_{ds:g}km_L{Lstrike:g}.dat", "w")
fi = open(f"ic_bp4_{ds:g}km_L{Lstrike:g}.in",  "w")
nvw = nnuc = 0
for i in range(imax):
    for j in range(jmax):
        x2 = (i + 0.5 - imax/2.0)*ds      # along strike [km], centered on 0
        x3 = (j + 0.5 - jmax/2.0)*ds      # down dip [km], centered on 0 (-Wf..+Wf)
        # a-distribution (BP4 eq 13): VW at |x3|<=H and |x2|<=l/2, VS beyond a
        # transition of width h; symmetric in x3 about 0
        r = max(abs(x3) - H, abs(x2) - lVW/2.0)/h
        a = min(a0 + max(0.0, r)*(a_max - a0), a_max)
        if a < 0.5*(a0 + b0): nvw += 1
        tau = tau0; v = Vinit
        if abs(x2 - nuc_x2) < wi/2.0 and abs(x3 - nuc_x3) < wi/2.0:
            if nuc_seed == "velocity":
                # higher slip rate Vi with the background state left unchanged; the
                # consistent prestress follows from the friction law (BP5-style seed)
                v   = Vi
                tau = sig*a*np.arcsinh((Vi/Vinit)*np.sinh((tau0 - eta*Vinit)/(sig*a))) + eta*Vi
            else:
                tau = tau_bump*tau0
            nnuc += 1
        # interact geometry: x1=0 fault plane, x2->y, x3->z=-(x3+z_off) so the
        # whole symmetric fault sits below z=0 (see z_off note above); strike=0,
        # dip=90 (vertical)
        fg.write(f"0.0 {x2*1e3:.6e} {-(x3+z_off)*1e3:.6e} 0.0 90.0 {half:.1f} {half:.1f} 0\n")
        fr.write(f"{a:.6e} {b0:.6e}\n")
        fi.write(f"{tau*1e6:.8e} {v:.6e}\n")     # MPa -> Pa
fg.close(); fr.close(); fi.close()
print(f"BP4 ds={ds:g} km, Lstrike={Lstrike:g} km : {imax}x{jmax} = {imax*jmax} cells")
seed_desc = f"velocity seed Vi={Vi:g} m/s (immediate)" if nuc_seed == "velocity" else f"prestress bump {tau_bump:g}x (slow)"
print(f"  VW cells={nvw} nucleation cells={nnuc} ; uniform tau0={tau0:.4f} MPa ; {seed_desc}")
print(f"  x3 in [-{Wf:g},{Wf:g}] km (symmetric) ; nucleation center ({nuc_x2:g},{nuc_x3:g}) km")
if z_off == 0.0:
    print(f"  fault at true position, symmetric about z=0 (needs the full-space z fix in eval_okada/dc3d.F)")
else:
    print(f"  fault shifted to depth center z_off={z_off:g} km (fallback for kernels rejecting z>0; BP4 x3 = depth - {z_off:g})")
print(f"  run with -full_space -dc {dc:g} -sigma_init {sig*1e6:g} -shear_modulus {G*1e9:g} -s_wave_speed {cs*1e3:g}")
