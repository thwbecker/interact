#!/usr/bin/env python3
"""
Generate SEAS BP5-QD (rectangular) input files for interact's rsf_solve.

Mirrors the HBI `bp5r` setup (https://strike.scec.org/cvws/seas/download/SEAS_BP5.pdf).
Vertical strike-slip fault, 100 km (along strike) x 40 km (depth), square cells of
side `ds` km.  A velocity-weakening (VW) rectangle is embedded in a velocity-
strengthening (VS) region, and a high-velocity nucleation patch seeds the first event.

Outputs files in interact's native formats:
  geom_bp5_<ds>km.in           x y z strike dip L W group  (L,W are HALF-lengths [m])
  rsf_bp5_<ds>km.dat            a  b                        (per patch, geometry order)
  ic_bp5_<ds>km.in              tau[Pa]  v[m/s]             (per patch, for -rsf_ic_file)
  ic_bp5_ssvinit_<ds>km.in      tau[Pa]  v[m/s]             (alternative IC, see below)
  dc_bp5_<ds>km.in              per-cell D_c [m]

The two IC files differ only in the initial state implied inside the nucleation
patch (rsf_solve recovers psi by inverting the flow law at the IC's tau and v):

  ic_bp5_...        BP5 spec convention: state at steady state for Vpl everywhere,
                    with the patch forced to slide at Vnuc.  In the patch the state
                    therefore starts b*ln(Vnuc/Vpl) = 17.2 b above steady state
                    (Omega = Vnuc/Vpl = 3e7).  The aging law relaxes such a lag at a
                    bounded rate and nucleates gently over seconds; the slip law
                    relaxes it at rate ~ v/dc, which is also unproblematic.  This
                    file is intended for the aging law (per the BP5 spec) and has
                    been fine for the slip law in the configurations tested here.

  ic_bp5_ssvinit_...  state at steady state for the LOCAL initial velocity, i.e.
                    psi_ss(Vnuc) in the patch, identical to the other file outside
                    it.  Intended for the PRZ law, whose state rate grows like
                    Omega^2 above steady state, so that the spec IC collapses an
                    explicit integrator's step to ~1e-7 s and turns the intended
                    gentle nucleation into an instantaneous strength drop of
                    b*ln(Vnuc/Vpl) = 0.52.  With this IC the patch instead starts
                    creeping at Vnuc on velocity-weakening friction and still
                    seeds the first event.  Note that a PRZ (or slip law) BP5
                    variant is off-spec either way; BP5 defines the aging law only,
                    so the IC is a modeling choice.  Aging and slip runs can also
                    use this file, but results will then differ from the BP5 spec.

Both IC conventions include the radiation-damping term eta*v in tau, as in HBI.
In the patch this leaves the recovered psi about a*eta*Vnuc/(sig*a) = eta*Vnuc/sig
= 5.6e-3 above psi_ss(Vnuc), i.e. a residual Omega of about exp(5.6e-3/b) = 1.2,
which appears harmless for all three laws but is noted here for completeness.

Usage:  python3 make_bp5.py [ds_km]      (default ds = 1.0 km  -> 100x40 = 4000 cells)
"""
import sys, numpy as np

ds = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0    # cell size [km]
Lx, Ld = 100.0, 40.0                                     # strike length, depth [km]
imax, jmax = int(round(Lx/ds)), int(round(Ld/ds))

# --- BP5-QD friction / medium parameters (SEAS / HBI conventions) ---
a0, a_max = 0.004, 0.04   # direct-effect a in VW core / VS exterior
b0  = 0.03                # evolution parameter b (uniform)
dc0 = 0.14                # characteristic slip distance D_c [m] (VW/VS bulk)
dc_nuc = 0.13             # reduced D_RS inside the nucleation patch (SEAS BP5):
                          #   smaller h*, pins re-nucleation to this patch
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
fs = open(f"ic_bp5_ssvinit_{ds:g}km.in", "w")
fd = open(f"dc_bp5_{ds:g}km.in",  "w")
nvw = nnuc = 0
for i in range(imax):
    for j in range(jmax):
        x   = (i + 0.5 - imax/2)*ds      # along-strike coordinate [km], centered
        dep = (j + 0.5)*ds               # depth [km], 0 at the free surface
        # a-distribution: VW core |x|<=30, |dep-10|<=6 ; linear 2 km taper to VS
        r = max(abs(dep-10)-6, abs(x)-30)/2.0
        a = min(a0 + r*(a_max-a0), a_max); a = max(a, a0)
        v = Vpl
        dc = dc0
        if abs(x+24) < 6 and abs(dep-10) < 6:     # 12x12 km nucleation patch
            v = Vnuc; dc = dc_nuc; nnuc += 1
        if a < 0.01: nvw += 1
        # steady-state-at-Vpl initial stress evaluated at the local velocity v
        # (independent of D_RS; D_RS enters only the state evolution rate)
        tau = sig*a*np.arcsinh(0.5*v/V0*np.exp((f0 + b0*np.log(V0/Vpl))/a)) + eta*v
        # alternative: stress such that the recovered state is at steady state
        # for the LOCAL initial velocity (psi_ss(v)); identical to tau outside
        # the nucleation patch, where v = Vpl.  intended for the PRZ law; see
        # the docstring for the rationale and caveats
        tau_ssv = sig*a*np.arcsinh(0.5*v/V0*np.exp((f0 + b0*np.log(V0/v))/a)) + eta*v
        # interact geometry: strike=0 => strike dir +y, dip=90 => normal +x, dip +z(up)
        #   so BP5 along-strike -> interact y, depth -> z=-dep, fault plane at x=0
        fg.write(f"0.0 {x*1e3:.6e} {-dep*1e3:.6e} 0.0 90.0 {half:.1f} {half:.1f} 0\n")
        fr.write(f"{a:.6e} {b0:.6e}\n")
        fi.write(f"{tau*1e6:.8e} {v:.6e}\n")          # MPa -> Pa
        fs.write(f"{tau_ssv*1e6:.8e} {v:.6e}\n")      # MPa -> Pa
        fd.write(f"{dc:.6e}\n")                        # per-cell D_c [m]
fg.close(); fr.close(); fi.close(); fs.close(); fd.close()
print(f"ds={ds:g} km : {imax}x{jmax} = {imax*jmax} cells ; VW={nvw} nucleation={nnuc}")
