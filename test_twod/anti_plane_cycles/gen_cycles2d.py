#!/usr/bin/env python3
"""
gen_cycles2d.py: generate single- or multi-fault 2-D antiplane
earthquake-cycle setups for rsf_solve (vertical shear zones, depth-only
property variation, out-of-plane slip), for the experiment suite in this
directory (driver: run_cycles2d; docs: README_cycles2d.md).

Two friction profiles:

  bp1   the SEAS benchmark BP1-QD profile (Erickson & Jiang 2018;
        Erickson et al., SRL 2020): a = 0.010 above 15 km, linear to
        amax = 0.025 over 15-18 km, 0.025 to the fault bottom;
        b = 0.015; sigma = 50 MPa; Dc = 8 mm; V0 = 1e-6; f0 = 0.6;
        Vp = 1e-9 m/s; G = 32.04 GPa; cs = 3464 m/s; BP1's exact
        uniform initial condition (prestress at the amax steady state
        of Vinit = 1e-9, so the velocity-weakening region starts above
        steady state and the sequence spins up on its own).  With
        D_km = 40 this IS BP1-QD (in the backslip formulation the
        steady creep below 40 km contributes zero stress, so the
        40-km rate-and-state fault is the entire model).  Community
        sanity values: events about every 100 yr (78 yr between later
        events at fine resolution), ~3 m surface slip per event,
        nucleation near 12-13 km, arrest near 18 km.

  bp5   the BP5-flavored profile of gen_bp2d.py (VS cap to 4 km,
        2-km transitions, VW core, VS below; a 0.004/0.04, b 0.03,
        sigma 25 MPa, Dc 0.14 m default -> pass -dc; G = 32.04 GPa),
        with a small nucleation seed on fault 0 so the first event
        starts immediately.  Robust workhorse for tM/Trec and D/H
        scans (elastic recurrence ~127 yr at D=20, seed excluded).

Multi-fault: nfault parallel vertical faults spaced spacing_km apart;
each fault is its own GROUP, so with -rsf_monitor_by_group 1 and
-rsf_catalog you get per-fault catalogs rsf_catalog.gNNN.dat and
monitors, and the field frames tmp_rsf/rsf_{vel,tau}.gNNN.*.bin are per
fault as well.

usage: gen_cycles2d.py profile ds_km D_km nfault spacing_km [prefix]
       e.g.  gen_cycles2d.py bp1 0.1  40   1      10        cyc
             gen_cycles2d.py bp5 0.5  20   8      10        cyc
writes <prefix>_{geom.in,rsf.dat,ic.in} and prints the rsf_solve
options that must accompany them (-dc, -sigma_init or sigma file, -v0,
-f0, -vpl, -shear_modulus, -s_wave_speed).
"""
import sys
import numpy as np

if len(sys.argv) < 6:
    sys.stderr.write(__doc__)
    sys.exit(1)
profile = sys.argv[1]
ds      = float(sys.argv[2])   # cell size [km]
D       = float(sys.argv[3])   # fault depth extent [km]
nfault  = int(sys.argv[4])
spacing = float(sys.argv[5])   # fault spacing [km]
prefix  = sys.argv[6] if len(sys.argv) > 6 else "cyc"
assert profile in ("bp1", "bp5")

if profile == "bp1":
    G, cs = 32.04e9, 3464.0
    sig0, f0, V0, Vp, dc = 50.0e6, 0.6, 1e-6, 1e-9, 0.008
    a0, amax, b0 = 0.010, 0.025, 0.015
    H_vw, h_tr = 15.0, 3.0
    Vinit = 1e-9
else:
    G, cs = 32.04e9, 3464.0
    sig0, f0, V0, Vp, dc = 25.0e6, 0.6, 1e-6, 1e-9, 0.14
    a0, amax, b0 = 0.004, 0.04, 0.03
    h_s, h_t = 4.0, 2.0            # VS cap, transition widths [km]
eta = G/(2.0*cs)

nz = int(round(D/ds))
half = ds*1e3/2.0
fg = open(f"{prefix}_geom.in", "w")
fr = open(f"{prefix}_rsf.dat", "w")
fi = open(f"{prefix}_ic.in", "w")

if profile == "bp1":
    # BP1 exact IC: uniform prestress at the amax steady state of Vinit
    tau0 = sig0*amax*np.arcsinh(Vinit/(2*V0)*np.exp((f0 + b0*np.log(V0/Vinit))/amax)) + eta*Vinit
for kf in range(nfault):
    x = kf*spacing*1e3
    for j in range(nz):
        z = (j+0.5)*ds
        if profile == "bp1":
            if z < H_vw:            a = a0
            elif z < H_vw + h_tr:   a = a0 + (amax-a0)*(z-H_vw)/h_tr
            else:                   a = amax
            b = b0
            v, tau = Vinit, tau0
        else:
            # BP5-flavor: VS cap [0,h_s], taper, VW core to D-h_t-h_s?, use
            # gen_bp2d convention: VW core from h_s+h_t to D-h_t-2 (VS toe)
            zvw1, zvw2 = h_s + h_t, D - h_t - 2.0
            if z < h_s:               a = amax
            elif z < zvw1:            a = amax + (a0-amax)*(z-h_s)/h_t
            elif z < zvw2:            a = a0
            elif z < zvw2 + h_t:      a = a0 + (amax-a0)*(z-zvw2)/h_t
            else:                     a = amax
            b = b0
            # nucleation seed on fault 0 near the VW center
            zc = 0.5*(zvw1+zvw2)
            v = 3e-2 if (kf == 0 and abs(z-zc) <= 1.0) else Vp
            tau = sig0*(f0 + (a-b)*np.log(Vp/V0)) + eta*v
        fg.write(f"{x:.6e} {-z*1e3:.6e} 0 0 -90 {half:.6e} 0 {kf}\n")
        fr.write(f"{a:.6f} {b:.6f}\n")
        fi.write(f"{tau:.8e} {v:.6e}\n")
fg.close(); fr.close(); fi.close()

Lb = G*dc/(b0*sig0)
hstar = (2.0/np.pi)*G*dc/((b0-a0)*sig0)
print(f"gen_cycles2d: {profile}, {nfault} fault(s) x {nz} cells, ds {ds} km, D {D} km, "
      f"spacing {spacing} km")
print(f"  resolution: Lb = {Lb:.0f} m (need ds <~ Lb/3 = {Lb/3:.0f} m), "
      f"h* = {hstar:.0f} m; cell = {ds*1e3:.0f} m")
print(f"  pass: -dc {dc} -sigma_init {sig0:.3e} -f0 {f0} -v0 {V0} -vpl {Vp} "
      f"-shear_modulus {G:.4e} -s_wave_speed {cs}")
