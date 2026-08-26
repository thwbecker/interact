#!/usr/bin/env python3
"""
gen_mn19.py: generate the Miyake and Noda (EPS 2019) replica for rsf_solve.

Their model: rate-weakening patch of width 2R = 1 km (smoothed boxcar in
b, a uniform) on a planar antiplane fault in a UNIFORM Maxwell-viscoelastic
INFINITE medium; sigma = 100 MPa, mu = 30 GPa, cs = 3000 m/s, f0 = 0.6 at
the reference velocity Vpl = 1e-9 m/s, aging law, a = 0.02, a/b = 0.8 in
the patch; loading: far-field stress with the rate-strengthening fault
creeping at Vpl (equivalent to backslip at Vpl).  Their relaxation time
t_c equals -ve_tmaxwell (mode 1, single exact pole); the loading timescale
t_load = 2 a sigma R/(mu Vpl) = 6.67e7 s.

Approximations vs the paper: our fault strip is buried 30 km deep in a
half-space (free surface 52R away, negligible) and is FINITE (length
16R with welded exterior) instead of 4.096R-periodic; and rsf_solve is
quasi-dynamic while their SBIEM is fully dynamic (they show QD shifts
things quantitatively).  Comparison is therefore of regime and scaling.

usage: gen_mn19.py [ds_m] [Rc_over_R] [len_R] [prefix] [npatch]
       defaults    6       0.4         16     mn19     1
(Rc_over_R sets L via Rubin-Ampuero: Rc = (1/pi)(b/(b-a))^2 mu L/(b sig);
pass the printed -dc to rsf_solve.  npatch > 1 places identical patches
every 4.096R, approximating their periodic boundaries; lenR is then
overridden to npatch*4.096.  Patch cells are group 0, creep cells 1.)
"""
import sys
import numpy as np

ds     = float(sys.argv[1]) if len(sys.argv) > 1 else 6.0
RcR    = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
lenR   = float(sys.argv[3]) if len(sys.argv) > 3 else 16.0
prefix = sys.argv[4]        if len(sys.argv) > 4 else "mn19"
npatch = int(sys.argv[5])   if len(sys.argv) > 5 else 1
if npatch > 1:
    lenR = npatch*4.096

G    = 30.0e9; cs = 3000.0; eta = G/(2.0*cs)
sig  = 100.0e6; f0 = 0.6; Vpl = 1.0e-9
R    = 500.0
a0   = 0.02; bmax = a0/0.8            # a/b = 0.8 in the patch
L    = RcR*R*np.pi*(bmax*sig/G)*((bmax-a0)/bmax)**2
zdep = 30.0e3                          # patch center depth [m]

# smoothed boxcar for b (Noda and Lapusta 2010 style): transition over 0.1R.
# Their Fig. 2: b = 0 OUTSIDE the patch (a-b = +0.02); floored at 0.001
# because the psi-based engine divides by b (state term then negligible)
b_out = 0.001
wtr   = 0.1*R
n  = int(round(lenR*R/ds))
half = ds/2.0
fg = open(f"{prefix}_geom.in","w"); fr = open(f"{prefix}_rsf.dat","w")
fi = open(f"{prefix}_ic.in","w")
for j in range(n):
    x = (j+0.5)*ds - lenR*R/2.0        # along-depth coordinate rel. strip center
    z = zdep + x
    # smooth boxcar: 1 inside |x|<R-wtr, 0 outside |x|>R+wtr, cos taper
    # between; with npatch > 1 measured from the nearest patch center
    # (patch centers every 4.096R)
    if npatch > 1:
        k = np.clip(np.round((x/(R*4.096)) + (npatch-1)/2.0), 0, npatch-1)
        ax = abs(x - (k - (npatch-1)/2.0)*4.096*R)
    else:
        ax = abs(x)
    if ax < R - wtr: w = 1.0
    elif ax > R + wtr: w = 0.0
    else: w = 0.5*(1.0+np.cos(np.pi*(ax-(R-wtr))/(2*wtr)))
    b = b_out + (bmax-b_out)*w
    grp = 0 if w > 0.5 else 1
    # M&N Eq (23) initial perturbation: ln(Vpl th/L) = ln(Vpl/Vini)(1+0.01x/2R),
    # Vini = 0.999 Vpl; theta above steady state -> v below Vpl at steady tau
    eps = np.log(1.0/0.999)*(1.0 + 0.01*x/(2*R))
    v = Vpl*np.exp(-(b/a0)*eps)
    tau = sig*f0 + eta*Vpl                          # steady-state stress at Vpl
    fg.write(f"0 {-z:.6e} 0 0 -90 {half:.6e} 0 {grp}\n")
    fr.write(f"{a0:.6f} {b:.6f}\n")
    fi.write(f"{tau:.8e} {v:.6e}\n")
fg.close(); fr.close(); fi.close()
Lb = G*L/(bmax*sig)
tload = 2*a0*sig*R/(G*Vpl)
print(f"gen_mn19: {n} cells ds {ds} m, Rc/R {RcR}, L(dc) {L:.6e} m, "
      f"Lb {Lb:.1f} m (need ds < Lb/3 = {Lb/3:.1f}), t_load {tload/3.15576e7:.3f} yr")
print(f"  pass: -dc {L:.6e} -sigma_init 100e6 -f0 0.6 -v0 1e-9 -vpl 1e-9 "
      f"-shear_modulus 30e9 -s_wave_speed 3000 -state_law 1")
