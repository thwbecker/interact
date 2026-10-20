#!/bin/sh
# SEAS BP4-QD with rsf_solve: vertical planar fault in a whole space, rate and
# state (aging law), fixed rake (strike slip). The rate-and-state region is
# |x3| <= Wf = 40 km; along strike the fault is truncated at |x2| <= X2MAX
# (default 50 km, 17 km of velocity-strengthening fault beyond the transition
# zone); outside the modelled patches the backslip formulation lets the fault
# creep at Vp. usage, from the repository root:
#   fully_dynamic/c_port/run_bp4_qd.sh [dx_m] [tend_yr] [x2max_km] [extra rsf_solve options]
# e.g. run_bp4_qd.sh 500 1000 50 "-use_hmatrix 1"   (HACApK) for the suggested 500 m
dx=${1-2000}
tend=${2-1000}
x2max=${3-50}
extra=${4-}
bin=`pwd`/bin/rsf_solve
[ -x $bin ] || { echo "run from the repository root"; exit 1; }
d=bp4_qd_dx$dx
mkdir -p $d && cd $d
python3 - $dx $x2max << 'PY'
import sys, numpy as np
dx=float(sys.argv[1]); x2max=float(sys.argv[2])*1e3
rho, cs = 2670.0, 3464.0; mu = rho*cs**2; eta = mu/(2*cs)
a0, amax, b0, sig, L, Vp, Vinit, V0, f0 = 0.0065, 0.025, 0.013, 50e6, 0.04, 1e-9, 1e-9, 1e-6, 0.6
H, l, h, Wf, wi, hi = 15e3, 60e3, 3e3, 40e3, 12e3, 1.5e3
n2 = int(round(2*x2max/dx)); n3 = int(round(2*Wf/dx))
x2 = -x2max + (np.arange(n2) + 0.5)*dx; x3 = -Wf + (np.arange(n3) + 0.5)*dx
tau0 = sig*amax*np.arcsinh(Vinit/(2*V0)*np.exp((f0 + b0*np.log(V0/Vinit))/amax)) + eta*Vinit
zoff = -(Wf + 10e3)                      # interact z of x3 = 0 (all patches below z = 0; whole space)
with open("geom.in","w") as g, open("rsf.dat","w") as r, open("ic.in","w") as ic:
    for i in range(n2):
        for j in range(n3):
            X2, X3 = x2[i], x3[j]
            rr = max(abs(X3) - H, abs(X2) - l/2)/h
            a = a0 if rr <= 0 else (amax if rr >= 1 else a0 + rr*(amax - a0))
            nuc = (-l/2 + hi <= X2 <= -l/2 + hi + wi) and (-H + hi <= X3 <= -H + hi + wi)
            g.write(f"{X2:.6e} 0.000000e+00 {zoff - X3:.6e} 90 90 {dx/2:.6e} {dx/2:.6e} 0\n")
            r.write(f"{a:.6f} {b0:.6f}\n")
            ic.write(f"{(1.1 if nuc else 1.0)*tau0:.8e} {Vinit:.6e}\n")
st = [(-36.0, 0.0), (-22.5, -7.5), (-16.5, -12.0), (-16.5, 0.0), (-16.5, 12.0), (0.0, -21.0), (0.0, -12.0), (0.0, 0.0),
      (0.0, 12.0), (0.0, 21.0), (16.5, -12.0), (16.5, 0.0), (16.5, 12.0), (36.0, 0.0)]
with open("stations.in","w") as f:
    for s2, s3 in st:
        i = int(np.argmin(np.abs(x2 - s2*1e3))); j = int(np.argmin(np.abs(x3 - s3*1e3)))
        f.write(f"strk{int(round(s2*10)):+04d}dp{int(round(s3*10)):+04d} {i*n3 + j}\n")
print(f"BP4: {n2} x {n3} = {n2*n3} squares of {dx:g} m, tau0 {tau0/1e6:.4f} MPa")
PY
$bin -geom_file geom.in -rsf_file rsf.dat -rsf_ic_file ic.in -rsf_stations stations.in \
  -dc 0.04 -sigma_init 50e6 -v0 1e-6 -f0 0.6 -vpl 1e-9 -shear_modulus 3.2038e10 -s_wave_speed 3464 \
  -full_space 1 -rsf_slip_mode 0 -rtol 1e-6 -stop_time_yr $tend -track_events 1 -vel_event 1e-3 \
  -dt_monitor_yr 10 $extra > run.log 2>&1
echo "events (onset, yr):"; awk '$3==1 {printf "%.2f ", $2} END {print ""}' rsf_events.dat
