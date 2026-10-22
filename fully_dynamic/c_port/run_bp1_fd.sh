#!/bin/sh
# SEAS BP1 with rsf_solve, QD and FD, at a given element size.
# usage: fully_dynamic/c_port/run_bp1_fd.sh [dz_m] [tend_yr] [lmax] [extra FD options]
#   e.g. run_bp1_fd.sh 25 520 12000 "-fd_nl 1" for second-order marching (FD run only)
#   dz 50 (default) or 25; tend 520; lmax 6000 for 50 m, 12000 for 25 m
# Writes geom_bp1.in, rsf_bp1.dat, ic_bp1.in in the current directory,
# runs QD into bp1_qd_dz<dz>/ and FD into bp1_fd_dz<dz>/.
# Prototype reference (fully_dynamic/README.md): QD events 196.36, 271.98,
# 349.82, 428.12, 506.48 yr (25 m); FD 196.36, 302.11, 420.81 yr (25 m),
# 196.37, 301.80, 420.54 (50 m); FD mean slip per event 1.556, 1.303, 1.459 m.
dz=${1-50}
tend=${2-520}
lmax=${3-6000}
extra=${4-}
bin=`pwd`/bin/rsf_solve
[ -x $bin ] || { echo "run from the repository root (bin/rsf_solve not found)"; exit 1; }
python3 - $dz << 'PY'
import sys, numpy as np
dz=float(sys.argv[1]); Wf=40e3; n=int(Wf/dz); H=15e3; h=3e3; a0=0.010; amax=0.025; b=0.015
with open("geom_bp1.in","w") as g, open("rsf_bp1.dat","w") as r, open("ic_bp1.in","w") as ic:
    for i in range(n):
        z=(i+0.5)*dz
        a = a0 if z<H else (a0+(amax-a0)*(z-H)/h if z<H+h else amax)
        g.write(f"0.000000e+00 {-z:.6e} 0 0 -90 {dz/2:.6e} 0 0\n")
        r.write(f"{a:.6f} {b:.6f}\n")
        ic.write(f"{26.5461e6:.8e} 1.000000e-09\n")
# stations at 2.5, 7.5, 12.5 km: SEAS names and the element indices
i25=int(2.5e3/dz-0.5); i75=int(7.5e3/dz-0.5); i125=int(12.5e3/dz-0.5)
open("stations_bp1.in","w").write(f"dp025 {i25}\ndp075 {i75}\ndp125 {i125}\n")
open("stations_bp1.idx","w").write(f"{i25},{i75},{i125}\n")
print("wrote", n, "elements; station indices", i25, i75, i125)
PY
st=`cat stations_bp1.idx`
common="-geom_file ../geom_bp1.in -rsf_file ../rsf_bp1.dat -rsf_ic_file ../ic_bp1.in -rsf_stations ../stations_bp1.in \
  -dc 0.008 -sigma_init 50e6 -v0 1e-6 -f0 0.6 -vpl 1e-9 -shear_modulus 3.204e10 -s_wave_speed 3464 \
  -rtol 1e-6 -stop_time_yr $tend -track_events 1 -vel_event 1e-3 -dt_monitor_yr 50"
mkdir -p bp1_qd_dz$dz bp1_fd_dz$dz
( cd bp1_qd_dz$dz && $bin $common > run.log 2>&1 )
echo "QD events (yr):"; awk '$3==1 {printf "%.3f ", $2} END {print ""}' bp1_qd_dz$dz/rsf_events.dat
( cd bp1_fd_dz$dz && $bin $common -fd_mode 1 -fd_lmax $lmax -fd_stations $st $extra > run.log 2>&1 )
echo "FD phases:"; grep "rsf_fd: phase" bp1_fd_dz$dz/run.log
