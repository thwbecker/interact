#!/bin/sh
# SEAS BP4 with rsf_solve, using the repository's BP4 generator
# (seas_bp_tests/bp4/make_bp4.py: whole space, symmetric fault |x3| <= 40 km,
# velocity seed Vi = 0.03 m/s in the nucleation square with the background
# state unchanged, as in the published comparisons of Jiang et al. 2022, so the
# first event starts immediately and event times line up with theirs; fixed
# rake). QD by default; add "-fd_mode 1" to the extra options for the fully
# dynamic variant. usage, from the repository root:
#   fully_dynamic/c_port/run_bp4_qd.sh [ds_km] [Lstrike_km] [tend_yr] [extra rsf_solve options]
# e.g. run_bp4_qd.sh 0.5 120 1000 "-use_hmatrix 1"
ds=${1-1}
L=${2-120}
tend=${3-1000}
extra=${4-}
root=`pwd`
bin=$root/bin/rsf_solve
[ -x $bin ] || { echo "run from the repository root"; exit 1; }
d=bp4_ds${ds}_L${L}
mkdir -p $d && cd $d
python3 $root/seas_bp_tests/bp4/make_bp4.py $ds $L
res=`printf '%gkm' $ds`; tag=L`printf '%g' $L`
$bin -geom_file geom_bp4_${res}_${tag}.in -rsf_file rsf_bp4_${res}_${tag}.dat -rsf_ic_file ic_bp4_${res}_${tag}.in \
  -full_space 1 -shear_modulus 3.204e10 -s_wave_speed 3464 -f0 0.6 -dc 0.04 -vpl 1e-9 -v0 1e-6 -sigma_init 50e6 \
  -rtol 1e-4 -ts_rk_type 3bs -ts_max_steps 2000000 -stop_time_yr $tend -track_events 1 -vel_event 1e-3 $extra > run.log 2>&1
echo "events (onset, yr):"; awk '$3==1 {printf "%.2f ", $2} END {print ""}' rsf_events.dat
