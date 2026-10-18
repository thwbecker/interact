#!/bin/sh
# Stepover sequences with rsf_solve, QD and FD: two parallel in-plane
# segments (first 24 km), step W (km; +W puts segment 2 at +y, the
# dilational side for interact's left-lateral positive slip), overlap (km),
# velocity-strengthening margin (km), 100 m elements.
# usage: fully_dynamic/c_port/run_step_fd.sh W overlap margin [tend_yr] [qd|fd|both] [L2_km] [dba2]
#   L2_km  length of segment 2 (default 24), dba2  b - a on its weakening part
#   (default 0.005, as segment 1); unequal values let the receiver's cycle
#   phase vary from event to event.
# Output directory step_<mode>_W<W>_ov<ov>_vs<vs>[_L<L2>_ba<dba2>]/ with
# rsf_events.dat, fltst_s1.dat, fltst_s2.dat (segment centres) and the phase
# files; summarize_rsf_step.py turns it into the jump table.
W=${1-1.0}
ov=${2-2.0}
vs=${3-0.5}
tend=${4-700}
mode=${5-both}
L2=${6-24}
dba2=${7-0.005}
bin=`pwd`/bin/rsf_solve
gen=`pwd`/src/testing/rsf_solve/gen_step.py
here=`cd \`dirname $0\` && pwd`
[ -x $bin ] || { echo "run from the repository root"; exit 1; }
suffix=""
if [ "$L2" != "24" ] || [ "$dba2" != "0.005" ]; then suffix="_L${L2}_ba${dba2}"; fi
for m in qd fd; do
  [ "$mode" = "both" ] || [ "$mode" = "$m" ] || continue
  d=step_${m}_W${W}_ov${ov}_vs${vs}${suffix}
  mkdir -p $d && cd $d
  python3 $gen $W $ov $vs $L2 $dba2 > gen.log
  n1=`awk '{print $1}' segments.txt`
  n2=`awk '{print $2}' segments.txt`
  printf "s1 %d\ns2 %d\n" `expr $n1 / 2` `expr $n1 + $n2 / 2` > stations.in
  opts="-geom_file geom.in -rsf_file rsf.dat -rsf_ic_file ic.in -rsf_stations stations.in \
    -dc 0.016 -sigma_init 50e6 -v0 1e-6 -f0 0.6 -vpl 1e-9 -shear_modulus 3.24e10 -s_wave_speed 3464 \
    -full_space 1 -calc_sigma_dot 1 -rtol 1e-6 -stop_time_yr $tend -track_events 1 -vel_event 1e-3 -dt_monitor_yr 50"
  if [ "$m" = "fd" ]; then opts="$opts -fd_mode 1"; fi
  $bin $opts > run.log 2>&1
  cd ..
  python3 $here/summarize_rsf_step.py $d
done
