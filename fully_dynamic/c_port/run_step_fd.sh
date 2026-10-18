#!/bin/sh
# Stepover sequences with rsf_solve, QD and FD: two parallel 24 km in-plane
# segments, step W (km, sign: +W puts segment 2 at +y; with interact's
# left-lateral positive slip +W is the dilational side), overlap (km),
# velocity-strengthening margin (km), 100 m elements.
# usage: fully_dynamic/c_port/run_step_fd.sh W overlap margin [tend_yr] [mode: qd|fd|both]
# Output directory step_<mode>_W<W>_ov<ov>_vs<vs>/ with rsf_events.dat and the
# two station files fltst_s1.dat, fltst_s2.dat (centres of the segments),
# which fully_dynamic/c_port/summarize_rsf_step.py turns into the jump table.
W=${1-1.0}
ov=${2-2.0}
vs=${3-0.5}
tend=${4-700}
mode=${5-both}
bin=`pwd`/bin/rsf_solve
gen=`pwd`/src/testing/rsf_solve/gen_step.py
[ -x $bin ] || { echo "run from the repository root"; exit 1; }
for m in qd fd; do
  [ "$mode" = "both" ] || [ "$mode" = "$m" ] || continue
  d=step_${m}_W${W}_ov${ov}_vs${vs}
  mkdir -p $d && cd $d
  python3 $gen $W $ov $vs > gen.log
  printf "s1 120\ns2 360\n" > stations.in
  opts="-geom_file geom.in -rsf_file rsf.dat -rsf_ic_file ic.in -rsf_stations stations.in \
    -dc 0.016 -sigma_init 50e6 -v0 1e-6 -f0 0.6 -vpl 1e-9 -shear_modulus 3.24e10 -s_wave_speed 3464 \
    -full_space 1 -calc_sigma_dot 1 -rtol 1e-6 -stop_time_yr $tend -track_events 1 -vel_event 1e-3 -dt_monitor_yr 50"
  if [ "$m" = "fd" ]; then opts="$opts -fd_mode 1"; fi
  $bin $opts > run.log 2>&1
  cd ..
  python3 `dirname $0`/summarize_rsf_step.py $d
done
