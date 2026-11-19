#!/bin/bash
#
# 2-D in-plane megathrust test with the fully dynamic extension (fd_plan.md M7, M8): a Cascadia
# section from test_twod/cascadia resampled to elements of about h m (coarsen_2d.py), run
#   q  quasi-dynamic alone (reference)
#   s  setup only (tables, with the reflected waves of the free surface) with MPI on np ranks and with OpenMP on np threads,
#      to compare the times of the table build
#   m  QD with dynamic phases, MPI on np ranks (one thread each)
#   o  QD with dynamic phases, OpenMP on np threads (one rank)
#   c  compare: events of q, m and o; MPI against OpenMP (the station files must be identical, the
#      tables and the history sums are the same for any split over ranks or threads); wall times
# Results go to outdir/summary.txt. Run from the repository root after make bin/rsf_solve.
#
# usage: run_2d_cascadia.sh [parts] [h] [np] [mpi command] [years] [section dir] [outdir] [options] [stations]
#   defaults: "q s m o c" 500 16 mpirun 700 test_twod/cascadia/PD13_nosigma 2d_cascadia_out "" auto
#   section dir: any directory in the format of test_twod/cascadia/PD13_nosigma (geom_*.in and the
#     rsf_*.in files), e.g. one written by make_2d_subduction.py (idealized megathrust and splays)
#   options: further rsf_solve options for all runs (QD and FD alike), e.g. "-calc_sigma_dot 1" for
#     normal-stress coupling, which a QD against FD comparison of a dipping fault should have
#   stations: FD station elements (-fd_stations), auto: a fifth and a half of the megathrust
# The reflected half-plane tables take FD_HPDWN_* from the environment (FD_HPDWN_Q=2 is adequate
# for faults whose tops stay about a kilometre below the surface, fd_plan.md M10), e.g.
#   OPENBLAS_NUM_THREADS=1 FD_HPDWN_Q=2 fully_dynamic/subduction2d/run_2d_cascadia.sh "q o c" 1000 32 mpirun 700 \
#     test_twod/cascadia/PD13_nosigma 2d_pd13_1km_sd "-calc_sigma_dot 1"
#
parts=${1-"q s m o c"}
h=${2-500}
np=${3-16}
mpi=${4-mpirun}
years=${5-700}
sec=${6-test_twod/cascadia/PD13_nosigma}
out=${7-2d_cascadia_out}
xopt=${8-""}
stations=${9-auto}
root=`pwd`
bin=$root/bin/rsf_solve
mkdir -p $out
out=`cd $out && pwd`
sum=$out/summary.txt
echo "run_2d_cascadia `date`, `hostname`: parts '$parts', h $h m, $np ranks/threads, $years yr, $sec, options '$xopt'" >> $sum
env | grep '^FD_HPDWN' | sed 's/^/  environment: /' >> $sum
# the resampled section
if [ ! -s $out/geom/rsf_ab.in ]; then
  case $sec in
    /*) secp=$sec ;;
    *) secp=$root/$sec ;;
  esac
  python3 $root/fully_dynamic/subduction2d/coarsen_2d.py $secp $out/geom $h | tee -a $sum
fi
geom=`cd $out/geom && ls geom_*.in | head -1`
base="-geom_file ../geom/$geom -rsf_file ../geom/rsf_ab.in -rsf_dc_file ../geom/rsf_dc.in -rsf_sigma_file ../geom/rsf_sigma.in \
 -rsf_vpl_file ../geom/rsf_vpl.in -rsf_ic_file ../geom/rsf_ic.in -state_law 1 -shear_modulus 3.0e10 -s_wave_speed 3464 \
 -f0 0.6 -v0 1e-6 -vpl 1.2e-9 -vel_event 1e-3 -track_events 1 -ts_rk_type 3bs -rtol 1e-5 $xopt"
if [ "$stations" = "auto" ]; then
  stations=`awk '$8 == 0 {n++} END {a = int(n/5); b = int(n/2); if(a < 1) a = 1; print a "," b}' $out/geom/$geom`
fi
fdopt="-fd_mode 1 -fd_stations $stations"
# run name, launcher (mpi ranks, 0 for none), threads, extra options
run () {
  name=$1
  nr=$2
  nt=$3
  shift 3
  rm -rf $out/$name
  mkdir -p $out/$name
  cd $out/$name
  s0=`date +%s`
  if [ $nr -gt 0 ]; then
    OMP_NUM_THREADS=$nt $mpi -np $nr $bin $base "$@" > run.log 2>&1
  else
    OMP_NUM_THREADS=$nt $bin $base "$@" > run.log 2>&1
  fi
  rc=$?
  echo "$name: rc $rc, wall $(( `date +%s` - s0 )) s" >> $sum
  grep 'split by receiver rows\|OpenMP threads\|dense pairs\|static gate\|half plane' run.log | sed "s/^/  $name: /" >> $sum
  grep 'rsf_fd: phase' run.log | cut -c1-170 | sed "s/^/  $name: /" >> $sum
  grep 'PETSC ERROR\|Segmentation\|Abort' run.log | head -3 | sed "s/^/  $name: /" >> $sum
  cd $root
}
for p in $parts; do
  case $p in
    q) run qd 0 1 -stop_time_yr $years ;;
    s) run setup_mpi $np 1 $fdopt -stop_time_yr 1e-9
       run setup_omp 0 $np $fdopt -stop_time_yr 1e-9 ;;
    m) run fd_mpi $np 1 $fdopt -stop_time_yr $years ;;
    o) run fd_omp 0 $np $fdopt -stop_time_yr $years ;;
    c) python3 $root/fully_dynamic/subduction2d/compare_runs_2d.py $out | sed 's/^/c: /' >> $sum
       python3 - $out << 'PYEOF' >> $sum
import sys, os, glob, filecmp
import numpy as np
out = sys.argv[1]
fm = sorted(glob.glob(os.path.join(out, "fd_mpi", "fd_phase_*_st*.dat")))
fo = sorted(glob.glob(os.path.join(out, "fd_omp", "fd_phase_*_st*.dat")))
if fm and fo:
    same = [os.path.basename(a) for a in fm if os.path.exists(os.path.join(out, "fd_omp", os.path.basename(a)))
            and filecmp.cmp(a, os.path.join(out, "fd_omp", os.path.basename(a)), shallow=False)]
    print(f"c: MPI against OpenMP: {len(same)} of {len(fm)} station files identical"
          + ("" if len(same) == len(fm) == len(fo) else " (FAIL: the results should not depend on ranks or threads)"))
PYEOF
       ;;
  esac
done
echo "done `date`" >> $sum
