#!/bin/bash
#
# 2-D in-plane megathrust test with the fully dynamic extension (fd_plan.md M7, M8): a Cascadia
# section from test_twod/cascadia resampled to elements of about h m (coarsen_2d.py), run
#   q  quasi-dynamic alone (reference)
#   s  setup only (tables, virtual free surface) with MPI on np ranks and with OpenMP on np threads,
#      to compare the times of the table build
#   m  QD with dynamic phases, MPI on np ranks (one thread each)
#   o  QD with dynamic phases, OpenMP on np threads (one rank)
#   c  compare: events of q, m and o; MPI against OpenMP (the station files must be identical, the
#      tables and the history sums are the same for any split over ranks or threads); wall times
# Results go to outdir/summary.txt. Run from the repository root after make bin/rsf_solve.
#
# usage: run_2d_cascadia.sh [parts] [h] [np] [mpi command] [years] [section dir] [outdir]
#   defaults: "q s m o c" 500 16 mpirun 700 test_twod/cascadia/PD13_nosigma 2d_cascadia_out
#
parts=${1-"q s m o c"}
h=${2-500}
np=${3-16}
mpi=${4-mpirun}
years=${5-700}
sec=${6-test_twod/cascadia/PD13_nosigma}
out=${7-2d_cascadia_out}
root=`pwd`
bin=$root/bin/rsf_solve
mkdir -p $out
out=`cd $out && pwd`
sum=$out/summary.txt
echo "run_2d_cascadia `date`, `hostname`: parts '$parts', h $h m, $np ranks/threads, $years yr, $sec" >> $sum
# the resampled section
if [ ! -s $out/geom/rsf_ab.in ]; then
  python3 $root/fully_dynamic/coarsen_2d.py $root/$sec $out/geom $h | tee -a $sum
fi
geom=`cd $out/geom && ls geom_*.in | head -1`
base="-geom_file ../geom/$geom -rsf_file ../geom/rsf_ab.in -rsf_dc_file ../geom/rsf_dc.in -rsf_sigma_file ../geom/rsf_sigma.in \
 -rsf_vpl_file ../geom/rsf_vpl.in -rsf_ic_file ../geom/rsf_ic.in -state_law 1 -shear_modulus 3.0e10 -s_wave_speed 3464 \
 -f0 0.6 -v0 1e-6 -vpl 1.2e-9 -vel_event 1e-3 -track_events 1 -ts_rk_type 3bs -rtol 1e-5"
fdopt="-fd_mode 1 -fd_surf_margin 10000 -fd_stations 40,200"
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
  grep 'split by receiver rows\|OpenMP threads\|dense pairs\|static gate\|virtual free surface:' run.log | sed "s/^/  $name: /" >> $sum
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
    c) python3 - $out << 'PYEOF' >> $sum
import sys, os, glob, filecmp
import numpy as np
out = sys.argv[1]
def events(d):
    f = os.path.join(out, d, "rsf_events.dat")
    if not os.path.exists(f):
        return None
    a = np.loadtxt(f, comments="#", ndmin=2)
    if a.size == 0:
        return []
    on = a[a[:, 2] == 1]
    off = a[a[:, 2] == -1]
    return [(on[i, 1], off[i, 4] if i < len(off) else np.nan) for i in range(len(on))]
for d in ("qd", "fd_mpi", "fd_omp"):
    e = events(d)
    if e is None:
        print(f"c: {d}: no events file")
        continue
    print(f"c: {d}: {len(e)} events, onsets (yr) / mean slip at arrest (m): " + ", ".join(f"{t:.2f} / {s:.3f}" for t, s in e))
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
