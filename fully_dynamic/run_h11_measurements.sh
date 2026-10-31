#!/bin/bash
# Measurements for fd_plan.md H11, to run on a multi-core machine from the repository root
# after "make bin/fd3d_test bin/rsf_solve". Every part writes its output to $out/<part>.log and a
# one-line summary to $out/summary.txt; send back the directory (or at least summary.txt).
#
# usage: fully_dynamic/run_h11_measurements.sh [parts] [max ranks] [mpirun command] [output dir]
#   parts: any of t m s g b (default "t m s g b"):
#     t  regression tests (run_tests.sh, about 15 minutes on one core)
#     m  memory of the partitioned TPV14 cross tables (D1) at 1000, 750, 500 and 250 m (table
#        build only; 250 m needs about 30 to 40 GB, skipped when the run fails)
#     s  MPI speed-up: TPV14 at 500 m for 2 s on 1, 2, 4, ... up to [max ranks] ranks, with a
#        check that the final states are identical on all rank counts
#     g  graded surface strips: TPV14 at 750 and 500 m (whole run, 10 s) and the main fault
#        alone at 500 m against the mirror plus opening layer
#     b  the band preconditioned strip solve against the dense LU where it is chosen (TPV14,
#        1 km faults, strips of 333 m graded to +-15 km, 3 s; several GB)
#     r  resolution ladder against the SCEC stations: TPV14 and TPV15 (mirror) at each cell size of
#        [r ladder] over 12 s on [max ranks] ranks, with the on-fault stations of the SCEC files in
#        [scec dir] (named tpv14_faultst020dp000.txt, tpv14_branchst020dp000.txt,
#        tpv15_faultst-020dp000.txt, ...: tpv<14|15>_<fault|branch>st<s>dp<z>, s and z in 100 m;
#        plain text or the HTML saved from the SCEC site); 250 m needs about 70 GB in all
#   max ranks: largest number of MPI ranks (default 16)
#   mpirun command: default "mpirun"
#   output dir: default h11_out
#   scec dir: directory with SCEC station files (optional, for part r)
#   r ladder: cell sizes of part r (default "1000 750 500 375 250"; 200 or 150 m need much more
#     memory, see the m part for the trend)
parts=${1-"t m s g b"}
nmax=${2-16}
mpi=${3-mpirun}
out=${4-h11_out}
scec=${5-}
ladder=${6-"1000 750 500 375 250"}
bin=`pwd`/bin/fd3d_test
d=`pwd`/src/testing/rsf_solve
mkdir -p $out
out=`cd $out && pwd`
sum=$out/summary.txt
echo "H11 measurements `date`, `hostname`, `nproc` cores" >> $sum

peak(){                               # run a command, report wall time and peak memory if possible
  local s=`date +%s`
  if [ -x /usr/bin/time ]; then
    /usr/bin/time -v -o time.err "$@"
    grep 'Maximum resident' time.err | awk '{printf "peak %.0f MB, ", $6/1024}'
  else
    "$@"
  fi
  echo "wall $(( `date +%s` - s )) s"
}

for p in $parts; do
  case $p in
  t)
    ( src/testing/rsf_solve/run_tests.sh > $out/t.log 2>&1 )
    echo "t: run_tests: `grep -c PASS $out/t.log` PASS, `grep -c FAIL $out/t.log` FAIL" >> $sum
    ;;
  m)
    mkdir -p $out/m && cd $out/m
    for dx in 1000 750 500 250; do
      r=`peak env FD3D_DP_DIAG=1 $bin tpv14 $dx 14 0.05 > m_$dx.log 2>&1; tail -1 m_$dx.log`
      echo "m: tpv14 $dx m: `grep 'cross tables' m_$dx.log` ($r)" >> $sum
      grep 'fd3d_multi: pair' m_$dx.log >> $sum
    done
    cd - > /dev/null
    ;;
  s)
    mkdir -p $out/s && cd $out/s
    n=1
    while [ $n -le $nmax ]; do
      mkdir -p n$n && cd n$n
      s0=`date +%s`
      env FD3D_DUMP=dump.txt $mpi -np $n $bin tpv14 500 14 2 > run.log 2>&1
      echo "s: tpv14 500 m, 2 s, $n ranks: wall $(( `date +%s` - s0 )) s, `grep 'cross tables' run.log`" >> $sum
      cd ..
      if [ $n -gt 1 ]; then
        if cmp -s n1/dump.txt n$n/dump.txt; then echo "s: $n ranks identical to 1 rank" >> $sum
        else echo "s: $n ranks DIFFER from 1 rank" >> $sum; fi
      fi
      n=$(( 2*n ))
    done
    cd - > /dev/null
    ;;
  g)
    mkdir -p $out/g && cd $out/g
    for dx in 750 500; do
      mkdir -p tpv14s_$dx && cd tpv14s_$dx
      r=`peak env FD3D_TPV14_STATIONS=0:2000:0,1:2000:0 FD3D_TPV14S_GRADE=1.3 FD_STRIPS_DIAG=1 $bin tpv14s $dx 14 10 0 1 90 15000 > run.log 2>&1; tail -1 run.log`
      echo "g: `grep '^tpv14s: dx' run.log | tail -1` ($r)" >> $sum
      mkdir -p ../tpv14_$dx && cd ../tpv14_$dx
      env FD3D_TPV14_STATIONS=0:2000:0,1:2000:0 $bin tpv14 $dx 14 10 > run.log 2>&1
      echo "g: mirror: `tail -1 run.log`" >> $sum
      python3 - ../tpv14_$dx/tpv14_rt_%s.dat ../tpv14s_$dx/tpv14s_rt_%s.dat >> $sum << 'PYEOF'
import sys, math
for part in ('main', 'branch'):
    a = [tuple(map(float, l.split())) for l in open(sys.argv[1] % part)]
    b = [tuple(map(float, l.split())) for l in open(sys.argv[2] % part)]
    d = [y[2] - x[2] for x, y in zip(a, b) if x[2] >= 0 and y[2] >= 0]
    top = [y[2] - x[2] for x, y in zip(a, b) if x[2] >= 0 and y[2] >= 0 and x[1] < 3000]
    print(f"g:   {part}: strips - mirror rupture time rms {math.sqrt(sum(t*t for t in d)/max(len(d), 1)):.3f} s, top 3 km mean {sum(top)/max(len(top), 1):+.3f} s")
PYEOF
      cd ..
    done
    FD3D_TPV14S_GRADE=1.3 $bin tpv14one 500 10 1 15000 > tpv14one_500.log 2>&1
    grep 'against' tpv14one_500.log | sed 's/^/g: /' >> $sum
    cd - > /dev/null
    ;;
  b)
    mkdir -p $out/b && cd $out/b
    for sv in bj lu; do
      r=`peak env FD_STRIPS_SOLVER=$sv FD_STRIPS_DIAG=1 FD3D_TPV14S_GRADE=1.3 $bin tpv14s 1000 14 3 0 3 90 15000 > $sv.log 2>&1; tail -1 $sv.log`
      echo "b: solver $sv: `grep 'GMRES iterations' $sv.log` ($r)" >> $sum
    done
    cd - > /dev/null
    ;;
  r)
    # TPV14 and TPV15 (mirror) at each cell size of the ladder on nmax ranks, with the on-fault
    # stations of the SCEC files in $scec: names tpv<14|15>_<fault|branch>st<s>dp<z>* with s and z
    # in 100 m (e.g. tpv14_faultst020dp000.txt: main fault, 2 km along strike, surface;
    # tpv15_faultst-020dp000.txt: main fault, -2 km); without files the stations of the
    # benchmark files used so far
    if [ -n "$scec" ]; then
      scec=`cd $scec && pwd`
    fi
    for w in 14 15; do
      st=""
      for f in `ls $scec/tpv${w}_*st*dp* 2> /dev/null`; do
        nm=`basename $f`
        k=0
        case $nm in
          tpv${w}_branch*) k=1 ;;
        esac
        sz=`echo $nm | sed 's/^tpv[0-9]*_[a-z]*st\(-\{0,1\}[0-9]*\)dp\([0-9]*\).*/\1 \2/' | awk '{printf "%d:%d", $1*100, $2*100}'`
        st="$st${st:+,}$k:$sz"
      done
      if [ -z "$st" ]; then
        if [ $w -eq 14 ]; then st="0:2000:0,1:2000:0"; else st="0:-2000:0"; fi
      fi
      mkdir -p $out/r/tpv$w && cd $out/r/tpv$w
      for dx in $ladder; do
        mkdir -p dx$dx && cd dx$dx
        s0=`date +%s`
        FD3D_TPV14_STATIONS=$st $mpi -np $nmax $bin tpv14 $dx $w 12 > run.log 2>&1
        echo "r: tpv$w, $dx m, $nmax ranks: `grep "^tpv$w: dx" run.log | tail -1 | cut -c1-120` (wall $(( `date +%s` - s0 )) s)" >> $sum
        for f in `ls $scec/tpv${w}_*st*dp* 2> /dev/null`; do
          nm=`basename $f`
          k=0
          case $nm in
            tpv${w}_branch*) k=1 ;;
          esac
          sz=`echo $nm | sed 's/^tpv[0-9]*_[a-z]*st\(-\{0,1\}[0-9]*\)dp\([0-9]*\).*/\1 \2/' | awk '{printf "%d_%d", $1*100, $2*100}'`
          if [ -s tpv${w}_st_${k}_${sz}.dat ]; then
            python3 $d/tpv14_stcmp.py "$f" tpv${w}_st_${k}_${sz}.dat | sed "s/^/r: tpv$w $dx m: /" >> $sum
          fi
        done
        cd ..
      done
      cd - > /dev/null
    done
    ;;
  esac
done
echo "done `date`" >> $sum
