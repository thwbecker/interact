#!/bin/sh
# Tests of the fully dynamic extension of rsf_solve (src/rsf_solve/fd*.c, rsf_fd.c; plan and
# results in fully_dynamic/docs/fd_plan.md). Run from the repository root after
#   make rsf_testing bin/rsf_solve
# usage: src/testing/rsf_solve/run_tests.sh [tiers] [MPI launch command]
#   tiers, comma separated (default unit,ruptures):
#     unit     the C modules against the prototype and exact references (about 3 minutes)
#     ruptures 2-D and 3-D ruptures of the test drivers fd_dynamic and fd3d_test against the
#              prototype and earlier values (about 10 minutes)
#     rsf      rsf_solve end to end (e2e/*.sh: static gates, cycles, 2-D half plane, 3-D
#              planes and surfaces; about 30 minutes on one core)
#     mpi      identical results on 1, 2, 3 and 4 MPI ranks (3-D marcher, dense couplings,
#              2-D faults; 30 minutes or more on one core)
#     all      all of them
#   MPI launch command: e.g. "mpirun", "mpiexec", "srun"; default the first of mpirun,
#   mpirun --oversubscribe, mpirun --allow-run-as-root --oversubscribe that starts 4 ranks
# Optional libraries for some checks (skipped without them): TDD_LIB, Crouch and Starfield's
# 2-D kernels for the Python checkers (default /tmp/tdd/libtdd.so; build:
#   gfortran -O2 -fPIC -shared -cpp -Isrc/includes -DUSE_DOUBLE_PRECISION src/green/tdd_coeff.F -o /tmp/tdd/libtdd.so)
# Prints one line per check (PASS, FAIL or SKIP, wall seconds, the check's own last line) and a
# summary; returns non-zero if any check failed. The prototype tests are separate:
# fully_dynamic/proto/tests/run_proto_tests.sh.
d=src/testing/rsf_solve
tmp=tmp_fd_dynamic
tiers=${1-unit,ruptures}
mpi=${2-}
npass=0
nfail=0
nskip=0
failed=""
t_start=`date +%s`
mkdir -p $tmp

# result name status detail: status 0 PASS, 1 FAIL, 2 SKIP
result() {
  case $2 in
    0) npass=$((npass + 1)); tag="PASS";;
    2) nskip=$((nskip + 1)); tag="SKIP";;
    *) nfail=$((nfail + 1)); tag="FAIL"; failed="$failed
  $1";;
  esac
  printf "%-4s %-58s %s\n" "$tag" "$1" "$3"
}

# check name command: runs the command (a string, in $tmp) and reports its exit status with the
# wall time and the last line of its output; a last line containing FAIL counts as failed too
check() {
  name=$1
  shift
  t0=`date +%s`
  ( cd $tmp && eval "$*" ) > $tmp/check.out 2>&1
  st=$?
  last=`grep -v '^ *$' $tmp/check.out | tail -1 | cut -c1-90`
  echo "$last" | grep -q FAIL && st=1
  result "$name" $st "($((`date +%s` - t0)) s) $last"
}

# skip name reason
skip() {
  result "$1" 2 "($2)"
}

want() {
  case ",$tiers," in
    *",$1,"*|*",all,"*) return 0;;
  esac
  return 1
}

have_python() {
  command -v python3 > /dev/null 2>&1
}

have_tdd() {
  have_python && [ -s ${TDD_LIB-/tmp/tdd/libtdd.so} ]
}

# fsrup_check name sense top_ref mid_ref bot_ref fault_top [env]: the 2-D half-plane rupture of
# fd_dynamic fsrup, |slip| top / middle / bottom within 1 cm of the reference values
fsrup_check() {
  name=$1
  t0=`date +%s`
  ( cd $tmp && eval "$7 ../bin/fd_dynamic fsrup $2 1 3000 200 $6 3 30" ) > $tmp/fsrup.out 2>/dev/null
  msg=`awk -v t=$3 -v m=$4 -v b=$5 '{for(i=1;i<=NF;i++) if($i=="bottom"){top=$(i+1); mid=$(i+3); bot=$(i+5)}}
       END {if(top == ""){print "no result"; exit 1}
            d=0; e=top-t; if(e<0)e=-e; if(e>d)d=e; e=mid-m; if(e<0)e=-e; if(e>d)d=e; e=bot-b; if(e<0)e=-e; if(e>d)d=e;
            print "slip " top " / " mid " / " bot " m, largest difference " d " m"; if(d>0.01) exit 1}' $tmp/fsrup.out`
  result "$name" $? "($((`date +%s` - t0)) s) $msg"
}

find_mpi() {
  if [ -z "$mpi" ] && command -v mpirun > /dev/null 2>&1; then
    for m in "mpirun" "mpirun --oversubscribe" "mpirun --allow-run-as-root --oversubscribe"; do
      if $m -np 4 true > /dev/null 2>&1; then
        mpi=$m
        break
      fi
    done
  fi
}

# ---------------------------------------------------------------------------------------------
if want unit; then
  echo "== unit: C modules against the prototype and exact references"
  check "FFT" "../bin/test_fd_fft"
  check "segments (polylines, element geometry)" "../bin/test_fd_segments"
  for f in $d/ref/ref_inplane_*.txt $d/ref/ref_antiplane_*.txt; do
    check "2-D tables `basename $f .txt | sed 's/ref_//'`" "../bin/test_fd_tables ../$f"
  done
  # the coupled rate-and-state solve of the dynamic phase, slip in both directions (2026-11-18)
  check "rate and state solve, Newton" "../bin/test_fd_rsf"
  check "rate and state solve, elementwise fallback" "TFR_MAXIT=1 ../bin/test_fd_rsf"
  if have_python; then
    # reflected half-plane stresses by DWN (fd_plan.md M10) against the prototype
    check "half plane DWN against the prototype" "python3 ../$d/py/test_fd_hpdwn.py ../bin/test_fd_hpdwn"
    # non-parallel 2-D pairs (fd_plan.md N1): dense tables for a kinked pair
    for a in 30 -45; do
      for m in 0 1; do
        check "dense 2-D pair, angle $a, mode $m" "FD_CT_OFF=1 ../bin/test_fd_dense $a $m dense_${a}_${m}.txt > /dev/null && python3 ../$d/py/test_fd_dense.py $a $m dense_${a}_${m}.txt"
      done
    done
  else
    skip "Python checkers (DWN, dense pairs)" "no python3"
  fi
  if have_tdd; then
    # batch reflected tables: DWN (M10) and Cagniard-de Hoop (M11) against the prototypes and
    # Crouch and Starfield
    check "half plane DWN batch tables" "../bin/test_fd_hpdwn_batch 4 250 500 40 > hpb.txt && python3 ../$d/py/test_fd_hpdwn_batch.py hpb.txt"
    check "half plane Cagniard-de Hoop batch tables" "../bin/test_fd_hpdwn_batch 4 250 500 40 cdh > hpc.txt && python3 ../$d/py/test_fd_hpcdh.py hpc.txt"
  else
    skip "half plane batch tables (DWN, Cagniard-de Hoop)" "no TDD_LIB"
  fi
  # 3-D kernels, surface layer and surface stability, quasi-2-D limit, plane pairs
  check "3-D kernels" "../bin/fd3d_test kernels ../$d/ref/ref_fd3d_kernels.txt"
  check "3-D surface layer" "../bin/fd3d_test layer ../$d/ref/ref_fd3d_layer.txt"
  check "3-D surface layer (FFT)" "FD3D_LAYER_FFT=1 ../bin/fd3d_test layer ../$d/ref/ref_fd3d_layer.txt"
  check "3-D virtual surface" "../bin/fd3d_test surf ../$d/ref/ref_fd3d_surf.txt"
  check "3-D quasi-2-D limit" "../bin/fd3d_test q2d ../$d/ref/ref_fd3d_q2d.txt"
  check "3-D parallel plane pair" "../bin/fd3d_test pair"
  check "3-D parallel plane pair (1)" "../bin/fd3d_test pair 1"
  # several 3-D planes of any orientation (G1): a fault split into two planes reproduces the
  # single plane (upright, dipping and turned), strike and dip slip; the dense coupling equals
  # the Toeplitz pair for parallel planes
  for c in "0 90 2 1" "35 -30 2 1" "0 90 2 2" "35 -30 2 2"; do
    check "3-D split planes ($c)" "../bin/fd3d_test multi $c"
  done
  check "3-D dense coupling against the Toeplitz pair" "../bin/fd3d_test multipair"
  # H6: the generic global-frame element response against the one-fault surface code, rotated
  check "3-D H6 frames" "../bin/fd3d_test h6frames"
  check "3-D H6 strips" "../bin/fd3d_test h6strips"
fi

# ---------------------------------------------------------------------------------------------
if want ruptures; then
  echo "== ruptures: test drivers against the prototype and earlier values"
  # TPV205-2D at 100 m against the prototype station series, direct and partitioned convolution
  for B in "" 16; do
    check "TPV205-2D (convolution ${B:-direct})" "FD_CONV_B=$B ../bin/fd_dynamic tpv205 100 12 ../$d/ref/ref_tpv205_dx100.txt"
  done
  # Harris and Day case D, two segments with normal-stress coupling, both step senses
  for W in 500 -500; do
    ref=$d/ref/ref_hd_D_W`echo $W | sed 's/-/m/'`_dx100.txt
    check "Harris and Day D (W $W)" "FD_CONV_B=16 ../bin/fd_dynamic hd $W 100 8 ../$ref"
  done
  # second-order marching (Stepper.enable_nl)
  check "Harris and Day D (W 500), second order" "FD_NL=1 ../bin/fd_dynamic hd 500 100 8 ../$d/ref/ref_hd_nl_W500_dx100.txt"
  check "TPV205-2D, second order" "FD_NL=1 ../bin/fd_dynamic tpv205 100 12 ../$d/ref/ref_tpv205_nl_dx100.txt"
  # 2-D in-plane half plane with exact kernels (M10, M11): spontaneous rupture on a 3 km fault
  # of 200 m dipping 30 degrees, buried (top 400 m) or reaching the surface, thrust and normal.
  # DWN tables (default): within 1 cm of the values of 2026-11-05 at FD_HPDWN_Q = 8 (the default
  # q = 4 is 6 mm lower at the top of the thrust reaching the surface).
  # Cagniard-de Hoop tables (FD_HP_METHOD=cdh): within 1 cm of the values of 2026-11-19; up to
  # 15 mm from DWN (bottom of the thrust reaching the surface), where DWN smooths the fronts.
  fsrup_check "half plane rupture DWN, buried thrust" -1 0.673 1.322 0.663 400 ""
  fsrup_check "half plane rupture DWN, buried normal" 1 0.229 0.916 0.563 400 ""
  fsrup_check "half plane rupture DWN, thrust to the surface" -1 2.305 2.087 0.898 0 ""
  fsrup_check "half plane rupture DWN, normal to the surface" 1 1.562 1.555 0.773 0 ""
  fsrup_check "half plane rupture CdH, buried thrust" -1 0.672 1.323 0.657 400 "FD_HP_METHOD=cdh"
  fsrup_check "half plane rupture CdH, buried normal" 1 0.229 0.916 0.563 400 "FD_HP_METHOD=cdh"
  fsrup_check "half plane rupture CdH, thrust to the surface" -1 2.300 2.090 0.883 0 "FD_HP_METHOD=cdh"
  fsrup_check "half plane rupture CdH, normal to the surface" 1 1.561 1.555 0.770 0 "FD_HP_METHOD=cdh"
  # 3-D: TPV3 (slip weakening) and TPV101 (rate and state, 6 s) at 500 m against the prototype
  check "3-D TPV3" "../bin/fd3d_test tpv3 500 ../$d/ref/ref_tpv3_dx500.txt"
  check "3-D TPV101" "../bin/fd3d_test tpv101 500 ../$d/ref/ref_tpv101_dx500.txt 6"
  # stability of the 3-D virtual surface per wavenumber (several minutes)
  check "3-D virtual surface stability" "../bin/fd3d_test surfstab 4"
  # 3-D mirror-image free surface against the explicit doubled fault
  check "3-D mirror surface" "../bin/fd3d_test mirror 500 6"
  check "3-D mirror surface, second order" "FD3D_NL=1 ../bin/fd3d_test mirror 500 6"
  # 3-D stepover smoke test: case-B strength, 500 m step, 500 m elements: segment 1 must rupture
  check "3-D stepover (segment 1 triggered)" "../bin/fd3d_test stepover -500 500 12 0.25 0.645 > so.out; grep -q 'segment 1 triggered' so.out && grep -q 'ruptured 2[0-9].[0-9] km (max slip 1.1' so.out && echo 'segment 1 triggered: PASS'"
  # several faults with one virtual surface (G3b): a dipping fault split along dip reproduces it
  for c in 1 2; do
    check "3-D several faults, one surface ($c)" "../bin/fd3d_test msurf $c 2"
  done
  # H6 check (1), short: two faults sharing a strike, partitioned strips against the several-fault
  # surface (G3b) on the same cells, prescribed slip (the full check: fd3d_test h6step 12)
  check "3-D H6 strips against the several-fault surface" "../bin/fd3d_test h6step 8 0"
  if have_python; then
    # D1: partitioned dense lag tables against all lags, TPV14 at 1 km over 10 s, exact dumps
    # within 1e-7 relative; H8: band block-Jacobi strip solve against dense LU within 1e-10
    check "3-D TPV14 partitioned tables against all lags" "FD3D_DP=0 FD3D_DUMP=dp0.txt ../bin/fd3d_test tpv14 1000 14 10 > /dev/null 2>&1 && FD3D_DP=1 FD3D_DUMP=dp1.txt ../bin/fd3d_test tpv14 1000 14 10 > /dev/null 2>&1 && python3 ../$d/py/fd3d_dumpdiff.py dp0.txt dp1.txt 1e-7"
    check "3-D strips, band preconditioned solve against LU" "FD_STRIPS_SOLVER=lu FD3D_DUMP=sv0.txt ../bin/fd3d_test h6step 12 2 > /dev/null 2>&1 && FD_STRIPS_SOLVER=bj FD3D_DUMP=sv1.txt ../bin/fd3d_test h6step 12 2 > /dev/null 2>&1 && python3 ../$d/py/fd3d_dumpdiff.py sv0.txt sv1.txt 1e-10"
  else
    skip "3-D dump comparisons (TPV14, strips)" "no python3"
  fi
fi

# ---------------------------------------------------------------------------------------------
if want rsf; then
  echo "== rsf: rsf_solve end to end"
  check "QD/FD sequence against the prototype" "cd .. && $d/e2e/test_fd_sequence.sh"
  check "2-D half plane static gate (DWN)" "cd .. && $d/e2e/test_fd_halfplane.sh"
  check "2-D half plane static gate (Cagniard-de Hoop)" "cd .. && FD_HP_METHOD=cdh $d/e2e/test_fd_halfplane.sh"
  check "3-D static gate" "cd .. && $d/e2e/test_fd3d_gate.sh"
  check "3-D parallel planes" "cd .. && $d/e2e/test_fd3d_pair.sh"
  check "3-D several planes" "cd .. && $d/e2e/test_fd3d_multi.sh"
  check "3-D several faults, one surface" "cd .. && $d/e2e/test_fd3d_msurf.sh"
  check "3-D dipping thrust, virtual surface" "cd .. && $d/e2e/test_fd3d_dip.sh"
fi

# ---------------------------------------------------------------------------------------------
if want mpi; then
  echo "== mpi: identical results on several ranks"
  find_mpi
  if [ -n "$mpi" ]; then
    echo "MPI launch command: $mpi"
    mkdir -p $tmp/mpi1 $tmp/mpi3
    # the 3-D marcher, bit-identical on 1 and 3 ranks (S2); the virtual surface split by
    # wavenumber (S3): TPV10 at 1 km for 3 s
    check "3-D marcher on 1 and 3 ranks" "cd mpi1 && FD3D_NL=1 ../../bin/fd3d_test tpv101 500 - 4 > out 2>/dev/null; cd ../mpi3 && FD3D_NL=1 $mpi -np 3 ../../bin/fd3d_test tpv101 500 - 4 > out 2>/dev/null; cd .. && cmp mpi1/tpv101_rt.dat mpi3/tpv101_rt.dat && echo identical"
    check "3-D virtual surface on 1 and 3 ranks" "cd mpi1 && ../../bin/fd3d_test tpv10 1000 30000 0 2 tv > /dev/null 2>&1; cd ../mpi3 && $mpi -np 3 ../../bin/fd3d_test tpv10 1000 30000 0 2 tv > /dev/null 2>&1; cd .. && cmp mpi1/tv_st12000.dat mpi3/tv_st12000.dat && cmp mpi1/tv_st00000.dat mpi3/tv_st00000.dat && echo identical"
    # dense couplings split by receiver rows (H7): several planes (30 degree bend, images),
    # several faults with one surface, strips of different strikes, TPV14 at 1 km; exact dumps
    # of the final state on 1, 2 and 4 ranks identical
    for c in "multi 30 90 2 1 1" "msurf 1 2" "h6step 8 2" "tpv14 1000 14 2"; do
      t=`echo $c | cut -d' ' -f1`
      check "3-D dense couplings, $t, on 1, 2 and 4 ranks" "for n in 1 2 4; do mkdir -p h7_$t/n\$n; ( cd h7_$t/n\$n && rm -f dump.txt && FD3D_DUMP=dump.txt $mpi -np \$n ../../../bin/fd3d_test $c > out 2>/dev/null ); done; [ -s h7_$t/n1/dump.txt ] && cmp h7_$t/n1/dump.txt h7_$t/n2/dump.txt && cmp h7_$t/n1/dump.txt h7_$t/n4/dump.txt && echo identical"
    done
    # rsf_solve: 2-D faults (M8, M10; curved fault and half plane) and 3-D (H7)
    check "rsf_solve 2-D faults on 1, 2 and 4 ranks" "cd .. && $d/e2e/test_fd_mpi2d.sh \"$mpi\""
    check "rsf_solve 2-D faults on 1, 2 and 4 ranks (CdH)" "cd .. && FD_HP_METHOD=cdh $d/e2e/test_fd_mpi2d.sh \"$mpi\""
    check "rsf_solve 3-D on 1, 2 and 4 ranks" "cd .. && $d/e2e/test_fd3d_mpi.sh \"$mpi\""
  else
    skip "MPI checks" "no MPI launch command that starts 4 ranks"
  fi
fi

echo "== $npass passed, $nfail failed, $nskip skipped ($((`date +%s` - t_start)) s; tiers $tiers)"
if [ $nfail -gt 0 ]; then
  echo "failed:$failed"
  exit 1
fi
exit 0
