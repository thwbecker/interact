#!/bin/sh
# unit tests of the fully dynamic extension; run from the repository root after make rsf_testing
d=src/testing/rsf_solve
fail=0
bin/test_fd_fft || fail=1
bin/test_fd_segments || fail=1
for f in $d/ref/ref_inplane_*.txt $d/ref/ref_antiplane_*.txt; do bin/test_fd_tables $f || fail=1; done
mkdir -p tmp_fd_dynamic
# TPV205-2D at 100 m against the prototype station series, direct and partitioned convolution (a few s each)
for B in "" 16; do
  ( cd tmp_fd_dynamic && FD_CONV_B=$B ../bin/fd_dynamic tpv205 100 12 ../$d/ref/ref_tpv205_dx100.txt > fd_dynamic.out 2>/dev/null ); st=$?
  tail -1 tmp_fd_dynamic/fd_dynamic.out | sed "s/:/ (partition ${B:-direct}):/"
  [ $st -eq 0 ] || fail=1
done
# Harris and Day case D, two segments with normal-stress coupling, both step senses, against prototype traces
for W in 500 -500; do
  ref=$d/ref/ref_hd_D_W`echo $W | sed 's/-/m/'`_dx100.txt
  ( cd tmp_fd_dynamic && FD_CONV_B=16 ../bin/fd_dynamic hd $W 100 8 ../$ref > hd.out 2>/dev/null ); st=$?
  tail -1 tmp_fd_dynamic/hd.out | sed "s/:/ (W $W):/"
  [ $st -eq 0 ] || fail=1
done
# 2-D second-order marching with normal-stress coupling: Harris-Day W = 500 m at 100 m against the prototype
( cd tmp_fd_dynamic && FD_NL=1 ../bin/fd_dynamic hd 500 100 8 ../$d/ref/ref_hd_nl_W500_dx100.txt > nlhd.out 2>&1 ); st=$?
echo "(second-order marching) $(tail -1 tmp_fd_dynamic/nlhd.out | sed 's/:/ (W 500):/')"; [ $st -eq 0 ] || fail=1
# 2-D second-order marching: TPV205 at 100 m against the prototype (Stepper.enable_nl)
( cd tmp_fd_dynamic && FD_NL=1 ../bin/fd_dynamic tpv205 100 12 ../$d/ref/ref_tpv205_nl_dx100.txt > nl.out 2>&1 ); st=$?
echo "(second-order marching) $(tail -1 tmp_fd_dynamic/nl.out)"; [ $st -eq 0 ] || fail=1
# 3-D: TPV3 (slip weakening) and TPV101 (rate and state, 6 s) at 500 m against the prototype (about 15 s)
( cd tmp_fd_dynamic && ../bin/fd3d_test tpv3 500 ../$d/ref/ref_tpv3_dx500.txt > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && ../bin/fd3d_test tpv101 500 ../$d/ref/ref_tpv101_dx500.txt 6 > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
# 3-D mirror-image free surface against the explicit doubled fault (a few s)
( cd tmp_fd_dynamic && ../bin/fd3d_test mirror 500 6 > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && ../bin/fd3d_test kernels ../$d/ref/ref_fd3d_kernels.txt > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && ../bin/fd3d_test layer ../$d/ref/ref_fd3d_layer.txt > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && FD3D_LAYER_FFT=1 ../bin/fd3d_test layer ../$d/ref/ref_fd3d_layer.txt > fd3d.out 2>/dev/null ); st=$?
echo "(FFT layer) $(tail -1 tmp_fd_dynamic/fd3d.out)"; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && ../bin/fd3d_test surf ../$d/ref/ref_fd3d_surf.txt > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && ../bin/fd3d_test surfstab 4 > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && ../bin/fd3d_test q2d ../$d/ref/ref_fd3d_q2d.txt > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && ../bin/fd3d_test pair > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
( cd tmp_fd_dynamic && ../bin/fd3d_test pair 1 > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
# 3-D stepover smoke test: case-B strength, 500 m step, 500 m elements: segment 1 must rupture
( cd tmp_fd_dynamic && ../bin/fd3d_test stepover -500 500 12 0.25 0.645 > fd3d.out 2>/dev/null )
if grep -q "segment 1 triggered" tmp_fd_dynamic/fd3d.out && grep -q "ruptured 2[0-9].[0-9] km (max slip 1.1" tmp_fd_dynamic/fd3d.out; then echo "fd3d_test stepover: PASS"
else echo "fd3d_test stepover: FAIL"; cat tmp_fd_dynamic/fd3d.out; fail=1; fi
# several 3-D planes of any orientation (fd_plan.md G1): a fault split into two planes must
# reproduce the single plane (upright, and rigidly rotated: dipping and turned), strike and
# dip slip; the dense coupling must equal the Toeplitz pair for parallel planes
for c in "0 90 2 1" "35 -30 2 1" "0 90 2 2" "35 -30 2 2"; do
  ( cd tmp_fd_dynamic && ../bin/fd3d_test multi $c > fd3d.out 2>/dev/null ); st=$?
  tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
done
( cd tmp_fd_dynamic && ../bin/fd3d_test multipair > fd3d.out 2>/dev/null ); st=$?
tail -1 tmp_fd_dynamic/fd3d.out; [ $st -eq 0 ] || fail=1
# non-parallel 2-D pairs (fd_plan.md N1): dense tables against the prototype for a kinked pair
if command -v python3 > /dev/null 2>&1; then
  for a in 30 -45; do for m in 0 1; do
    ( cd tmp_fd_dynamic && ../bin/test_fd_dense $a $m dense_${a}_${m}.txt > /dev/null 2>&1 && \
      python3 ../fully_dynamic/test_fd_dense.py $a $m dense_${a}_${m}.txt > fd3d.out 2>&1 ); st=$?
    tail -1 tmp_fd_dynamic/fd3d.out; grep -q PASS tmp_fd_dynamic/fd3d.out || fail=1
  done; done
fi
# MPI (fd_plan.md S2): the 3-D marcher must give bit-identical results on 1 and 3 ranks
if command -v mpirun > /dev/null 2>&1; then
  mkdir -p tmp_fd_dynamic/mpi1 tmp_fd_dynamic/mpi3
  ( cd tmp_fd_dynamic/mpi1 && FD3D_NL=1 ../../bin/fd3d_test tpv101 500 - 4 > out 2>/dev/null )
  ( cd tmp_fd_dynamic/mpi3 && FD3D_NL=1 mpirun --oversubscribe -np 3 ../../bin/fd3d_test tpv101 500 - 4 > out 2>/dev/null \
      || FD3D_NL=1 mpirun --allow-run-as-root --oversubscribe -np 3 ../../bin/fd3d_test tpv101 500 - 4 > out 2>/dev/null )
  # the virtual surface (tables split by wavenumber, S3): TPV10 at 1 km for 3 s
  ( cd tmp_fd_dynamic/mpi1 && ../../bin/fd3d_test tpv10 1000 30000 0 2 tv > /dev/null 2>&1 )
  ( cd tmp_fd_dynamic/mpi3 && ( mpirun --oversubscribe -np 3 ../../bin/fd3d_test tpv10 1000 30000 0 2 tv > /dev/null 2>&1 || \
      mpirun --allow-run-as-root --oversubscribe -np 3 ../../bin/fd3d_test tpv10 1000 30000 0 2 tv > /dev/null 2>&1 ) )
  if cmp -s tmp_fd_dynamic/mpi1/tv_st12000.dat tmp_fd_dynamic/mpi3/tv_st12000.dat && cmp -s tmp_fd_dynamic/mpi1/tv_st00000.dat tmp_fd_dynamic/mpi3/tv_st00000.dat; then echo "fd3d_test MPI surface (1 and 3 ranks bit-identical): PASS"
  else echo "fd3d_test MPI surface (1 and 3 ranks bit-identical): FAIL"; fail=1; fi
  if cmp -s tmp_fd_dynamic/mpi1/tpv101_rt.dat tmp_fd_dynamic/mpi3/tpv101_rt.dat; then echo "fd3d_test MPI (1 and 3 ranks bit-identical): PASS"
  else echo "fd3d_test MPI (1 and 3 ranks bit-identical): FAIL"; fail=1; fi
else
  echo "fd3d_test MPI: skipped (no mpirun)"
fi
( cd tmp_fd_dynamic && FD3D_NL=1 ../bin/fd3d_test mirror 500 6 > fd3d.out 2>/dev/null ); st=$?
echo "(second-order marching) $(tail -1 tmp_fd_dynamic/fd3d.out)"; [ $st -eq 0 ] || fail=1
[ $fail -eq 0 ] && echo "rsf_solve fd tests: all PASS" || echo "rsf_solve fd tests: FAIL"
exit $fail
