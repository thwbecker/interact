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
( cd tmp_fd_dynamic && FD3D_NL=1 ../bin/fd3d_test mirror 500 6 > fd3d.out 2>/dev/null ); st=$?
echo "(second-order marching) $(tail -1 tmp_fd_dynamic/fd3d.out)"; [ $st -eq 0 ] || fail=1
[ $fail -eq 0 ] && echo "rsf_solve fd tests: all PASS" || echo "rsf_solve fd tests: FAIL"
exit $fail
