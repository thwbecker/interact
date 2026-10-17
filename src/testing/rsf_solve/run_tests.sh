#!/bin/sh
# unit tests of the fully dynamic extension; run from the repository root after make rsf_testing
d=src/testing/rsf_solve
fail=0
bin/test_fd_fft || fail=1
bin/test_fd_segments || fail=1
for f in $d/ref/ref_inplane_*.txt $d/ref/ref_antiplane_*.txt; do bin/test_fd_tables $f || fail=1; done
# TPV205-2D at 100 m against the prototype station series (about 2 s)
mkdir -p tmp_fd_dynamic
( cd tmp_fd_dynamic && ../bin/fd_dynamic tpv205 100 12 ../$d/ref/ref_tpv205_dx100.txt > fd_dynamic.out 2>/dev/null ); st=$?
tail -2 tmp_fd_dynamic/fd_dynamic.out
[ $st -eq 0 ] || fail=1
[ $fail -eq 0 ] && echo "rsf_solve fd tests: all PASS" || echo "rsf_solve fd tests: FAIL"
exit $fail
