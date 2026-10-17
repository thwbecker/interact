#!/bin/sh
# unit tests of the fully dynamic extension; run from the repository root after make rsf_testing
d=src/testing/rsf_solve
fail=0
bin/test_fd_fft || fail=1
bin/test_fd_segments || fail=1
for f in $d/ref/ref_*.txt; do bin/test_fd_tables $f || fail=1; done
[ $fail -eq 0 ] && echo "rsf_solve fd tests: all PASS" || echo "rsf_solve fd tests: FAIL"
exit $fail
