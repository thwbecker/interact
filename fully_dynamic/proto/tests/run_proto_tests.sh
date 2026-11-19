#!/bin/sh
# Tests of the Python prototypes (fully_dynamic/proto). Files with pytest functions run under
# pytest, the others as scripts (they print PASS or FAIL and return non-zero on failure).
# usage: fully_dynamic/proto/tests/run_proto_tests.sh [quick|all] (default quick: without the
# sequence and stepover regressions test_cycle.py and test_hd.py, about 2 minutes each)
# Needs python3 with numpy, scipy, sympy and pytest; test_stepper.py's TPV205 case needs the
# untracked reference data in fully_dynamic/proto/ref (see ../README.md) and is skipped without it;
# the half-plane checks need the tdd library (TDD_LIB, see src/testing/rsf_solve/run_tests.sh).
tier=${1-quick}
dir=`dirname $0`
fail=0
for f in $dir/test_*.py; do
  b=`basename $f`
  if [ "$tier" = "quick" ] && { [ "$b" = "test_cycle.py" ] || [ "$b" = "test_hd.py" ]; }; then
    echo "$b: skipped (quick)"
    continue
  fi
  case $b in
    test_fs3d_*)
      if [ ! -s ${DC3D_LIB-/tmp/okd/libdc3d.so} ]; then
        echo "$b: skipped (no Okada library DC3D_LIB, default /tmp/okd/libdc3d.so; build: see fully_dynamic/proto/fs3d_static.py)"
        continue
      fi;;
    test_halfplane_*)
      if [ ! -s ${TDD_LIB-/tmp/tdd/libtdd.so} ]; then
        echo "$b: skipped (no tdd library TDD_LIB, see src/testing/rsf_solve/run_tests.sh)"
        continue
      fi;;
  esac
  t0=`date +%s`
  if grep -q '^def test_' $f; then
    python3 -m pytest -q -p no:cacheprovider $f > /tmp/proto_test.out 2>&1
  else
    python3 $f > /tmp/proto_test.out 2>&1
  fi
  st=$?
  t1=`date +%s`
  if [ $st -eq 0 ]; then
    echo "$b: PASS ($((t1 - t0)) s) $(tail -1 /tmp/proto_test.out | cut -c1-80)"
  else
    echo "$b: FAIL ($((t1 - t0)) s)"
    tail -5 /tmp/proto_test.out | sed 's/^/    /'
    fail=1
  fi
done
[ $fail -eq 0 ] && echo "prototype tests: all PASS" || echo "prototype tests: FAIL"
exit $fail
