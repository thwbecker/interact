#!/bin/bash
#
# run_tdbench.sh
#
# Builds the triangular-dislocation kernels tdstresshs.f90 / tddisphs.f90 plus
# tdbench_driver.F90 under several optimisation regimes, times each, and
# reports the speedup and the numerical deviation from the strict -O2 build.
#
# The point of the experiment: the value-preserving common-subexpression work
# already in the sources is something -O2 does on its own, so the interesting
# question is the reciprocal/reassociation transforms a strict compiler will
# NOT do (x/rb -> x*(1/rb), etc.). -freciprocal-math is the compiler doing that
# reduction automatically; -ffast-math adds reassociation on top. The accuracy
# columns say whether either is safe for these artefact-prone kernels.
#
# All parameters are hardcoded plain assignments below (no environment-variable
# overrides). Run from the directory that holds the three source files.

# ============================ CONFIG =======================================
FC=gfortran
CC=gcc
SRC=.                       # dir with tdstresshs.f90, tddisphs.f90, tdbench_driver.F90
INC=../includes             # interact headers (precision_double.h, properties.h)
WDIR=tdbench_build          # scratch / output dir
REPS=3                      # timing repetitions per regime; the minimum is kept
COMMON="-cpp -ffree-line-length-none -DUSE_DOUBLE_PRECISION -I$INC"

# regime name : extra flags (strict O2 is the accuracy + timing reference)
REG_NAMES="O2strict O3 O3recip O3fast"
flags_for(){
  case "$1" in
    O2strict) echo "-O2" ;;
    O3)       echo "-O3 -funroll-loops" ;;
    O3recip)  echo "-O3 -funroll-loops -freciprocal-math -fno-math-errno -fno-trapping-math" ;;
    O3fast)   echo "-O3 -funroll-loops -ffast-math" ;;
    *) echo "ERROR unknown regime $1" >&2; exit 1 ;;
  esac
}
# ===========================================================================

mkdir -p "$WDIR"
for f in "$SRC/tdstresshs.f90" "$SRC/tddisphs.f90" "$SRC/tdbench_driver.F90" "$SRC/mysincos_shim.c" "$INC/precision_double.h"; do
  [ -e "$f" ] || { echo "ERROR: missing $f (check CONFIG)"; exit 1; }
done

printf "%-9s %10s %10s %10s %9s %14s %14s %14s\n" \
       regime stress_s disp_s total_s speedup max_abs_dev max_rel_dev rms_dev

REF="$WDIR/res_O2strict.dat"
declare -A TOTAL
ref_total=""

for reg in $REG_NAMES; do
  ff=$(flags_for "$reg")
  # compile from the source dir (so -I$INC stays correct), objects into $WDIR
  $FC $COMMON $ff -c "$SRC/tdstresshs.f90"     -o "$WDIR/tdstresshs.o" && \
  $FC $COMMON $ff -c "$SRC/tddisphs.f90"       -o "$WDIR/tddisphs.o"   && \
  $FC $COMMON $ff -c "$SRC/tdbench_driver.F90" -o "$WDIR/driver.o"     && \
  $CC -O2 -c "$SRC/mysincos_shim.c" -o "$WDIR/mysincos_shim.o"          && \
  $FC $ff "$WDIR/driver.o" "$WDIR/tdstresshs.o" "$WDIR/tddisphs.o" "$WDIR/mysincos_shim.o" -lm -o "$WDIR/bench" \
    || { echo "ERROR: build failed for regime $reg"; exit 1; }

  # run REPS times, keep the minimum total wall time; save the result dump once
  best=""
  st=""; dt=""; ck=""
  for k in $(seq 1 $REPS); do
    out=$( cd "$WDIR" && ./bench )
    s=$(echo "$out" | awk '/STRESS_TIME_S/{print $2}')
    d=$(echo "$out" | awk '/DISP_TIME_S/{print $2}')
    tot=$(awk -v a="$s" -v b="$d" 'BEGIN{printf "%.5f", a+b}')
    if [ -z "$best" ] || awk -v t="$tot" -v b="$best" 'BEGIN{exit !(t<b)}'; then
      best="$tot"; st="$s"; dt="$d"
      ck=$(echo "$out" | awk '/CHECKSUM/{print $2}')
    fi
  done
  cp "$WDIR/td_results.dat" "$WDIR/res_${reg}.dat"
  TOTAL[$reg]="$best"
  [ "$reg" = "O2strict" ] && ref_total="$best"

  # deviation vs the strict reference
  if [ "$reg" = "O2strict" ]; then
    devs="0 0 0"
  else
    devs=$(python3 - "$REF" "$WDIR/res_${reg}.dat" <<'PY'
import sys
ref=[float(x) for line in open(sys.argv[1]) for x in line.split()]
tst=[float(x) for line in open(sys.argv[2]) for x in line.split()]
mxa=mxr=0.0; ss=0.0; n=0
for a,b in zip(ref,tst):
    d=abs(a-b); mxa=max(mxa,d)
    if abs(a)>0: mxr=max(mxr,d/abs(a))
    ss+=d*d; n+=1
rms=(ss/n)**0.5 if n else 0.0
print(f"{mxa:.6e} {mxr:.6e} {rms:.6e}")
PY
)
  fi
  set -- $devs
  spd=$(awk -v r="$ref_total" -v t="$best" 'BEGIN{ if(t>0) printf "%.3f", r/t; else printf "NA"}')
  printf "%-9s %10s %10s %10s %9s %14s %14s %14s\n" "$reg" "$st" "$dt" "$best" "$spd" "$1" "$2" "$3"
done

echo ""
echo "checksums should agree to round-off across regimes (else a result changed):"
for reg in $REG_NAMES; do
  ( cd "$WDIR" && ./bench >/dev/null 2>&1 )
done
echo "results dumps: $WDIR/res_*.dat   (15 cols: stress 1-6, strain 7-12, disp 13-15)"
