#!/bin/bash
#
# bigwham_vs_okada_test.sh
#
# Correctness / stress test for the BigWham full-space backend (use_hmatrix 5)
# against the native full-space Okada operator (the dense reference that
# compress_interaction_matrix builds internally with -full_space 1).
#
# It makes faults with makefault, assembles the interaction matrix through both
# BigWham and Okada, and compares the two operators. The comparison metric is
# compress's own forward test, |b - b_h| / |b| for b = A_okada x and
# b_h = A_bigwham x on a random slip vector x. With a sign or local-frame error
# this is O(0.1) to O(2); when the two operators agree it is at the level of the
# requested ACA tolerance (and at machine precision when no compression happens).
# A complementary inverse (GMRES) test is also reported on one case.
#
# BigWham is an infinite-medium kernel, so every comparison here is run with
# -full_space 1; comparing BigWham to the Okada half-space operator would not be
# apples to apples and is never done.
#
# Usage (positional arguments, all optional, with defaults):
#   ./bigwham_vs_okada_test.sh [strike] [dip] [n] [m] [eps_aca] [nrandom]
#
#   strike   patch strike, deg clockwise from N        (default 30)
#   dip      patch dip,    deg down from horizontal     (default 60)
#   n        patches along strike for the planar mesh   (default 8)
#   m        patches along dip    for the planar mesh   (default 6)
#   eps_aca  BigWham ACA tolerance                       (default 1e-4)
#   nrandom  matvecs used for the timing loop            (default 100)
#
# Examples:
#   ./bigwham_vs_okada_test.sh                 # strike 30, dip 60, 8x6 mesh
#   ./bigwham_vs_okada_test.sh 0   90          # vertical N-S fault
#   ./bigwham_vs_okada_test.sh 170 30 12 10    # oblique strike, shallow dip, bigger
#
set -u

strike=${1:-30}
dip=${2:-60}
n=${3:-8}
m=${4:-6}
eps_aca=${5:-1e-4}
nrandom=${6:-10}

# ---- configuration (edit these paths if you run from elsewhere) -------------
BIN=${BIN:-../bin/compress_interaction_matrix}
MAKEFAULT=${MAKEFAULT:-../bin/makefault}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# accuracy threshold for the PASS/FAIL flag: a real frame/sign bug shows up as
# O(0.1)..O(2), far above this, while a correct operator agrees to the ACA
# tolerance or better. Floor at 1e-9 for the no-compression (tiny N) cases.
pass_tol=$(awk -v e="$eps_aca" 'BEGIN{t=50*e; if(t<1e-9)t=1e-9; printf "%.3e",t}')

for exe in "$BIN" "$MAKEFAULT"; do
    if [ ! -x "$exe" ]; then
        echo "error: cannot find/execute $exe"
        echo "       build with: make bin/compress_interaction_matrix bin/makefault"
        exit 1
    fi
done

echo "================================================================================"
echo " BigWham (full-space) vs native full-space Okada"
echo "   strike=$strike dip=$dip   planar mesh n x m = $n x $m   eps_aca=$eps_aca"
echo "   PASS if |b-b_h|/|b| < $pass_tol"
echo "================================================================================"
printf "%-46s %6s %-12s %-9s %s\n" "case" "N" "relerr" "compress" "verdict"
printf -- "--------------------------------------------------------------------------------\n"

# run one comparison: $1 = geom file, $2 = label, $3 = extra compress opts
run_cmp(){
    local geom="$1" label="$2" extra="${3:-}" out relerr comp N verdict
    N=$(grep -cve '^[[:space:]]*$' "$geom")
    out=$("$BIN" -geom_file "$geom" -use_hmatrix 5 -full_space 1 \
                 -bigwham_eps_aca "$eps_aca" -nrandom "$nrandom" $extra 2>&1)
    relerr=$(echo "$out" | grep "b-b_h" | awk '{print $NF}')
    comp=$(echo "$out" | grep -oE "compression ratio [0-9.eE+-]+" | head -1 | awk '{print $3}')
    if [ -z "$relerr" ]; then
        printf "%-46s %6s %-12s %-9s %s\n" "$label" "$N" "NO-OUTPUT" "${comp:--}" "FAIL"
        echo "$out" | grep -iE "signal|error|not compiled" | head -2 | sed 's/^/      /'
        return
    fi
    verdict=$(awk -v r="$relerr" -v t="$pass_tol" 'BEGIN{print (r+0<=t+0)?"PASS":"FAIL"}')
    printf "%-46s %6s %-12s %-9s %s\n" "$label" "$N" "$relerr" "${comp:--}" "$verdict"
}

# ---------------------------------------------------------------------------
# 1. single isolated patch at the requested orientation (self term, no
#    compression: this isolates the local-frame mapping for this strike/dip)
"$MAKEFAULT" -strike "$strike" -dip "$dip" -n 1 -m 1 -z -5 > "$WORK/single.in" 2>/dev/null
run_cmp "$WORK/single.in" "single patch (strike $strike, dip $dip)"

# 2. planar fault, uniform orientation, n x m patches (adds the cross terms
#    between equally oriented patches)
"$MAKEFAULT" -strike "$strike" -dip "$dip" -n "$n" -m "$m" -z -5 > "$WORK/plane.in" 2>/dev/null
run_cmp "$WORK/plane.in" "planar fault (uniform $strike/$dip)"

# 3. per-patch randomized orientation around the requested strike/dip: every
#    patch gets a different local frame, the strongest test of the corner
#    ordering -> local-frame mapping across many orientations in one mesh
"$MAKEFAULT" -strike "$strike" -dip "$dip" -n "$n" -m "$m" -z -5 \
             -srand 25 -drand 25 -seed 1 > "$WORK/rand.in" 2>/dev/null
run_cmp "$WORK/rand.in" "randomized orientations (+-25 deg)"

# 4. fixed nine-orientation panel: strike {0,60,170} x dip {0,30,90}, spaced out
{
  k=0
  for s in 0 60 170; do for d in 0 30 90; do
      x=$(awk -v k="$k" 'BEGIN{printf "%.3f", 3.0*int(k/3)}')
      y=$(awk -v k="$k" 'BEGIN{printf "%.3f", 3.0*(k%3)}')
      printf "%s %s -5.0 %s %s 0.5 0.5 0\n" "$x" "$y" "$s" "$d"
      k=$((k+1))
  done; done
} > "$WORK/panel9.in"
run_cmp "$WORK/panel9.in" "9-orientation panel (cross terms)"

# 5. full-space depth invariance: the same isolated patch placed shallow and
#    deep must give the same operator (the free surface is suppressed). This
#    also re-checks BigWham vs Okada at each depth.
"$MAKEFAULT" -strike "$strike" -dip "$dip" -n 1 -m 1 -z -1   > "$WORK/shallow.in" 2>/dev/null
"$MAKEFAULT" -strike "$strike" -dip "$dip" -n 1 -m 1 -z -200 > "$WORK/deep.in"    2>/dev/null
run_cmp "$WORK/shallow.in" "single patch at z=-1   (full-space)"
run_cmp "$WORK/deep.in"    "single patch at z=-200 (full-space)"

printf -- "--------------------------------------------------------------------------------\n"

# ---------------------------------------------------------------------------
# 6. ACA tolerance sweep: the error should track eps_aca. ACA only compresses
#    once the problem is large enough, so this uses a dedicated larger mesh
#    (compression ratio near or below 1 at small N just means the H format is
#    not yet smaller than dense; the accuracy trend is the point here).
echo ""
sn=30; sm=20
"$MAKEFAULT" -strike "$strike" -dip "$dip" -n "$sn" -m "$sm" -z -5 > "$WORK/sweep.in" 2>/dev/null
sN=$(grep -cve '^[[:space:]]*$' "$WORK/sweep.in")
echo "ACA tolerance sweep on a $sn x $sm = $sN patch mesh (error should track eps_aca):"
printf "%-14s %-12s %-10s\n" "eps_aca" "relerr" "compression"
for tol in 1e-2 1e-4 1e-6; do
    out=$("$BIN" -geom_file "$WORK/sweep.in" -use_hmatrix 5 -full_space 1 \
                 -bigwham_eps_aca "$tol" -nrandom 1 2>&1)
    r=$(echo "$out" | grep "b-b_h" | awk '{print $NF}')
    c=$(echo "$out" | grep -oE "compression ratio [0-9.eE+-]+" | head -1 | awk '{print $3}')
    printf "%-14s %-12s %-10s\n" "$tol" "${r:-FAIL}" "${c:--}"
done

# ---------------------------------------------------------------------------
# 7. inverse (GMRES) test on the planar mesh: complementary check that A^{-1}
#    agrees, not just the forward operator. The relerr line is printed by the
#    tool only when nrandom = 0, so this run uses that.
echo ""
echo "Inverse (GMRES) test on the planar mesh (|x-x_h|/|x|, KSP-tolerance level):"
"$BIN" -geom_file "$WORK/plane.in" -use_hmatrix 5 -full_space 1 \
       -bigwham_eps_aca "$eps_aca" -nrandom 0 -test_forward 0 2>&1 \
    | grep -E "x-x_h" | sed 's/^/   /' | head -1

# ---------------------------------------------------------------------------
# 8. half-space guard: without -full_space the tool should warn that the
#    comparison is not apples to apples
echo ""
echo "Half-space guard (running use_hmatrix 5 without -full_space should warn):"
"$BIN" -geom_file "$WORK/single.in" -use_hmatrix 5 -nrandom 1 2>&1 \
    | grep -i "full-space only" | sed 's/^/   /' | head -1

echo ""
echo "done."
