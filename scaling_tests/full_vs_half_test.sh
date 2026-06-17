#!/bin/bash
#
# full_vs_half_test.sh
#
# Compare the interact Okada operator in HALF-space (standard, default) and
# FULL-space (infinite medium, -full_space 1) on a slightly bigger problem,
# using compress_interaction_matrix. This is meant as the staging ground for
# adding a native full-space H-matrix backend (BigWham, use_hmatrix 5): once
# that backend exists, it is full-space only and is validated against the
# native full-space Okada operator produced here with -full_space 1.
#
# What each run reports (compress_interaction_matrix is self-checking): it
# builds a dense reference with the SAME kernel/mode, then the requested
# H-matrix backend, and prints the relative forward error |b-b_h|/|b| of
# b = A x, the assembly time, and dense/H matvec timings.
#
# Two things this script highlights:
#   1) HALF vs FULL physics: the dense forward norm |b| (for a fixed unit-slip
#      pattern) differs between the two modes; near a free surface the
#      difference is large, which is exactly why a full-space backend must be
#      compared against a full-space reference, not the half-space one.
#   2) per backend, in each mode, the H-vs-dense self-consistency |b-b_h|/|b|.
#      In FULL mode this is "backend vs native full-space Okada (dense)", which
#      is the apples-to-apples metric the BigWham backend will be judged on.
#
# Backends (compress_interaction_matrix -use_hmatrix N):
#   0 dense (reference, always built internally)
#   1 HTOOL    (MPI)            2 H2OPUS (serial here; symmetric-only)
#   3 HACApK   (MPI)            4 hmmvp  (MPI)
#   5 BigWham  (full-space only; NOT YET COMPILED -> guarded, see ENABLE_BIGWHAM)
#
# Results are config-specific (this geometry, accuracy settings, machine);
# other settings may differ. Treat the numbers as relative comparisons within
# a single run, not absolute statements about any backend.
#
# usage:
#   ./full_vs_half_test.sh                 # defaults below
#   NFAULT=60 MFAULT=40 ./full_vs_half_test.sh
#   CORES=8 ./full_vs_half_test.sh         # single core count
#   ENABLE_BIGWHAM=1 ./full_vs_half_test.sh   # once use_hmatrix 5 is built in
#
set -u

# ---------------------------------------------------------------- configuration
BIN=${BIN:-../bin/compress_interaction_matrix}
MAKEFAULT=${MAKEFAULT:-../bin/makefault}
MPIRUN=${MPIRUN:-mpirun}
GEOM=${GEOM:-geom_fvh.in}
# "slightly bigger": 2 x NFAULT x MFAULT patches (two parallel fault planes).
# 2 x 40 x 25 = 2000 patches by default; bump for a real scaling test.
NFAULT=${NFAULT:-40}
MFAULT=${MFAULT:-25}
CORES=${CORES:-1}                 # space separated list, e.g. "1 4 16"; MPI ranks
NRANDOM=${NRANDOM:-50}            # matvecs timed per run
OUT=${OUT:-full_vs_half.dat}
CSV=${CSV:-full_vs_half.csv}

# roughly comparable accuracy across backends (same spirit as hmat_scaling_test.sh)
HTOOL_OPTS=${HTOOL_OPTS:-"-mat_htool_epsilon 3e-5 -mat_htool_eta 10"}
HACAPK_OPTS=${HACAPK_OPTS:-"-hacapk_ztol 1e-1"}
HMMVP_OPTS=${HMMVP_OPTS:-"-hmmvp_tol 1e-7"}
BIGWHAM_OPTS=${BIGWHAM_OPTS:-"-bigwham_eps_aca 1e-4 -bigwham_eta 3 -bigwham_leaf 32"}  # placeholder names

# active H-matrix backends to sweep (dense reference is always built by the tool).
# 2 (H2OPUS) omitted by default: serial + symmetric-only here. Add if wanted.
BACKENDS=${BACKENDS:-"1 3 4"}
# BigWham is full-space-only and not yet compiled; enable once use_hmatrix 5 works.
ENABLE_BIGWHAM=${ENABLE_BIGWHAM:-0}

opts_for(){ # backend -> accuracy options
    case "$1" in
        1) echo "$HTOOL_OPTS" ;;
        2) echo "" ;;
        3) echo "$HACAPK_OPTS" ;;
        4) echo "$HMMVP_OPTS" ;;
        5) echo "$BIGWHAM_OPTS" ;;
        *) echo "" ;;
    esac
}
name_for(){ case "$1" in 1) echo htool;; 2) echo h2opus;; 3) echo hacapk;; 4) echo hmmvp;; 5) echo bigwham;; *) echo "be$1";; esac; }

# ------------------------------------------------------------------- geometry
if [ ! -s "$GEOM" ]; then
    if [ ! -x "$MAKEFAULT" ]; then
        echo "$0: need $MAKEFAULT (build with: make bin/makefault)"; exit 1
    fi
    echo "$0: generating $GEOM with 2 x $NFAULT x $MFAULT patches"
    "$MAKEFAULT" -n "$NFAULT" -m "$MFAULT"                > "$GEOM"
    "$MAKEFAULT" -n "$NFAULT" -m "$MFAULT" -x 1 -grp 1   >> "$GEOM"
fi
N=$(grep -cve '^[[:space:]]*$' "$GEOM")
echo "$0: N = $N patches, cores = [$CORES], nrandom = $NRANDOM"

# ------------------------------------------------------------------ run helper
# runs one (mode, backend, cores) and echoes:  err asm mv_dense mv_H absb
run_one(){ # $1 mpiprefix  $2 use_hmatrix  $3 opts  $4 full_space(0/1)
    local out err asm td th absb fsflag
    fsflag=""; [ "$4" -eq 1 ] && fsflag="-full_space 1"
    out=$($1 "$BIN" -geom_file "$GEOM" -use_hmatrix "$2" $3 $fsflag -nrandom "$NRANDOM" 2>&1)
    err=$(echo  "$out" | grep "b-b_h"             | awk '{printf "%.2e",$NF}')
    absb=$(echo "$out" | grep -oE "\|b\| = +[0-9.eE+]+" | head -1 | awk '{print $3}')
    asm=$(echo  "$out" | grep "H matrix assembly" | awk '{printf "%.2f",$(NF-1)}')
    td=$(echo   "$out" | grep "dense solves"      | awk '{print $4}' | tr -d s)
    th=$(echo   "$out" | grep -E "H-matrix\(|Hmatrix\(" | awk '{print $4}' | tr -d s)
    echo "${err:-FAIL} ${asm:--} ${td:--} ${th:--} ${absb:--}"
}

: > "$OUT"
echo "mode,backend,np,N,relerr_H_vs_dense,asm_s,mv_dense_s,mv_H_s,absb" > "$CSV"

printf "%-5s %-8s %4s %6s  %-10s %9s %11s %9s  %-16s\n" \
       mode backend np N relerr asm_s mv_dense_s mv_H_s "|b|(dense fwd)" | tee -a "$OUT"
printf -- '------------------------------------------------------------------------------------------\n' | tee -a "$OUT"

ABSB_HALF=""; ABSB_FULL=""

for np in $CORES; do
    pre=""; [ "$np" -gt 1 ] && pre="$MPIRUN -np $np"
    for mode in half full; do
        fs=0; [ "$mode" = full ] && fs=1
        # build the active backend list for this mode
        list="$BACKENDS"
        if [ "$mode" = full ] && [ "$ENABLE_BIGWHAM" -eq 1 ]; then list="$list 5"; fi
        for be in $list; do
            # BigWham is full-space only: never run it in half mode
            if [ "$be" -eq 5 ] && [ "$mode" = half ]; then continue; fi
            read err asm td th absb < <(run_one "$pre" "$be" "$(opts_for "$be")" "$fs")
            printf "%-5s %-8s %4d %6d  %-10s %9s %11s %9s  %-16s\n" \
                   "$mode" "$(name_for "$be")" "$np" "$N" "$err" "$asm" "$td" "$th" "$absb" | tee -a "$OUT"
            echo "$mode,$(name_for "$be"),$np,$N,$err,$asm,$td,$th,$absb" >> "$CSV"
            # remember the dense forward norm once per mode (any backend run prints it)
            if [ "$mode" = half ] && [ -z "$ABSB_HALF" ] && [ "$absb" != "-" ]; then ABSB_HALF="$absb"; fi
            if [ "$mode" = full ] && [ -z "$ABSB_FULL" ] && [ "$absb" != "-" ]; then ABSB_FULL="$absb"; fi
        done
    done
done

# ------------------------------------------------------ half vs full physics
echo "" | tee -a "$OUT"
if [ -n "$ABSB_HALF" ] && [ -n "$ABSB_FULL" ]; then
    rel=$(awk -v h="$ABSB_HALF" -v f="$ABSB_FULL" 'BEGIN{ if(h+0!=0) printf "%.3f", (f-h)/h; else print "NA"}')
    echo "HALF vs FULL (dense forward |b| for the same unit-slip pattern):" | tee -a "$OUT"
    echo "   |b|_half = $ABSB_HALF   |b|_full = $ABSB_FULL   (|b|_full-|b|_half)/|b|_half = $rel" | tee -a "$OUT"
    echo "   -> nonzero relative difference is the free-surface contribution that the" | tee -a "$OUT"
    echo "      full-space backend (BigWham) must NOT be expected to reproduce." | tee -a "$OUT"
fi

echo "" | tee -a "$OUT"
echo "NOTE: when use_hmatrix 5 (BigWham) is compiled in, run with ENABLE_BIGWHAM=1." | tee -a "$OUT"
echo "      In FULL mode the 'relerr' column for bigwham is then exactly the" | tee -a "$OUT"
echo "      apples-to-apples error between BigWham (full-space) and the native" | tee -a "$OUT"
echo "      full-space Okada dense reference. BigWham is never run in HALF mode." | tee -a "$OUT"
echo "wrote $OUT and $CSV"
