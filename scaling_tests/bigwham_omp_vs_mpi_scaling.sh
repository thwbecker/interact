#!/bin/bash
#
# bigwham_omp_vs_mpi_scaling.sh
#
# Full-space assembly and matvec scaling, comparing:
#   - BigWham  (use_hmatrix 5): OpenMP-parallel, full-space native kernel
#   - HTOOL    (use_hmatrix 1): MPI-parallel,    Okada full-space (-full_space 1)
#   - HACApK   (use_hmatrix 3): MPI-parallel,    Okada full-space
#   - HMMVP    (use_hmatrix 4): MPI-parallel,    Okada full-space
#
# BigWham exposes two independent thread knobs and they must be set separately:
#   - assembly is driven by the ambient OMP_NUM_THREADS
#   - the matvec is driven by -bigwham_nthreads (its internal n_openMP_threads)
# BigWham's threaded matvec uses a private per-thread accumulator plus a final
# reduction, so on a small or well-compressed operator it stops scaling and can
# even slow down once the thread count is high (the reduction and the per-thread
# 3N allocation overtake the useful work). Tying the two knobs together hides
# this, so this script measures them in two separate phases:
#
#   PHASE 1  assembly scaling: vary the assembly parallelism P
#            (OMP threads for BigWham, MPI ranks for the Okada backends).
#            BigWham's matvec is pinned to 1 thread here so it does not pollute
#            the wall time; we only read the assembly number.
#
#   PHASE 2  matvec scaling: BigWham is assembled once with many OMP threads
#            (fast), then its matvec is timed while sweeping -bigwham_nthreads.
#            The Okada backends are run under P MPI ranks and their matvec is
#            timed. This is the number that matters for matvec-bound cycle runs.
#
# The Okada full-space dense reference is built every run for the relerr column,
# so each backend's relerr is against the same native full-space Okada operator
# (apples to apples, BigWham is also full-space). That dense reference is O(N^2)
# in time and memory, so it, not the H assembly, caps how large N can go here.
#
# Timings mean something only with at least P real cores. Results are specific to
# this geometry, these tolerances, and this machine; a different BigWham build or
# a different HACApK / HTOOL / hmmvp version may behave differently.
#
# Usage (positional arguments, all optional):
#   ./bigwham_omp_vs_mpi_scaling.sh [n] [m] [procs] [nrandom]
#
#   n        patches along strike                       (default 30)
#   m        patches along dip                          (default 20)
#   procs    quoted list of parallelism levels          (default "1 2 4")
#   nrandom  matvecs in the timing loop (phase 2)        (default 50)
#
# Examples:
#   ./bigwham_omp_vs_mpi_scaling.sh                       # laptop, quick
#   ./bigwham_omp_vs_mpi_scaling.sh 60 40 "1 2 4 8 16 32" 100   # big node
#
set -u

n=${1:-30}
m=${2:-20}
procs=${3:-"1 2 4"}
nrandom=${4:-50}

# ---- configuration (edit for your machine) ----------------------------------
BIN=${BIN:-../bin/compress_interaction_matrix}
MAKEFAULT=${MAKEFAULT:-../bin/makefault}
# launcher for the MPI (Okada) backends. On a laptop with few cores use
# "mpirun --oversubscribe"; on a cluster use your scheduler launcher
# (ibrun, srun, ...). BigWham always runs on one rank and does not use this.
MPIRUN=${MPIRUN:-mpirun}
STRIKE=${STRIKE:-30}
DIP=${DIP:-70}
# fixed OMP thread count used to assemble BigWham during the phase-2 matvec
# sweep (we want assembly fast and out of the way there). Default: the largest
# entry in procs.
ASM_THREADS=${ASM_THREADS:-$(for p in $procs; do echo "$p"; done | sort -n | tail -1)}
# roughly comparable accuracy across backends; tune per backend as you like.
HTOOL_OPTS=${HTOOL_OPTS:-"-mat_htool_epsilon 1e-4 -mat_htool_eta 10"}
HACAPK_OPTS=${HACAPK_OPTS:-"-hacapk_ztol 1e-2"}
HMMVP_OPTS=${HMMVP_OPTS:-"-hmmvp_tol 1e-4"}
BIGWHAM_OPTS=${BIGWHAM_OPTS:-"-bigwham_eps_aca 1e-4 -bigwham_eta 3 -bigwham_leaf 32"}

GEOM="$(mktemp -d)/scale.in"
OUT=${OUT:-bigwham_omp_vs_mpi.dat}
CSV=${CSV:-bigwham_omp_vs_mpi.csv}
trap 'rm -rf "$(dirname "$GEOM")"' EXIT

for exe in "$BIN" "$MAKEFAULT"; do
    if [ ! -x "$exe" ]; then
        echo "error: cannot find/execute $exe"
        echo "       build with: make bin/compress_interaction_matrix bin/makefault"
        exit 1
    fi
done

"$MAKEFAULT" -strike "$STRIKE" -dip "$DIP" -n "$n" -m "$m" -z -5 > "$GEOM" 2>/dev/null
N=$(grep -cve '^[[:space:]]*$' "$GEOM")

echo "================================================================================"
echo " OpenMP BigWham vs MPI (HTOOL/HACApK/hmmvp), full-space, N = $N patches"
echo "   strike=$STRIKE dip=$DIP   P = [$procs]   nrandom=$nrandom"
echo "   phase 1 = assembly scaling, phase 2 = matvec scaling"
echo "   BigWham assembly threads = OMP_NUM_THREADS, BigWham matvec threads = -bigwham_nthreads"
echo "================================================================================"

name_for(){ case "$1" in 1) echo htool;; 3) echo hacapk;; 4) echo hmmvp;; 5) echo bigwham;; *) echo "be$1";; esac; }
opts_for(){ case "$1" in 1) echo "$HTOOL_OPTS";; 3) echo "$HACAPK_OPTS";; 4) echo "$HMMVP_OPTS";; 5) echo "$BIGWHAM_OPTS";; *) echo "";; esac; }
par_for(){  case "$1" in 5) echo omp;; *) echo mpi;; esac; }

# pull one number out of a compress run
asm_of(){  echo "$1" | grep "H matrix assembly took" | awk '{print $(NF-1)}'; }
mv_of(){   echo "$1" | grep "H-matrix(" | grep "solves" | grep -v inverse | awk '{print $4}' | tr -d s; }
err_of(){  echo "$1" | grep "b-b_h" | awk '{print $NF}'; }
comp_of(){ echo "$1" | grep -oE "compression ratio [0-9.eE+-]+" | head -1 | awk '{print $3}'; }
percall(){ awk -v t="$1" -v k="$2" 'BEGIN{ if(t=="" || t=="-"){print "-"} else printf "%.4f", 1000.0*t/k }'; }

# run a backend; $1 backend  $2 P(matvec/ranks)  $3 omp_threads(bigwham asm)  $4 nrandom
run_be(){
    local be="$1" P="$2" ompP="$3" nr="$4" pre=""
    if [ "$be" -eq 5 ]; then
        OMP_NUM_THREADS="$ompP" OPENBLAS_NUM_THREADS=1 \
            "$BIN" -geom_file "$GEOM" -use_hmatrix 5 -full_space 1 \
                   -bigwham_nthreads "$P" $(opts_for 5) -nrandom "$nr" 2>&1
    else
        [ "$P" -gt 1 ] && pre="$MPIRUN -np $P"
        OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
            $pre "$BIN" -geom_file "$GEOM" -use_hmatrix "$be" -full_space 1 \
                   $(opts_for "$be") -nrandom "$nr" 2>&1
    fi
}

: > "$OUT"
echo "phase,P,backend,par,N,asm_s,mv_ms_per_call,relerr,compression" > "$CSV"

# ---- PHASE 1: assembly scaling ---------------------------------------------
echo "" | tee -a "$OUT"
echo "PHASE 1  assembly scaling (read asm_s; P = OMP threads for bigwham, MPI ranks for the rest)" | tee -a "$OUT"
printf "%-4s %-8s %-4s %7s %-12s %-13s %-11s\n" P backend par N asm_s relerr compression | tee -a "$OUT"
printf -- "--------------------------------------------------------------------------------\n" | tee -a "$OUT"
for P in $procs; do
    for be in 5 1 3 4; do
        # bigwham: assembly uses OMP=P; pin its matvec to 1 thread + nrandom 1 so
        # the (anti-scaling) matvec does not waste wall time in this phase.
        if [ "$be" -eq 5 ]; then out=$(run_be 5 1 "$P" 1); else out=$(run_be "$be" "$P" 1 1); fi
        a=$(asm_of "$out"); e=$(err_of "$out"); c=$(comp_of "$out")
        printf "%-4s %-8s %-4s %7s %-12s %-13s %-11s\n" \
               "$P" "$(name_for "$be")" "$(par_for "$be")" "$N" "${a:--}" "${e:-FAIL}" "${c:--}" | tee -a "$OUT"
        echo "1,$P,$(name_for "$be"),$(par_for "$be"),$N,${a:--},-,${e:-FAIL},${c:--}" >> "$CSV"
    done
    echo "" | tee -a "$OUT"
done

# ---- PHASE 2: matvec scaling -----------------------------------------------
echo "PHASE 2  matvec scaling (read mv_ms_per_call; bigwham matvec threads = -bigwham_nthreads=P," | tee -a "$OUT"
echo "         bigwham assembled with OMP_NUM_THREADS=$ASM_THREADS; Okada backends under P MPI ranks)" | tee -a "$OUT"
printf "%-4s %-8s %-4s %7s %-15s %-13s\n" P backend par N mv_ms_per_call relerr | tee -a "$OUT"
printf -- "--------------------------------------------------------------------------------\n" | tee -a "$OUT"
for P in $procs; do
    for be in 5 1 3 4; do
        if [ "$be" -eq 5 ]; then out=$(run_be 5 "$P" "$ASM_THREADS" "$nrandom"); else out=$(run_be "$be" "$P" 1 "$nrandom"); fi
        mv=$(mv_of "$out"); e=$(err_of "$out"); pc=$(percall "${mv:--}" "$nrandom")
        printf "%-4s %-8s %-4s %7s %-15s %-13s\n" \
               "$P" "$(name_for "$be")" "$(par_for "$be")" "$N" "$pc" "${e:-FAIL}" | tee -a "$OUT"
        echo "2,$P,$(name_for "$be"),$(par_for "$be"),$N,-,$pc,${e:-FAIL},-" >> "$CSV"
    done
    echo "" | tee -a "$OUT"
done

echo "wrote $OUT and $CSV"
echo ""
echo "reading the numbers:"
echo "  - phase 1 asm_s: BigWham scales with OMP threads, the Okada backends with MPI ranks."
echo "    BigWham assembles the full 3N x 3N traction operator, the others the N x N"
echo "    strike-strike block, so compare slope vs P within each group, not absolute time."
echo "  - phase 2 mv_ms_per_call: watch for BigWham's matvec flattening or worsening as"
echo "    -bigwham_nthreads grows; the minimum is the thread count to pin for matvec-bound"
echo "    (rsf_solve) runs. The Okada-backend matvec scales with MPI ranks."
echo "  - relerr is each backend vs the native full-space Okada dense; the compression"
echo "    column uses backend-specific conventions, so do not compare it across backends."
