#!/bin/bash
#
# bigwham_omp_vs_mpi_scaling.sh
#
# Assembly (and matvec) scaling on a full-space problem, comparing:
#   - HTOOL        (use_hmatrix 1): MPI-parallel,    Okada full-space (-full_space 1)
#   - HACApK       (use_hmatrix 3): MPI-parallel,    Okada full-space
#   - HMMVP        (use_hmatrix 4): MPI-parallel,    Okada full-space
#   - BigWham      (use_hmatrix 5): OpenMP-parallel, full-space native kernel
#
# For each parallelism level P the script runs BigWham on a single MPI rank with
# P OpenMP threads, and the Okada backends on P MPI ranks with one thread each.
# So P means "OpenMP threads" for BigWham and "MPI ranks" for the others; the
# column header reflects that. The Okada full-space dense matrix is built by the
# tool on every run, so the relerr column is, for every backend, the error of
# that backend against the same native full-space Okada reference (apples to
# apples, since BigWham is also full-space).
#
# This reports wall-clock assembly time, forward matvec time (summed over the
# timing loop), accuracy, and the reported compression ratio. Timings are only
# meaningful on a machine with at least P real cores; on an oversubscribed box
# (e.g. a laptop with fewer cores than P, or a single-core container) the
# numbers still print but should not be read as scaling. Results are specific to
# this geometry, these accuracy settings, and this machine; a different BigWham
# build or a different HACApK / HTOOL / hmmvp version may behave differently.
#
# Usage (positional arguments, all optional):
#   ./bigwham_omp_vs_mpi_scaling.sh [n] [m] [procs] [nrandom]
#
#   n        patches along strike                       (default 30)
#   m        patches along dip                          (default 20)
#   procs    quoted list of parallelism levels          (default "1 2 4")
#   nrandom  matvecs in the timing loop                 (default 20)
#
# Examples:
#   ./bigwham_omp_vs_mpi_scaling.sh
#   ./bigwham_omp_vs_mpi_scaling.sh 60 40 "1 4 16 48" 100
#
set -u

n=${1:-30}
m=${2:-20}
procs=${3:-"1 2 4"}
nrandom=${4:-20}

# ---- configuration (edit for your machine) ----------------------------------
BIN=${BIN:-../bin/compress_interaction_matrix}
MAKEFAULT=${MAKEFAULT:-../bin/makefault}
# launcher for the MPI (Okada) backends. On a laptop with few cores use
# "mpirun --oversubscribe"; on a cluster use your scheduler launcher (ibrun,
# srun, ...). BigWham always runs on one rank, so it does not use this.
MPIRUN=${MPIRUN:-mpirun}
STRIKE=${STRIKE:-30}
DIP=${DIP:-70}
# roughly comparable accuracy across backends (same spirit as the existing
# hmat scaling tests); tune per backend as you like.
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

# ---- geometry: one planar full-space fault, n x m patches -------------------
"$MAKEFAULT" -strike "$STRIKE" -dip "$DIP" -n "$n" -m "$m" -z -5 > "$GEOM" 2>/dev/null
N=$(grep -cve '^[[:space:]]*$' "$GEOM")
echo "================================================================================"
echo " OpenMP BigWham vs MPI (HTOOL/HACApK/hmmvp), full-space, N = $N patches"
echo "   strike=$STRIKE dip=$DIP   parallelism P = [$procs]   nrandom=$nrandom"
echo "   P = OpenMP threads for BigWham, MPI ranks for the Okada backends"
echo "================================================================================"

opts_for(){ case "$1" in
    1) echo "$HTOOL_OPTS";; 3) echo "$HACAPK_OPTS";;
    4) echo "$HMMVP_OPTS";; 5) echo "$BIGWHAM_OPTS";; *) echo "";; esac; }
name_for(){ case "$1" in 1) echo htool;; 3) echo hacapk;; 4) echo hmmvp;; 5) echo bigwham;; *) echo "be$1";; esac; }
par_for(){  case "$1" in 5) echo omp;; *) echo mpi;; esac; }

# run one backend at parallelism P; echoes: asm_s mv_s relerr compression
run_one(){ # $1 backend  $2 P
    local be="$1" P="$2" pre="" out asm mv err comp
    if [ "$be" -eq 5 ]; then
        # BigWham: single rank, P OpenMP threads
        out=$(OMP_NUM_THREADS="$P" OPENBLAS_NUM_THREADS="$P" \
              "$BIN" -geom_file "$GEOM" -use_hmatrix 5 -full_space 1 \
                     -bigwham_nthreads "$P" $(opts_for 5) -nrandom "$nrandom" 2>&1)
    else
        # Okada backend: P MPI ranks, one thread each
        [ "$P" -gt 1 ] && pre="$MPIRUN -np $P"
        out=$(OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
              $pre "$BIN" -geom_file "$GEOM" -use_hmatrix "$be" -full_space 1 \
                     $(opts_for "$be") -nrandom "$nrandom" 2>&1)
    fi
    asm=$(echo "$out" | grep "H matrix assembly took" | awk '{print $(NF-1)}')
    mv=$(echo  "$out" | grep "H-matrix(" | grep "solves" | grep -v inverse | awk '{print $4}' | tr -d s)
    err=$(echo "$out" | grep "b-b_h" | awk '{print $NF}')
    comp=$(echo "$out" | grep -oE "compression ratio [0-9.eE+-]+" | head -1 | awk '{print $3}')
    echo "${asm:--} ${mv:--} ${err:-FAIL} ${comp:--}"
}

: > "$OUT"
printf "%-4s %-8s %-4s %7s %-12s %-12s %-13s %-11s\n" \
       P backend par N asm_s mv_s relerr compression | tee -a "$OUT"
printf -- "--------------------------------------------------------------------------------\n" | tee -a "$OUT"
echo "P,backend,par,N,asm_s,mv_s,relerr,compression" > "$CSV"

for P in $procs; do
    for be in 5 1 3 4; do
        read asm mv err comp < <(run_one "$be" "$P")
        printf "%-4s %-8s %-4s %7s %-12s %-12s %-13s %-11s\n" \
               "$P" "$(name_for "$be")" "$(par_for "$be")" "$N" "$asm" "$mv" "$err" "$comp" | tee -a "$OUT"
        echo "$P,$(name_for "$be"),$(par_for "$be"),$N,$asm,$mv,$err,$comp" >> "$CSV"
    done
    echo "" | tee -a "$OUT"
done

echo "wrote $OUT and $CSV"
echo ""
echo "reminder: assembly/matvec timings are meaningful only with >= P real cores;"
echo "          the relerr column is each backend vs the native full-space Okada dense."
