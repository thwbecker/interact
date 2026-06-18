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
# IMPORTANT about BigWham threading (confirmed in the BigWham source):
#   Every BigWham parallel region (block assembly in hmat.cc:187/217 and the
#   matvec in hmat.cc:278) uses num_threads(n_openMP_threads_), i.e. the value
#   passed as -bigwham_nthreads. So that one flag is BigWham's thread count for
#   BOTH assembly and matvec; OMP_NUM_THREADS does not drive those regions. The
#   one exception is the kernel-entry generator (bie_matrix_generator.h:103),
#   which has no num_threads clause and uses the ambient OMP_NUM_THREADS. If
#   OMP_NUM_THREADS is larger than -bigwham_nthreads, that single region
#   oversubscribes against the rest and assembly time explodes. So this script
#   always sets OMP_NUM_THREADS == -bigwham_nthreads == P for BigWham.
#
# Because BigWham uses one thread count for both phases, you cannot assemble with
# many threads and matvec with few in a single build. This script therefore
# reports assembly and matvec in one table per P, so the tradeoff is explicit:
# BigWham assembly speeds up with P while its matvec (private per-thread
# accumulator plus a reduction) flattens and then slows down. For a
# matvec-bound run (rsf_solve) you would pick a low -bigwham_nthreads; for a
# one-shot assemble-and-discard you would pick a high one.
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
#   nrandom  matvecs in the timing loop                 (default 50)
#
# Examples:
#   ./bigwham_omp_vs_mpi_scaling.sh                          # laptop, quick
#   ./bigwham_omp_vs_mpi_scaling.sh 60 40 "1 2 4 8 16 32 48" 100   # big node
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
# roughly comparable accuracy across backends; tune per backend as you like.
HTOOL_OPTS=${HTOOL_OPTS:-"-mat_htool_epsilon 1e-4 -mat_htool_eta 10"}
HACAPK_OPTS=${HACAPK_OPTS:-"-hacapk_ztol 1e-2"}
HMMVP_OPTS=${HMMVP_OPTS:-"-hmmvp_tol 1e-4"}
BIGWHAM_OPTS=${BIGWHAM_OPTS:-"-bigwham_eps_aca 1e-4 -bigwham_eta 3 -bigwham_leaf 32"}
#BIGWHAM_OPTS=${BIGWHAM_OPTS:-"-bigwham_eps_aca 1e-3 -bigwham_eta 5 -bigwham_leaf 32"}

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
#"$MAKEFAULT" -y 2 -strike "$STRIKE" -dip "$DIP" -n "$n" -m "$m" -z -5 >> "$GEOM" 2>/dev/null
N=$(grep -cve '^[[:space:]]*$' "$GEOM")

echo "================================================================================"
echo " OpenMP BigWham vs MPI (HTOOL/HACApK/hmmvp), full-space, N = $N patches"
echo "   strike=$STRIKE dip=$DIP   P = [$procs]   nrandom=$nrandom"
echo "   P = OMP threads (= -bigwham_nthreads) for BigWham, MPI ranks for the Okada backends"
echo "================================================================================"

name_for(){ case "$1" in 1) echo htool;; 3) echo hacapk;; 4) echo hmmvp;; 5) echo bigwham;; *) echo "be$1";; esac; }
opts_for(){ case "$1" in 1) echo "$HTOOL_OPTS";; 3) echo "$HACAPK_OPTS";; 4) echo "$HMMVP_OPTS";; 5) echo "$BIGWHAM_OPTS";; *) echo "";; esac; }
par_for(){  case "$1" in 5) echo omp;; *) echo mpi;; esac; }

asm_of(){  echo "$1" | grep "H matrix assembly took" | awk '{print $(NF-1)}'; }
mv_of(){   echo "$1" | grep "H-matrix(" | grep "solves" | grep -v inverse | awk '{print $4}' | tr -d s; }
err_of(){  echo "$1" | grep "b-b_h" | awk '{print $NF}'; }
comp_of(){ echo "$1" | grep -oE "compression ratio [0-9.eE+-]+" | head -1 | awk '{print $3}'; }
percall(){ awk -v t="$1" -v k="$2" 'BEGIN{ if(t=="" || t=="-"){print "-"} else printf "%.4f", 1000.0*t/k }'; }

# run a backend at parallelism P; echoes the full compress output
run_be(){
    local be="$1" P="$2" pre=""
    if [ "$be" -eq 5 ]; then
        # BigWham: one rank, P threads for assembly AND matvec, BLAS single-
        # threaded (no nesting), and OMP_NUM_THREADS == P so the generator
        # region does not oversubscribe against the num_threads(P) regions.
        OMP_NUM_THREADS="$P" OPENBLAS_NUM_THREADS=1 \
            "$BIN" -geom_file "$GEOM" -use_hmatrix 5 -full_space 1 \
                   -bigwham_nthreads "$P" $(opts_for 5) -nrandom "$nrandom" 2>&1
    else
        [ "$P" -gt 1 ] && pre="$MPIRUN -np $P"
        OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
            $pre "$BIN" -geom_file "$GEOM" -use_hmatrix "$be" -full_space 1 \
                   $(opts_for "$be") -nrandom "$nrandom" 2>&1
    fi
}

: > "$OUT"
printf "%-4s %-8s %-4s %7s %-12s %-15s %-13s %-11s\n" \
       P backend par N asm_s mv_ms_per_call relerr compression | tee -a "$OUT"
printf -- "--------------------------------------------------------------------------------\n" | tee -a "$OUT"
echo "P,backend,par,N,asm_s,mv_ms_per_call,relerr,compression" > "$CSV"

for P in $procs; do
    for be in 5 1 3 4; do
        out=$(run_be "$be" "$P")
        a=$(asm_of "$out"); mv=$(mv_of "$out"); e=$(err_of "$out"); c=$(comp_of "$out")
        pc=$(percall "${mv:--}" "$nrandom")
        printf "%-4s %-8s %-4s %7s %-12s %-15s %-13s %-11s\n" \
               "$P" "$(name_for "$be")" "$(par_for "$be")" "$N" "${a:--}" "$pc" "${e:-FAIL}" "${c:--}" | tee -a "$OUT"
        echo "$P,$(name_for "$be"),$(par_for "$be"),$N,${a:--},$pc,${e:-FAIL},${c:--}" >> "$CSV"
    done
    echo "" | tee -a "$OUT"
done

echo "wrote $OUT and $CSV"
echo ""
echo "reading the numbers:"
echo "  - BigWham uses ONE thread count (-bigwham_nthreads, kept == OMP_NUM_THREADS == P"
echo "    here) for both assembly and matvec. Expect asm_s to fall with P while"
echo "    mv_ms_per_call flattens then rises: that is the assembly-vs-matvec tradeoff."
echo "    For matvec-bound (rsf_solve) runs, read mv_ms_per_call and pick the P at its"
echo "    minimum, accepting the slower one-time assembly there."
echo "  - the Okada backends scale assembly and matvec with MPI ranks."
echo "  - BigWham assembles the full 3N x 3N traction operator, the others the N x N"
echo "    strike-strike block, so compare slope vs P within each group, not absolute time."
echo "  - relerr is each backend vs the native full-space Okada dense; the compression"
echo "    column uses backend-specific conventions, so do not compare it across backends."
