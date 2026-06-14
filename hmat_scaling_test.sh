#!/bin/bash
#
# benchmark H matrix backends of compress_interaction_matrix vs core count
#
#   use_hmatrix 1 (HTOOL)  - MPI parallel (mpirun -np P)
#   use_hmatrix 3 (HACApK) - MPI parallel (mpirun -np P)
#   use_hmatrix 4 (hmmvp)  - MPI parallel (mpirun -np P): distributed assembly
#                            to a temp file + MpiHmat matvec (root-gather)
#
# every run rebuilds the (MPI-distributed) dense reference and reports the
# relative error of b = A x, the H assembly time, and timings for NRANDOM
# dense and H matrix-vector products, so each line is self-checking.
#
# tolerances below are the matched ~1e-5 error class settings from
# compress_interaction_matrix.md; adjust as needed.
#
# notes:
# - HACApK matvecs gather x to all ranks and use a ring exchange:
#   expect sublinear MPI scaling, that is part of what this measures
# - hmmvp now runs under MPI (was OpenMP); assembly is distributed across
#   ranks and the matvec gathers x to the root, computes distributed, and
#   scatters y back (expect a root-funnel cost like HACApK)
# - set MPIRUN="mpirun --oversubscribe" etc. if needed
#
BIN=${BIN:-compress_interaction_matrix}
MPIRUN=${MPIRUN:-$PETSC_DIR/build/bin/mpirun}
GEOM=${GEOM:-geom.in}
NFAULT=${NFAULT:-120}            # makefault -n; 120 x 60 x 2 = 14400 patches
MFAULT=${MFAULT:-60}
CORES=${CORES:-"1 2 4 8 16 24 48"}
NRANDOM=${NRANDOM:-100}
HTOOL_OPTS=${HTOOL_OPTS:-"-mat_htool_epsilon 3e-5 -mat_htool_eta 10"}
HACAPK_OPTS=${HACAPK_OPTS:-"-hacapk_ztol 1e-4"}
HMMVP_OPTS=${HMMVP_OPTS:-"-hmmvp_tol 1e-6"}
OUT=${OUT:-hmat_scaling.dat}

# geometry
if [ ! -s $GEOM ]; then
    echo "$0: generating $GEOM with 2 x $NFAULT x $MFAULT patches"
    makefault -n $NFAULT -m $MFAULT > $GEOM
    makefault -n $NFAULT -m $MFAULT -x 1 -grp 1 >> $GEOM
fi
N=$(wc -l < $GEOM)

run_one(){ # $1 command prefix, $2 use_hmatrix, $3 options, $4 label, $5 cores
    local out err ha td th
    out=$($1 $BIN -geom_file $GEOM -use_hmatrix $2 $3 -nrandom $NRANDOM 2>&1)
    err=$(echo "$out" | grep "b-b_h" | awk '{printf "%.2e",$NF}')
    ha=$(echo  "$out" | grep "H matrix assembly" | awk '{printf "%.1f",$(NF-1)}')
    td=$(echo  "$out" | grep "dense solves" | awk '{print $4}' | tr -d s)
    th=$(echo  "$out" | grep -E "Hmatrix\(|H-matrix\(" | awk '{print $4}' | tr -d s)
    printf "%-8s %5d %5d  err %-9s  asm %8s s  mv_dense %8s s  mv_H %8s s\n" \
	   "$4" $N $5 "${err:-FAIL}" "${ha:--}" "${td:--}" "${th:--}" | tee -a $OUT
}

echo "# N=$N nrandom=$NRANDOM cores: $CORES  ($(date))" | tee $OUT
# thread hygiene, applies to ALL runs:
# - OMP_NUM_THREADS=1 keeps MPI ranks from oversubscribing; hmmvp
#   raises its own pool internally via omp_set_num_threads(nthreads)
#   from -hmmvp_nthreads, so it does NOT need the env var
# - OPENBLAS/MKL pins keep BLAS single threaded: PETSc's dense MatMult
#   and hmmvp's per-block dgemv otherwise nest BLAS threads inside
#   block-level parallelism and thrash catastrophically
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

unset OMP_PLACES; export OMP_PROC_BIND=false
# this was a test for HMMVP earlier, did not do much/hurt
#export OMP_PROC_BIND=close
#export OMP_PLACES=cores 

for P in $CORES; do
    run_one "$MPIRUN --bind-to core --map-by core -np $P"  1 "$HTOOL_OPTS"  "HTOOL"  $P
done
for P in $CORES; do
    run_one "$MPIRUN --bind-to core --map-by core -np $P"  3 "$HACAPK_OPTS" "HACApK" $P
done
for P in $CORES; do
    run_one "$MPIRUN --bind-to core --map-by core -np $P"  4 "$HMMVP_OPTS" "hmmvp" $P
done
echo "$0: done, results in $OUT"
