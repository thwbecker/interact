#!/bin/bash
#
# Diagnose the Ax = b GMRES stagnation on the UCERF interaction operator:
# run a small matrix of KSP configurations on ONE cheap H-matrix operator
# (default: HMMVP global MREM at 1e-4, about 27 ms per apply) with the
# true residual monitored per iteration (summarized in the output), so
# one pass shows (a) at
# what residual level each method stalls and (b) whether unrestarted
# GMRES or BiCGstab breaks the stall. Jacobi is kept on throughout; it
# is harmless, and for this near-uniform mesh essentially a scalar
# scaling (which is why it did not change the iteration counts).
#
# Configurations:
#   gmres30    : the sweep's setting (restart 30), baseline for the stall
#   gmres2000  : effectively unrestarted GMRES; restart-length memory is
#                2000 vectors x N/ncore x 8 bytes per rank (about 88 MB
#                per rank at N = 265464 on 24 ranks), the standard cure
#                for restart stagnation
#   bcgs       : BiCGstab, no restart concept, 2 applies per iteration
#   gmres30_ns : restart 30 WITHOUT preconditioning, to confirm the
#                jacobi-is-a-no-op reading (should match gmres30)
#
# Run from the directory holding the geometry. Parameters are plain
# assignments; no environment variables.

ncore=24
bin=../bin/compress_interaction_matrix
geom=geom.in
maxit=10000
rtol=1e-6
backend="-use_hmatrix 4 -hmmvp_inorm 3 -hmmvp_tol 1e-4"

common="-geom_file $geom -make_matrix_externally -skip_dense -nrandom 0 -nsolve 1"

run_one () {     # label extra-ksp-flags...
    label=$1; shift
    log=log.diag.$label
    if ! grep -q "hmat_solve backend" $log 2> /dev/null ; then
	mpirun -np $ncore $bin $common $backend \
	       -ksp_rtol $rtol -ksp_max_it $maxit \
	       -ksp_monitor_true_residual -ksp_converged_reason \
	       "$@" &> $log
    fi
    printf "%-12s : " $label
    gawk '/hmat_solve backend/ {
            for(i=1;i<=NF;i++){
              if($i=="its_total")   its=$(i+1)
              if($i=="per_solve_s") ss=$(i+1)
              if($i=="reason")      rs=$(i+1)
            }
            printf "its %6i  time %9.1f s  reason %3i  ", its, ss, rs }' $log
    # first, an early, a middle, and the last true-residual samples
    grep "true resid norm" $log | gawk 'NR==1{f=$0} NR==11{e=$0} {l=$0; n++}
         END{ if(n>0){
                split(f,a," "); split(l,b," ");
                printf "true resid: it %s %.1e -> it %s %.1e", a[1], a[10], b[1], b[10]
              } else printf "no residual monitor output (see log)"
            }'
    echo
}

echo "diagnosing on: $backend (rtol $rtol, maxit $maxit)"
run_one gmres30    -pc_type jacobi -ksp_type gmres -ksp_gmres_restart 30
run_one gmres2000  -pc_type jacobi -ksp_type gmres -ksp_gmres_restart 2000
run_one bcgs       -pc_type jacobi -ksp_type bcgs
run_one gmres30_ns -pc_type none   -ksp_type gmres -ksp_gmres_restart 30
