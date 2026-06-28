#!/bin/bash
#
# Characterize each backend on the fundamental primitive, the accuracy and cost
# of b = A x, independent of the BP5 application. Three things are swept per
# backend, norm mode, and tolerance:
#
#   generic_err   : relative matvec error on a random x, |A_h x - A x|/|A x|,
#                   averaged over several random draws. This estimates the
#                   operator error in a typical direction and is THE fundamental
#                   accuracy axis. (Requires the generic-accuracy probe added to
#                   compress_interaction_matrix.c, which prints
#                   "random-x rel err: mean = ... max = ...".)
#   coherent_err  : the x=1 uniform-slip error, kept as the BP5-relevant
#                   secondary metric (the loading direction that sets timing).
#   stored_scalars, mbytes : memory cost (core-independent).
#   matvec_ms     : per-matvec time. NOTE this is core-count dependent; htool's
#                   matvec parallelizes well and hmmvp's does not, so run this on
#                   the target node count to get a meaningful speed ranking.
#
# The nominal eps is NOT a fair comparison axis: in matrix-relative mode the
# generic error is about eps, in block-local mode it is far below eps, and it
# differs by ~100x between packages at matched eps. Compare on generic_err
# (accuracy) against stored_scalars and matvec_ms (cost), not on eps.
#
# imode 1 = block-local (htool epsilon, hacapk absolute, hmmvp BREM)
# imode 3 = matrix-relative / global (hacapk and hmmvp only; HBI default)
#
# Two compress runs per config (the dense reference in pass 2 is the only cost).
# Run from bp5/ with the PETSc environment loaded, as for run_new_tests.
#
# Output: hmat_bAx.<res>.dat with columns
#   backend imode eps generic_err coherent_err stored_scalars mbytes matvec_ms

for res in 1km 0.5km 0.25km;do                # 1km | 0.5km | 0.25km
ncore=24
nrandom=100

bin=../bin/compress_interaction_matrix
geom=geom_bp5_${res}.in
out=hmat_bAx.${res}.dat
shear=3.204e10
swave=3464

htool_eps="1e-3 7e-4 5e-4 1e-4 1e-5"
hacapk_eps="1e-2 1e-3 1e-4 1e-5"
hmmvp_eps="1e-3 5e-4 1e-4 1e-5 1e-6"

storage_awk='
/hmat_storage backend/ && $0 !~ /backend dense/ {
  for(i=1;i<=NF;i++){ if($i=="m")N=$(i+1); if($i=="stored_scalars")s=$(i+1); if($i=="mbytes")mb=$(i+1) }
}
/compression ratio:/ && ratio=="" {ratio=$NF}
END{ if(s=="NA" && ratio>0){ s=int(N*N/ratio); mb=s*8/1048576 } print s, mb }'

# pass-2 extractor: generic (mean random-x) error, coherent (x=1) error, matvec ms
acc_awk='
/random-x rel err/ {ge=$7}
/\|b-b_h\|\/\|b\| =/ {ce=$NF}
/it took/ && /matrix/ { for(i=1;i<=NF;i++){ if($i ~ /^[0-9.]+s$/){t=$i; sub(/s$/,"",t)}; if($i=="for")c=$(i+1) } }
END{ printf "%s %s %.4f\n", ge, ce, (c>0?1000*t/c:0) }'

echo "# backend imode eps generic_err coherent_err stored_scalars mbytes matvec_ms   ($res, ncore=$ncore)" > $out

run_one(){   # $1 label  $2 imode  $3 eps  rest: backend flags
    local lbl="$1" imode="$2" eps="$3"; shift 3
    read s mb < <(mpirun -np $ncore $bin -geom_file $geom -make_matrix_externally -skip_dense \
                    -shear_modulus $shear -s_wave_speed $swave "$@" 2>&1 | awk "$storage_awk")
    read ge ce mv < <(mpirun -np $ncore $bin -geom_file $geom -nrandom $nrandom \
                    -shear_modulus $shear -s_wave_speed $swave "$@" 2>&1 | awk "$acc_awk")
    echo "$lbl $imode $eps $ge $ce $s $mb $mv" | tee -a $out
}

# htool: block-local only
for eps in $htool_eps; do
    run_one HTOOL 1 $eps -use_hmatrix 1 -mat_htool_eta 3 -mat_htool_epsilon $eps
done

# hacapk: block-local (imode 1) and matrix-relative (imode 3)
for eps in $hacapk_eps; do
    run_one HACApk 1 $eps -use_hmatrix 3 -hacapk_inorm 1 -hacapk_eta 2 -hacapk_ztol $eps
    run_one HACApk 3 $eps -use_hmatrix 3 -hacapk_inorm 3 -hacapk_eta 2 -hacapk_ztol $eps
done

# hmmvp: block-local BREM (imode 1) and matrix-relative MREM (imode 3)
for eps in $hmmvp_eps; do
    run_one HMMVP 1 $eps -use_hmatrix 4 -hmmvp_inorm 1 -hmmvp_tol $eps
    run_one HMMVP 3 $eps -use_hmatrix 4 -hmmvp_inorm 3 -hmmvp_tol $eps
done

echo "wrote $out"
echo "plot stored_scalars and matvec_ms (cost) vs generic_err (accuracy): the"
echo "lowest curve at a target generic_err is the best b = A x backend on that cost."
done
