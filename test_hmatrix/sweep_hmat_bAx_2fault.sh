#!/bin/bash
#
# Two-fault version of sweep_hmat_bAx.sh. Same accuracy-and-cost characterization
# of b = A x (generic error, coherent error, stored size, matvec time, both norm
# modes), but on the interaction operator of a controlled two-fault geometry
# (make_two_fault.py) rather than a single BP5 fault: two parallel 2:1 (L:W)
# faults, separated by 0.5 W in the fault-normal direction and offset along
# strike by 0.05 W, same strike. This gives an operator with two separated
# clusters and a coherent cross-fault far-field, exactly the structure a
# matrix-relative tolerance tends to discard.
#
# Resolution (jmax = down-dip cells) is set so the TOTAL patch count matches the
# single-fault 0.5km (N=16000) and 0.25km (N=64000) runs, so the matrix size and
# dense-reference cost are equivalent and the comparison is at matched problem
# size, just restructured into two faults. total N = 4*jmax^2:
#   jmax = 63  -> 126x63  per fault, 15876 total (~ single-fault 0.5km)
#   jmax = 126 -> 252x126 per fault, 63504 total (~ single-fault 0.25km)
# A coarser two-fault case is not representative of production compression.
#
# Output columns (same as the single-fault sweep):
#   backend imode eps generic_err coherent_err stored_scalars mbytes matvec_ms
#
# Run from test_hmatrix/ with the PETSc environment loaded, as for run_new_tests.

for jmax in 30 60 120 150;do                 # 63 (~16000, single 0.5km-equiv) | 126 (~64000, single 0.25km-equiv)
ncore=24
nrandom=100

bin=../bin/compress_interaction_matrix
shear=3.204e10
swave=3464

htool_eps=" 1e-3 5e-4 1e-4 1e-5"
hacapk_eps="1e-3 1e-4 1e-5 1e-6"
hmmvp_eps=" 1e-3 1e-4 1e-5 1e-6"   # block-local hmmvp over-resolves at tight
                                       # eps on large/structured operators; watch
                                       # storage at the 1e-5/1e-6 end (it can
                                       # approach a large fraction of the dense N^2)

# generate the controlled two-fault geometry (2:1 aspect, 0.5 W separation,
# 0.05 W along-strike offset, set in make_two_fault.py). only the geom is used;
# the rsf/ic/dc files it also writes are for rsf_solve and are unused here.
python3 make_two_fault.py $jmax
geom=geom_2f_j${jmax}.in
out=hmat_bAx_2f.j${jmax}.dat
log=run_2f_j${jmax}.log          # full raw compress output per run, for diagnostics

storage_awk='
/hmat_storage backend/ && $0 !~ /backend dense/ {
  for(i=1;i<=NF;i++){ if($i=="m")N=$(i+1); if($i=="stored_scalars")s=$(i+1); if($i=="mbytes")mb=$(i+1) }
}
/compression ratio:/ && ratio=="" {ratio=$NF}
END{ if(s=="NA" && ratio>0){ s=int(N*N/ratio); mb=s*8/1048576 } print s, mb }'

acc_awk='
/random-x rel err/ {ge=$7}
/\|b-b_h\|\/\|b\| =/ {ce=$NF}
/it took/ && /matrix/ { for(i=1;i<=NF;i++){ if($i ~ /^[0-9.]+s$/){t=$i; sub(/s$/,"",t)}; if($i=="for")c=$(i+1) } }
END{ printf "%s %s %.4f\n", ge, ce, (c>0?1000*t/c:0) }'

echo "# backend imode eps generic_err coherent_err stored_scalars mbytes matvec_ms   (two-fault jmax=$jmax, $geom, ncore=$ncore)" > $out
: > $log

run_one(){   # $1 label  $2 imode  $3 eps  rest: backend flags
    local lbl="$1" imode="$2" eps="$3"; shift 3
    local p1 p2
    echo "### $lbl imode$imode eps=$eps  flags: $*" >> $log
    # pass 1: storage (capture full output, tee to log, then parse)
    p1=$(mpirun -np $ncore $bin -geom_file $geom -make_matrix_externally -skip_dense \
            -shear_modulus $shear -s_wave_speed $swave "$@" 2>&1)
    printf '%s\n' "$p1" >> $log
    read s mb < <(printf '%s\n' "$p1" | awk "$storage_awk")
    # pass 2: generic + coherent error and matvec
    p2=$(mpirun -np $ncore $bin -geom_file $geom -nrandom $nrandom \
            -shear_modulus $shear -s_wave_speed $swave "$@" 2>&1)
    printf '%s\n' "$p2" >> $log
    read ge ce mv < <(printf '%s\n' "$p2" | awk "$acc_awk")
    # a real run fills all of generic, coherent, and storage; if any is empty the
    # compress call failed (degenerate ACA block, MPI/init error, OOM, ...). do not
    # write a blank row that races on silently: mark it FAILED and surface the cause.
    if [ -z "$ge" ] || [ -z "$ce" ] || [ -z "$s" ]; then
        echo "WARNING: $lbl imode$imode eps=$eps produced empty fields; compress failed, see $log" >&2
        printf '%s\n' "$p1" "$p2" | grep -iE 'error|stop|abort|fail|not initialized|MPI_|out of memory' | head -4 >&2
        echo "$lbl $imode $eps FAILED FAILED ${s:-NA} ${mb:-NA} ${mv:-NA}" | tee -a $out
    else
        echo "$lbl $imode $eps $ge $ce $s $mb $mv" | tee -a $out
    fi
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
echo "plot the two costs (mbytes, matvec_ms) against generic_err per backend/mode;"
echo "overlay on the single-fault hmat_bAx.<res>.dat at matched total N to see what"
echo "the two-fault structure (separated clusters, cross-fault far-field) changes."
done
