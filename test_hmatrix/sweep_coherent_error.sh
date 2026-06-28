#!/bin/bash
#
# Sweep the coherent loading error per backend, norm mode, and tolerance. The
# coherent error is the relative error of the interaction operator applied to a
# uniform-slip vector, |b-b_h|/|b|, i.e. the steady loading that sets the
# recurrence interval. It is computed in a single matvec against the dense
# reference, so it needs no cycle run, and it is the fair package-independent
# accuracy axis: the same nominal eps buys very different physical accuracy in
# each package and norm mode, so eps is not a usable common metric.
#
# This version sweeps BOTH norm modes so they can be compared directly:
#   imode 1 = block-local   (htool epsilon, hacapk absolute, hmmvp BREM)
#   imode 3 = matrix-relative / global  (hacapk and hmmvp only; HBI default)
# htool has only a block-local (block-relative) mode and is run once.
#
# Two compress runs per config (the only cost is the dense reference in pass 2,
# cheap on a real node):
#   pass 1 (-make_matrix_externally -skip_dense): storage report.
#   pass 2 (default path, builds the dense reference): x=1 coherent error and
#           the H-matrix matvec timing.
#
# Run from bp5/ with the PETSc environment loaded, as for run_new_tests.
#
# Output: hmat_coherent.<res>.dat with columns
#   backend imode eps coherent_err stored_scalars mbytes matvec_ms

res=1km                 # 1km | 0.5km | 0.25km
ncore=24
nrandom=100

bin=../bin/compress_interaction_matrix
geom=geom_bp5_${res}.in
out=hmat_coherent.${res}.dat
shear=3.204e10
swave=3464

htool_eps="1e-3 7e-4 5e-4 1e-4 1e-5"
hacapk_eps="1e-2 1e-3 1e-4 1e-5"
hmmvp_eps="1e-3 5e-4 1e-4 1e-5 1e-6"

# pass-1 extractor: storage (htool reports NA, recovered from compression ratio;
# N taken from the m field so it is resolution-agnostic)
storage_awk='
/hmat_storage backend/ && $0 !~ /backend dense/ {
  for(i=1;i<=NF;i++){ if($i=="m")N=$(i+1); if($i=="stored_scalars")s=$(i+1); if($i=="mbytes")mb=$(i+1) }
}
/compression ratio:/ && ratio=="" {ratio=$NF}
END{ if(s=="NA" && ratio>0){ s=int(N*N/ratio); mb=s*8/1048576 } print s, mb }'

# pass-2 extractor: coherent error and per-matvec ms (the /matrix/ guard picks
# the H-matrix solve-timing line, not the dense one)
cohmv_awk='
/\|b-b_h\|\/\|b\| =/ {ce=$NF}
/it took/ && /matrix/ { for(i=1;i<=NF;i++){ if($i ~ /^[0-9.]+s$/){t=$i; sub(/s$/,"",t)}; if($i=="for")c=$(i+1) } }
END{ printf "%s %.4f\n", ce, (c>0?1000*t/c:0) }'

echo "# backend imode eps coherent_err stored_scalars mbytes matvec_ms   ($res, ncore=$ncore)" > $out

run_one(){   # $1 label  $2 imode  $3 eps  rest: backend flags
    local lbl="$1" imode="$2" eps="$3"; shift 3
    read s mb < <(mpirun -np $ncore $bin -geom_file $geom -make_matrix_externally -skip_dense \
                    -shear_modulus $shear -s_wave_speed $swave "$@" 2>&1 | awk "$storage_awk")
    read ce mv < <(mpirun -np $ncore $bin -geom_file $geom -nrandom $nrandom \
                    -shear_modulus $shear -s_wave_speed $swave "$@" 2>&1 | awk "$cohmv_awk")
    echo "$lbl $imode $eps $ce $s $mb $mv" | tee -a $out
}

# htool: block-local only (epsilon is block-relative; no matrix-relative mode)
for eps in $htool_eps; do
    run_one HTOOL 1 $eps -use_hmatrix 1 -mat_htool_eta 3 -mat_htool_epsilon $eps
done

# hacapk: imode 1 (absolute, block-local) and imode 3 (matrix-relative, global; HBI default)
for eps in $hacapk_eps; do
    run_one HACApk 1 $eps -use_hmatrix 3 -hacapk_inorm 1 -hacapk_eta 2 -hacapk_ztol $eps
    run_one HACApk 3 $eps -use_hmatrix 3 -hacapk_inorm 3 -hacapk_eta 2 -hacapk_ztol $eps
done

# hmmvp: imode 1 (BREM, block-local) and imode 3 (MREM, matrix-relative, global; default)
for eps in $hmmvp_eps; do
    run_one HMMVP 1 $eps -use_hmatrix 4 -hmmvp_inorm 1 -hmmvp_tol $eps
    run_one HMMVP 3 $eps -use_hmatrix 4 -hmmvp_inorm 3 -hmmvp_tol $eps
done

echo "wrote $out"
echo "compare imode 1 vs imode 3 per backend: plot coherent_err vs eps, and cost"
echo "(matvec_ms or mbytes) vs coherent_err. block-local (imode 1) should sit well"
echo "below matrix-relative (imode 3) in coherent error at matched eps."
