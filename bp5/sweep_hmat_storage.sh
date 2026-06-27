#!/bin/bash
#
# Sweep H-matrix storage and a per-matvec cost estimate vs tolerance for the
# three backends at one BP5 resolution, using compress_interaction_matrix so
# it only builds the interaction matrix (no cycle integration). This is the
# cheap way to extend the storage-vs-tolerance comparison to 0.5km and 0.25km:
# assembly is fast on a real node but far too slow on a single core.
#
# The backend flags mirror the uniform block-local settings of run_new_tests
# (htool eta 3, hacapk inorm 1 eta 2, hmmvp inorm 1) so the points line up
# with the accuracy runs. It also adds an hmmvp global-MREM row (inorm 3) for
# contrast, since that is what run_new_tests currently uses for hmmvp.
#
# Assumes the PETSc environment is already loaded (modules), same as
# run_new_tests. Run from the bp5/ directory.
#
# Output: hmat_storage.<res>.dat with columns
#   backend  eps  stored_scalars  mbytes  matvec_ms
# For htool, stored_scalars is recovered from the printed compression ratio
# (htool reports NA on the storage line in this PETSc build). matvec_ms is a
# rank-0 CPU-time estimate over (nrandom+1) applies; for precise parallel
# matvec timing prefer the MatMult line from the -log_view of the rsf runs.

# parameters (plain assignments; edit here)
res=1km                                 # 1km | 0.5km | 0.25km
ncore=24                                # MPI ranks
nrandom=200                             # matvec timing applies (set 0 to skip, faster)

bin=../bin/compress_interaction_matrix
geom=geom_bp5_${res}.in
out=hmat_storage.${res}.dat

shear=3.204e10
swave=3464

# eps lists per backend, matching run_new_tests
htool_eps="1e-3 7e-4 5e-4 1e-4 1e-5"
hacapk_eps="1e-2 1e-3 1e-4 1e-5"
hmmvp_eps="1e-3 5e-4 1e-4 1e-5 1e-6"

# extractor: pulls stored scalars, Mbyte, and per-matvec ms from one run.
# N (for the htool ratio inversion) is read from the m field on the backend
# storage line, so this works at any resolution. The /matrix/ guard picks the
# H-matrix solve-timing line rather than the dense one.
extract='
/hmat_storage backend/ && $0 !~ /backend dense/ {
  for(i=1;i<=NF;i++){
    if($i=="backend")        b=$(i+1)
    if($i=="m")              N=$(i+1)
    if($i=="stored_scalars") s=$(i+1)
    if($i=="mbytes")         mb=$(i+1)
  }
}
/compression ratio:/ && ratio=="" {ratio=$NF}
/it took/ && /matrix/ {
  for(i=1;i<=NF;i++){
    if($i ~ /^[0-9.]+s$/){tt=$i; sub(/s$/,"",tt)}
    if($i=="for") cnt=$(i+1)
  }
}
END{
  if(s=="NA" && ratio>0){ s=int(N*N/ratio); mb=s*8/1048576 }
  printf "%s %s %s %.2f %.4f\n", B, EPS, s, mb, (cnt>0?1000*tt/cnt:0)
}'

echo "# backend eps stored_scalars mbytes matvec_ms   ($res, ncore=$ncore, nrandom=$nrandom)" > $out

for eps in $htool_eps; do
    mpirun -np $ncore $bin -geom_file $geom -make_matrix_externally -skip_dense -nrandom $nrandom \
           -shear_modulus $shear -s_wave_speed $swave \
           -use_hmatrix 1 -mat_htool_eta 3 -mat_htool_epsilon $eps 2>&1 \
        | awk -v B=HTOOL -v EPS=$eps "$extract" >> $out
done

for eps in $hacapk_eps; do
    mpirun -np $ncore $bin -geom_file $geom -make_matrix_externally -skip_dense -nrandom $nrandom \
           -shear_modulus $shear -s_wave_speed $swave \
           -use_hmatrix 3 -hacapk_inorm 1 -hacapk_eta 2 -hacapk_ztol $eps 2>&1 \
        | awk -v B=HACApk -v EPS=$eps "$extract" >> $out
done

for eps in $hmmvp_eps; do
    mpirun -np $ncore $bin -geom_file $geom -make_matrix_externally -skip_dense -nrandom $nrandom \
           -shear_modulus $shear -s_wave_speed $swave \
           -use_hmatrix 4 -hmmvp_inorm 1 -hmmvp_tol $eps 2>&1 \
        | awk -v B=HMMVP_BREM -v EPS=$eps "$extract" >> $out
    mpirun -np $ncore $bin -geom_file $geom -make_matrix_externally -skip_dense -nrandom $nrandom \
           -shear_modulus $shear -s_wave_speed $swave \
           -use_hmatrix 4 -hmmvp_inorm 3 -hmmvp_tol $eps 2>&1 \
        | awk -v B=HMMVP_MREM -v EPS=$eps "$extract" >> $out
done

echo "wrote $out"
