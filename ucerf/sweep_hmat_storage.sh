#!/bin/bash
#
# Sweep H-matrix storage and measured per-matvec cost vs tolerance for the
# H-matrix backends on a given geometry, using compress_interaction_matrix
# with -skip_dense (build and time the operator only, no dense reference,
# no cycle integration), and report speedup and memory reduction against a
# dense baseline.
#
# Dense baseline logic:
#   make_dense_reference=1 : run compress_interaction_matrix
#       -dense_reference_only ONCE (with ncore_dense ranks; each rank's
#       local block local_rows*n must stay below 2^31 entries, so at large
#       N this needs correspondingly many ranks; the tool errors out with
#       the required advice if not), and store the measured assembly time
#       and per-matvec time in dense_ref_file. Subsequent runs can then
#       set make_dense_reference=0.
#   make_dense_reference=0 : if dense_ref_file exists, read the measured
#       baseline from it; otherwise ESTIMATE the dense matvec as memory
#       bandwidth bound, N*N*8 bytes / (bw_gbs GB/s), and mark the speedup
#       column "est". Estimated numbers ignore parallel scaling and cache
#       effects and are only meant to set the scale.
#
# Requires the 2026-07 tool version (-skip_dense honored in the
# -make_matrix_externally path, hmat_matvec and dense_reference output
# lines, -dense_reference_only mode).
#
# Backend flags mirror the uniform block-local settings of run_new_tests
# (htool eta 3, hacapk inorm 1 eta 2, hmmvp inorm 1) plus an hmmvp
# global-MREM row (inorm 3) for contrast.
#
# Storage extraction per backend: HMMVP from its assembly line ("stored
# scalars"); HTOOL from the MatView compression ratio as N*N/ratio;
# HACApK prints no storage figure (NA).
#
# Run from the directory holding the geometry. All parameters are plain
# assignments below; edit them there (no environment variables).

ncore=48                 # MPI ranks for the H-matrix sweep runs
ncore_dense=48           # MPI ranks for the one dense reference run
nrandom=200              # matvec timing applies (0 skips the timing line!)
nsolve=3                 # Ax=b GMRES solves to time per point (0: off)
solver_pc=jacobi         # KSP preconditioner for the solve timing (jacobi
                         # or none). jacobi needs the 2026-07 binary, which
                         # provides MatGetDiagonal on the shell H operators
                         # via a diagonal cache filled at assembly; without
                         # convergence (its at the cap) solve times measure
                         # iteration throughput only
bin=../bin/compress_interaction_matrix
geom=geom.in
out=hmat_storage.dat
dense_ref_file=dense_reference.dat
make_dense_reference=1   # 1: (re)create dense_ref_file first, then sweep
bw_gbs=200               # node memory bandwidth [GB/s] for the estimate

htool_eps="1e-3 1e-4 1e-5 1e-6"
hacapk_eps="1e-3 1e-4 1e-5 1e-6"
hmmvp_eps="1e-3 1e-4 1e-5 1e-6"

# ---------------------------------------------------------------- dense
if [ $make_dense_reference -eq 1 ]; then
    mpirun -np $ncore_dense $bin -geom_file $geom -make_matrix_externally \
	   -use_hmatrix 0 -dense_reference_only -nrandom $nrandom -nsolve $nsolve \
	   -pc_type $solver_pc &> log.dense_ref
    gawk 'BEGIN{si="NA";ss="NA"}
          /dense_solve m/ {
            for(i=1;i<=NF;i++){
              if($i=="its_total")   sits=$(i+1)
              if($i=="nsolve")      ns=$(i+1)
              if($i=="per_solve_s") ss=$(i+1)
            }
            if(ns+0>0) si=sits/ns
          }
          /dense_reference m/ {
            for(i=1;i<=NF;i++){
              if($i=="m")             N=$(i+1)
              if($i=="assembly_s")    as=$(i+1)
              if($i=="per_matvec_ms") ms=$(i+1)
            }
            printf "# measured dense reference: ncore m assembly_s per_matvec_ms solve_its per_solve_s\n%i %i %s %s %s %s\n", NC, N, as, ms, si, ss
          }' NC=$ncore_dense log.dense_ref > $dense_ref_file
    if [ ! -s $dense_ref_file ]; then
	echo "dense reference run failed; see log.dense_ref" ; exit 1
    fi
    echo "wrote $dense_ref_file"
fi

npatch=`grep -cv '^#' $geom`
if [ -s $dense_ref_file ]; then
    read d_nc d_n d_as d_ms d_si d_ss << EOF
`grep -v '^#' $dense_ref_file`
EOF
    d_si=${d_si:-NA}
    d_ss=${d_ss:-NA}
    d_tag="meas_np$d_nc"
else
    # bandwidth-bound estimate: N*N*8 bytes per apply
    d_ms=`echo $npatch $bw_gbs | gawk '{printf "%.4f", $1*$1*8/($2*1e9)*1000}'`
    d_as="NA"
    d_si="NA"
    d_ss="NA"
    d_tag="est_${bw_gbs}GBs"
fi

# one row from one log: label, eps, stored scalars, MByte, ratio, matvec
# ms, speedup vs dense. the hmat_matvec line is the success marker.
extract='
/hmat_matvec backend/ {
  for(i=1;i<=NF;i++){
    if($i=="m")          N=$(i+1)
    if($i=="per_matvec") ms=$(i+1)
  }
  ok=1
}
/hmat_assembly backend/ {
  for(i=1;i<=NF;i++) if($i=="assembly_s") as=$(i+1)
}
/hmat_solve backend/ {
  for(i=1;i<=NF;i++){
    if($i=="its_total")   sits=$(i+1)
    if($i=="nsolve")      ns=$(i+1)
    if($i=="per_solve_s") ss=$(i+1)
  }
  if(ns+0>0) si=sits/ns
}
/hmat_storage backend/ {      # primary: printed by calc_petsc_Isn_matrices
  for(i=1;i<=NF;i++){         # for every backend (NA where unavailable)
    if($i=="stored_scalars" && $(i+1)!="NA") s=$(i+1)
    if($i=="dense_ratio"    && $(i+1)!="NA") r=$(i+1)
  }
}
/stored scalars/ {            # fallback: hmmvp assembly line
  if(s=="") for(i=1;i<=NF;i++) if($(i+1)=="stored") s=$i
}
/compression ratio/ { if(r=="") r=$NF }   # fallback: htool MatView
END{
  if(!ok){ printf "# FAILED %s %s (no hmat_matvec line; see log)\n", B, EPS; exit }
  if(s=="" && r!="" && r+0>0) s=int(N*N/r)
  if(s!="" && (r=="" || r+0<=0)) r=N*N/s
  mb=(s!="")?(s*8/1048576):("NA")
  sp=(ms+0>0)?sprintf("%.1f",DMS/ms):("NA")
  asp=(as!="" && DAS!="NA" && as+0>0)?sprintf("%.1f",DAS/as):("NA")
  ssp=(ss!="" && DSS!="NA" && ss+0>0)?sprintf("%.1f",DSS/ss):("NA")
  printf "%s %s %s %s %s %.4f %s %s %s %s %s %s %s\n", B, EPS, (s!=""?s:"NA"), mb, (r!=""?r:"NA"), ms, sp, (as!=""?as:"NA"), asp, (ss!=""?ss:"NA"), (si!=""?si:"NA"), ssp, DTAG
}'

echo "# dense baseline: per_matvec $d_ms ms assembly ${d_as} s per_solve ${d_ss} s solve_its ${d_si} [$d_tag]" > $out
echo "# backend eps stored_scalars mbytes compression_ratio matvec_ms matvec_speedup assembly_s assembly_speedup per_solve_s solve_its solve_speedup dense_tag  (geom=$geom npatch=$npatch ncore=$ncore nrandom=$nrandom nsolve=$nsolve)" >> $out

# resumable: if the log already holds a completed run (hmat_matvec marker),
# only re-extract; delete the log to force a rerun of that point
run_one () {                 # label logfile extra-flags...
    label=$1; log=$2; shift 2
    if ! grep -q "hmat_matvec backend" $log 2> /dev/null ; then
	mpirun -np $ncore $bin -geom_file $geom -make_matrix_externally \
	       -skip_dense -nrandom $nrandom -nsolve $nsolve \
	       -pc_type $solver_pc "$@" &> $log
    fi
    gawk -v B=$label -v EPS=$eps -v DMS=$d_ms -v DAS=$d_as -v DSS=$d_ss -v DTAG=$d_tag "$extract" $log >> $out
}

for eps in $htool_eps; do
    run_one HTOOL log.$eps.htools -use_hmatrix 1 -mat_htool_eta 3 -mat_htool_epsilon $eps
done

for eps in $hacapk_eps; do
    run_one HACApK log.$eps.hacapk -use_hmatrix 3 -hacapk_inorm 1 -hacapk_eta 2 -hacapk_ztol $eps
done

for eps in $hmmvp_eps; do
    run_one HMMVP log.$eps.hmmvp -use_hmatrix 4 -hmmvp_inorm 1 -hmmvp_tol $eps
done

for eps in $hmmvp_eps; do
    run_one HMMVP_MREM log.$eps.hmmvp_mrem -use_hmatrix 4 -hmmvp_inorm 3 -hmmvp_tol $eps
done

echo "wrote $out"
