#!/usr/bin/env bash
#
# rsf_solve_scaling_test.sh -- MPI scaling benchmark for rsf_solve's matrix
#   backends (dense, HTOOL/ACA, HACApK, hmmvp) on the SEAS BP5 case.
#
# Sweeps MPI rank counts x backends, runs rsf_solve with -log_view, and parses
# assembly time, mean matvec time, total wallclock, #steps, and H-matrix
# memory into a table (also written as CSV).
#
# DEFAULTS ARE THE MATCHED ~1e-6 ERROR-BAND SETTINGS derived empirically from
# the forward-operator sweep in compress_interaction_matrix.md (N=14400):
#
#     HTOOL  -mat_htool_epsilon 3e-5   -> ~6.6e-7
#     HACApK -hacapk_ztol       1e-1   -> ~2.2e-7   (NOTE: ztol is very
#            conservative for the smooth Okada kernel - 1e-1 is much looser
#            than nominal and is what keeps HACApK from running near-dense;
#            do NOT use ztol 1e-4 here, that is ~near-dense and not comparable)
#     hmmvp  -hmmvp_tol         1e-7   -> ~1.6e-6   (floored ~1e-6 at this N by
#            hmmvp's stochastic Frobenius-norm estimate; tightening buys little)
#
# so the four backends solve the SAME problem at comparable operator accuracy
# and the only thing that differs is the H-matrix machinery. RTOL below is the
# (separate) time-integration tolerance and is held fixed across backends.
#
# Usage:   ./rsf_solve_scaling_test.sh [RES] [stop_yr]   # e.g. ./rsf_solve_scaling_test.sh 1km 1000
# Knobs are POSITIONAL ($1,$2) or hardcoded in CONFIG; nothing is read from the
# environment (except $PETSC_DIR for the launcher), so the run is deterministic.
# ---------------------------------------------------------------------------
set -u

# ============================ CONFIG =======================================
RB=../bin/rsf_solve                 # path to rsf_solve binary
BP5=../bp5                          # dir with the BP5 input files
RES=${1-0.5km}                      # $1: resolution tag, e.g. 1km / 0.5km (default 0.5km)
NPLIST="1 2 4 8 16 24 48"       # MPI rank counts to test
BACKENDS="dense htool hacapk hmmvp"  # subset of {dense,htool,hacapk,hmmvp}

# expanded run time so the matvec-dominated phase is well sampled and the
# one-off assembly is amortised.
STOP_YR=${2-1000}                   # $2: run length in years (default 1000)

RTOL=1e-4                         # time-integration rtol (NOT H tol)
HEPS=3e-5                         # HTOOL  -mat_htool_epsilon  (~6.6e-7)
ZTOL=1e-1                         # HACApK -hacapk_ztol        (~2.2e-7)
HMMVP_TOL=1e-7               # hmmvp  -hmmvp_tol          (~1.6e-6)

MPIRUN=$PETSC_DIR/build/bin/mpirun   # set to "mpirun --oversubscribe"
                                           #   to test more ranks than cores
EXTRA_MPI="--bind-to core --map-by core"   # MPI pinning, as in the compress sweep

# --- thread hygiene, applies to ALL runs (mirrors hmat_scaling_test.sh) ---
#  single-thread BLAS: threaded BLAS nested inside block-level parallelism
#  thrashes; hmmvp raises its own pool via -hmmvp_nthreads if asked.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset OMP_PLACES; export OMP_PROC_BIND=false
# If PETSc shared libs aren't on the loader path, set them here:
# export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
# ===========================================================================

GEOM=$BP5/geom_bp5_${RES}.in
RSF=$BP5/rsf_bp5_${RES}.dat
IC=$BP5/ic_bp5_${RES}.in
DC=$BP5/dc_bp5_${RES}.in
for f in "$RB" "$GEOM" "$RSF" "$IC" "$DC"; do
  [ -e "$f" ] || { echo "ERROR: missing $f (check CONFIG)"; exit 1; }
done

COMMON="-geom_file $GEOM -rsf_file $RSF -rsf_ic_file $IC -rsf_dc_file $DC \
 -shear_modulus 3.204e10 -s_wave_speed 3464 -f0 0.6 -dc 0.14 -vpl 1e-9 \
 -v0 1e-6 -sigma_init 25e6 -rtol $RTOL -stop_time_yr $STOP_YR \
 -ts_max_steps 3000000 -print_interval_yr 1e9 -log_view"

flags_for(){  # backend -> rsf_solve flags
  case "$1" in
    dense)  echo "-use_hmatrix 0" ;;
    htool)  echo "-use_hmatrix 1 -mat_htool_epsilon $HEPS" ;;  # ACA is rsf_solve default
    hacapk) echo "-use_hmatrix 3 -hacapk_ztol $ZTOL" ;;
    hmmvp)  echo "-use_hmatrix 4 -hmmvp_tol $HMMVP_TOL" ;;
    *) echo "ERROR: unknown backend $1" >&2; exit 1 ;;
  esac
}

echo "# rsf_solve scaling: RES=$RES stop_yr=$STOP_YR rtol=$RTOL  ($(date))"
echo "# matched ~1e-6 band: htool eps=$HEPS  hacapk ztol=$ZTOL  hmmvp tol=$HMMVP_TOL"
CSV=rsf_solve_scaling_${RES}.csv
echo "backend,np,total_s,assembly_s,step_s,nsteps,matvec_ms,mem_MB" > "$CSV"
printf "%-8s %4s %10s %11s %9s %8s %11s %9s\n" \
       backend np total_s assembly_s step_s nsteps matvec_ms mem_MB

for be in $BACKENDS; do
  for np in $NPLIST; do
    # hmmvp assembly is master/worker: np=1 does no parallel work (and has
    # been seen to fail intermittently on some hosts holding the full
    # undistributed dense reference); it is effectively serial until np>=2.
    log=run_${be}_np${np}.log
    $MPIRUN $EXTRA_MPI -np "$np" "$RB" $COMMON $(flags_for "$be") > "$log" 2>&1
    # --- parse -log_view ---
    total=$(awk '/Time \(sec\):/{print $3; exit}' "$log")
    read mm_c mm_t < <(awk '/^MatMult /{print $2, $4; exit}' "$log")
    read ts_c ts_t < <(awk '/^TSStep /{print $2, $4; exit}' "$log")
    [ -z "${ts_t:-}" ] && ts_t=0
    [ -z "${mm_c:-}" ] && { mm_c=0; mm_t=0; }
    asm=$(awk -v t="$total" -v s="$ts_t" 'BEGIN{printf "%.2f", t-s}')
    mv=$(awk -v t="$mm_t" -v c="$mm_c" 'BEGIN{printf (c>0)?"%.3f":"NA", (c>0)?t/c*1000:0}')
    # --- memory (N read from rsf_solve's 'NNNN patches' banner -> works at any RES) ---
    N=$(awk '/patches/{for(i=1;i<=NF;i++) if($i ~ /^patches/){print $(i-1); exit}}' "$log")
    case "$be" in
      hacapk) mem=$(awk -F= '/Memory of the H-matrix/{gsub(/[^0-9.]/,"",$2);printf "%.1f",$2; exit}' "$log") ;;
      htool)  cr=$(awk -F: '/compression ratio:/{print $2+0; exit}' "$log")
              mem=$(awk -v c="$cr" -v n="$N" 'BEGIN{printf (c>0)?"%.1f":"NA",(c>0)?n*n*8/1048576/c:0}') ;;
      hmmvp)  # 'hmmvp(MPI) 4000 by 4000, 3336045 stored scalars, ...' -> field before 'stored'
              mem=$(awk '/stored scalars/{for(i=1;i<=NF;i++) if($i=="stored"){printf "%.1f",$(i-1)*8/1048576; exit}}' "$log") ;;
      dense)  mem=$(awk -v n="$N" 'BEGIN{printf "%.1f", n*n*8/1048576}') ;;
    esac
    printf "%-8s %4s %10s %11s %9s %8s %11s %9s\n" \
           "$be" "$np" "${total:-NA}" "$asm" "${ts_t:-NA}" "${ts_c:-NA}" "$mv" "${mem:-NA}"
    echo "$be,$np,${total:-NA},$asm,${ts_t:-NA},${ts_c:-NA},$mv,${mem:-NA}" >> "$CSV"
  done
done
echo ""
echo "wrote $CSV  (and per-run logs run_<backend>_np<n>.log)"
echo "matvec_ms is per-MatMult mean (steady-state cost); compare across backends"
echo "at the np you actually run. See compress_interaction_matrix.md for the"
echo "matched-band rationale and the forward-operator scaling these settings came from."
