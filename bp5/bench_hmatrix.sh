#!/usr/bin/env bash
#
# bench_hmatrix.sh -- MPI scaling benchmark for rsf_solve's H-matrix backends
#                     (dense, HTOOL/ACA, HACApK) on the SEAS BP5 1 km case.
#
# Sweeps a list of MPI rank counts x backends, runs rsf_solve with -log_view,
# and parses assembly time, mean matvec time, total wallclock, #steps, and
# H-matrix memory into a table (also written as CSV).
#
# Usage:   ./bench_hmatrix.sh [stop_yr]      # default stop_yr below
# Edit the CONFIG block for your paths / core list / inputs.
# ---------------------------------------------------------------------------
set -u

# ============================ CONFIG =======================================
RB=${RB:-../bin/rsf_solve}                 # path to rsf_solve binary
BP5=${BP5:-../bp5}                          # dir with the BP5 input files
#RES=${RES:-1km}                            # resolution tag: 1km (4000) or 2km (1000)
RES=${RES:-0.5km}                            # resolution tag: 1km (4000) or 2km (1000)
NPLIST=${NPLIST:-"1 2 4 8 16 24 48"}          # MPI rank counts to test
BACKENDS=${BACKENDS:-"dense htool hacapk"} # subset of {dense,htool,hacapk}
STOP_YR=${1:-${STOP_YR:-60}}              # short run: assembly + steady matvec timing
RTOL=${RTOL:-1e-4}
ZTOL=${ZTOL:-1e-4}                         # HACApK -hacapk_ztol
HEPS=${HEPS:-1e-4}                         # HTOOL  -mat_htool_epsilon
MPIRUN=${MPIRUN:-$PETSC_DIR/build/bin/mpirun}                   # set to "mpirun --oversubscribe" to test
                                           #   more ranks than physical cores
EXTRA_MPI=${EXTRA_MPI:-}                   # e.g. "--bind-to core --map-by core"
# If PETSc shared libs aren't on the loader path, set them here:
# export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# ===========================================================================

GEOM=$BP5/geom_bp5_${RES}.in
RSF=$BP5/rsf_bp5_${RES}.dat
IC=$BP5/ic_bp5_${RES}.in
DC=$BP5/dc_bp5_${RES}.in
for f in "$RB" "$GEOM" "$RSF" "$IC" "$DC"; do
  [ -e "$f" ] || { echo "ERROR: missing $f (check CONFIG)"; exit 1; }
done

NCELL=$(grep -cve '^[[:space:]]*$' "$GEOM")   # patches = non-blank lines in geom file

COMMON="-geom_file $GEOM -rsf_file $RSF -rsf_ic_file $IC -rsf_dc_file $DC \
 -shear_modulus 3.204e10 -s_wave_speed 3464 -f0 0.6 -dc 0.14 -vpl 1e-9 \
 -v0 1e-6 -sigma_init 25e6 -rtol $RTOL -stop_time_yr $STOP_YR \
 -ts_max_steps 300000 -print_interval_yr 1e9 -log_view"

flags_for(){  # backend -> rsf_solve flags
  case "$1" in
    dense)  echo "-use_hmatrix 0" ;;
    htool)  echo "-use_hmatrix 1 -mat_htool_epsilon $HEPS" ;;  # ACA is rsf_solve default
    hacapk) echo "-use_hmatrix 3 -hacapk_ztol $ZTOL" ;;
    *) echo "ERROR: unknown backend $1" >&2; exit 1 ;;
  esac
}

CSV=bench_hmatrix_${RES}.csv
echo "backend,np,total_s,assembly_s,step_s,nsteps,matvec_ms,mem_MB" > "$CSV"
printf "%-8s %4s %10s %11s %9s %8s %11s %9s\n" \
       backend np total_s assembly_s step_s nsteps matvec_ms mem_MB

DENSE_MB=""   # filled from the np=1 dense run for HTOOL %-of-dense if desired
for be in $BACKENDS; do
  for np in $NPLIST; do
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
    # --- memory ---
    case "$be" in
      hacapk) mem=$(awk -F= '/Memory of the H-matrix/{gsub(/[^0-9.]/,"",$2);printf "%.1f",$2; exit}' "$log") ;;
      htool)  cr=$(awk -F: '/compression ratio:/{print $2+0; exit}' "$log")
              mem=$(awk -v c="$cr" -v n="$NCELL" 'BEGIN{d=n*n*8/1048576; printf (c>0)?"%.1f":"NA",(c>0)?d/c:0}') ;;
      dense)  mem=$(awk -v n="$NCELL" 'BEGIN{printf "%.1f", n*n*8/1048576}') ;;
    esac
    printf "%-8s %4s %10s %11s %9s %8s %11s %9s\n" \
           "$be" "$np" "${total:-NA}" "$asm" "${ts_t:-NA}" "${ts_c:-NA}" "$mv" "${mem:-NA}"
    echo "$be,$np,${total:-NA},$asm,${ts_t:-NA},${ts_c:-NA},$mv,${mem:-NA}" >> "$CSV"
  done
done
echo ""
echo "wrote $CSV  (and per-run logs run_<backend>_np<n>.log)"
