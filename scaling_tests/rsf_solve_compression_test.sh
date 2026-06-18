#!/usr/bin/env bash
#
# rsf_solve_compression_test.sh -- compression-vs-speed sweep for rsf_solve's
#   H-matrix backends (HTOOL/ACA, HACApK, hmmvp) on the single SEAS BP5 fault.
#
# Companion to rsf_solve_scaling_test.sh. That script fixes the tolerances at a
# matched ~1e-6 band and sweeps MPI ranks (parallel scaling). This one fixes the
# rank count and sweeps a compression knob, so the axis is compression vs speed.
# Two knobs are selectable via SWEEP, and the resolution can be a list so the
# compression-versus-N behaviour can be read across sizes.
#
#   SWEEP=tol  (default): vary each backend's accuracy tolerance
#                htool  -mat_htool_epsilon
#                hacapk -hacapk_ztol
#                hmmvp  -hmmvp_tol
#   SWEEP=eta            : hold tolerance at the matched ~1e-6 band and vary the
#                admissibility eta (the tree-shape knob), larger eta = looser
#                admissibility = more blocks compressed
#                htool  -mat_htool_eta
#                hmmvp  -hmmvp_eta
#                hacapk has NO eta knob, so it appears once at its matched ztol.
#
# Knobs that are NOT reachable here, so they are not swept:
#   - leaf / minimum cluster size: this PETSc build ignores
#     -mat_htool_minclustersize (reported as an unused option), -mat_htool_
#     maxblocksize does not change compression, and hmmvp and hacapk expose no
#     leaf option. Sweeping leaf size would need a change to interact's C option
#     parsing in petsc_interact.c (and possibly the hmmvp/HACApK shims).
#
# Notes on reading the result:
#   - The single BP5 fault is deterministic, so compression, memory, assembly,
#     and per-matvec cost are clean functions of the knob; nsteps is the
#     dynamics-neutrality proxy (flat nsteps => the operator change did not
#     perturb the integration).
#   - This measures compression and speed, NOT forward operator accuracy. For
#     the tolerance sweep the tol->error mapping is in
#     compress_interaction_matrix.md. For the eta sweep there is no such map, so
#     to choose eta on accuracy grounds the cleaner instrument is
#     compress_interaction_matrix, which reports the forward error b-b_h
#     directly; here a large eta that quietly hurts accuracy would only show up
#     if it also moved nsteps.
#   - Tolerances are each package's own knob and are not comparable across
#     packages at equal nominal value (see compress_interaction_matrix.md).
#
# Usage:   ./rsf_solve_compression_test.sh [RES] [np] [stop_yr]
#   RES      resolution tag, or a space-separated list, e.g. "2km 1km 0.5km"
#            (default 1km). A list runs each in turn for the compression-vs-N view.
#   np       fixed MPI rank count for every run  (default 8)
#   stop_yr  run length [yr] (default 50; compression, memory, assembly, and
#            per-matvec cost do not need a full event, so a short interseismic
#            run suffices; only total_s scales with this)
#   SWEEP    environment knob: tol (default) or eta
# ---------------------------------------------------------------------------
set -u

# ============================ CONFIG =======================================
wdir=${HOST:-compr}; mkdir -p "$wdir"; cd "$wdir" || exit 1
RB=../../bin/rsf_solve
BP5=../../bp5
RESLIST=${1:-1km}                   # one tag or a list, e.g. "2km 1km 0.5km"
NP=${2:-8}
STOP_YR=${3:-50}
SWEEP=${SWEEP:-tol}                 # tol | eta

BACKENDS="htool hacapk hmmvp"

# --- SWEEP=tol lists (the compression knob for each) ---
HEPS_LIST="1e-2 1e-3 1e-4 1e-5 1e-6"          # htool  -mat_htool_epsilon
HZTOL_LIST="1e-1 3e-2 1e-2 3e-3 1e-3"         # hacapk -hacapk_ztol
HMTOL_LIST="1e-3 1e-4 1e-5 1e-6 1e-7"         # hmmvp  -hmmvp_tol

# --- SWEEP=eta lists (tolerance held at the matched ~1e-6 band) ---
HETA_LIST="2 5 10 25 50 100"                  # htool  -mat_htool_eta
HMETA_LIST="2 3 5 8"                          # hmmvp  -hmmvp_eta
HEPS_MATCH=3e-5                               # htool  epsilon at matched band
HZTOL_MATCH=1e-1                              # hacapk ztol    at matched band
HMTOL_MATCH=1e-7                              # hmmvp  tol     at matched band

DENSE_REF=${DENSE_REF:-auto}        # dense baseline row? yes/no/auto
                                    #  (auto = yes for 1km/2km, no for 0.5km)
RTOL=1e-4                           # time-integration rtol (NOT an H tol)
MPIRUN=mpirun                       # set "mpirun --oversubscribe" to exceed cores
EXTRA_MPI="--bind-to core --map-by core"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset OMP_PLACES; export OMP_PROC_BIND=false
# export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
# ===========================================================================

[ "$RB" -nt /dev/null ] 2>/dev/null || [ -e "$RB" ] || { echo "ERROR: missing $RB"; exit 1; }
case "$SWEEP" in tol|eta) ;; *) echo "ERROR: SWEEP must be tol or eta"; exit 1 ;; esac
KNOBCOL=$([ "$SWEEP" = eta ] && echo eta || echo tol)

# list of knob values for a backend under the current SWEEP
knob_list(){
  if [ "$SWEEP" = tol ]; then
    case "$1" in htool) echo "$HEPS_LIST";; hacapk) echo "$HZTOL_LIST";; hmmvp) echo "$HMTOL_LIST";; esac
  else
    case "$1" in htool) echo "$HETA_LIST";; hmmvp) echo "$HMETA_LIST";; hacapk) echo "na";; esac
  fi
}
# rsf_solve flags for a backend at a given knob value under the current SWEEP
flags_for(){  # $1 backend  $2 value
  if [ "$SWEEP" = tol ]; then
    case "$1" in
      htool)  echo "-use_hmatrix 1 -mat_htool_epsilon $2" ;;
      hacapk) echo "-use_hmatrix 3 -hacapk_ztol $2" ;;
      hmmvp)  echo "-use_hmatrix 4 -hmmvp_tol $2" ;;
    esac
  else
    case "$1" in
      htool)  echo "-use_hmatrix 1 -mat_htool_epsilon $HEPS_MATCH -mat_htool_eta $2" ;;
      hmmvp)  echo "-use_hmatrix 4 -hmmvp_tol $HMTOL_MATCH -hmmvp_eta $2" ;;
      hacapk) echo "-use_hmatrix 3 -hacapk_ztol $HZTOL_MATCH" ;;
    esac
  fi
}

COMMON=""   # set per resolution inside the loop
run_and_parse(){  # $1=backend  $2=knob-label  $3=flags
  local be="$1" lab="$2" flags="$3"
  local log="run_${RES}_${be}_${SWEEP}${lab}.log"
  $MPIRUN $EXTRA_MPI -np "$NP" "$RB" $COMMON $flags > "$log" 2>&1
  local total mm_c mm_t ts_c ts_t asm mv N denseMB mem compr stored
  total=$(awk '/Time \(sec\):/{print $3; exit}' "$log")
  read mm_c mm_t < <(awk '/^MatMult /{print $2, $4; exit}' "$log")
  read ts_c ts_t < <(awk '/^TSStep /{print $2, $4; exit}' "$log")
  [ -z "${ts_t:-}" ] && ts_t=0; [ -z "${ts_c:-}" ] && ts_c=NA
  [ -z "${mm_c:-}" ] && { mm_c=0; mm_t=0; }
  asm=$(awk -v t="${total:-0}" -v s="$ts_t" 'BEGIN{printf "%.2f", t-s}')
  mv=$(awk -v t="$mm_t" -v c="$mm_c" 'BEGIN{printf (c>0)?"%.3f":"NA",(c>0)?t/c*1000:0}')
  N=$(awk '/patches/{for(i=1;i<=NF;i++) if($i ~ /^patches/){print $(i-1); exit}}' "$log")
  denseMB=$(awk -v n="${N:-0}" 'BEGIN{printf "%.1f", n*n*8/1048576}')
  case "$be" in
    dense)  mem=$denseMB ;;
    hacapk) mem=$(awk -F= '/Memory of the H-matrix/{gsub(/[^0-9.]/,"",$2);printf "%.1f",$2; exit}' "$log") ;;
    htool)  cr=$(awk -F: '/compression ratio:/{print $2+0; exit}' "$log")
            mem=$(awk -v c="$cr" -v d="$denseMB" 'BEGIN{printf (c>0)?"%.1f":"NA",(c>0)?d/c:0}') ;;
    hmmvp)  mem=$(awk '/stored scalars/{for(i=1;i<=NF;i++) if($i=="stored"){printf "%.1f",$(i-1)*8/1048576; exit}}' "$log") ;;
  esac
  if [ -n "${mem:-}" ] && [ "$mem" != NA ] && [ "${denseMB%.*}" -gt 0 ] 2>/dev/null; then
    compr=$(awk -v d="$denseMB" -v m="$mem" 'BEGIN{printf (m>0)?"%.1f":"NA",(m>0)?d/m:0}')
    stored=$(awk -v d="$denseMB" -v m="$mem" 'BEGIN{printf (d>0)?"%.2f":"NA",(d>0)?100*m/d:0}')
  else compr=NA; stored=NA; fi
  printf "%-7s %6s %8s %8s %9s %11s %10s %10s %7s\n" \
         "$be" "$lab" "$compr" "$stored" "${mem:-NA}" "$asm" "$mv" "${total:-NA}" "$ts_c"
  echo "$be,$lab,$compr,$stored,${mem:-NA},$asm,$mv,${total:-NA},$ts_c" >> "$CSV"
}

echo "# rsf_solve compression vs speed: SWEEP=$SWEEP RES=\"$RESLIST\" np=$NP stop_yr=$STOP_YR rtol=$RTOL  ($(date))"
[ "$SWEEP" = eta ] && echo "# eta sweep holds tol at matched band: htool eps=$HEPS_MATCH, hacapk ztol=$HZTOL_MATCH, hmmvp tol=$HMTOL_MATCH (hacapk has no eta knob)"
echo "# compr_x = dense / H-matrix memory ; stored% = 100 / compr_x ; matvec_ms = mean per MatMult"

for RES in $RESLIST; do
  GEOM=$BP5/geom_bp5_${RES}.in; RSF=$BP5/rsf_bp5_${RES}.dat
  IC=$BP5/ic_bp5_${RES}.in;     DC=$BP5/dc_bp5_${RES}.in
  for f in "$GEOM" "$RSF" "$IC" "$DC"; do
    [ -e "$f" ] || { echo "ERROR: missing $f (skipping RES=$RES)"; continue 2; }
  done
  dref=$DENSE_REF; [ "$dref" = auto ] && { [ "$RES" = "0.5km" ] && dref=no || dref=yes; }
  COMMON="-geom_file $GEOM -rsf_file $RSF -rsf_ic_file $IC -rsf_dc_file $DC \
 -shear_modulus 3.204e10 -s_wave_speed 3464 -f0 0.6 -dc 0.14 -vpl 1e-9 \
 -v0 1e-6 -sigma_init 25e6 -rtol $RTOL -stop_time_yr $STOP_YR \
 -ts_max_steps 3000000 -print_interval_yr 1e9 -log_view"

  CSV=rsf_solve_compression_${SWEEP}_${RES}.csv
  echo "backend,$KNOBCOL,compr_x,stored_pct,mem_MB,assembly_s,matvec_ms,total_s,nsteps" > "$CSV"
  echo ""
  echo "## RES=$RES  (SWEEP=$SWEEP, np=$NP)"
  printf "%-7s %6s %8s %8s %9s %11s %10s %10s %7s\n" \
         backend "$KNOBCOL" compr_x stored% mem_MB assembly_s matvec_ms total_s nsteps
  [ "$dref" = yes ] && run_and_parse dense "-" "-use_hmatrix 0"
  for be in $BACKENDS; do
    for v in $(knob_list "$be"); do run_and_parse "$be" "$v" "$(flags_for "$be" "$v")"; done
  done
  echo "wrote $CSV"
  [ -n "${HOST:-}" ] && cp "$CSV" "../rsf_solve_compression_${SWEEP}_${RES}.$HOST.csv"
done

echo ""
echo "reading it: within a backend, looser tol or larger eta gives more"
echo "compression (larger compr_x, less memory). Across the RES list, compr_x at"
echo "fixed accuracy is the compression-versus-N curve. matvec_ms is the steady"
echo "per-apply cost that dominates a long cycle; assembly_s is paid once."
cd ..
