#!/usr/bin/env bash
#
# rsf_solve_compression_test.sh -- compression-vs-speed sweep for rsf_solve's
#   H-matrix backends (HTOOL/ACA, HACApK, hmmvp) on the single SEAS BP5 fault.
#
# This is the companion to rsf_solve_scaling_test.sh.  That script fixes the
# tolerances at a matched ~1e-6 band and sweeps MPI ranks to compare PARALLEL
# SCALING.  This script instead fixes the rank count and sweeps each backend's
# accuracy tolerance, so the axis is COMPRESSION vs SPEED: how aggressively each
# package compresses the operator as its tolerance is loosened, and what that
# buys (and costs) in memory, assembly time, and per-matvec time.
#
# The single BP5 fault is deterministic and reproducible, so unlike the
# two-fault case there is no chaotic-recurrence confound; compression, memory,
# assembly, and the steady per-MatMult cost are clean functions of the
# tolerance.
#
# Tolerance is each package's own knob and they are NOT comparable across
# packages at equal nominal value: the same number means different operator
# accuracy for each (see compress_interaction_matrix.md for the tol->error
# mapping at N=14400).  In particular HACApK ztol is conservative for the smooth
# Okada kernel: ztol ~1e-1 already gives ~1e-6 accuracy, and tightening toward
# 1e-4 drives it near-dense (little compression).  Read each backend's curve on
# its own; cross-backend comparison should be done at matched accuracy, which is
# what rsf_solve_scaling_test.sh is for.
#
# Usage:   ./rsf_solve_compression_test.sh [RES] [np] [stop_yr]
#   RES      resolution tag: 2km / 1km / 0.5km   (default 1km)
#   np       fixed MPI rank count for every run  (default 8)
#   stop_yr  run length                          (default 50; compression,
#            memory, assembly, and per-matvec cost do not need a full event, so
#            a short interseismic run suffices and is much faster; only total_s
#            scales with this).  The slow runs in the sweep are the TIGHT-tol
#            ones (hacapk 3e-3/1e-3, htool 1e-6, hmmvp 1e-7): they barely
#            compress (near-dense), so trim them if you only want the compressed
#            regime.
# Knobs are positional or hardcoded in CONFIG; nothing else is read from the
# environment (except PETSC for the launcher) so the run is deterministic.
# ---------------------------------------------------------------------------
set -u

# ============================ CONFIG =======================================
wdir=${HOST:-compr}                 # working subdirectory
mkdir -p "$wdir"; cd "$wdir" || exit 1
RB=../../bin/rsf_solve              # path to rsf_solve binary
BP5=../../bp5                       # dir with the BP5 input files
RES=${1:-1km}                       # resolution tag
NP=${2:-8}                          # fixed MPI rank count
STOP_YR=${3:-50}                    # run length [yr]; 50 is plenty since

# per-backend tolerance sweep (the compression knob for each).  Loose -> tight,
# i.e. most -> least compression.
BACKENDS="htool hacapk hmmvp"       # subset of {htool,hacapk,hmmvp}
HEPS_LIST="1e-2 1e-3 1e-4 1e-5 1e-6"          # HTOOL  -mat_htool_epsilon
HZTOL_LIST="1e-1 3e-2 1e-2 3e-3 1e-3"         # HACApK -hacapk_ztol (1e-1 ~matched
                                              #  band; tighter -> near-dense)
HMTOL_LIST="1e-3 1e-4 1e-5 1e-6 1e-7"         # hmmvp  -hmmvp_tol (floors ~1e-6 at
                                              #  moderate N; tighter buys little)
DENSE_REF=${DENSE_REF:-auto}        # run a dense baseline row? yes/no/auto
                                    #  (auto = yes for 1km/2km, no for 0.5km
                                    #   where the dense N^2 operator is heavy)

RTOL=1e-4                           # time-integration rtol (NOT an H tol)
MPIRUN=mpirun                       # set "mpirun --oversubscribe" to exceed cores
EXTRA_MPI="--bind-to core --map-by core"      # MPI pinning, as in the sibling script

# thread hygiene (single-thread BLAS; hmmvp raises its own pool only if asked)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset OMP_PLACES; export OMP_PROC_BIND=false
# export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
# ===========================================================================

GEOM=$BP5/geom_bp5_${RES}.in; RSF=$BP5/rsf_bp5_${RES}.dat
IC=$BP5/ic_bp5_${RES}.in;     DC=$BP5/dc_bp5_${RES}.in
for f in "$RB" "$GEOM" "$RSF" "$IC" "$DC"; do
  [ -e "$f" ] || { echo "ERROR: missing $f (check CONFIG)"; exit 1; }
done
[ "$DENSE_REF" = auto ] && { [ "$RES" = "0.5km" ] && DENSE_REF=no || DENSE_REF=yes; }

COMMON="-geom_file $GEOM -rsf_file $RSF -rsf_ic_file $IC -rsf_dc_file $DC \
 -shear_modulus 3.204e10 -s_wave_speed 3464 -f0 0.6 -dc 0.14 -vpl 1e-9 \
 -v0 1e-6 -sigma_init 25e6 -rtol $RTOL -stop_time_yr $STOP_YR \
 -ts_max_steps 3000000 -print_interval_yr 1e9 -log_view"

flags_for(){  # backend + tol -> rsf_solve flags
  case "$1" in
    dense)  echo "-use_hmatrix 0" ;;
    htool)  echo "-use_hmatrix 1 -mat_htool_epsilon $2" ;;
    hacapk) echo "-use_hmatrix 3 -hacapk_ztol $2" ;;
    hmmvp)  echo "-use_hmatrix 4 -hmmvp_tol $2" ;;
  esac
}

echo "# rsf_solve compression vs speed: RES=$RES np=$NP stop_yr=$STOP_YR rtol=$RTOL  ($(date))"
echo "# compr_x = dense / H-matrix memory ; stored% = 100 / compr_x ; matvec_ms = mean per MatMult"
CSV=rsf_solve_compression_${RES}.csv
echo "backend,tol,compr_x,stored_pct,mem_MB,assembly_s,matvec_ms,total_s,nsteps" > "$CSV"
printf "%-7s %6s %8s %8s %9s %11s %10s %10s %7s\n" \
       backend tol compr_x stored% mem_MB assembly_s matvec_ms total_s nsteps

run_and_parse(){  # $1=backend  $2=tol(label "-")  $3=flags
  local be="$1" tol="$2" flags="$3"
  local log="run_${be}_${tol}.log"
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
         "$be" "$tol" "${compr}" "${stored}" "${mem:-NA}" "$asm" "$mv" "${total:-NA}" "${ts_c}"
  echo "$be,$tol,$compr,$stored,${mem:-NA},$asm,$mv,${total:-NA},${ts_c}" >> "$CSV"
}

[ "$DENSE_REF" = yes ] && run_and_parse dense "-" "$(flags_for dense)"
for be in $BACKENDS; do
  case "$be" in
    htool)  list="$HEPS_LIST" ;;
    hacapk) list="$HZTOL_LIST" ;;
    hmmvp)  list="$HMTOL_LIST" ;;
  esac
  for tol in $list; do run_and_parse "$be" "$tol" "$(flags_for "$be" "$tol")"; done
done

echo ""
echo "wrote $CSV  (and per-run logs run_<backend>_<tol>.log)"
echo "reading it: within a backend, looser tol -> larger compr_x (less memory) and"
echo "usually faster matvec_ms, at the cost of operator accuracy.  The efficient"
echo "setting is the loosest tol still inside your target accuracy band (see"
echo "compress_interaction_matrix.md for tol->error).  matvec_ms is the steady"
echo "per-apply cost that dominates a long cycle; assembly_s is paid once."
[ -n "${HOST:-}" ] && cp "$CSV" "../rsf_solve_compression_${RES}.$HOST.csv"
cd ..
