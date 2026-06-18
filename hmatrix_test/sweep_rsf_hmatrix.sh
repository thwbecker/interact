#!/usr/bin/env bash
#
# sweep_rsf_hmatrix.sh -- two-fault rsf_solve accuracy-vs-speed sweep for the
#   HTOOLS, HACAPK, and HMMVP H-matrix backends (BigWham deliberately excluded).
#
# Builds an interacting two-fault BP5-style problem (make_two_fault.py), runs a
# DENSE reference once, then runs each H-matrix backend across a range of its
# accuracy tolerance.  For every run it records speed (assembly, matvec, total,
# step count from -log_view) and accuracy versus the dense solve (max and RMS
# deviation of the log10(max|V|) trace, and event-timing differences, via
# rsf_accuracy.py).  Results go to a table on screen and a CSV.
#
# Usage:   ./sweep_rsf_hmatrix.sh [ds_km] [stop_yr] [procs]
#   ds_km    cell size for both faults        (default 2.0  -> 2000 patches)
#   stop_yr  simulated time                    (default 260  -> past 1st event)
#   procs    parallelism per run               (default 4)
# Edit the CONFIG block for paths, separation, tolerance lists, and MPIRUN.
#
# Notes on parallelism: the MPI backends (dense, HTOOLS, HACAPK) run on `procs`
# MPI ranks; HMMVP, which is OpenMP, runs on 1 rank with `procs` threads.  This
# gives each backend procs-way parallelism in its native mode.  On a single-core
# box set procs=1 (or set MPIRUN to oversubscribe) for a functional check; the
# timings are only meaningful with real cores.
# ---------------------------------------------------------------------------
set -u

# ============================ CONFIG =======================================
HERE=$(cd "$(dirname "$0")" && pwd)
RB=${RB:-$HERE/../bin/rsf_solve}            # path to rsf_solve binary
GEN=${GEN:-$HERE/make_two_fault.py}
ACC=${ACC:-$HERE/rsf_accuracy.py}

DS=${1:-${DS:-2.0}}                         # cell size [km]
STOP_YR=${2:-${STOP_YR:-260}}               # simulated time [yr]
PROCS=${3:-${PROCS:-4}}                     # parallelism per run
SEP=${SEP:-4.0}                             # fault-normal separation [km]
LX=${LX:-100}                               # along-strike length [km]
LD=${LD:-40}                                # depth extent [km]

RTOL=${RTOL:-1e-4}                          # ODE relative tolerance
DT_MON=${DT_MON:-2.0}                       # monitor cadence [yr] (also caps step)
VEL_EVENT=${VEL_EVENT:-1e-3}                # event slip-rate threshold [m/s]
MAXSTEPS=${MAXSTEPS:-400000}
RUN_TIMEOUT=${RUN_TIMEOUT:-0}               # per-run wall limit [s]; 0 = no limit.
                                            #  Guards against a backend whose matvec
                                            #  is slow in a given environment (e.g.
                                            #  HMMVP streams its operator from disk,
                                            #  so its throughput depends on the
                                            #  filesystem) stalling the whole sweep.

# tolerance lists per backend (the primary accuracy knob for each)
HEPS_LIST=${HEPS_LIST:-"1e-2 1e-3 1e-4 1e-6"}    # HTOOLS -mat_htool_epsilon
HZTOL_LIST=${HZTOL_LIST:-"1e-2 1e-3 1e-4 1e-6"}  # HACAPK -hacapk_ztol
HMTOL_LIST=${HMTOL_LIST:-"1e-2 1e-3 1e-5"}       # HMMVP  -hmmvp_tol
BACKENDS=${BACKENDS:-"dense htool hacapk hmmvp"} # subset to run

MPIRUN=${MPIRUN:-mpirun}                    # e.g. ibrun / srun on a cluster,
                                            #  "mpirun --oversubscribe --allow-run-as-root" locally
# export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
# ===========================================================================

command -v python3 >/dev/null || { echo "ERROR: need python3"; exit 1; }
[ -x "$RB" ] || { echo "ERROR: rsf_solve not found/executable at $RB"; exit 1; }

# match the generator's %g-formatted filename tag (e.g. 4.0 -> 4, 0.5 -> 0.5)
DSG=$(printf '%g' "$DS"); SEPG=$(printf '%g' "$SEP")
TAG="${DSG}km_${SEPG}km"
OUT="$HERE/out_2f_${DS}km_sep${SEP}km_${STOP_YR}yr_np${PROCS}"
mkdir -p "$OUT"; cd "$OUT" || exit 1

# --- inputs ---------------------------------------------------------------
GEOM="$OUT/geom_2f_${TAG}.in"
RSF="$OUT/rsf_2f_${TAG}.dat"
IC="$OUT/ic_2f_${TAG}.in"
DC="$OUT/dc_2f_${TAG}.in"
if [ ! -s "$GEOM" ]; then
  echo "generating two-fault inputs (ds=$DS sep=$SEP Lx=$LX Ld=$LD) ..."
  python3 "$GEN" "$DS" "$SEP" "$LX" "$LD" || { echo "ERROR: generator failed"; exit 1; }
fi
NCELL=$(grep -cve '^[[:space:]]*$' "$GEOM" 2>/dev/null)
[ -s "$GEOM" ] && [ -n "$NCELL" ] || { echo "ERROR: geometry $GEOM not found after generation"; exit 1; }
echo "two-fault problem: $NCELL patches, stop=$STOP_YR yr, rtol=$RTOL, dt_mon=$DT_MON yr, procs=$PROCS"

COMMON="-geom_file $GEOM -rsf_file $RSF -rsf_ic_file $IC -rsf_dc_file $DC \
 -shear_modulus 3.204e10 -s_wave_speed 3464 -f0 0.6 -dc 0.14 -vpl 1e-9 \
 -v0 1e-6 -sigma_init 25e6 -rtol $RTOL -stop_time_yr $STOP_YR \
 -dt_monitor_yr $DT_MON -vel_event $VEL_EVENT -track_events \
 -ts_max_steps $MAXSTEPS -print_interval_yr 1e9 -log_view"

# backend -> rsf_solve type flag + per-run launch parallelism
type_flag(){ case "$1" in
  dense) echo "-use_hmatrix 0" ;; htool) echo "-use_hmatrix 1" ;;
  hacapk) echo "-use_hmatrix 3" ;; hmmvp) echo "-use_hmatrix 4" ;;
esac; }

run_one(){  # $1=cfgname  $2=backend  $3=tol-flags
  local name="$1" be="$2" tf="$3" dir="$OUT/$1"
  mkdir -p "$dir"; ( cd "$dir" || exit 1
    local np=$PROCS omp=1 extra="" to=""
    if [ "$be" = hmmvp ]; then np=1; omp=$PROCS; extra="-hmmvp_nthreads $PROCS"; fi
    [ "${RUN_TIMEOUT:-0}" != 0 ] && to="timeout ${RUN_TIMEOUT}"
    OMP_NUM_THREADS=$omp OPENBLAS_NUM_THREADS=1 \
      $to $MPIRUN -np "$np" "$RB" $COMMON $(type_flag "$be") $tf $extra > run.log 2>&1
  )
}

# --- parse -log_view from a run.log into "total asm matvec_ms nsteps mem" ---
parse_speed(){
  local log="$1/run.log" be="$2"
  [ -s "$log" ] || { echo "NA NA NA NA NA"; return; }
  local total mm_c mm_t ts_c ts_t asm mv mem
  total=$(awk '/Time \(sec\):/{print $3; exit}' "$log")
  read mm_c mm_t < <(awk '/^MatMult /{print $2, $4; exit}' "$log")
  read ts_c ts_t < <(awk '/^TSStep /{print $2, $4; exit}' "$log")
  [ -z "${total:-}" ] && { echo "NA NA NA NA NA"; return; }
  [ -z "${ts_t:-}" ] && ts_t=0; [ -z "${ts_c:-}" ] && ts_c=NA
  [ -z "${mm_c:-}" ] && { mm_c=0; mm_t=0; }
  asm=$(awk -v t="$total" -v s="$ts_t" 'BEGIN{printf "%.2f", t-s}')
  mv=$(awk -v t="$mm_t" -v c="$mm_c" 'BEGIN{printf (c>0)?"%.3f":"NA",(c>0)?t/c*1000:0}')
  case "$be" in
    dense)  mem=$(awk -v n="$NCELL" 'BEGIN{printf "%.1f", n*n*8/1048576}') ;;
    hacapk) mem=$(awk -F= '/Memory of the H-matrix/{gsub(/[^0-9.]/,"",$2);printf "%.1f",$2; exit}' "$log") ;;
    htool)  cr=$(awk -F: '/compression ratio:/{print $2+0; exit}' "$log")
            mem=$(awk -v c="$cr" -v n="$NCELL" 'BEGIN{d=n*n*8/1048576; printf (c>0)?"%.1f":"NA",(c>0)?d/c:0}') ;;
    *)      mem=NA ;;
  esac
  echo "${total} ${asm} ${mv} ${ts_c} ${mem:-NA}"
}

CSV="$OUT/sweep_summary.csv"
echo "backend,setting,total_s,assembly_s,matvec_ms,nsteps,mem_MB,max_dlog10V,rms_dlog10V,d_first_event_yr,d_recurrence_yr,n_events" > "$CSV"
hdr(){ printf "%-7s %-9s %9s %10s %9s %7s %8s %11s %11s %12s %12s %4s\n" \
  backend setting total_s asm_s mv_ms nstep mem_MB maxdV rmsdV dEv1_yr dRec_yr nEv; }
row(){ printf "%-7s %-9s %9s %10s %9s %7s %8s %11s %11s %12s %12s %4s\n" "$@"; }

echo ""; echo "=== running DENSE reference ==="
run_one dense dense ""
read D_total D_asm D_mv D_steps D_mem < <(parse_speed "$OUT/dense" dense)
echo "dense: total=${D_total}s assembly=${D_asm}s matvec=${D_mv}ms steps=${D_steps} mem=${D_mem}MB"
DREF="$OUT/dense"
# dense self-row (zero deviation by construction)
declare -a ROWS
ROWS+=("dense - $D_total $D_asm $D_mv $D_steps $D_mem 0 0 0 0 -")

# --- H-matrix backends across their tolerance knob -------------------------
sweep_backend(){  # $1=backend  $2=tolflag-name  $3="list"
  local be="$1" flagname="$2"; shift 2
  for tol in $@; do
    local name="${be}_${tol}"
    echo ""; echo "=== $be  $flagname=$tol ==="
    run_one "$name" "$be" "$flagname $tol"
    read total asm mv steps mem < <(parse_speed "$OUT/$name" "$be")
    local accline; accline=$(python3 "$ACC" "$DREF" "$OUT/$name" 2>>"$OUT/accuracy.log")
    # accline: maxdev=.. rmsdev=.. dt1=.. drec=.. nev_ref=.. nev_run=.. rec_ref=..
    eval "$accline"
    echo "  speed: total=${total}s asm=${asm}s mv=${mv}ms steps=${steps} mem=${mem}MB"
    ROWS+=("$be $tol $total $asm $mv $steps $mem $maxdev $rmsdev $dt1 $drec $nev_run")
  done
}

for be in $BACKENDS; do
  case "$be" in
    dense)  : ;;  # already done
    htool)  sweep_backend htool  -mat_htool_epsilon $HEPS_LIST ;;
    hacapk) sweep_backend hacapk -hacapk_ztol      $HZTOL_LIST ;;
    hmmvp)  sweep_backend hmmvp  -hmmvp_tol         $HMTOL_LIST ;;
  esac
done

echo ""; echo "================================================================================"
echo " two-fault rsf_solve sweep: accuracy (vs dense) and speed   N=$NCELL  stop=$STOP_YR yr"
echo "   dense reference recurrence basis; deviations are run minus dense, log10 units / yr"
echo "================================================================================"
hdr
for r in "${ROWS[@]}"; do row $r; echo "$r" | tr ' ' ',' >> "$CSV"; done
# (CSV header was written above with named columns; per-row append keeps order)
echo ""
echo "wrote $CSV"
echo "per-run logs + rsf_monitor.dat / rsf_events.dat under $OUT/<config>/"
echo ""
echo "reading the table: at matched tolerance the deviation columns (maxdV, rmsdV,"
echo "dEv1_yr, dRec_yr) say how close each backend tracks the dense sequence; the"
echo "speed columns (asm_s, mv_ms, total_s) say what it costs.  The best setting per"
echo "package is the loosest tolerance whose deviation is still in the noise of the"
echo "ODE rtol (here ${RTOL}); tightening past that point usually buys accuracy you"
echo "cannot see while still costing assembly time.  Confirm that holds for this"
echo "interacting two-fault case rather than assuming it from the single-fault note."
