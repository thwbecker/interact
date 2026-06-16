#!/usr/bin/env bash
#
# hbi_bp5_scaling_test.sh -- MPI scaling benchmark for HBI (sozawa94/hbi) on the
#   SEAS BP5 rectangular case, as a companion to rsf_solve_scaling_test.sh so the
#   two codes can be compared head-to-head on the same problem.
#
# HBI is HACApK-only (lattice H-matrix). This runs it "as configured from GitHub"
# plus the one-line HACApK bounds fix bp5/hbi_bp5_fix.patch
#   (m_HACApK_use.f90: do il=max(ndnr_s,1),ndnr_e  -- guards a rank that owns zero
#    near-field rows in the param(61)==3 diagonal-scaling path).
#
# Sweeps MPI rank counts, runs HBI on bp5r.in (imax*jmax = N elements), and parses
# HBI's stdout for assembly, total matvec, and total wall time + step count into a
# table/CSV laid out like the rsf_solve one.
#
# ===================  TWO THINGS TO GET RIGHT FOR A FAIR COMPARISON  ==========
#
# (1) H-MATRIX TOLERANCE.  HBI's ACA tolerance is eps_h, HARDCODED at
#     main_LH.f90:145  (eps_r=1d-4; eps_h=1d-4).  It is NOT a .in keyword.
#     For the smooth BP5 Okada kernel, eps_h=1d-4 is ~NEAR-DENSE (cf. interact's
#     HACApK ztol 1e-4 -> 82% of dense, err ~2e-10): HBI then does far more matvec
#     work than its accuracy needs, exactly as rsf_solve-HACApK did before we
#     loosened ztol.  So:
#       - to reproduce rsf_solve's matched ~1e-6 band (HACApK ztol 1e-1 -> 2.2e-7),
#         set  eps_h = 1d-1  in main_LH.f90:145 and recompile;
#       - to compare "out of the box" at HBI's default accuracy, leave eps_h=1d-4
#         AND run rsf_solve-HACApK at -hacapk_ztol 1e-4 so both are near-dense.
#     Either way, compare HBI only against an rsf_solve run at the SAME accuracy.
#     Record which you used in HBI_EPS_H below (label only; it cannot set eps_h).
#
# (2) PHYSICAL PARAMETERS must match rsf_solve's BP5 inputs.  In bp5r.in set the
#     rigidity / S-wave speed to rsf_solve's (-shear_modulus 3.204e10 Pa,
#     -s_wave_speed 3464 m/s; BP5 spec mu=32.04 GPa, Vs=3.464 km/s).  HBI >=2026.1.0
#     reads rigidity/Vs/poisson from the .in file; older builds may hardcode them.
#     f0=0.6, b=0.03, dc=0.14, sigma=25 MPa, a from bp5r_param.dat already match BP5.
# =============================================================================
#
# Usage:   ./hbi_bp5_scaling_test.sh [RES] [tmax_yr] [eps_h_label]
#          e.g.  ./hbi_bp5_scaling_test.sh 1km 1000 1d-1
# Knobs are POSITIONAL ($1,$2,$3) or hardcoded in CONFIG; nothing is read from the
# environment (except $PETSC_DIR for the launcher), so the run is deterministic.
# Needs ./examples/bp5r_<RES>.in and bp5r_<RES>_param.dat (make with bp5r_param_gen.py).
#
# eps_h ($3) is the HACApK ACA tolerance; it is injected as an .in keyword (read at
# main_LH.f90:171, before matrix generation), so it overrides the hardcoded default
# with no recompile. eps_h=1d-1 gives ~12.3 MB at N=4000, matching rsf-hacapk's 12.1 MB
# (stock 1d-4 is near-dense, ~25.3 MB). Build step needs no eps_h argument anymore.
#
# IMPORTANT: HBI's lattice H-matrix needs a SQUARE MPI process grid (npgl==npgt),
# so the rank count MUST be a perfect square (1,4,9,16,25,36,...). Non-squares abort
# with 'Process grid must have square shape!'. NPLIST is set accordingly; the loop
# also skips any non-square np defensively.
# ---------------------------------------------------------------------------
set -u

# ============================ CONFIG =======================================
RES=${1-1km}                                # $1: resolution tag: 1km (N=4000) or 0.5km (N=16000)
TMAX_YR=${2-1000}                           # $2: run length in years (default 1000; align with rsf_solve)
HBI_EPS_H=${3-1d-1}                          # $3: ACA tolerance, injected into the .in (real knob; 1d-1 ~= rsf ztol 1e-1)

HBI="./hbi/lhbiem"                          # path to the patched HBI binary (Makefile TARGET=lhbiem)
INTPL="./examples/bp5r_${RES}.in"       # BP5 rectangular .in for this resolution
PARAM="./examples/bp5r_${RES}_param.dat" # matching per-element params (rows == imax*jmax)
NPLIST="1 4 9 16 25 36"                     # MPI ranks: MUST be perfect squares (see note below)
NSTEP=2000000                               # high step cap so tmax is the stop criterion
INTERVAL=100000                             # large so output/checkpoint I/O does not contaminate timing

MPIRUN=$PETSC_DIR/build/bin/mpirun          # launcher ($PETSC_DIR is the one toolchain env input)
EXTRA_MPI="--bind-to core --map-by core"    # MPI pinning, as in the other sweeps

# HBI is hybrid MPI/OpenMP; pin one thread per rank for a pure-MPI sweep (to line
# up with the rsf_solve sweep). To explore HBI's OpenMP axis, edit to >1 and drop ranks.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset OMP_PLACES; export OMP_PROC_BIND=false
# ===========================================================================

for f in "$HBI" "$INTPL" "$PARAM"; do
  [ -e "$f" ] || { echo "ERROR: missing $f (check CONFIG)"; exit 1; }
done

# element count from the template (imax * jmax) for the memory/normalisation columns
IMAX=$(awk '/^imax/{print $2; exit}' "$INTPL"); JMAX=$(awk '/^jmax/{print $2; exit}' "$INTPL")
N=$(( IMAX * JMAX ))

# build a per-run .in: copy template, override tmax / nstep / interval, keep the
# parameter_file pointing at $PARAM (HBI resolves it relative to its run dir).
RUNIN=hbi_bp5_run.in
# eps_h IS a real knob: read_inputfile() (main_LH.f90:171) runs before HACApK_generate
# (493), so an 'eps_h' line in the .in overrides the hardcoded default (line 145) WITHOUT
# recompiling. We inject HBI_EPS_H here (replacing any existing eps_h line, else appending).
make_runin(){
  awk -v tmax="$TMAX_YR" -v nstep="$NSTEP" -v intv="$INTERVAL" -v pf="$PARAM" -v eh="$HBI_EPS_H" '
    /^tmax /     {print "tmax", tmax, "!set by scaling script"; next}
    /^nstep /    {print "nstep", nstep, "!set by scaling script"; next}
    /^interval / {print "interval", intv, "!set by scaling script"; next}
    /^parameter_file /{print "parameter_file \"" pf "\""; next}
    /^eps_h /    {print "eps_h", eh, "!set by scaling script"; seen=1; next}
    {print}
    END{ if(!seen) print "eps_h", eh, "!set by scaling script" }
  ' "$INTPL" > "$RUNIN"
}
make_runin

echo "# HBI BP5 scaling: RES=$RES N=$N (imax=$IMAX jmax=$JMAX) tmax=$TMAX_YR yr eps_h=$HBI_EPS_H  ($(date))"
echo "# NOTE eps_h set via .in keyword (overrides hardcoded default; no recompile needed)"
CSV=hbi_bp5_scaling_${RES}.csv
echo "code,np,N,eps_h,total_s,assembly_s,matvec_s,nsteps,matvec_ms_per_step,mem_MB" > "$CSV"
printf "%-5s %4s %6s %7s %9s %11s %10s %8s %16s %8s\n" \
       code np N eps_h total_s assembly_s matvec_s nsteps mv_ms/step mem_MB

for np in $NPLIST; do
  # HBI's lattice H-matrix requires a SQUARE process grid (npgl==npgt), so np must
  # be a perfect square; non-squares abort with 'Process grid must have square shape!'.
  sq=$(awk -v n="$np" 'BEGIN{s=int(sqrt(n)+0.5); print (s*s==n)?1:0}')
  if [ "$sq" != "1" ]; then
    echo "  skip np=$np (not a perfect square; HBI lattice needs a square process grid)"
    continue
  fi
  log=hbi_run_np${np}.log
  $MPIRUN $EXTRA_MPI -np "$np" "$HBI" "$RUNIN" > "$log" 2>&1
  # --- parse HBI stdout (rank-0 prints) ---
  #  'Finished all initial processing, time (s)=  T'   -> setup+assembly
  #  ' time(s)  Ttot  timer  timeH'  -> total = field 2, matvec(timeH) = field 4
  #  exit line '... at time step=  K' -> step count
  #  (the 'time for matvec(s)' line prints sum(st_ctl%time), which is UNINITIALISED
  #   in this HBI build - 1e-309/NaN - so we take timeH from the 'time(s)' line.)
  asm=$(awk -F= '/Finished all initial processing/{v=$2; gsub(/[^0-9.eE+-]/,"",v); print v; exit}' "$log")
  total=$(awk '/^[[:space:]]*time\(s\)/{print $2; exit}' "$log")
  mv_tot=$(awk '/^[[:space:]]*time\(s\)/{print $4; exit}' "$log")
  nsteps=$(awk '/at time step=/{v=$NF} END{print v}' "$log")
  # per-step matvec (proxy; HBI does 6 RKCK stage matvecs per accepted step -
  # divide by 6*(accepted+rejected) for a true per-matvec, see notes below)
  mvms=$(awk -v t="${mv_tot:-}" -v n="${nsteps:-}" 'BEGIN{if(t!=""&&n+0>0)printf "%.4f",t/n*1000; else print "NA"}')
  # memory: HBI/HACApK may print an H-matrix memory line; best-effort, else NA
  mem=$(awk -F= 'tolower($0)~/memory.*h-?matrix|h-?matrix.*memory/{v=$2; gsub(/[^0-9.]/,"",v); print v; exit}' "$log")
  printf "%-5s %4s %6s %7s %9s %11s %10s %8s %16s %8s\n" \
         hbi "$np" "$N" "$HBI_EPS_H" "${total:-NA}" "${asm:-NA}" "${mv_tot:-NA}" \
         "${nsteps:-NA}" "$mvms" "${mem:-NA}"
  echo "hbi,$np,$N,$HBI_EPS_H,${total:-NA},${asm:-NA},${mv_tot:-NA},${nsteps:-NA},$mvms,${mem:-NA}" >> "$CSV"
done

echo ""
echo "wrote $CSV  (per-run logs hbi_run_np<n>.log, run input $RUNIN)"
cat <<'EOF'

CROSS-CODE COMPARISON NOTES
 - HBI's matvec scaling vs cores reads directly off matvec_s: nsteps is
   deterministic across rank counts (same adaptive RKCK sequence), so the total
   matvec time is an honest per-core comparison without normalising.
 - For rsf_solve-vs-HBI you want PER-MATVEC time (isolates the H-engine from the
   ODE solver, which differs between codes):
     rsf_solve : exact, from -log_view  (MatMult time / MatMult count)
     HBI       : matvec_s / (6 * (accepted+rejected steps)).  HBI does 6 matvecs
                 per RKCK step; rejected steps also cost 6. matvec_ms_per_step
                 above divides by accepted steps only, so it OVERSTATES per-matvec
                 by the rejection fraction and the factor 6 - use it only to
                 compare HBI-to-HBI, not directly to rsf_solve's per-matvec.
   The directly-comparable bottom line is total wallclock for the SAME BP5 problem
   to the SAME tmax, at matched accuracy (see header note 1) - that folds in both
   the H-engine and the solver, which is what an end user actually pays.
 - Match accuracy first (header note 1): comparing HBI@eps_h=1d-4 (near-dense)
   against rsf_solve@ztol=1e-1 would be apples-to-oranges.
EOF
