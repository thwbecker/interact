#!/bin/bash
#
# rsf_solve_assembly_diag.sh
#
# Localises the large-N htool "assembly" cost in rsf_solve by putting three
# numbers side by side at each MPI rank count:
#
#   1. rsf_Is_build_wall   the wall time of the calc_petsc_Isn_matrices() call
#                          for the shear (Is) matrix, from the PetscBarrier +
#                          PetscTime bracket added in rsf_solve.c.
#   2. rsf_MatAsmEnd       the MatAssemblyEnd event time from rsf's own
#                          -log_view (should equal the bracket; if the bracket
#                          is small but this is large, the cost is deferred).
#   3. cim_assembly        the standalone compress_interaction_matrix htool
#                          assembly at the SAME geometry and the SAME htool
#                          options, built H-only via -skip_dense.
#
# If rsf_Is_build_wall tracks cim_assembly, the rsf build path is fine. If
# rsf_Is_build_wall (and rsf_MatAsmEnd) are much larger than cim_assembly at
# the same N and np, the gap is in rsf_solve's build path, not in htool.
#
# All parameters are hardcoded as plain assignments below (no environment-
# variable overrides). Edit them in place.
#
# Requires: rsf_solve and compress_interaction_matrix rebuilt from the
# instrumented sources, and the BP5 input files for the chosen resolution.

# ============================ CONFIG =======================================
wdir=assembly_diag                  # working subdir (created under scaling_tests)
mkdir -p "$wdir"
cd "$wdir" || exit 1

RB=../../bin/rsf_solve              # instrumented rsf_solve
CIM=../../bin/compress_interaction_matrix   # standalone (needs -skip_dense)
BP5=../../bp5                       # dir holding the BP5 input files

RES=0.25km                         # resolution tag (0.25km ~ N=64000)
NPLIST="8 16 24 48"                # MPI rank counts (np=1,2 are slow at 0.25km)

STOP_YR=0.01                       # short: we only need the one-off build
RTOL=1e-4                          # time-integration rtol (NOT the H tol)
HEPS=3e-5                          # htool -mat_htool_epsilon (eta + compressor
                                   # come from set_htools_defaults, exactly as
                                   # rsf_solve uses them)

MPIRUN=mpirun
EXTRA_MPI="--bind-to core --map-by core"

# thread hygiene: pin to 1 so threaded BLAS / block parallelism do not
# confound the timings (matches the other scaling scripts).
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
# If PETSc shared libs are not already on the loader path, set it here, e.g.:
# export LD_LIBRARY_PATH=/path/to/petsc/arch/lib:$LD_LIBRARY_PATH
# ===========================================================================

GEOM=$BP5/geom_bp5_${RES}.in
RSF=$BP5/rsf_bp5_${RES}.dat
IC=$BP5/ic_bp5_${RES}.in
DC=$BP5/dc_bp5_${RES}.in
for f in "$RB" "$CIM" "$GEOM" "$RSF" "$IC" "$DC"; do
  [ -e "$f" ] || { echo "ERROR: missing $f (check CONFIG)"; exit 1; }
done

COMMON="-geom_file $GEOM -rsf_file $RSF -rsf_ic_file $IC -rsf_dc_file $DC \
 -shear_modulus 3.204e10 -s_wave_speed 3464 -f0 0.6 -dc 0.14 -vpl 1e-9 \
 -v0 1e-6 -sigma_init 25e6 -rtol $RTOL -stop_time_yr $STOP_YR \
 -ts_max_steps 5 -print_interval_yr 1e9 -log_view -options_left"

CSV=rsf_assembly_diag_${RES}.csv
echo "# rsf_solve assembly diagnostic: RES=$RES heps=$HEPS  ($(date))" | tee "$CSV.hdr"
echo "np,rsf_Is_build_s,rsf_In_build_s,rsf_MatAsmEnd_s,rsf_MatMult_s,rsf_MatView_s,cim_assembly_s,cim_MatAsmEnd_s" > "$CSV"
printf "%4s %14s %14s %14s %12s %12s %13s %14s\n" \
       np rsf_Is_build rsf_In_build rsf_MatAsmEnd rsf_MatMult rsf_MatView cim_assembly cim_MatAsmEnd

# pull a -log_view event's max-time column (field 4) by exact event name
logview_time(){  # $1 logfile  $2 event name
  awk -v e="$2" '$1==e {print $4; exit}' "$1"
}

for np in $NPLIST; do
  rlog=rsf_${RES}_np${np}.log
  clog=cim_${RES}_np${np}.log

  # --- rsf_solve, htool backend, instrumented build ---
  $MPIRUN $EXTRA_MPI -np "$np" "$RB" $COMMON -use_hmatrix 1 -mat_htool_epsilon "$HEPS" \
      > "$rlog" 2>&1
  ris=$(awk '/Is \(shear\).*build wall/{print $(NF-1); exit}' "$rlog")
  rin=$(awk '/In \(normal\).*build wall/{print $(NF-1); exit}' "$rlog")
  rae=$(logview_time "$rlog" MatAssemblyEnd)
  rmm=$(logview_time "$rlog" MatMult)
  rmv=$(logview_time "$rlog" MatView)

  # --- standalone compress, SAME geometry + htool options, H only ---
  $MPIRUN $EXTRA_MPI -np "$np" "$CIM" -geom_file "$GEOM" -use_hmatrix 1 \
      -mat_htool_epsilon "$HEPS" -nrandom 0 -skip_dense -log_view \
      > "$clog" 2>&1
  cas=$(awk '/H matrix assembly took/{print $(NF-1); exit}' "$clog")
  cae=$(logview_time "$clog" MatAssemblyEnd)

  # print NA for any blank parse, without using default-substitution syntax
  for v in ris rin rae rmm rmv cas cae; do
    [ -z "$(eval echo \$$v)" ] && eval "$v=NA"
  done

  printf "%4s %14s %14s %14s %12s %12s %13s %14s\n" \
         "$np" "$ris" "$rin" "$rae" "$rmm" "$rmv" "$cas" "$cae"
  echo "$np,$ris,$rin,$rae,$rmm,$rmv,$cas,$cae" >> "$CSV"
done

echo ""
echo "logs: $wdir/rsf_${RES}_np*.log  $wdir/cim_${RES}_np*.log"
echo "csv:  $wdir/$CSV"
echo ""
echo "Read it like this:"
echo "  rsf_Is_build ~ rsf_MatAsmEnd      -> the build call IS the assembly cost"
echo "  rsf_Is_build >> cim_assembly      -> the gap is in rsf's build path, not htool"
echo "  rsf_Is_build ~ cim_assembly       -> rsf build is fine; look at MatView / MatMult"
