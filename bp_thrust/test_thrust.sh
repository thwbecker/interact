#!/bin/bash
#
# Task 2 end-to-end test: dip-slip on a listric thrust, exercising and validating
# the normal-stress path (-calc_sigma_dot) via -rsf_slip_mode 1.
#
# Run from the thrust/ directory:  bash test_thrust.sh
# Assumes PETSC_DIR / PETSC_ARCH are set in the environment (as for any PETSc
# build of interact) and that bin/rsf_solve has been rebuilt after applying the
# Task 2 changes.  Generates the geometry first with make_thrust.py.
#
# All parameters are hardcoded below as plain assignments (no environment-driven
# overrides, no fallbacks).  The run uses the dense backend on one process.

# ---- parameters (edit here) ------------------------------------------------
binary=../bin/rsf_solve       # interact rsf_solve, relative to thrust/
shear_modulus=3.204e10        # G [Pa]  (Ozawa: 32.04 GPa)
s_wave_speed=3464             # c_s [m/s]
f0=0.6                        # reference friction
dc=0.02                       # D_c [m]  (coarse at ds=1 km, see make_thrust.py)
vpl=1e-9                      # loading rate [m/s]
v0=1e-6                       # reference velocity [m/s]
sigma_init=58e6              # uniform initial normal stress [Pa]
rtol=1e-4                     # ODE relative tolerance
stop_time_yr=0.0004           # short: sigma spread appears almost immediately once slip starts
print_interval_yr=0.0002      # stats cadence [yr]
slip_mode=1                   # 0 strike, 1 dip (thrust)
workdir=test_thrust_run       # output directory
# ---------------------------------------------------------------------------

# generate the geometry and inputs
python3 make_thrust.py || { echo "make_thrust.py failed"; exit 1; }

for f in $binary geom_thrust.in rsf_thrust.dat ic_thrust.in sigma_thrust.in; do
    if [ ! -e $f ]; then echo "test_thrust.sh: missing $f"; exit 1; fi
done

rm -rf $workdir
mkdir -p $workdir
cd $workdir

inputs="-geom_file ../geom_thrust.in -rsf_file ../rsf_thrust.dat -rsf_ic_file ../ic_thrust.in"
common="$inputs -shear_modulus $shear_modulus -s_wave_speed $s_wave_speed -f0 $f0 -dc $dc -vpl $vpl -v0 $v0 -sigma_init $sigma_init -rtol $rtol -stop_time_yr $stop_time_yr -print_interval_yr $print_interval_yr -rsf_slip_mode $slip_mode -use_hmatrix 0"

echo ""
echo "==== A) calc_sigma_dot ON: normal stress should evolve ===="
../$binary $common -calc_sigma_dot > on.log 2>&1
grep -E 'Is .shear|In .normal' on.log
echo "sigma spread over time (max_sigma, min_sigma, spread), sigma0=$sigma_init:"
grep -v '^#' rsf_monitor.dat | awk 'END{init="'"$sigma_init"'"}
  {printf "  t=%11.6f yr  max_sig=%.5e  min_sig=%.5e  spread=%.4e\n",$3,$8,$9,$8-$9}' | tail -4
mv rsf_monitor.dat monitor_on.dat

echo ""
echo "==== B) calc_sigma_dot OFF (control): normal stress must stay pinned ===="
../$binary $common > off.log 2>&1
grep -E 'Is .shear|In .normal' off.log
echo "distinct (max_sigma,min_sigma) across the whole OFF run:"
grep -v '^#' rsf_monitor.dat | awk '{printf "  %.6e %.6e\n",$8,$9}' | sort -u
echo "(a single line equal to $sigma_init $sigma_init means the In path is correctly inactive)"

echo ""
echo "==== C) per-cell sigma file (-rsf_sigma_file): depth-dependent sigma0 ===="
../$binary $inputs -rsf_sigma_file ../sigma_thrust.in -shear_modulus $shear_modulus -s_wave_speed $s_wave_speed -f0 $f0 -dc $dc -vpl $vpl -v0 $v0 -sigma_init $sigma_init -rtol $rtol -stop_time_yr 0.0001 -print_interval_yr 0.00005 -rsf_slip_mode $slip_mode -calc_sigma_dot -use_hmatrix 0 > sig.log 2>&1
grep -E 'read per-cell initial normal stress' sig.log
echo "t=0 sigma range (should span the file, not be uniform $sigma_init):"
grep -v '^#' rsf_monitor.dat | head -1 | awk '{printf "  max_sig=%.4e  min_sig=%.4e\n",$8,$9}'

cd ..
echo ""
echo "test_thrust.sh: done. Outputs in $workdir/"
echo "Interpretation: A shows sigma moving away from uniform once dip slip starts"
echo "(clamping and unclamping near the bend and free surface); B shows it pinned"
echo "when the flag is off; C shows the per-cell sigma0 input is honored."
