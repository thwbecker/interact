#!/bin/bash
#
# Set up and run the Ozawa et al. (2023) nonplanar thrust with rsf_solve, as a
# mesh-convergence sweep. For each resolution in the list below it generates the
# geometry with make_thrust.py (ozawa case) and runs a dip-slip, normal-stress
# evolving earthquake-cycle simulation.
#
# Run from the bp_thrust/ directory:  bash run_ozawa.sh
#
# All model parameters are hardcoded as plain assignments below. There is no
# environment-variable indirection. OMP_NUM_THREADS is fixed to 1 because the
# parallelism here is MPI over patches, not threads. Adjust ncore, the
# resolution list, the backend, and stop_time_yr for your machine and for how
# many cycles you want.

export OMP_NUM_THREADS=1

# ---- run configuration (edit here) -----------------------------------------
binary=../../bin/rsf_solve       # rsf_solve, relative to each per-resolution run dir
generator=../make_thrust.py      # geometry generator, relative to each run dir
ncore=24                         # MPI ranks
resolutions="2.0 1.0 0.5"        # cell sizes [km] for the convergence sweep; extend downward

# backend: 0 dense, 3 HACApK, 4 hmmvp. Dense is fine at 2 km; use an H-matrix
# backend for the finer meshes. This mirrors the BP5 create_run default.
hmat=4
hmat_flag="-use_hmatrix 4 -hmmvp_tol 1e-4"

# elastic and rate-and-state parameters (Ozawa et al. 2023)
shear_modulus=3.204e10           # mu [Pa]  (32.04 GPa)
s_wave_speed=3464                # c_s [m/s]
f0=0.6                           # reference friction
dc=0.02                          # D_c [m]
vpl=1e-9                         # plate rate [m/s]
v0=1e-6                          # reference velocity [m/s]
sigma_init=58e6                  # uniform initial normal stress [Pa]
                                 # (initial shear 100 MPa comes from ic_thrust.in)

# normal-stress limiter (matches HBI limitsigma). On a strongly slipping bending
# thrust the normal stress can be driven toward zero or negative, which is
# unphysical and numerically fragile, so clamp it. Adjust the floor to taste.
min_sigma=1e6                    # floor [Pa]
max_sigma=300e6                  # ceiling [Pa]

# time stepping and output
rtol=1e-4
ts_rk_type=3bs
stop_time_yr=3000                # long enough for several cycles past the initial transient
print_interval_yr=0.02           # stats / catalog cadence [yr]
field_step_interval=250          # slip-rate field frame cadence [accepted steps]
rupture_vth=0.1                  # rupture-front threshold [m/s] for the catalog / rupture time
# ---------------------------------------------------------------------------

for res in $resolutions; do
    wdir=ozawa_$res\km
    echo "==== resolution $res km -> $wdir ===="
    rm -rf $wdir
    mkdir -p $wdir
    cd $wdir

    # generate the ozawa-case inputs at this resolution
    python3 $generator $res ozawa
    if [ ! -s geom_thrust.in ]; then
        echo "run_ozawa.sh: generation failed at $res km"; cd ..; continue
    fi

    lfile=log.ozawa.$res.dat
    date > $lfile
    mpirun -np $ncore $binary \
        -geom_file geom_thrust.in -rsf_file rsf_thrust.dat -rsf_ic_file ic_thrust.in \
        -shear_modulus $shear_modulus -s_wave_speed $s_wave_speed \
        -f0 $f0 -dc $dc -vpl $vpl -v0 $v0 -sigma_init $sigma_init \
        -rsf_slip_mode 1 -calc_sigma_dot \
        -limit_sigma -min_sigma $min_sigma -max_sigma $max_sigma \
        $hmat_flag -rtol $rtol -ts_rk_type $ts_rk_type \
        -stop_time_yr $stop_time_yr -print_interval_yr $print_interval_yr \
        -field_step_interval $field_step_interval \
        -rsf_catalog -rsf_rupture_time -rupture_vth $rupture_vth \
        -ts_max_steps 20000000 #>> $lfile 2>> $lfile
    exit
    date >> $lfile

    echo "  done. event catalog:"
    grep -v '^#' rsf_catalog.dat 2>/dev/null | \
        awk '{printf "    ev%s onset=%.4f yr dur=%.1f s n_rup=%s Mw=%.3f\n",$1,$2,$4,$5,$13}'
    cd ..
done

echo ""
echo "run_ozawa.sh: convergence sweep done."
echo "Compare the stabilized (post-transient) event sequence across the"
echo "ozawa_*km directories: recurrence, ruptured area, and Mw should converge"
echo "as the cell size decreases. The first event in each run is an"
echo "initialization transient from the uniform 100 MPa overstress and should be"
echo "discarded. Grid the rupture-time field (rsf_rupture_time.dat) and the"
echo "slip-rate frames (tmp_rsf/) to inspect where events nucleate and how the"
echo "normal stress changes near the bend and the free surface."
