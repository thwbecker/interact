#!/bin/bash
#
# Run the SEAS BP4-QD whole-space benchmark with interact's rsf_solve, modeled
# on the BP5 create_run script.  All parameters are hardcoded here as plain
# assignments; there is no environment-variable indirection.
#
# BP4's frictional domain is infinite along strike, so the along-strike domain
# size is a numerical knob.  The benchmark asks for results at two resolutions
# and two domain sizes to show the solution is resolution and domain-size
# independent, so the config list below sweeps both.  Coarse cases can run
# dense; the fine cases need the hmmvp H-matrix backend on many cores.
#
# usage: ./run_bp4.sh

# --- fixed physical parameters (BP4 Table 1) ---
shear_modulus=3.204e10          # G [Pa]  (= rho cs^2, rho=2670, cs=3464)
s_wave_speed=3464               # cs [m/s]
sigma_init=50e6                 # effective normal stress [Pa]
dc=0.04                         # L = D_c [m] (uniform)
f0=0.6
v0=1e-6
vpl=1e-9
rupture_vth=1e-1                # BP4 rupture-front threshold V_th = 0.1 m/s

# --- run controls ---
ncore=24
tmax=1000                       # BP4 final time t_f [yr]
rtol=1e-4
hmmvp_tol=1e-4
field_step_interval=250
print_interval_yr=0.05

# --- backend: 0 dense, 4 hmmvp block-local.  dense only for the coarse cases. ---
# each entry is "ds_km Lstrike_km backend"
configs=(
    "2.0 100 0"
    "1.0 100 4"
    "0.5 100 4"
    "0.5 150 4"
    "0.25 100 4"
    "0.25 150 4"
)

export OMP_NUM_THREADS=1

for cfg in "${configs[@]}"; do
    set -- $cfg
    ds=$1; Lstrike=$2; hmat=$3
    # format ds and Lstrike the same way make_bp4.py does (python {:g}), so the
    # filenames match: 2.0 -> 2km, 0.5 -> 0.5km, 100 -> L100
    res=$(printf '%gkm' "$ds")
    lstr=$(printf '%g' "$Lstrike")
    tag="L${lstr}"

    # generate the input files for this resolution and domain size
    python3 make_bp4.py $ds $Lstrike

    wdir=bp4.$res.$tag
    mkdir -p $wdir
    cd $wdir

    lfile=log.bp4.$res.$tag.$hmat.dat
    if [ ! -s $lfile ]; then
        echo "$0: in `pwd`, writing to $lfile"
        rm -rf tmp_rsf/*
        date > $lfile
        mpirun -np $ncore ../../bin/rsf_solve \
            -geom_file   ../"geom_bp4_"$res"_"$tag".in" -rsf_file ../"rsf_bp4_"$res"_"$tag".dat" \
            -rsf_ic_file ../"ic_bp4_"$res"_"$tag".in" \
            -full_space -use_hmatrix $hmat -hmmvp_tol $hmmvp_tol \
            -shear_modulus $shear_modulus -s_wave_speed $s_wave_speed \
            -f0 $f0 -dc $dc -vpl $vpl -v0 $v0 -sigma_init $sigma_init \
            -rtol $rtol -ts_rk_type 3bs -stop_time_yr $tmax -ts_max_steps 2000000 \
            -field_step_interval $field_step_interval -print_interval_yr $print_interval_yr \
            -rsf_catalog -rupture_vth $rupture_vth -log_view 2>> $lfile >> $lfile
        date >> $lfile
    else
        echo "$0: found $lfile in $wdir, skipping"
    fi
    cd ..
done
