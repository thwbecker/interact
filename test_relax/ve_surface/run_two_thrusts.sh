#!/bin/bash
#
# run_two_thrusts.sh: surface displacements through time for two thrust
# faults in an elastic layer over a Maxwell half-space, without gravity
# and with the interface-buoyancy gravity that rsf_solve's -ve_mode 3
# kernels use, then GMT figures.
#
#   fault A: dip 30, surface breaking, bottom at 0.5 H   (Rundle 1982 geometry)
#   fault B: dip 60, surface breaking, bottom at 0.866 H (SEAS BP3 geometry,
#            34.6 km in the 40 km plate of run_bp3ve_demo)
#
# Plate and half-space as in run_bp3ve_demo / bp3_ve_kernels.py: H = 40
# km, mu = 32.04 GPa, nu = 0.25, rho 2800 / 3300, g = 9.81, half-space
# Maxwell in shear with constant bulk modulus.  Times in Maxwell times
# tM = eta/mu (NOT the Savage-Prescott 2 eta/mu).  Uniform slip of 1 m;
# displacements in units of slip.
#
# usage: run_two_thrusts.sh [times] [H_km] [xmax_H] [prefix]
#        defaults         "0.1 0.5 1 2 10 100"  40  50  tt
#
# Needs ve_thrust_relax.py next to this script and GMT 6 for the plots
# (plot_ve_thrust_gmt.sh).  About 20 s per case, four cases.
#
times=${1-"0.1 0.5 1 2 10 100"}
H=${2-40}
xmax=${3-50}
pre=${4-tt}
sdir=`dirname "$0"`
for case in A B;do
    case $case in
	A) dip=30; ext=0.5 ;;
	B) dip=60; ext=0.866 ;;
    esac
    for g in 0 1;do
	out=${pre}_${case}_g$g
	[ -f $out.npz ] && { echo "$0: $out exists, skipping"; continue; }
	python3 $sdir/ve_thrust_relax.py --H $H --mu1 32.04 --nu 0.25 --rho1 2800 --rho2 3300 --g 9.81 \
	    --dip $dip --extent $ext --taper uniform --gravity $g --times $times \
	    --xmax $xmax --nx 2001 --out $out --no-plot > $out.log 2>&1 || { echo "$0: $out failed, see $out.log"; exit 1; }
	echo "$0: $out done"
    done
done
$sdir/plot_ve_thrust_gmt.sh $pre "$times" $H
