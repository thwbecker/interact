#!/bin/bash
#
# run_thrust_relax.sh: surface displacements through time for one thrust
# fault in an elastic layer over a Maxwell half-space, without gravity
# and with the interface-buoyancy gravity that rsf_solve's -ve_mode 3
# kernels use, then GMT figures (plot_ve_thrust_gmt.sh).
#
# Plate and half-space as in run_bp3ve_demo / bp3_ve_kernels.py: H = 40
# km, mu = 32.04 GPa, nu = 0.25, rho 2800 / 3300, g = 9.81, half-space
# Maxwell in shear with constant bulk modulus.  Times in Maxwell times
# tM = eta/mu (NOT the Savage-Prescott 2 eta/mu).  Uniform slip of 1 m;
# displacements in units of slip; x in units of H, trace at 0, fault
# dips towards +x, hanging wall at x > 0.
#
# usage: run_thrust_relax.sh [dip_deg] [extent] [times] [H_km] [xmax_H] [prefix]
#        defaults            60        0.866    "0.1 0.5 1 2 10 100"  40  50  thrust
#
# examples:
#   SEAS BP3 geometry (default): dip 60, surface breaking to 34.6 km in
#     the 40 km plate, i.e. extent 0.866; the fault of run_bp3ve_demo
#       run_thrust_relax.sh 60 0.866
#   Rundle (1982) Figs. 2/3 geometry: dip 30, down-dip width W = H,
#     bottom at 0.5 H (his tau_a = 2 tM, so his 5 and 45 tau_a are
#     times 10 and 90)
#       run_thrust_relax.sh 30 0.5 "0.1 0.5 1 2 10 90"
#
# Needs ve_thrust_relax.py next to this script and GMT 6 for the plots.
# About 20 s per gravity setting.  Run from a scratch directory; output
# <prefix>_g{0,1}.{npz,txt,log}, <prefix>_g{0,1}_profiles.txt and the
# PDFs <prefix>_near.pdf, <prefix>_far.pdf.
#
dip=${1-60}
ext=${2-0.866}
times=${3-"0.1 0.5 1 2 10 100"}
H=${4-40}
xmax=${5-50}
pre=${6-thrust}
sdir=`dirname "$0"`
for g in 0 1;do
    out=${pre}_g$g
    [ -f $out.npz ] && { echo "$0: $out exists, skipping"; continue; }
    python3 $sdir/ve_thrust_relax.py --H $H --mu1 32.04 --nu 0.25 --rho1 2800 --rho2 3300 --g 9.81 \
	--dip $dip --extent $ext --taper uniform --gravity $g --times $times \
	--xmax $xmax --nx 2001 --out $out --no-plot > $out.log 2>&1 || { echo "$0: $out failed, see $out.log"; exit 1; }
    echo "$0: $out done"
done
$sdir/plot_ve_thrust_gmt.sh $pre "$times" $H "dip $dip, bottom $ext H"
