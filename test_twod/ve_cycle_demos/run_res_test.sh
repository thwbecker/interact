#!/bin/bash
#
# run_res_test.sh: coarse-state resolution test for the BP3 VE sweep.
#
# Runs the sigma-coupled, gravity-on case (ve_g1) on both branches at
# a second cell size and the two Maxwell times given, to the same stop
# time as the 25 m production runs, into res_<ds>_g1 and
# res_<ds>_g1_normal.  bp3ve_compare.py then sets the coarse state
# (period, members, sign and size of the recurrence change, t_settle)
# against the 25 m directories.
#
# usage: run_res_test.sh [ds_km] [tmlist_yr] [stop_yr] [ncore] [ref_thrust] [ref_normal] [dip]
#        defaults        0.05     "50 250"    15000     16      ve_g1        ve_g1_normal 60
#
# Run from test_twod/ve_cycle_demos.  Cost at 50 m is roughly a
# quarter of 25 m per run (800 patches); at 12.5 m roughly eight
# times 25 m.  The kernel cache is per directory, so the two branches
# each pay one kernel generation.
#
ds=${1-0.05}
tms=${2-"50 250"}
stop=${3-15000}
ncore=${4-16}
refT=${5-ve_g1}
refN=${6-ve_g1_normal}
dip=${7-60}
sdir=`dirname "$0"`
tag=`echo $ds | sed 's/^0*//; s/\.//'`	# 0.05 -> 05, 0.0125 -> 0125
dT=res_${tag}_g1; dN=res_${tag}_g1_normal
[ -x ./run_bp3ve_demo ] || { echo "$0: run from test_twod/ve_cycle_demos"; exit 1; }
if [ "$dip" != "60" ];then
    # needs the dip argument of the patched run_bp3ve_demo (see
    # run_bp3ve_demo_dip.patch); the unpatched script ignores it and
    # would silently run dip 60
    grep -q '^dip=' run_bp3ve_demo || { echo "$0: run_bp3ve_demo has no dip argument"; exit 1; }
fi
for tm in $tms;do
    ./run_bp3ve_demo $dT "$tm" $stop 1 $ds 1 "" $ncore 1  $dip || exit 1
    ./run_bp3ve_demo $dN "$tm" $stop 1 $ds 1 "" $ncore -1 $dip || exit 1
done
echo "=== coarse-state comparison against $refT and $refN"
python3 $sdir/bp3ve_compare.py $refT $dT
python3 $sdir/bp3ve_compare.py $refN $dN
