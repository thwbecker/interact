#!/bin/bash
#
# test_restart_ve.sh: restart fidelity of viscoelastic cycle runs.
#
# Runs the same case twice through run_bp3ve_demo: once continuously
# to stop_yr, once to split_yr and then extended to stop_yr from the
# checkpoint (the extension path run_bp3ve_demo uses when called again
# with a longer stop time).  Both the elastic reference and the
# viscoelastic case go through the restart.  check_restart.py then
# compares the two catalogs event by event.
#
# usage: test_restart_ve.sh [ds_km] [tM_yr] [stop_yr] [ncore] [branch] [tol_yr] [split_yr]
#        defaults          0.1      50      1500      4       1        0.05     1200
#
# The split point is late by default (three cycles before the end) so
# that the comparison tests the restart itself and not the growth of
# solver-tolerance differences over many cycles; a cycle this close to
# a bifurcation can amplify rtol-level differences within a few
# hundred years.  If the test fails, rerun with a later split before
# concluding anything, and compare two from-scratch runs as a control.
#
# Run from test_twod/ve_cycle_demos (next to run_bp3ve_demo).  Output
# directories rst_cont_<tag> and rst_split_<tag>; remove them to rerun.
#
ds=${1-0.1}
tm=${2-50}
stop=${3-1500}
ncore=${4-4}
branch=${5-1}
tol=${6-0.05}
half=${7-1200}
sdir=`dirname "$0"`
tag=ds${ds}_tm${tm}_b${branch}
cont=rst_cont_$tag
split=rst_split_$tag
[ -x ./run_bp3ve_demo ] || { echo "$0: run from test_twod/ve_cycle_demos"; exit 1; }
echo "=== continuous run to $stop yr -> $cont"
./run_bp3ve_demo $cont "$tm" $stop 1 $ds 1 "" $ncore $branch || exit 1
echo "=== split run: to $half yr, then extended to $stop yr -> $split"
# the split run takes the continuous run's kernel file, so that the two
# runs see byte-identical kernels and only the restart differs (two
# independent kernel generations differ at roundoff, see READING notes)
mkdir -p $split
cp $cont/bp3_tm${tm}_g1n1.prony $split/ || exit 1
./run_bp3ve_demo $split "$tm" $half 1 $ds 1 "" $ncore $branch || exit 1
./run_bp3ve_demo $split "$tm" $stop 1 $ds 1 "" $ncore $branch || exit 1
echo "=== comparison"
for t in el tm$tm;do
    grep -c "^# restarted" $split/cat_$t.dat | \
	sed "s/^/$split\/cat_$t.dat: restarted lines: /"
    grep -h "restarting from" $split/run_$t.log | tail -1
done
python3 $sdir/check_restart.py $cont $split $tm $tol
