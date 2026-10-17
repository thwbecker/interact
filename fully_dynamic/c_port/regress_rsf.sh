#!/bin/sh
# regress_rsf.sh: bit-identity and timing harness for rsf_solve
#
# Builds a reference rsf_solve from a named commit (git worktree) and the
# working tree, runs a fixed case set with both, compares md5 sums of every
# output file and reports timing ratios. Run from the repository root.
#
# usage: fully_dynamic/c_port/regress_rsf.sh [ref_commit] [nrep]
#   ref_commit  commit of the reference binary (default: HEAD)
#   nrep        timing repetitions per case (default 3)
# environment: REFBIN=<path> uses a prebuilt reference binary instead of
#   building ref_commit; EXTRA_CASES="name|gen|opts;..." adds cases
#
# The case set: 2-D antiplane rate-and-state cycles (gen_bp2d.py from
# test_twod/anti_plane), elastic and with -ve_mode 1 and 2, short enough to
# finish in minutes. 3-D cases can be added through EXTRA_CASES once an
# example with a quick runtime is chosen (interact_examples).
ref=${1-HEAD}
nrep=${2-3}
root=`pwd`
wt=/tmp/rsf_regress_ref_`id -u`
log=$root/regress_rsf.log
fail=0
date > $log
echo "reference commit: $ref" >> $log

# 1. reference binary: REFBIN=<path> uses a prebuilt one (e.g. a copy of
#    bin/rsf_solve made before applying a patch); otherwise the commit is
#    built in a git worktree. Built, untracked external libraries of the
#    main tree (hmat/: HACApK, hmmvp) are linked into the worktree, since
#    a fresh checkout does not contain them.
if [ -n "$REFBIN" ]; then
  refbin=$REFBIN
  [ -x "$refbin" ] || { echo "REFBIN=$refbin is not executable"; exit 1; }
  echo "reference binary: $refbin (prebuilt)" >> $log
else
  if [ ! -f $wt/.commit ] || [ "`cat $wt/.commit`" != "`git rev-parse $ref`" ]; then
    rm -rf $wt; git worktree prune
    git worktree add -f $wt $ref > /dev/null 2>&1 || { echo "cannot create worktree for $ref"; exit 1; }
    for d in hmat; do
      if [ -d $root/$d ]; then rm -rf $wt/$d; ln -s $root/$d $wt/$d; fi
    done
    ( cd $wt && make bin/rsf_solve > make_ref.log 2>&1 ) || { echo "reference build failed, see $wt/make_ref.log"; exit 1; }
    git rev-parse $ref > $wt/.commit
  fi
  refbin=$wt/bin/rsf_solve
fi
# 2. current binary
make bin/rsf_solve > make_cur.log 2>&1 || { echo "current build failed, see make_cur.log"; exit 1; }
curbin=$root/bin/rsf_solve

# 3. cases: name | generator command | rsf_solve options
common="-dc 0.14 -sigma_init 25e6 -rtol 1e-4 -ts_rk_type 3bs -dt_monitor_yr 2 -track_events -stop_time_yr 3000"
gen="python3 $root/test_twod/anti_plane/gen_bp2d.py 0.5 1 20 20 6"
cases="ap_elastic|$gen|-geom_file geom_bp2d.in -rsf_file rsf_bp2d.dat -rsf_ic_file ic_bp2d.in $common
ap_ve2|$gen|-geom_file geom_bp2d.in -rsf_file rsf_bp2d.dat -rsf_ic_file ic_bp2d.in $common -ve_mode 2 -ve_tmaxwell_yr 100 -ve_plate_h 40e3
ap_ve1|$gen|-geom_file geom_bp2d.in -rsf_file rsf_bp2d.dat -rsf_ic_file ic_bp2d.in $common -ve_mode 1 -ve_tmaxwell_yr 100"
if [ -n "$EXTRA_CASES" ]; then
  cases="$cases
`echo "$EXTRA_CASES" | tr ';' '\n'`"
fi

run_case() {                      # $1 binary  $2 dir  $3 generator  $4 options
  mkdir -p $2 && cd $2 && rm -f *.dat *.log
  eval $3 > gen.log 2>&1
  start=`date +%s.%N`
  $1 $4 > run.log 2>&1
  end=`date +%s.%N`
  cd $root
  echo "$end - $start" | bc
}

echo "$cases" | while IFS='|' read name g opts; do
  [ -z "$name" ] && continue
  tref=""; tcur=""
  i=1
  while [ $i -le $nrep ]; do
    t1=`run_case $refbin tmp_regress/${name}_ref "$g" "$opts"`
    t2=`run_case $curbin tmp_regress/${name}_cur "$g" "$opts"`
    tref="$tref $t1"; tcur="$tcur $t2"
    i=`expr $i + 1`
  done
  mref=`echo $tref | tr ' ' '\n' | sort -g | head -1`
  mcur=`echo $tcur | tr ' ' '\n' | sort -g | head -1`
  ratio=`echo "scale=4; $mcur / $mref" | bc`
  # identity: every output file of the reference must exist and match
  ident=PASS
  for f in tmp_regress/${name}_ref/*.dat; do
    b=`basename $f`
    if [ ! -f tmp_regress/${name}_cur/$b ]; then ident="FAIL (missing $b)"; break; fi
    s1=`md5sum < $f`; s2=`md5sum < tmp_regress/${name}_cur/$b`
    if [ "$s1" != "$s2" ]; then ident="FAIL ($b differs)"; break; fi
  done
  speed=PASS
  if [ `echo "$ratio > 1.02 || $ratio < 0.98" | bc` -eq 1 ]; then speed="WARN"; fi
  printf "%-12s identity %-24s time ref %7.2fs cur %7.2fs ratio %6.4f %s\n" "$name" "$ident" "$mref" "$mcur" "$ratio" "$speed" | tee -a $log
  case "$ident" in FAIL*) echo FAIL > tmp_regress/.fail;; esac
done
if [ -f tmp_regress/.fail ]; then echo "regress_rsf: FAIL (see $log)"; rm -f tmp_regress/.fail; exit 1; fi
echo "regress_rsf: identity PASS for all cases (timing in $log)"
