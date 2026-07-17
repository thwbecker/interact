#!/bin/bash
#
# Per-event onset-time misfit vs the dense reference, as a function of event
# number, for one backend at one resolution. This is the diagnostic that
# separates two different things that a single overlaid "10th event" plot
# conflates:
#
#   1. is the H-matrix operator accurate (small per-event misfit), and
#   2. is the BP5 sequence a stable limit cycle (misfit bounded or decaying
#      with event number) or chaotic (misfit growing with event number).
#
# If the misfit stays flat or decays as the event number increases, the
# perturbation from the H-matrix approximation is not amplifying: all eps
# runs phase-lock onto the dense cycle, which is exactly why late events
# (the 10th) overlap. If it grows, the system is chaotic and tighter eps
# only delays divergence; in that case the eps curves separate, ordered by
# tolerance, at late events. It also flags event-count mismatches, since a
# missing or extra small event shifts the index and makes a fixed "10th
# event" compare physically different events (an alignment artifact, not
# dynamics).
#
# Run from bp5/. Output columns: eps event_no dense_yr run_yr misfit_yr

res=1km                 # 1km | 0.5km | 0.25km
hmat=3                  # 1 htool, 3 hacapk, 4 hmmvp
dense=ntest.0/rsf_events.${res}.0.inf.dat
out=event_misfit.${res}.${hmat}.dat

nd=$(awk '$3==1{c++} END{print c+0}' $dense)
echo "# eps event_no dense_yr run_yr misfit_yr   ($res hmat=$hmat, dense has $nd onsets)" > $out

for f in ntest.${hmat}/rsf_events.${res}.${hmat}.*.dat; do
    eps=$(echo $f | awk -F. '{print $(NF-1)}')
    ne=$(awk '$3==1{c++} END{print c+0}' $f)
    if [ "$ne" != "$nd" ]; then
        echo "# NOTE eps=$eps has $ne onsets vs dense $nd: index alignment unreliable, treat as count mismatch" >> $out
    fi
    # align onset times by event index (paste), emit per-event absolute misfit
    paste <(awk '$3==1{print $2}' $dense) <(awk '$3==1{print $2}' $f) \
        | awk -v eps=$eps 'NF==2{ n++; d=$2-$1; if(d<0)d=-d; printf "%s %d %.4f %.4f %.4f\n", eps, n, $1, $2, d }' >> $out
done

echo "wrote $out"
echo "plot misfit_yr (col 5) vs event_no (col 2), one line per eps:"
echo "  growing with event_no  -> chaotic (eps ordering visible at late events)"
echo "  flat or decaying        -> stable limit cycle (late events overlap; method is faithful)"
