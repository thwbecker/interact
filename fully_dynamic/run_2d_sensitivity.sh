#!/bin/bash
#
# Sensitivity of a 2-D QD + FD run to the numerical parameters of the fully dynamic phases
# (fd_plan.md M9): one resampled section (outdir/geom of run_2d_cascadia.sh, or made here with
# coarsen_2d.py), run up to years with each variant, compared with the base run:
#   base    defaults (-fd_tol 1e-2, -fd_tol_surf 1e-3, -fd_vsw 1e-3, -fd_vback 1e-4)
#   tsurf4  -fd_tol_surf 1e-4      (windows of the pairs with surface elements)
#   tol3    -fd_tol 1e-3           (windows of the fault pairs)
#   vsw4    -fd_vsw 1e-4 -fd_vback 1e-5 (switch earlier; vback must stay below vsw)
#   vback5  -fd_vback 1e-5         (hand back to QD later)
#   hnear4  -fd_surf_h_near h/4    (surface elements a quarter of the fault elements near tops)
# For each: event onsets, the phases (steps, peak slip rate, mean slip, wall time), the mean slip
# at the end (from rsf_monitor.dat), and for the variants with the same phase boundaries the
# largest and rms difference of the phase slip per element against the base (fd_phase_1_slip.dat).
# Results in outdir/sensitivity.txt. Run from the repository root after make bin/rsf_solve.
#
# usage: run_2d_sensitivity.sh [h] [np] [mpi command] [years] [outdir] [calc_sigma_dot 0|1] [variants] [section dir]
#   defaults: 500 16 mpirun 30 2d_sensitivity_out 1 "base tsurf4 tol3 vsw4 vback5 hnear4"
#             test_twod/cascadia/PD13_nosigma
#
h=${1-500}
np=${2-16}
mpi=${3-mpirun}
years=${4-30}
out=${5-2d_sensitivity_out}
sd=${6-1}
variants=${7-"base tsurf4 tol3 vsw4 vback5 hnear4"}
sec=${8-test_twod/cascadia/PD13_nosigma}
root=`pwd`
bin=$root/bin/rsf_solve
mkdir -p $out
out=`cd $out && pwd`
sum=$out/sensitivity.txt
echo "run_2d_sensitivity `date`, `hostname`: h $h m, $np ranks, $years yr, calc_sigma_dot $sd, variants '$variants'" >> $sum
if [ ! -s $out/geom/rsf_ab.in ]; then
  python3 $root/fully_dynamic/coarsen_2d.py $root/$sec $out/geom $h | tee -a $sum
fi
python3 $root/fully_dynamic/check_2d_resolution.py $out/geom >> $sum
geom=`cd $out/geom && ls geom_*.in | head -1`
hmin=`awk 'NR == 1 || 2*$6 < m {m = 2*$6} END {print m}' $out/geom/$geom`
base="-geom_file ../geom/$geom -rsf_file ../geom/rsf_ab.in -rsf_dc_file ../geom/rsf_dc.in -rsf_sigma_file ../geom/rsf_sigma.in \
 -rsf_vpl_file ../geom/rsf_vpl.in -rsf_ic_file ../geom/rsf_ic.in -state_law 1 -shear_modulus 3.0e10 -s_wave_speed 3464 \
 -f0 0.6 -v0 1e-6 -vpl 1.2e-9 -vel_event 1e-3 -track_events 1 -ts_rk_type 3bs -rtol 1e-5 -calc_sigma_dot $sd \
 -stop_time_yr $years -dt_monitor_yr 0.01 -fd_mode 1 -fd_surf_margin 10000 -fd_stations 10,20"
for v in $variants; do
  case $v in
    base) opt="" ;;
    tsurf4) opt="-fd_tol_surf 1e-4" ;;
    tol3) opt="-fd_tol 1e-3" ;;
    vsw4) opt="-fd_vsw 1e-4 -fd_vback 1e-5" ;;
    vback5) opt="-fd_vback 1e-5" ;;
    hnear4) opt="-fd_surf_h_near `awk -v h=$hmin 'BEGIN {print h/4.0}'`" ;;
    *) echo "unknown variant $v" >> $sum
       continue ;;
  esac
  rm -rf $out/$v
  mkdir -p $out/$v
  cd $out/$v
  s0=`date +%s`
  OMP_NUM_THREADS=1 $mpi -np $np $bin $base $opt > run.log 2>&1
  echo "$v ($opt): rc $?, wall $(( `date +%s` - s0 )) s" >> $sum
  grep 'rsf_fd: phase' run.log | cut -c1-170 | sed "s/^/  $v: /" >> $sum
  grep 'PETSC ERROR\|Segmentation\|Abort' run.log | head -3 | sed "s/^/  $v: /" >> $sum
  cd $root
done
python3 - $out "$variants" $years << 'PYEOF' >> $sum
import sys, os
import numpy as np
out, variants, years = sys.argv[1], sys.argv[2].split(), float(sys.argv[3])
def load(v):
    r = {}
    d = os.path.join(out, v)
    f = os.path.join(d, "rsf_events.dat")
    if os.path.exists(f):
        a = np.loadtxt(f, comments="#", ndmin=2)
        r["onsets"] = a[a[:, 2] == 1, 1] if a.size else np.array([])
    f = os.path.join(d, "rsf_monitor.dat")
    if os.path.exists(f):
        m = np.loadtxt(f, comments="#", ndmin=2)
        if len(m):
            r["slip_end"] = np.interp(min(years, m[-1, 2]), m[:, 2], m[:, 5])
            r["t_end"] = m[-1, 2]
    f = os.path.join(d, "fd_phase_1_slip.dat")
    if os.path.exists(f):
        s = np.loadtxt(f, ndmin=2)
        r["phase_slip"] = s[s[:, 1] >= 0, 2]          # fault elements (surface rows have index -1)
    return r
R = {v: load(v) for v in variants}
b = R["base"] if "base" in R else load("base")      # an earlier base run in outdir serves as well
for v in variants:
    r = R[v]
    line = f"c: {v}: onsets (yr) " + " ".join(f"{t:.4f}" for t in r.get("onsets", []))
    if "slip_end" in r:
        line += f"; mean slip at {min(years, r['t_end']):.2f} yr {r['slip_end']:.5f} m"
        if v != "base" and "slip_end" in b:
            line += f" ({(r['slip_end'] - b['slip_end'])/max(abs(b['slip_end']), 1e-30)*100:+.3f} % of base)"
    if v in ("tsurf4", "tol3", "hnear4") and "phase_slip" in r and "phase_slip" in b and len(r["phase_slip"]) == len(b["phase_slip"]):
        dd = r["phase_slip"] - b["phase_slip"]
        sc = np.abs(b["phase_slip"]).max()
        line += f"; phase 1 slip against base: largest {np.abs(dd).max()/sc*100:.3f} %, rms {np.sqrt(np.mean(dd**2))/sc*100:.3f} % of the largest slip"
    print(line)
PYEOF
echo "done `date`" >> $sum
