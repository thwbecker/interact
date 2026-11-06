#!/bin/bash
# rsf_solve -fd_mode 1 for 2-D in-plane faults on several MPI ranks (fd_plan.md M8, M10): a curved
# fault (polyline, dense tables; event at 11 yr) and a 30 degree thrust in a half plane (exact
# kernels; event at 24 yr), each with one dynamic event. The tables of 1, 2 and 4 ranks must be
# identical (-fd_table_hash); the station series are identical with reference BLAS and must agree
# to 1e-3 otherwise (see below). Run from the repository root after make bin/rsf_solve; a few
# minutes for the curved fault, longer for the half plane (its reflected tables are costly).
# usage: test_fd_mpi2d.sh [MPI launch command] (default: the first of mpirun, mpirun --oversubscribe,
#   mpirun --allow-run-as-root --oversubscribe that starts 4 ranks)
mpi=${1-}
if [ -z "$mpi" ]; then
  for m in "mpirun" "mpirun --oversubscribe" "mpirun --allow-run-as-root --oversubscribe"; do
    if $m -np 4 true > /dev/null 2>&1; then
      mpi=$m
      break
    fi
  done
fi
d=tmp_fd_mpi2d
rm -rf $d
mkdir -p $d
cd $d || exit 1
python3 - << 'PYEOF'
import numpy as np
def write(name, c, strike, dx, a, b):
    n = len(c)
    V0, f0, Vpl, sig = 1e-6, 0.6, 1e-9, 50e6
    psi = f0 + b*np.log(V0/Vpl)
    Vi = np.full(n, Vpl)
    Vi[n//2 - 2:n//2 + 3] *= 1.5
    tau0 = a*sig*np.arcsinh(Vi/(2*V0)*np.exp(psi/a))
    with open(f"geom_{name}.in", "w") as g, open(f"rsf_{name}.dat", "w") as r, open(f"ic_{name}.in", "w") as ic:
        for i in range(n):
            g.write(f"{c[i, 0]:.9e} {c[i, 1]:.9e} 0 {strike[i]:.9f} 90 {dx/2:.6e} 0 0\n")
            r.write(f"{a[i]:.6f} {b[i]:.6f}\n")
            ic.write(f"{tau0[i]:.8e} {Vi[i]:.6e}\n")
# curved fault: arc of 30 km radius, 60 elements of 200 m, full plane
R, dx, n = 30e3, 200.0, 60
th = np.radians(-11.5) + (np.arange(n) + 0.5)*dx/R
c = np.stack([R*np.sin(th), -R*(1 - np.cos(th))], axis=1)
d = np.gradient(c, axis=0)
strike = 90.0 - np.degrees(np.arctan2(d[:, 1], d[:, 0]))
k = np.abs(np.arange(n) - n/2) < 0.35*n
write("arc", c, strike, dx, np.where(k, 0.010, 0.020), np.full(n, 0.015))
# thrust dipping 30 degrees, top 1 km deep, 40 elements of 200 m, half plane
n = 40
s = np.array([np.cos(np.radians(30.0)), -np.sin(np.radians(30.0))])
c = np.array([0.0, -1000.0])[None, :] + ((np.arange(n) + 0.5)*dx)[:, None]*s[None, :]
k = np.abs(np.arange(n) - n/2) < 0.35*n
write("hp", c, np.full(n, 90.0 - np.degrees(np.arctan2(s[1], s[0]))), dx, np.where(k, 0.010, 0.020), np.full(n, 0.015))
PYEOF
common="-dc 0.016 -sigma_init 50e6 -v0 1e-6 -f0 0.6 -vpl 1e-9 -shear_modulus 3.24e10 -s_wave_speed 3464 -rtol 1e-6 -track_events 1 -vel_event 1e-3 -fd_mode 1 -fd_stations 10,20 -fd_table_hash"
fail=0
for c in "arc:-full_space 1 -stop_time_yr 13" "hp:-stop_time_yr 26"; do
  name=${c%%:*}
  opt=${c#*:}
  for n in 1 2 4; do
    mkdir -p $name/n$n
    ( cd $name/n$n && OMP_NUM_THREADS=1 $mpi -np $n ../../../bin/rsf_solve -geom_file ../../geom_$name.in -rsf_file ../../rsf_$name.dat \
        -rsf_ic_file ../../ic_$name.in $common $opt > run.log 2>&1 )
  done
  ok=1
  # the tables must be identical (an order-independent hash of their bits, -fd_table_hash)
  h1=`grep 'table hash' $name/n1/run.log | awk '{print $NF}'`
  for n in 2 4; do
    hn=`grep 'table hash' $name/n$n/run.log | awk '{print $NF}'`
    if [ -z "$h1" ] || [ "$h1" != "$hn" ]; then
      ok=0
      echo "test_fd_mpi2d $name: table hash on $n ranks '$hn' differs from 1 rank '$h1'"
    fi
  done
  nst=`ls $name/n1/fd_phase_*_st*.dat 2> /dev/null | wc -l`
  if [ $nst -eq 0 ]; then
    ok=0
    echo "test_fd_mpi2d $name: no station files on 1 rank; end of $name/n1/run.log:"
    tail -5 $name/n1/run.log
  fi
  # the station series: identical with reference BLAS; an optimised BLAS sums the QD part's MatMult
  # in an order that depends on the rows per rank, and the runs then differ in the last bits,
  # amplified over the QD years: accepted up to 1e-3 relative (4e-5 seen on a 64-core machine)
  same=1
  for n in 2 4; do
    for f in $name/n1/fd_phase_*_st*.dat; do
      [ -s $f ] || continue
      g=$name/n$n/`basename $f`
      if [ ! -s $g ]; then
        ok=0
        echo "test_fd_mpi2d $name: $g missing; end of $name/n$n/run.log:"
        tail -5 $name/n$n/run.log
      elif ! cmp -s $f $g; then
        same=0
        r=`paste $f $g | awk '{h = NF/2; for(i=1;i<=h;i++){d=$i-$(i+h); if(d<0)d=-d; a=$i; if(a<0)a=-a; if(d>dm[i])dm[i]=d; if(a>am[i])am[i]=a}}
          END {w=0; for(i in dm){r=(am[i]>0)?dm[i]/am[i]:dm[i]; if(r>w)w=r} printf "%.2e", w}'`
        echo "test_fd_mpi2d $name: `basename $f` on $n ranks differs from 1 rank by $r (largest relative)"
        awk -v r=$r 'BEGIN{exit !(r > 1e-3)}' && ok=0
      fi
    done
  done
  if [ $ok -eq 1 ]; then
    if [ $same -eq 1 ]; then echo "test_fd_mpi2d $name: PASS (tables and runs identical on 1, 2 and 4 ranks, $nst station files)"
    else echo "test_fd_mpi2d $name: PASS (tables identical on 1, 2 and 4 ranks; runs within 1e-3)"; fi
  else echo "test_fd_mpi2d $name: FAIL"; fail=1; fi
done
exit $fail
