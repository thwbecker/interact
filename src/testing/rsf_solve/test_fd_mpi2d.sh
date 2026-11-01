#!/bin/bash
# rsf_solve -fd_mode 1 for 2-D in-plane faults on several MPI ranks (fd_plan.md M8): a curved
# fault (polyline, dense tables; event at 11 yr) and a 30 degree thrust in a half plane (virtual
# free surface; event at 24 yr), each with one dynamic event; the station series of 1, 2 and 4 ranks must be identical byte for
# byte. Run from the repository root after make bin/rsf_solve; a few minutes.
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
common="-dc 0.016 -sigma_init 50e6 -v0 1e-6 -f0 0.6 -vpl 1e-9 -shear_modulus 3.24e10 -s_wave_speed 3464 -rtol 1e-6 -track_events 1 -vel_event 1e-3 -fd_mode 1 -fd_stations 10,20"
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
  nst=`ls $name/n1/fd_phase_*_st*.dat 2> /dev/null | wc -l`
  [ $nst -gt 0 ] || ok=0
  for n in 2 4; do
    for f in $name/n1/fd_phase_*_st*.dat; do
      cmp -s $f $name/n$n/`basename $f` || ok=0
    done
  done
  if [ $ok -eq 1 ]; then echo "test_fd_mpi2d $name: PASS (1, 2 and 4 ranks identical, $nst station files)"
  else echo "test_fd_mpi2d $name: FAIL"; fail=1; fi
done
exit $fail
