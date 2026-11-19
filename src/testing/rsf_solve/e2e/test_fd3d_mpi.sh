#!/bin/bash
# rsf_solve -fd_mode 1 on 1, 2 and 4 MPI ranks (fd_plan.md H7): the dense couplings, their static
# gate and the strip solve are split by receiver rows, so the station series must be identical
# byte for byte for any number of ranks. Three 3-D paths, rate and state with a seed patch:
#  multi : two vertical strike-slip segments with a 30 degree bend, half space by images (G3a),
#          -calc_sigma_dot (several-plane cross tables, junction first lag)
#  msurf : a thrust whose lower half dips 50 degrees, one virtual surface (G3b)
#  strips: two surface-breaking faults dipping 60 degrees, strikes 90 and 120 degrees (H6)
# Run from the repository root after make bin/rsf_solve; about 15 minutes on one core.
# Usage: test_fd3d_mpi.sh [mpirun command] [list of rank counts]
mpi=${1-"mpirun --oversubscribe"}
ranks=${2-"1 2 4"}
d=tmp_fd3d_mpi
rm -rf $d
mkdir -p $d
cd $d || exit 1
python3 - << 'PYEOF'
import numpy as np
dx = 1000.0
def rsf_ic(name, n, seed):
    with open(f"rsf_{name}.dat", "w") as r, open(f"ic_{name}.in", "w") as ic:
        for k in range(n):
            a, b = (0.010, 0.015)
            r.write(f"{a:.6e} {b:.6e}\n")
            v = 0.03 if k in seed else 1e-9
            ic.write(f"{50e6*(a*np.arcsinh(v/(2e-6)*np.exp((0.6+b*np.log(1e-6/1e-9))/a))):.8e} {v:.6e}\n")
# multi: two vertical segments of 8 x 6 squares from the surface, the second turned by 30 degrees
L, D = 8, 6
rows = []
for seg in (0, 1):
    phi = 0.0 if seg == 0 else 30.0
    sd = np.array([np.sin(np.radians(phi)), np.cos(np.radians(phi)), 0.0])
    start = np.zeros(3) if seg == 0 else np.array([0.0, L*dx, 0.0])
    for i in range(L):
        for j in range(D):
            c = start + (i+0.5)*dx*sd + np.array([0, 0, -(j+0.5)*dx])
            rows.append(f"{c[0]:.9e} {c[1]:.9e} {c[2]:.9e} {phi:.9f} 90.0 {dx/2:.1f} {dx/2:.1f} {seg}\n")
open("geom_multi.in", "w").writelines(rows)
rsf_ic("multi", len(rows), [5*D+2, 5*D+3, 6*D+2, 6*D+3])
# msurf: thrust of 8 x 8 squares, top 1 km deep, upper half dipping 30, lower half 50 degrees
nx, nz, top = 8, 8, 1000.0
rows = []
P = np.array([0.0, 0.0, -top])
for j in range(nz):
    dd = 30.0 if j < 4 else 50.0
    v = np.array([0.0, -np.cos(np.radians(dd)), -np.sin(np.radians(dd))])
    if j == 0:
        c0 = P + 0.5*dx*v
    elif j == 4:
        vp = np.array([0.0, -np.cos(np.radians(30.0)), -np.sin(np.radians(30.0))])
        c0 = c0 + 0.5*dx*vp + 0.5*dx*v
    else:
        c0 = c0 + dx*v
    for i in range(nx):
        c = c0 + np.array([(i+0.5)*dx, 0, 0])
        rows.append(f"{c[0]:.9e} {c[1]:.9e} {c[2]:.9e} 90.0 {dd:.9f} {dx/2:.1f} {dx/2:.1f} 0\n")
open("geom_msurf.in", "w").writelines(rows)
rsf_ic("msurf", len(rows), [3*nx+3, 3*nx+4, 4*nx+3, 4*nx+4])
# strips: two faults of 4 x 3 squares dipping 60 degrees from the surface, strikes 90 and 120
# degrees, the second starting 4 km beyond the end of the first along x
rows = []
for f, (phi, x0) in enumerate(((90.0, np.array([0.0, 0.0, 0.0])), (120.0, np.array([8000.0, 0.0, 0.0])))):
    s = np.array([np.sin(np.radians(phi)), np.cos(np.radians(phi)), 0.0])
    h = np.array([np.cos(np.radians(phi)), -np.sin(np.radians(phi)), 0.0])
    dd = np.radians(60.0)
    v = np.cos(dd)*h + np.array([0, 0, -np.sin(dd)])
    for j in range(3):
        for i in range(4):
            c = x0 + (i+0.5)*dx*s + (j+0.5)*dx*v
            rows.append(f"{c[0]:.9e} {c[1]:.9e} {c[2]:.9e} {phi:.9f} 60.0 {dx/2:.1f} {dx/2:.1f} {f}\n")
open("geom_strips.in", "w").writelines(rows)
rsf_ic("strips", len(rows), [4+1, 4+2])
PYEOF
common="-sigma_init 50e6 -stop_time_yr 0.00003 -track_events 1 -vel_event 1e-3 -fd_mode 1"
opt_multi="-dc 0.02 -calc_sigma_dot 1 -full_space 0 -fd_stations 27,75"
opt_msurf="-dc 0.04 -rsf_slip_mode 1 -fd_surface_grade 0 -fd_surface_ext 2 -fd_surface_extend 0 -calc_sigma_dot 1 -fd_stations 27,43"
opt_strips="-dc 0.04 -rsf_slip_mode 1 -fd_surface_m 3 -fd_surface_band 2000 -fd_stations 5,17"
fail=0
for c in multi msurf strips; do
  eval o=\$opt_$c
  for n in $ranks; do
    mkdir -p $c/n$n
    ( cd $c/n$n && $mpi -np $n ../../../bin/rsf_solve -geom_file ../../geom_$c.in -rsf_file ../../rsf_$c.dat \
        -rsf_ic_file ../../ic_$c.in $common $o > run.log 2>&1 )
  done
  r="test_fd3d_mpi $c:"
  first=`echo $ranks | awk '{print $1}'`
  nst=`ls $c/n$first/fd_phase_*_st*.dat 2> /dev/null | wc -l`
  ok=1
  [ $nst -gt 0 ] || ok=0
  for n in $ranks; do
    for f in $c/n$first/fd_phase_*_st*.dat; do
      cmp -s $f $c/n$n/`basename $f` || ok=0
    done
  done
  ph=`grep -c "3-D phase" $c/n$first/run.log`
  if [ $ok -eq 1 ]; then r="$r PASS"; else r="$r FAIL"; fail=1; fi
  echo "$r (ranks $ranks; $nst station files, $ph dynamic phases on $first rank(s))"
done
exit $fail
