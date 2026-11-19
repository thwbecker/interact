#!/bin/sh
# rsf_solve -fd_mode 1 on two parallel vertical strike-slip segments (a 3-D stepover, 16 x 8
# squares of 1 km each, 2 km apart, 4 km overlap), rate and state, -calc_sigma_dot:
# whole space: the static gates of the tables against interact's Is (within and across the
# planes) and of the cross normal traction against -In must pass at 1e-8; half space (mirror):
# a dynamic phase on both planes runs and gives the same station series on 1 and 3 MPI ranks
# (if mpirun is available). About 30 s. Run from the repository root after make bin/rsf_solve.
d=tmp_fd3d_pair
rm -rf $d; mkdir -p $d/fs $d/r1 $d/r3
python3 - << 'PYEOF'
import numpy as np
dx=1000.0; L=16; D=8; ov=4; W=2000.0
with open("tmp_fd3d_pair/geom.in","w") as g, open("tmp_fd3d_pair/rsf.dat","w") as r, open("tmp_fd3d_pair/ic.in","w") as ic:
    for seg in (0,1):
        x = 0.0 if seg==0 else W
        y0 = 0.0 if seg==0 else (L-ov)*dx
        for i in range(L):
            for j in range(D):
                g.write(f"{x:.6e} {y0+(i+0.5)*dx:.6e} {-(j+0.5)*dx:.6e} 0.0 90.0 {dx/2:.1f} {dx/2:.1f} {seg}\n")
                edge=(i in (0,L-1)) or (j==D-1)
                a,b=(0.025,0.015) if edge else (0.010,0.015)
                r.write(f"{a:.6e} {b:.6e}\n")
                v=0.03 if ((seg==0) and (i in (7,8)) and (j in (3,4))) else 1e-9
                ic.write(f"{50e6*(a*np.arcsinh(v/(2e-6)*np.exp((0.6+b*np.log(1e-6/1e-9))/a))):.8e} {v:.6e}\n")
PYEOF
O="-geom_file ../geom.in -rsf_file ../rsf.dat -rsf_ic_file ../ic.in -sigma_init 50e6 -dc 0.02 -track_events 1 -vel_event 1e-3 -fd_mode 1 -fd_back_eps 0.05 -calc_sigma_dot 1"
( cd $d/fs && ../../bin/rsf_solve $O -full_space 1 -stop_time_yr 0.00001 > run.log 2>&1 )
g=`grep "static gate, two planes" $d/fs/run.log | sed 's/.*Is \([0-9.e+-]*\) within planes, \([0-9.e+-]*\) across, normal traction across from -In \([0-9.e+-]*\);.*/\1 \2 \3/'`
ok=`echo "$g" | awk '{ print ($1 == $1+0 && $1 < 1e-8 && $2 < 1e-8 && $3 < 1e-8) ? 1 : 0 }'`
( cd $d/r1 && ../../bin/rsf_solve $O -stop_time_yr 0.0001 -fd_stations 60,180 > run.log 2>&1 )
grep -q "3-D phase 1 (two planes)" $d/r1/run.log || ok=0
mpi="skipped"
if command -v mpirun > /dev/null 2>&1; then
  ( cd $d/r3 && ( mpirun --oversubscribe -np 3 ../../bin/rsf_solve $O -stop_time_yr 0.0001 -fd_stations 60,180 > run.log 2>&1 || \
      mpirun --allow-run-as-root --oversubscribe -np 3 ../../bin/rsf_solve $O -stop_time_yr 0.0001 -fd_stations 60,180 > run.log 2>&1 ) )
  if cmp -s $d/r1/fd_phase_1_st60.dat $d/r3/fd_phase_1_st60.dat && cmp -s $d/r1/fd_phase_1_st180.dat $d/r3/fd_phase_1_st180.dat; then mpi="1 and 3 ranks identical"; else mpi="1 and 3 ranks DIFFER"; ok=0; fi
fi
echo "test_fd3d_pair: `[ $ok = 1 ] && echo PASS || echo FAIL` (whole-space gates: within $(echo $g | cut -d' ' -f1), across $(echo $g | cut -d' ' -f2), normal $(echo $g | cut -d' ' -f3); $mpi)"
[ $ok = 1 ]
