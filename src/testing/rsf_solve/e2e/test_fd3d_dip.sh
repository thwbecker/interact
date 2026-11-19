#!/bin/sh
# rsf_solve -fd_mode 1 on a small dipping thrust in a half space (full virtual
# free surface, normal-stress coupling with -calc_sigma_dot): 8 x 4 squares of
# 1 km dipping 30 degrees, top at 500 m, rate and state with a seeded patch.
# Checks that the setup chooses the virtual surface, that the first dynamic
# phase stays stable (finite, peak slip rate below 10 m/s) and that thrust slip
# is promoted relative to the full space (peak slip rate and mean slip larger).
# About 2 to 4 minutes. Run from the repository root after make bin/rsf_solve.
d=tmp_fd3d_dip
rm -rf $d; mkdir -p $d/hs $d/fs
python3 - << 'PYEOF'
import numpy as np
dx=1000.0; nxs=8; nzs=4; dip=30.0; top=500.0; d=np.radians(dip)
with open("tmp_fd3d_dip/geom.in","w") as g, open("tmp_fd3d_dip/rsf.dat","w") as r, open("tmp_fd3d_dip/ic.in","w") as ic:
    for i in range(nxs):
        for j in range(nzs):
            y=(i+0.5-nxs/2)*dx; x=(j+0.5)*dx*np.cos(d); z=-top-(j+0.5)*dx*np.sin(d)
            g.write(f"{x:.6e} {y:.6e} {z:.6e} 0.0 {dip:.1f} {dx/2:.1f} {dx/2:.1f} 0\n")
            edge=(i in (0,nxs-1)) or (j in (0,nzs-1))
            a,b=(0.025,0.015) if edge else (0.010,0.015)
            r.write(f"{a:.6e} {b:.6e}\n")
            v=0.03 if ((i in (3,4)) and (j in (1,2))) else 1e-9
            ic.write(f"{50e6*(a*np.arcsinh(v/(2e-6)*np.exp((0.6+b*np.log(1e-6/1e-9))/a))):.8e} {v:.6e}\n")
PYEOF
O="-geom_file ../geom.in -rsf_file ../rsf.dat -rsf_ic_file ../ic.in -sigma_init 50e6 -dc 0.02 -rsf_slip_mode 1 -stop_time_yr 0.00002 -track_events 1 -vel_event 1e-3 -fd_mode 1 -fd_back_eps 0.05 -calc_sigma_dot 1 -fd_surface_ext 2"
( cd $d/hs && ../../bin/rsf_solve $O > run.log 2>&1 )
( cd $d/fs && ../../bin/rsf_solve $O -full_space 1 > run.log 2>&1 )
mkdir -p $d/ap
( cd $d/ap && ../../bin/rsf_solve $O -fd_hs_static 1 > run.log 2>&1 )
grep -q "virtual free surface" $d/hs/run.log || { echo "test_fd3d_dip: FAIL (virtual surface not selected)"; exit 1; }
hs=`grep "3-D phase 1:" $d/hs/run.log | sed 's/.*peak slip rate \([0-9.e+-]*\) m\/s, mean slip \([0-9.e+-]*\) m.*/\1 \2/'`
fs=`grep "3-D phase 1:" $d/fs/run.log | sed 's/.*peak slip rate \([0-9.e+-]*\) m\/s, mean slip \([0-9.e+-]*\) m.*/\1 \2/'`
ap=`grep "3-D phase 1 (1 planes):" $d/ap/run.log | sed 's/.*peak slip rate \([0-9.e+-]*\) m\/s, mean slip \([0-9.e+-]*\) m.*/\1 \2/'`
# -fd_hs_static (whole-space dynamics, exact half-space statics, fd_plan.md G4) must be within
# 5 percent of the virtual surface for this buried thrust
echo "$hs $fs $ap" | awk '{ ok = ($1 == $1+0) && ($1 < 10.0) && ($1 > $3) && ($2 > $4) && ($5 == $5+0) && ($5 > 0) && (($1 - $5)/$1 < 0.05) && (($5 - $1)/$1 < 0.05) && (($2 - $6)/$2 < 0.05) && (($6 - $2)/$2 < 0.05);
  printf("test_fd3d_dip: %s (half space: peak %.3f m/s, mean slip %.4f m; full space: %.3f m/s, %.4f m; -fd_hs_static: %.3f m/s, %.4f m)\n", ok ? "PASS" : "FAIL", $1, $2, $3, $4, $5, $6); exit !ok }'
