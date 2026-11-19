#!/bin/sh
# 3-D static gate: a vertical strike-slip fault of 16 x 8 squares (500 m) in a
# whole space; rsf_solve -fd_mode -1 builds the 3-D tables and compares their
# static limit with the assembled Is. Run from the repository root.
d=tmp_fd3d_gate
mkdir -p $d && cd $d
python3 - << 'PY'
dx=500.0; nx,nz=16,8
with open("geom.in","w") as g, open("rsf.dat","w") as r:
    for i in range(nx):
        for j in range(nz):
            x=-4000+(i+0.5)*dx; z=-1000-(j+0.5)*dx
            g.write(f"{x:.6e} 0.000000e+00 {z:.6e} 90 90 {dx/2:.6e} {dx/2:.6e} 0\n"); r.write("0.010 0.015\n")
PY
../bin/rsf_solve -geom_file geom.in -rsf_file rsf.dat -dc 0.016 -sigma_init 50e6 -full_space 1 \
  -stop_time_yr 0.001 -fd_mode -1 > run.log 2>&1
grep "3-D static gate" run.log
awk '/3-D static gate/ {for(i=1;i<=NF;i++) if($i=="Is") dev=$(i+1)} END {ok=(dev!="" && dev+0 < 1e-6); print "test_fd3d_gate: " (ok ? "PASS" : "FAIL") " (deviation " dev ")"; exit ok?0:1}' run.log
