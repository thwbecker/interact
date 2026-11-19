#!/bin/sh
# rsf_solve -fd_mode 1 with several dipping faults sharing the strike in a half space, one
# virtual surface (fd_plan.md G3b): a buried thrust (dip 30 degrees, 8 x 8 squares of 1 km,
# top 1 km deep, rate and state, seed patch), (a) straight: the one-fault surface path,
# (b) its lower half turned to 30.001 degrees: two faults on the several-fault surface path.
# (b) must reproduce (a) at two stations to the kink's own effect (2e-4 of the peak rate);
# (c) the lower half dipping 50 degrees with -calc_sigma_dot must run. About 3 minutes.
d=tmp_fd3d_msurf
rm -rf $d; mkdir -p $d && cd $d
python3 - << 'PYEOF'
import numpy as np
dx=1000.0; nx=8; nz=8; top=1000.0
def geo(tag, d1, d2):
    with open(f"geom_{tag}.in","w") as g:
        P=np.array([0.0, 0.0, -top])
        for j in range(nz):
            d = d1 if j < 4 else d2
            dd = np.radians(d); v = np.array([0.0, -np.cos(dd), -np.sin(dd)])
            if j == 0: c0 = P + 0.5*dx*v
            elif j < 4: c0 = c0 + dx*v
            elif j == 4:
                vp = np.array([0.0, -np.cos(np.radians(d1)), -np.sin(np.radians(d1))]); c0 = c0 + 0.5*dx*vp + 0.5*dx*v
            else: c0 = c0 + dx*v
            for i in range(nx):
                c = c0 + np.array([(i+0.5)*dx, 0, 0])
                g.write(f"{c[0]:.9e} {c[1]:.9e} {c[2]:.9e} 90.0 {d:.9f} {dx/2:.1f} {dx/2:.1f} 0\n")
geo("S", 30.0, 30.0); geo("K", 30.0, 30.001); geo("K50", 30.0, 50.0)
with open("rsf.dat","w") as r, open("ic.in","w") as ic:
    for j in range(nz):
        for i in range(nx):
            edge = i in (0,nx-1) or j in (0,nz-1)
            a,b=(0.025,0.015) if edge else (0.010,0.015)
            r.write(f"{a:.6e} {b:.6e}\n")
            v=0.03 if (i in (3,4) and j in (3,4)) else 1e-9
            ic.write(f"{50e6*(a*np.arcsinh(v/(2e-6)*np.exp((0.6+b*np.log(1e-6/1e-9))/a))):.8e} {v:.6e}\n")
PYEOF
O="-rsf_file ../rsf.dat -rsf_ic_file ../ic.in -sigma_init 50e6 -dc 0.04 -stop_time_yr 0.00003 -track_events 1 -vel_event 1e-3 -fd_mode 1 -rsf_slip_mode 1 -fd_stations 27,35 -fd_surface_grade 0 -fd_surface_ext 2 -fd_surface_extend 0"
for t in S K; do mkdir -p $t; ( cd $t && ../../bin/rsf_solve -geom_file ../geom_$t.in $O > run.log 2>&1 ); done
mkdir -p K50; ( cd K50 && ../../bin/rsf_solve -geom_file ../geom_K50.in $O -calc_sigma_dot 1 > run.log 2>&1 )
python3 - << 'PYEOF'
import numpy as np
ok = True; e = 0.0; s = 0.0
for st in (27, 35):
    try:
        a = np.loadtxt(f"S/fd_phase_1_st{st}.dat"); b = np.loadtxt(f"K/fd_phase_1_st{st}.dat")
    except Exception:
        ok = False; continue
    n = min(len(a), len(b)); e = max(e, np.abs(a[:n,2]-b[:n,2]).max()); s = max(s, np.abs(a[:,2]).max())
    if len(a) != len(b): ok = False
if not (s > 0 and e <= 2e-4*s): ok = False
ran = "3-D phase 1 (2 planes)" in open("K50/run.log").read()
if not ran: ok = False
print("test_fd3d_msurf:", "PASS" if ok else "FAIL", f"(0.001 degree kink against the straight fault {e:.2e} of {s:.4f} m/s; dip-50 lower half {'ran' if ran else 'did not run'})")
raise SystemExit(0 if ok else 1)
PYEOF
