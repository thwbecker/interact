#!/bin/sh
# rsf_solve -fd_mode 1 with several 3-D planes of different orientation (fd_plan.md G2, G3a):
# a vertical strike-slip fault of two 12 x 8 segments of 1 km squares, the second turned by
# 0.001 or 30 degrees about the vertical through the junction, rate and state, a seed
# patch, -calc_sigma_dot. Whole space: the static gates of the dense cross tables against
# interact's Is and In must pass at 1e-8 (30 degree bend), and the 0.001 degree bend (two planes, the
# current-step coupling at the junction kept) must reproduce the straight fault (one
# plane) at two stations to 1e-8 of the peak rate; half space (images): the same identity.
# The 30 degree bend must run. About 40 s. Run from the repository root.
d=tmp_fd3d_multi
rm -rf $d; mkdir -p $d && cd $d
python3 - << 'PYEOF'
import numpy as np
dx=1000.0; L=12; D=8
for tag, ang, ztop in (("w0",0.0,1000.0),("w001",0.001,1000.0),("w30",30.0,1000.0),("h0",0.0,0.0),("h001",0.001,0.0)):
    with open(f"geom_{tag}.in","w") as g:
        for seg in (0,1):
            phi = 0.0 if seg == 0 else ang
            sd = np.array([np.sin(np.radians(phi)), np.cos(np.radians(phi)), 0.0])
            start = np.zeros(3) if seg == 0 else np.array([0.0, L*dx, 0.0])
            for i in range(L):
                for j in range(D):
                    c = start + (i+0.5)*dx*sd + np.array([0,0,-(j+0.5)*dx-ztop])
                    g.write(f"{c[0]:.9e} {c[1]:.9e} {c[2]:.9e} {phi:.9f} 90.0 {dx/2:.1f} {dx/2:.1f} {seg}\n")
with open("rsf.dat","w") as r, open("ic.in","w") as ic:
    for seg in (0,1):
        for i in range(L):
            for j in range(D):
                edge=(seg==0 and i==0) or (seg==1 and i==L-1) or j in (0,D-1)
                a,b=(0.025,0.015) if edge else (0.010,0.015)
                r.write(f"{a:.6e} {b:.6e}\n")
                v=0.03 if (seg==0 and i in (5,6) and j in (3,4)) else 1e-9
                ic.write(f"{50e6*(a*np.arcsinh(v/(2e-6)*np.exp((0.6+b*np.log(1e-6/1e-9))/a))):.8e} {v:.6e}\n")
PYEOF
O="-rsf_file ../rsf.dat -rsf_ic_file ../ic.in -sigma_init 50e6 -dc 0.02 -stop_time_yr 0.00003 -track_events 1 -vel_event 1e-3 -fd_mode 1 -calc_sigma_dot 1 -fd_stations 44,140"
for tag in w0 w001 w30 h0 h001; do
  mkdir -p $tag
  fs=1; case $tag in h*) fs=0;; esac
  ( cd $tag && ../../bin/rsf_solve -geom_file ../geom_$tag.in $O -full_space $fs > run.log 2>&1 )
done
python3 - << 'PYEOF'
import re, numpy as np
ok = True; msg = []
g = open("w30/run.log").read()
m = re.search(r"from Is ([\d.e+-]+) within planes, ([\d.e+-]+) across, normal across from -In ([\d.e+-]+)", g)
if not m or max(float(m.group(i)) for i in (1,2,3)) > 1e-8: ok = False
msg.append("gates " + (" ".join(m.groups()) if m else "missing"))
for a, b, name in (("w0","w001","whole space"),("h0","h001","half space")):
    e = 0.0; s = 0.0
    for st in (44, 140):
        try:
            x = np.loadtxt(f"{a}/fd_phase_1_st{st}.dat"); y = np.loadtxt(f"{b}/fd_phase_1_st{st}.dat")
        except Exception:
            ok = False; continue
        n = min(len(x), len(y)); e = max(e, np.abs(x[:n,2]-y[:n,2]).max()); s = max(s, np.abs(x[:,2]).max())
        if len(x) != len(y): ok = False
    if not (s > 0 and e <= 1e-8*s): ok = False
    msg.append(f"{name} bend 0.001 deg against straight {e:.2e} of {s:.3f} m/s")
if "3-D phase 1 (2 planes)" not in open("w30/run.log").read(): ok = False; msg.append("30 degree bend did not run")
else: msg.append("30 degree bend ran")
print("test_fd3d_multi:", "PASS" if ok else "FAIL", "(" + "; ".join(msg) + ")")
raise SystemExit(0 if ok else 1)
PYEOF
