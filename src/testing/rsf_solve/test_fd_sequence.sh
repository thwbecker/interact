#!/bin/sh
# rsf_solve with -fd_mode 1 on the prototype's single-fault test (30 km in-plane
# fault, 16 km velocity-weakening patch, 100 m), two events: checks the first
# dynamic phase against the Python prototype (2013 steps over 1370 s, peak
# 2.43 m/s, mean slip 0.984 m) and the recurrence (82.1 yr). About 1 minute.
# Run from the repository root after make bin/rsf_solve.
d=tmp_fd_sequence
mkdir -p $d && cd $d
python3 - << 'PY'
import numpy as np
dx=100.0; L=30e3; n=int(L/dx); x=-L/2+(np.arange(n)+0.5)*dx
a=np.where(np.abs(x)<8e3,0.010,0.020); b=np.full(n,0.015); V0=1e-6; f0=0.6; Vpl=1e-9; sig=50e6
psi_ss=f0+b*np.log(V0/Vpl); Vinit=Vpl*np.where(np.abs(x)<2e3,1.5,1.0)
tau0=a*sig*np.arcsinh(Vinit/(2*V0)*np.exp(psi_ss/a))
with open("geom.in","w") as g, open("rsf.dat","w") as r, open("ic.in","w") as ic:
    for i in range(n):
        g.write(f"{x[i]:.6e} 0.000000e+00 0 90 90 {dx/2:.6e} 0 0\n"); r.write(f"{a[i]:.6f} {b[i]:.6f}\n"); ic.write(f"{tau0[i]:.8e} {Vinit[i]:.6e}\n")
PY
../bin/rsf_solve -geom_file geom.in -rsf_file rsf.dat -rsf_ic_file ic.in -dc 0.016 -sigma_init 50e6 \
  -v0 1e-6 -f0 0.6 -vpl 1e-9 -shear_modulus 3.24e10 -s_wave_speed 3464 -full_space 1 -rtol 1e-6 \
  -stop_time_yr 90 -track_events 1 -vel_event 1e-3 -dt_monitor_yr 50 -fd_mode 1 > run.log 2>&1
grep "rsf_fd: phase" run.log
awk '/rsf_fd: phase 1:/ {for(i=1;i<=NF;i++){if($i=="phase") ; if($i=="steps") st=$(i-1); if($i=="slip" && $(i-1)=="mean") ms=$(i+1); if($i=="rate") pk=$(i+1)}; t1=$5}
     /rsf_fd: phase 2:/ {t2=$5}
     END {ok=1;
          if(st<2000||st>2030){print "FAIL: phase 1 steps",st; ok=0}
          if(ms<0.975||ms>0.995){print "FAIL: phase 1 mean slip",ms; ok=0}
          if(pk<2.35||pk>2.52){print "FAIL: phase 1 peak slip rate",pk; ok=0}
          if(t2-t1<81.5||t2-t1>82.8){print "FAIL: recurrence",t2-t1; ok=0}
          if(ok) print "test_fd_sequence: PASS (phase 1 " st " steps, mean slip " ms " m, peak " pk " m/s, recurrence " t2-t1 " yr)"; else print "test_fd_sequence: FAIL"; exit ok?0:1}' run.log
