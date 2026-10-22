"""Harris and Day case D with the prototype, per-step trace for the C test
(src/testing/rsf_solve/fd_dynamic.c, case hd). usage: dump_hd.py W_m dx tend out"""
import sys, numpy as np, shells as sh, stepper as st
W=float(sys.argv[1]); dx=float(sys.argv[2]); tend=float(sys.argv[3]); out=sys.argv[4]
AL,BE,RHO=6000.,3464.,2700.; L=28e3; n=int(L/dx); overlap=5e3
f1=sh.Fault([-L/2,0.],[1.,0.],L,n,"f1"); f2=sh.Fault([L/2-overlap,W],[1.,0.],L,n,"f2")
dt=0.5*dx/BE; nsteps=int(round(tend/dt))
sy=sh.System([f1,f2],dt,AL,BE,RHO,tol=1e-2,lmax=nsteps+2,cap_ok=True)
x=np.concatenate([f1.centers[:,0],f2.centers[:,0]]); ft=np.full(sy.N,np.inf)
inner=(np.arange(sy.N)<n)&(np.abs(x)<2.5e3); ft[inner]=np.abs(x[inner])/(0.5*BE)
fr=st.SlipWeakening(0.75,0.51,0.10,force_time=ft,force_t0=0.2)
import os
S=st.Stepper(sy,fr,tau0=20e6,sigma0=33.3e6,normal_coupling=True)
if os.environ.get("FD_NL"): S.enable_nl()
pts=[n-1, n//2+50, n+int(5.5e3/dx), n+int(1e3/dx)]
with open(out,"w") as f:
    f.write(f"# hd W {W} dx {dx} n {n} pts {' '.join(map(str,pts))}: t, delta[pts], V[pts], sigma[pts]\n")
    for k in range(nsteps):
        S.step()
        f.write(f"{S.t:.10e} "+" ".join(f"{S.delta[i]:.10e}" for i in pts)+" "+" ".join(f"{S.V[i]:.10e}" for i in pts)+" "+" ".join(f"{S.sigma[i]:.10e}" for i in pts)+"\n")
print("done", nsteps, "steps; final max slip f2", S.delta[n:].max())
