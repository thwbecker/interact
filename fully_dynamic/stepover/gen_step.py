"""Two parallel in-plane segments for rsf_solve stepover sequences.
usage: gen_step.py W overlap margin [L2_km] [dba2]
  W        step width in km (+W: second segment at +y)
  overlap  overlap in km (negative: gap)
  margin   velocity-strengthening margin at each segment end in km
  L2_km    length of the second segment (default 24, same as the first)
  dba2     b - a on the velocity-weakening part of the second segment
           (default 0.005, same as the first)
Writes geom.in, rsf.dat, ic.in and segments.txt (n1 n2)."""
import numpy as np, sys
W=float(sys.argv[1])*1e3; ov=float(sys.argv[2])*1e3; vs_end=float(sys.argv[3])*1e3
L2=float(sys.argv[4])*1e3 if len(sys.argv) > 4 else 24e3
dba2=float(sys.argv[5]) if len(sys.argv) > 5 else 0.005
dx=100.0; L1=24e3; n1=int(L1/dx); n2=int(round(L2/dx))
x1=-L1+(np.arange(n1)+0.5)*dx; x2=-ov+(np.arange(n2)+0.5)*dx
xs=np.concatenate([x1,x2]); ys=np.concatenate([np.zeros(n1),np.full(n2,W)])
s_along=np.concatenate([x1+L1,x2+ov]); Lseg=np.concatenate([np.full(n1,L1),np.full(n2,L2)])
vs=(s_along<vs_end)|(s_along>Lseg-vs_end)
b=np.full(n1+n2,0.015)
a_vw=np.concatenate([np.full(n1,0.010),np.full(n2,0.015-dba2)])
a=np.where(vs,0.020,a_vw); V0=1e-6; f0=0.6; Vpl=1e-9; sig=50e6
psi_ss=f0+b*np.log(V0/Vpl); Vinit=np.full(n1+n2,Vpl); Vinit[:n1][np.abs(x1+L1/2)<2e3]*=1.5
tau0=a*sig*np.arcsinh(Vinit/(2*V0)*np.exp(psi_ss/a))
with open("geom.in","w") as g, open("rsf.dat","w") as r, open("ic.in","w") as ic:
    for i in range(n1+n2):
        g.write(f"{xs[i]:.6e} {ys[i]:.6e} 0 90 90 {dx/2:.6e} 0 0\n"); r.write(f"{a[i]:.6f} {b[i]:.6f}\n"); ic.write(f"{tau0[i]:.8e} {Vinit[i]:.6e}\n")
open("segments.txt","w").write(f"{n1} {n2}\n")
print("elements", n1+n2)
