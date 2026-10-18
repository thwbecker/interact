import numpy as np, sys
W=float(sys.argv[1])*1e3; ov=float(sys.argv[2])*1e3; vs_end=float(sys.argv[3])*1e3
dx=100.0; L1=L2=24e3; n1=int(L1/dx); n2=int(L2/dx)
x1=-L1+(np.arange(n1)+0.5)*dx; x2=-ov+(np.arange(n2)+0.5)*dx
xs=np.concatenate([x1,x2]); ys=np.concatenate([np.zeros(n1),np.full(n2,W)])
s_along=np.concatenate([x1+L1,x2+ov]); Lseg=np.concatenate([np.full(n1,L1),np.full(n2,L2)])
vs=(s_along<vs_end)|(s_along>Lseg-vs_end)
a=np.where(vs,0.020,0.010); b=np.full(n1+n2,0.015); V0=1e-6; f0=0.6; Vpl=1e-9; sig=50e6
psi_ss=f0+b*np.log(V0/Vpl); Vinit=np.full(n1+n2,Vpl); Vinit[:n1][np.abs(x1+L1/2)<2e3]*=1.5
tau0=a*sig*np.arcsinh(Vinit/(2*V0)*np.exp(psi_ss/a))
with open("geom.in","w") as g, open("rsf.dat","w") as r, open("ic.in","w") as ic:
    for i in range(n1+n2):
        g.write(f"{xs[i]:.6e} {ys[i]:.6e} 0 90 90 {dx/2:.6e} 0 0\n"); r.write(f"{a[i]:.6f} {b[i]:.6f}\n"); ic.write(f"{tau0[i]:.8e} {Vinit[i]:.6e}\n")
print("elements", n1+n2)
