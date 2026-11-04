"""Compare the batch reflected half-plane tables (bin/test_fd_hpdwn_batch output, fd_hpdwn.c) with
Crouch and Starfield's static half-plane solution (TDD_LIB, fs_static2d.py) and with the prototype
halfplane_dwn.py (fd_plan.md M10). Usage: python3 test_fd_hpdwn_batch.py output.txt
PASS: static within 1e-4, lag series within 2e-2 of the largest value (the batch evaluates the
reflected fronts with FD_HPDWN_Q times the lag Nyquist frequency, default 8)."""
import sys, os, numpy as np
os.environ.setdefault('TDD_LIB','/tmp/tdd/libtdd.so')
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import halfplane_dwn as hd, fs_static2d as fs
al,be,rho=6000.,3464.,2670.; fs.MU=rho*be**2; fs.NU=(al**2-2*be**2)/(2*(al**2-be**2))
L=open(sys.argv[1]).read().split('\n')
hdr=L[0].split(); n=int(hdr[2]); h=float(hdr[4]); top=float(hdr[6]); nt=int(hdr[8]); dt=float(hdr[10])
dip=np.radians(30); tx,tz=np.cos(dip),np.sin(dip)
X=[((i+0.5)*h*tx, top+(i+0.5)*h*tz) for i in range(n)]
K={}; S={}
for l in L[1:]:
    w=l.split()
    if not w: continue
    if w[0]=='K': K[(int(w[1]),int(w[2]))]=(float(w[3]),float(w[4]))
    else: S.setdefault((int(w[1]),int(w[2])),[]).append((float(w[4]),float(w[5])))
def project(sxx,sxz,szz):
    sr=np.array([tx,-tz]); nr=np.array([-sr[1],sr[0]]); sxy=-sxz; syy=szz
    trx=sxx*nr[0]+sxy*nr[1]; tr_y=sxy*nr[0]+syy*nr[1]
    return trx*sr[0]+tr_y*sr[1], trx*nr[0]+tr_y*nr[1]
# static: reflected = -(tdd half - tdd full), projected
wk=0; sk=0
for (i,j),(ks,kn) in K.items():
    xs,zs=X[j]; xr,zr=X[i]
    h1,_=fs.tdd(xr,-zr,xs,-zs,h/2,tx,-tz,1); f1,_=fs.tdd(xr,-zr,xs,-zs,h/2,tx,-tz,0)
    d=-(np.array(h1)-np.array(f1))
    rs,rn=project(d[0],-d[1],d[2])
    wk=max(wk,abs(ks-rs),abs(kn-rn)); sk=max(sk,abs(rs),abs(rn))
print('static reflected: largest |C - Crouch-Starfield| %.3e of %.3e (relative %.1e)'%(wk,sk,wk/sk))
# series for a few pairs against the prototype
ts=dt*np.arange(nt)
we=0; se=0
for (i,j) in ((0,0),(1,0),(3,1),(0,3)):
    xs,zs=X[j]; xr,zr=X[i]
    P=hd.ramp_stress(xr,zr,xs,zs,(tx,tz),h/2,ts[1:],al,be,rho,part='refl',T=2.5*ts.max(),N=2048,wI_T=18.0)
    ps,pn=project(P[:,0],P[:,1],P[:,2])
    c=np.array(S[(i,j)])[1:]
    we=max(we,np.abs(c[:,0]-ps).max(),np.abs(c[:,1]-pn).max()); se=max(se,np.abs(ps).max(),np.abs(pn).max())
print('reflected series: largest |C batch - prototype| %.3e of %.3e (relative %.1e)'%(we,se,we/se))
ok = wk/sk < 1e-4 and we/se < 2e-2
print('test_fd_hpdwn_batch.py: static %.1e, series %.1e: %s' % (wk/sk, we/se, 'PASS' if ok else 'FAIL'))
