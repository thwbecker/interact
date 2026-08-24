import numpy as np, subprocess
P="/home/claude/kaj/pom"
def run(mvec, t, stations, zr):
    with open("fd_sta.txt","w") as f:
        for x,y in stations: f.write(f"{x} {y}\n")
    subprocess.run([P,"maxwell","-m"]+[f"{v}" for v in mvec]+["-t",f"{t}","-tR","25","-H","20","-mulam","1","-zr",f"{zr}","-sta","fd_sta.txt","-o","fd.out"],check=True,capture_output=True)
    return np.atleast_2d(np.loadtxt("fd.out"))[:,2:]
h=0.02; lam=mu=1.0
worst=0.0
mss=[2,2,11,90,0,0,0,1,0,0]; mds=[2,2,10.5,30,0,0.866,0,0,1,0]
for mvec,mn in ((mss,"ss90"),(mds,"ds30")):
    for t in (10.0, 100.0, 500.0):
        for zr in (2.0, 8.0, 14.0):
            x,y = 12.0, 6.0
            C  = run(mvec,t,[(x,y),(x+h,y),(x-h,y),(x,y+h),(x,y-h)], zr)
            Cu = run(mvec,t,[(x,y)], zr-h); Cd = run(mvec,t,[(x,y)], zr+h)
            dx=(C[1,:3]-C[2,:3])/(2*h); dy=(C[3,:3]-C[4,:3])/(2*h); dzr=(Cd[0,:3]-Cu[0,:3])/(2*h)
            euu=-dzr[2]; eh=dx[0]+dy[1]
            szz=lam*(eh+euu)+2*mu*euu
            sxz=mu*(-dzr[0]+dx[2]); syz=mu*(-dzr[1]+dy[2])
            scl=max(abs(szz),abs(sxz),abs(syz),1e-18)
            dev=max(abs(C[0,3]-szz),abs(C[0,4]-sxz),abs(C[0,5]-syz))/scl
            worst=max(worst,dev)
print(f"FD matrix (2 mech x 3 t x 3 zr): worst rel dev = {worst:.1e}")
# free-surface limit after corrections
Cs = run(mss,100.0,[(10,0),(30,0)],0.0001)
print("tractions at zr->0 (post-correction):", " ".join(f"{v:.2e}" for v in Cs[:, 3:].flatten()))
# t->0
C0 = run(mss,1e-4,[(12,6)],5.0)
print("t->0 interior state:", " ".join(f"{v:.2e}" for v in C0[0]))
