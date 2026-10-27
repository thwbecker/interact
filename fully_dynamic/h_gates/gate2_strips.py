"""Gate 2b: exact wavefront strips from the nearest and farthest points of the source square,
and whether the response is exactly linear in t (static) after the S strip; then the storage per
source-receiver pair and component: F strips + 4 polynomial coefficients (I) + static slope (S),
against the present dense window."""
import sys, numpy as np
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), ".."))
import tada_quad as tq
cL, cT, rho, dx = 6000.0, 3464.0, 2670.0, 500.0
a = 0.5*dx; dt = 0.5*dx/cT
comps = ["s11","s22","s33","s12","s31","s32"]
def rminmax(X):
    # nearest and farthest point of the square |x1|,|x2| <= a in the plane x3 = 0
    n1 = np.clip(X[0], -a, a); n2 = np.clip(X[1], -a, a)
    rmin = np.sqrt((X[0]-n1)**2 + (X[1]-n2)**2 + X[2]**2)
    f1 = -a*np.sign(X[0]) if X[0] != 0 else a; f2 = -a*np.sign(X[1]) if X[1] != 0 else a
    rmax = np.sqrt((X[0]-f1)**2 + (X[1]-f2)**2 + X[2]**2)
    return rmin, rmax
rng = np.random.default_rng(2)
print("r[km]   F strip steps (P, S)   I steps   I residual(cubic)   after-S residual(linear)")
for r in (1.5e3, 3e3, 10e3, 30e3):
    for k in range(2):
        u = rng.normal(size=3); u /= np.linalg.norm(u); X = r*u
        rmin, rmax = rminmax(X)
        lP0, lP1 = int(np.floor(rmin/cL/dt)), int(np.ceil(rmax/cL/dt)) + 1
        lS0, lS1 = int(np.floor(rmin/cT/dt)), int(np.ceil(rmax/cT/dt)) + 1
        t = np.arange(1, lS1 + 40)*dt
        S = tq.rect_ramp_stress(np.full(t.size, X[0]), np.full(t.size, X[1]), X[2], t, a, a, cL, cT, rho, comp=1)
        Smax = max(np.abs(S[c]).max() for c in comps)
        mI = (t > lP1*dt) & (t < lS0*dt); mA = t > lS1*dt
        eI = eA = 0.0
        for c in comps:
            if mI.sum() >= 5:
                p = np.polyfit(t[mI]/t.max(), S[c][mI], 3); eI = max(eI, np.abs(np.polyval(p, t[mI]/t.max()) - S[c][mI]).max())
            p = np.polyfit(t[mA], S[c][mA], 1); eA = max(eA, np.abs(np.polyval(p, t[mA]) - S[c][mA]).max())
        print(f"{r/1e3:5.1f}   {lP1-lP0:3d}, {lS1-lS0:3d}            {max(lS0-lP1,0):4d}     {eI/Smax:9.1e}          {eA/Smax:9.1e}")
