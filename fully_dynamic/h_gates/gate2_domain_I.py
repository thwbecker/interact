"""Gate 2: does the ramp response of a square source factorise between the P and S fronts
(Ando 2016, Domain I)? For receivers at distance r in several directions, fit the response S(t)
(and the step increments W) on the interval between the end of the P strip and the start of the
S strip with low-degree polynomials in t, and report the relative residual."""
import sys, numpy as np
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), ".."))
import tada_quad as tq
cL, cT, rho, dx = 6000.0, 3464.0, 2670.0, 500.0
a = 0.5*dx; dt = 0.5*dx/cT; diag = np.sqrt(2)*a
comps = ["s11","s22","s33","s12","s31","s32"]
rng = np.random.default_rng(1)
print("r[km]  dir                      tP_end  tS_beg  nI(steps)  max rel residual for degree 1 / 2 / 3 / 4  (relative to max |S| over the full window)")
for r in (2e3, 5e3, 10e3, 30e3):
    for k in range(3):
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        X = r*u
        tP_end = (r + diag)/cL + 2*dt          # P strip: the source's arrivals plus a margin
        tS_beg = (r - diag)/cT - 2*dt          # S strip begins
        tfull = np.arange(1, int(((r + diag)/cT + 6*dt)/dt) + 1)*dt
        S = tq.rect_ramp_stress(np.full(tfull.size, X[0]), np.full(tfull.size, X[1]), X[2], tfull, a, a, cL, cT, rho, comp=1)
        Smax = max(np.abs(S[c]).max() for c in comps)
        m = (tfull > tP_end) & (tfull < tS_beg)
        if m.sum() < 6: continue
        res = []
        for deg in (1, 2, 3, 4):
            e = 0.0
            for c in comps:
                y = S[c][m]; tt = tfull[m]
                p = np.polyfit(tt/tt.max(), y, deg)
                e = max(e, np.abs(np.polyval(p, tt/tt.max()) - y).max())
            res.append(e/Smax)
        print(f"{r/1e3:5.0f}  ({u[0]:+.2f},{u[1]:+.2f},{u[2]:+.2f})  {tP_end:6.2f}  {tS_beg:6.2f}  {m.sum():6d}     " + "  ".join(f"{x:.1e}" for x in res))
