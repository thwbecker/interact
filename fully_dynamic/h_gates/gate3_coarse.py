"""Gate 3c: coarse-graining with non-uniform slip: 4 x 4 squares of 250 m with slip rates varying
linearly by +-50 percent (and randomly by +-50 percent) across a 1 km group, against one 1 km square
carrying the mean rate, at receivers 2 to 20 km away: largest difference of the traction history (cumulative increments)
over all lags, relative to the largest |W| of the fine group."""
import sys, numpy as np
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), ".."))
import tada_quad as tq
cL, cT, rho = 6000.0, 3464.0, 2670.0
dt = 0.5*500.0/cT; afine = 125.0; n = 4
rng = np.random.default_rng(7)
cells = [np.array([(i + 0.5)*2*afine - 500.0, (j + 0.5)*2*afine - 500.0, 0.0]) for i in range(n) for j in range(n)]
def W(X, c, a, t):
    x = X - c
    S = tq.rect_ramp_stress(np.full(t.size, x[0]), np.full(t.size, x[1]), x[2], t, a, a, cL, cT, rho, comp=1)
    return np.diff(np.array([S["s31"], S["s32"], S["s33"]]), axis=1)
for kind in ("linear", "random"):
    v = np.array([1.0 + 0.5*(c[0] + c[1])/750.0 for c in cells]) if kind == "linear" else 1.0 + 0.5*(2*rng.random(len(cells)) - 1)
    print(f"slip rates {kind} within the group (mean {v.mean():.2f}):")
    for r in (2e3, 3e3, 5e3, 10e3, 20e3):
        worst = 0.0
        for k in range(4):
            u = rng.normal(size=3); u /= np.linalg.norm(u); X = r*u
            t = np.arange(0, int((r + 1500)/cT/dt) + 8)*dt
            Wf = sum(vi*W(X, c, afine, t) for vi, c in zip(v, cells))
            Wc = v.mean()*W(X, np.zeros(3), 500.0, t)
            Cf, Cc = np.cumsum(Wf, axis=1), np.cumsum(Wc, axis=1); worst = max(worst, np.abs(Cf - Cc).max()/np.abs(Cf).max())
        print(f"   r {r/1e3:4.0f} km: {worst:.1e}")
