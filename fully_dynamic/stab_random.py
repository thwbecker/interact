"""Stability test after Noda, Sato and Kurihara (2020): traction-free in-plane
(mode II) fault of N elements with random initial traction; growth rate of
||V|| over the second half of the run (c_p t/dx from 160 to 320).
usage: stab_random.py [N] ; FD_NL=1 for the second-order scheme, FD_HT=0.5 (cS dt/dx)"""
import sys, os, numpy as np
import shells as sh, stepper as st
AL, BE, RHO = 6000.0, 3464.0, 2670.0
N = int(sys.argv[1]) if len(sys.argv) > 1 else 64
hT = float(os.environ.get("FD_HT", "0.5"))
dx = 100.0; dt = hT*dx/BE
nsteps = int(round(320*dx/(AL*dt)))
f = sh.Fault([-N*dx/2, 0.0], [1.0, 0.0], N*dx, N)
sysm = sh.System([f], dt, AL, BE, RHO, tol=1e-3, lmax=nsteps + 2, cap_ok=True)
rng = np.random.default_rng(1)
tau0 = rng.normal(0.0, 1e6, N)
class TractionFree:
    """tau = 0 exactly: solve H + A V = 0 (both signs of V allowed)"""
    def solve(self, H, K1, *args, **kw):
        import scipy.sparse as sps, scipy.sparse.linalg as spl
        A = K1.toarray() if sps.issparse(K1) else K1
        state = args[3] if len(args) > 3 else None
        return np.linalg.solve(A, -H), state
fr = TractionFree()
S = st.Stepper(sysm, fr, tau0=tau0, sigma0=1e6, normal_coupling=False)
if os.environ.get("FD_NL"): S.enable_nl()
t, nv = [], []
for m in range(nsteps):
    S.step(); t.append(S.t*AL/dx); nv.append(np.linalg.norm(S.V))
t = np.array(t); nv = np.log10(np.maximum(np.array(nv), 1e-300))
sel = (t >= 160) & (nv > -250)
if sel.sum() < 10: sel = (nv > -250) & (t > 20)
print("  log10 ||V|| at cp t/dx = 20, 80, 160, 240, 320:", " ".join(f"{np.interp(x, t, nv):.1f}" for x in (20, 80, 160, 240, 320)))
# envelope slope: fit to the running maxima over windows
w = max(1, int(len(t[sel])/20)); ts, vs = [], []
for i in range(0, sel.sum() - w + 1, w):
    j = np.argmax(nv[sel][i:i + w]); ts.append(t[sel][i + j]); vs.append(nv[sel][i + j])
g = np.polyfit(ts, vs, 1)[0]
print(f"{'NL' if os.environ.get('FD_NL') else 'first order'} hT {hT}: growth rate {g:+.4f} per cp t/dx (log10 ||V||: start {nv[0]:.2f}, end {nv[-1]:.2f})")
