"""Tests of the 3-D mirror + opening-layer prototype (fs3d_dyn.py):
1. causality: with a buried fault (top 500 m) the layer cannot affect the fault
   before 2 top / cP; until then the run equals the mirror-only run;
2. static limit: after a prescribed slip pulse the fault traction approaches
   the static traction of the final slip with the exact half-space operator
   (Okada, fs3d_static.py); the mirror alone is off by its static error;
3. stability of the surface layer without a filter (NL), and the first-order
   scheme for comparison."""
import sys, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import fs3d_dyn as F, fs3d_static as S

def okada_exact(f, top):
    S.MU = F.MU; S.NU = 0.25; S.LAM = 2*S.MU*S.NU/(1 - 2*S.NU)
    h = f.h; K = np.zeros((f.Nf, f.Nf))
    for j, src in enumerate(f.fault):
        for i, X in enumerate(f.fault):
            K[i, j] = -S.okada_stress(X[0] - src[0], 0.0, X[2], -(src[2] - h/2), 90.0, (-h/2, h/2), (0.0, h), (1.0, 0.0, 0.0), False)[0][0, 1]
    return K           # sign: Tada's fault traction convention is minus Okada's sigma_xy here

def run(top=500.0, T=3.5, trise=0.4, nl=True, L=45, hs=250.0):
    f = F.FS3D(Lx=2e3, W=1e3, top=top, h=500.0, hs=hs, ext=2e3, margin=1e3, L=L)
    g = F.FS3D(Lx=2e3, W=1e3, top=top, h=500.0, hs=hs, ext=2e3, margin=1e3, L=L, with_layer=False, verbose=False)
    tl, vs, m = F.kinematic(f, T=T, trise=trise, nl=nl)
    tm, _, mm = F.kinematic(g, T=T, trise=trise, nl=nl)
    t = (np.arange(len(tl)) + 1)*f.dt
    tc = 2*top/F.AL
    early = np.abs(tl[t < tc - f.dt] - tm[t < tc - f.dt]).max() if (t < tc - f.dt).any() else 0.0
    Kex = okada_exact(f, top)
    slip = m.delta[:f.Nf]
    stat_ex = Kex @ slip
    e_layer = np.abs(tl[-1] - stat_ex).max()/np.abs(stat_ex).max()
    e_mirror = np.abs(tm[-1] - stat_ex).max()/np.abs(stat_ex).max()
    print(f"{'NL' if nl else 'first order'} top {top:g} m, layer {hs:g} m: before {tc:.3f} s max |tau(layer) - tau(mirror)| = {early:.1e} Pa; "
          f"at {T:g} s |tau - exact static| / max: with layer {e_layer:.2e}, mirror only {e_mirror:.2e}; "
          f"surface opening rate max over the last 20 % {vs[int(0.8*len(vs)):].max():.2e} (peak {vs.max():.2e})")
    return early, e_layer, e_mirror, vs

if __name__ == "__main__":
    early, el, em, vs = run(nl=True)
    assert early < 1e-6 and el < 0.5*em, (early, el, em)
    run(nl=False)
    print("test_fs3d_dyn: PASS")
