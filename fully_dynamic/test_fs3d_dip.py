"""Tests of the dipping-fault prototype with a full virtual surface (fs3d_dip.py).
1. static: the fault's shear and normal traction operators through the surface
   (Schur complement of the late-time slopes) against Okada's half space
   (interact's dc3dmod.F, via fs3d_static.okada_stress); the full space for
   comparison."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fs3d_dip as D, fs3d_static as S

def okada_ops(f, dip, rake):
    S.MU = D.MU; S.NU = 0.25; S.LAM = 2*S.MU*S.NU/(1 - 2*S.NU)
    h = f.h; F = f.F; Nf = f.Nf
    Ksh = np.zeros((Nf, Nf)); Knh = np.zeros((Nf, Nf)); Ksf = np.zeros((Nf, Nf)); Knf = np.zeros((Nf, Nf))
    disl = (1.0, 0.0, 0.0) if rake == "strike" else (0.0, 1.0, 0.0)
    for j, c in enumerate(f.fault):
        ref = c - 0.5*h*F[1]                              # bottom edge centre
        for i, X in enumerate(f.fault):
            for full, Ks, Kn in ((False, Ksh, Knh), (True, Ksf, Knf)):
                st = S.okada_stress(X[0] - ref[0], X[1] - ref[1], X[2], -ref[2], dip, (-h/2, h/2), (0.0, h), disl, full)[0]
                Ks[i, j] = f.slip_dir @ st @ F[2]; Kn[i, j] = F[2] @ st @ F[2]
    return Ksh, Knh, Ksf, Knf

def static_check(dip=30.0, top=500.0, hs=250.0, ext=2e3, margin=1e3, rake="dip", Lx=2e3, W=1e3, h=500.0):
    f = D.FSDip(Lx=Lx, W=W, dip=dip, top=top, h=h, hs=hs, ext=ext, margin=margin, rake=rake, L=2)
    Nf = f.Nf; s = slice(Nf, None)
    X = np.linalg.solve(f.Ks[s, s], f.Ks[s, :Nf])
    Kss = f.Ks[:Nf, :Nf] - f.Ks[:Nf, s] @ X
    Kns = f.Kn[:, :Nf] - f.Kn[:, s] @ X
    Ksh, Knh, Ksf, Knf = okada_ops(f, dip, rake)
    sg = np.sign(np.sum(Kss*Ksh)); sn = np.sign(np.sum(Kns*Knh)) if np.abs(Knh).max() > 0 else 1.0
    e = lambda A, B: np.linalg.norm(A - B)/np.linalg.norm(B)
    print(f"dip {dip:g}, top {top:g} m, {rake} slip, surface {f.Ns} x 3 of {hs:g} m (+-{ext/1e3:g} km): "
          f"shear op: surface {e(sg*Kss, Ksh):.2e} (full space {e(sg*f.Ks[:Nf, :Nf], Ksh):.2e}, Okada full vs our full {e(sg*f.Ks[:Nf, :Nf], Ksf):.1e}); "
          f"normal op: surface {e(sn*Kns, Knh):.2e} (full space {e(sn*f.Kn[:, :Nf], Knh):.2e})")
    return e(sg*Kss, Ksh), e(sn*Kns, Knh)

def kinematic_check(dip=30.0, top=500.0, hs=250.0, ext=4e3, margin=2e3, rake="dip", T=3.5, L=45, nl=True):
    """causality, late-time static (shear and normal) against Okada, stability"""
    f = D.FSDip(dip=dip, top=top, hs=hs, ext=ext, margin=margin, rake=rake, L=L)
    g = D.FSDip(dip=dip, top=top, hs=hs, ext=ext, margin=margin, rake=rake, L=L, surface=False, verbose=False)
    tl, vs, m = D.kinematic(f, T=T, nl=nl)
    tg, _, mg = D.kinematic(g, T=T, nl=nl)
    t = (np.arange(len(tl)) + 1)*f.dt
    tc = 2*top/D.AL
    early = np.abs(tl[t < tc - f.dt] - tg[t < tc - f.dt]).max() if (t < tc - f.dt).any() else 0.0
    Ksh, Knh, _, _ = okada_ops(f, dip, rake)
    slip = m.delta[:f.Nf]
    sg = np.sign(np.sum(f.Ks[:f.Nf, :f.Nf]*okada_ops(g, dip, rake)[2]))
    ex_s = sg*(Ksh @ slip); ex_n = (Knh @ slip)
    sn_end = m.sn_hist[-1]
    sgn = np.sign(np.sum(sn_end*ex_n)) if np.abs(ex_n).max() > 0 else 1.0
    e_s = np.abs(tl[-1] - ex_s).max()/np.abs(ex_s).max(); e_n = np.abs(sgn*sn_end - ex_n).max()/np.abs(ex_n).max()
    print(f"kinematic {'NL' if nl else 'first order'} dip {dip:g} top {top:g} m: before {tc:.3f} s max |tau - tau(full space)| = {early:.1e} Pa; "
          f"at {T:g} s against Okada: shear {e_s:.2e} (full space {np.abs(tg[-1] - ex_s).max()/np.abs(ex_s).max():.2e}), normal {e_n:.2e}; "
          f"surface rate max over the last 20 % {vs[int(0.8*len(vs)):].max():.2e} (peak {vs.max():.2e})")
    return early, e_s, e_n, vs

if __name__ == "__main__":
    es, en = static_check(ext=4e3, margin=2e3)
    assert es < 1e-3 and en < 1e-2, (es, en)
    early, e_s, e_n, vs = kinematic_check()
    assert early < 1e-6 and e_s < 5e-3 and e_n < 3e-2 and vs[int(0.8*len(vs)):].max() < 0.2*vs.max(), (early, e_s, e_n)
    print("test_fs3d_dip: PASS")
