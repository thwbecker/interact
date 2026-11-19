"""Checks of the Cagniard-de Hoop prototype (halfplane_cdh.py, fd_plan.md M11):
1. reflected ramp responses of a 200 m segment dipping 30 degrees against halfplane_dwn.py for a
   buried and a shallow far-offset geometry (DWN at 8192 frequencies is itself accurate to about
   1e-3 near the fronts), PASS below 2e-3 of the largest value;
2. with TDD_LIB (fs_static2d.py): the slope at 20 and 40 t0 against Crouch and Starfield's reflected
   static (half minus full plane), PASS if it approaches it (below 2e-2 at 20 t0, below 5e-3 at 40 t0).
Usage: python3 check_halfplane_cdh.py"""
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import halfplane_cdh as hc
import halfplane_dwn as hd

st = (np.cos(np.radians(30)), np.sin(np.radians(30)))
ok = True
for (X0, zs0, zr) in ((700., 1500., 900.), (-3000., 600., 150.)):
    ts = np.array([0.6, 1.0, 1.6, 2.4])
    D = hd.ramp_stress(X0, zr, 0.0, zs0, st, 100.0, ts, hc.AL, hc.BE, hc.RHO, part="refl", T=8.0, N=8192, wI_T=18.0)
    C = np.array([hc.seg_total(t, X0, zr, 0.0, zs0, st, 100.0, nq=24) for t in ts])
    r = np.abs(C - D).max()/np.abs(D).max()
    ok = ok and r < 2e-3
    print(f"check_halfplane_cdh: X {X0:.0f} m, source {zs0:.0f} m, receiver {zr:.0f} m deep: against DWN {r:.1e}")
if os.path.exists(os.environ.get("TDD_LIB", "/tmp/tdd/libtdd.so")):
    os.environ.setdefault("TDD_LIB", "/tmp/tdd/libtdd.so")
    import fs_static2d as fs
    fs.MU = hc.RHO*hc.BE**2
    fs.NU = (hc.AL**2 - 2*hc.BE**2)/(2*(hc.AL**2 - hc.BE**2))
    X0, zs0, zr = 700., 1500., 900.
    h1, _ = fs.tdd(X0, -zr, 0.0, -zs0, 100.0, st[0], -st[1], 1)
    f1, _ = fs.tdd(X0, -zr, 0.0, -zs0, 100.0, st[0], -st[1], 0)
    K = -(np.array(h1) - np.array(f1))
    K = np.array([K[0], -K[1], K[2]])
    t0 = hc.first_arrival(X0, zr, 0.0, zs0, st, 100.0)
    for m, lim in ((20.0, 2e-2), (40.0, 5e-3)):
        T = m*t0
        d = (hc.seg_total(T*1.01, X0, zr, 0.0, zs0, st, 100.0) - hc.seg_total(T, X0, zr, 0.0, zs0, st, 100.0))/(0.01*T)
        r = np.abs(d - K).max()/np.abs(K).max()
        ok = ok and r < lim
        print(f"check_halfplane_cdh: slope at {m:.0f} t0 against the reflected static {r:.1e}")
print("check_halfplane_cdh: " + ("PASS" if ok else "FAIL"))
