"""Compare the C reflected half-plane stresses (bin/test_fd_hpdwn, fd_hpdwn.c) with the prototype
(halfplane_dwn.py, part "refl") for a dipping and a horizontal source (fd_plan.md M10).
Usage: python3 test_fd_hpdwn.py path/to/test_fd_hpdwn"""
import sys
import os
import subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import halfplane_dwn as hd

exe = sys.argv[1] if len(sys.argv) > 1 else "../bin/test_fd_hpdwn"
al, be, rho = 6000., 3464., 2670.
cases = [(0., 1500., np.cos(np.radians(30)), np.sin(np.radians(30)), 250., 600., 900.),
         (0., 800., 1., 0., 125., -400., 300.)]
ts = np.array([0.3, 0.6, 1.0, 1.5])
worst = 0.0
for (xs, zs, tx, tz, a, xr, zr) in cases:
    T, N, wIT = 4.0, 512, 18.0
    L = 2*al*T + 4*(abs(xr - xs) + a)
    M = int(L/(0.25*a)) + 1
    P = hd.ramp_stress(xr, zr, xs, zs, (tx, tz), a, ts, al, be, rho, part="refl", T=T, N=N, L=L, M=M, wI_T=wIT)
    inp = f"{xs} {zs} {tx} {tz} {a} {xr} {zr} {T} {N} {wIT} {L} {M} {len(ts)} " + " ".join(str(x) for x in ts)
    out = subprocess.run([exe], input=inp, capture_output=True, text=True).stdout.split("\n")
    C = np.array([[float(v) for v in line.split()] for line in out if line.strip()])
    worst = max(worst, np.abs(C - P).max()/np.abs(P).max())
print(f"test_fd_hpdwn.py: reflected half-plane stresses, C against the prototype: largest relative difference {worst:.1e}: {'PASS' if worst < 1e-10 else 'FAIL'}")
