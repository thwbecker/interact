"""Compare the dense tables of bin/test_fd_surf (fault plus virtual free-surface elements,
in-plane; fd_plan.md M7 step 2) with densesys.py of the prototype.
Usage: python3 test_fd_surf.py tables.txt [nf ns L]"""
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import densesys as ds
import inplane_open as io

# the C tables average the opening kernel over t -+ 1.7e-3 dt (densesys uses 1e-4 dt, which is
# ill-conditioned when a front is tangent to the source line at a sample time; fd_tables.c)
_seg_open = io.segment_ramp_stress_open


def _seg_open_shift(xr, yr, t, a, alpha, beta, rho, ngauss=12, near=0.5, tshift=0.0):
    dt0 = 0.5*250.0/3464.0
    return _seg_open(xr, yr, t, a, alpha, beta, rho, ngauss, near, 1.7e-3*dt0 if tshift > 0 else 0.0)


io.segment_ramp_stress_open = _seg_open_shift
# the prototype's static opening field uses fixed MU, NU; set them to the material's
import fs_static2d as fs
fs.MU = 2670.0*3464.0**2
fs.NU = (6000.0**2 - 2*3464.0**2)/(2*(6000.0**2 - 3464.0**2))

fn = sys.argv[1]
nf = int(sys.argv[2]) if len(sys.argv) > 2 else 4
ns = int(sys.argv[3]) if len(sys.argv) > 3 else 8
L = int(sys.argv[4]) if len(sys.argv) > 4 else 30
alpha, beta, rho = 6000.0, 3464.0, 2670.0
dt = 0.5*250.0/beta
dip = np.radians(30.0)
sf = np.array([np.cos(dip), -np.sin(dip)])
X = [np.array([0.0, -500.0]) + (j + 0.5)*500.0*sf for j in range(nf)]
X += [np.array([-1000.0 + (j + 0.5)*250.0, 0.0]) for j in range(ns)]
s = [sf]*nf + [np.array([1.0, 0.0])]*ns
a = [250.0]*nf + [125.0]*ns
kind = [0]*nf + [1]*ns
d = ds.DenseSys(np.array(X), np.array(s), np.array(a), np.array(kind), dt, alpha, beta, rho, L)


def col(seg, j):
    return d.idx[(j, 0)] if seg == 0 else (d.idx[(nf + j, 0)] if seg == 1 else d.idx[(nf + j, 1)])


def row(seg, i):
    return d.idx[(i, 0)] if seg == 0 else (d.idx[(nf + i, 0)] if seg == 1 else d.idx[(nf + i, 1)])


def elem(seg, i):
    return i if seg == 0 else nf + i


worst = 0.0
scale = 0.0
worstk = 0.0
scalek = 0.0
for line in open(fn):
    w = line.split()
    if w[0] == "K":
        A, B, i, j = map(int, w[1:5])
        ks, kn = float(w[5]), float(w[6])
        rs, cs = row(B, i), col(A, j)
        worstk = max(worstk, abs(ks - d.Ks[rs, cs]))
        scalek = max(scalek, abs(d.Ks[rs, cs]))
        if B == 0:
            worstk = max(worstk, abs(kn - d.Kn[elem(B, i), cs]))
    elif w[0] == "W":
        A, B, l, i, j = map(int, w[1:6])
        ws, wn = float(w[6]), float(w[7])
        rs, cs = row(B, i), col(A, j)
        worst = max(worst, abs(ws - d.W[l - 1, rs, cs]))
        scale = max(scale, abs(d.W[l - 1, rs, cs]))
        if B == 0:
            worst = max(worst, abs(wn - d.Wn[l - 1, elem(B, i), cs]))
ok = worst <= 1e-7*scale and worstk <= 1e-7*scalek
print(f"test_fd_surf.py: fault {nf} + surface {ns} elements, {L} lags: tables max diff {worst:.3e} of {scale:.3e}, "
      f"static max diff {worstk:.3e} of {scalek:.3e}: {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
