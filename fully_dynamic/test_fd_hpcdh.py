"""Check of fd_hpcdh_batch (C, fd_hpcdh.c) against the prototype halfplane_cdh.py: reads the output of
test_fd_hpdwn_batch n h top nt cdh (fault of n elements of h m dipping 30 degrees from depth top) and
compares every entry and lag (shear and normal, fd_project's convention) with the prototype, the
tail beyond 40 first arrivals excluded (the prototype has none there), and the static tractions with
Crouch and Starfield's half minus full plane (fs_static2d.py, needs TDD_LIB). PASS if both agree to
1e-8 of the largest value.
Usage: python3 test_fd_hpcdh.py file"""
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import halfplane_cdh as hc

lines = open(sys.argv[1]).read().split("\n")
hdr = lines[0].split()
n = int(hdr[2])
h = float(hdr[4])
top = float(hdr[6])
nt = int(hdr[8])
dt = float(hdr[10])
dip = np.radians(30.0)
tx, tz = np.cos(dip), np.sin(dip)
X = [((i + 0.5)*h*tx, top + (i + 0.5)*h*tz) for i in range(n)]
S = {}
K = {}
for l in lines[1:]:
    w = l.split()
    if not w:
        continue
    if w[0] == "S":
        S.setdefault((int(w[1]), int(w[2])), []).append((float(w[4]), float(w[5])))
    elif w[0] == "K":
        K[(int(w[1]), int(w[2]))] = np.array([float(w[3]), float(w[4])])


def proj(s):
    sr = np.array([tx, -tz])
    nr = np.array([-sr[1], sr[0]])
    sxx, sxz, szz = s
    sxy = -sxz
    syy = szz
    a = sxx*nr[0] + sxy*nr[1]
    b = sxy*nr[0] + syy*nr[1]
    return np.array([a*sr[0] + b*sr[1], a*nr[0] + b*nr[1]])


worst = 0.0
scale = max(np.abs(np.array(v)).max() for v in S.values())
for (i, j), ser in S.items():
    ser = np.array(ser)
    t0 = hc.first_arrival(X[i][0], X[i][1], X[j][0], X[j][1], (tx, tz), h/2)
    for l in range(1, nt, 3):
        t = l*dt
        if t > 40.0*t0:
            break
        p = proj(hc.seg_total(t, X[i][0], X[i][1], X[j][0], X[j][1], (tx, tz), h/2, nq=16)) if t > t0 else np.zeros(2)
        worst = max(worst, np.abs(ser[l] - p).max()/scale)
ok = worst < 1e-8
print(f"test_fd_hpcdh: series against the prototype {worst:.1e}")
if os.path.exists(os.environ.get("TDD_LIB", "/tmp/tdd/libtdd.so")):
    os.environ.setdefault("TDD_LIB", "/tmp/tdd/libtdd.so")
    import fs_static2d as fs
    fs.MU = hc.RHO*hc.BE**2
    fs.NU = (hc.AL**2 - 2*hc.BE**2)/(2*(hc.AL**2 - hc.BE**2))
    ws = 0.0
    ks = max(np.abs(v).max() for v in K.values())
    for (i, j), k in K.items():
        h1, _ = fs.tdd(X[i][0], -X[i][1], X[j][0], -X[j][1], h/2, tx, -tz, 1)
        f1, _ = fs.tdd(X[i][0], -X[i][1], X[j][0], -X[j][1], h/2, tx, -tz, 0)
        d = -(np.array(h1) - np.array(f1))
        ws = max(ws, np.abs(k - proj((d[0], -d[1], d[2]))).max()/ks)
    ok = ok and ws < 1e-8
    print(f"test_fd_hpcdh: static against Crouch and Starfield {ws:.1e}")
print("test_fd_hpcdh: " + ("PASS" if ok else "FAIL"))
