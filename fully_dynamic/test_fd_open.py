"""Compare the C opening kernel (bin/test_fd_open) with inplane_open.py at 200 random
receivers (off the segment, on it, near the source line) and times.
Usage: python3 test_fd_open.py path/to/test_fd_open"""
import sys
import os
import subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inplane_open as io

exe = sys.argv[1] if len(sys.argv) > 1 else "../bin/test_fd_open"
rng = np.random.default_rng(1)
a = 250.0
pts = []
for _ in range(200):
    kind = rng.integers(4)
    if kind == 0:
        xr, yr = rng.uniform(-3000, 3000), rng.uniform(-3000, 3000)
    elif kind == 1:
        xr, yr = rng.uniform(-a, a), rng.uniform(-0.4*a, 0.4*a)
    elif kind == 2:
        xr, yr = rng.uniform(-2000, 2000), rng.uniform(-400, 400)
    else:
        xr, yr = rng.uniform(-500, 500), rng.uniform(300, 600)
    pts.append((xr, yr, rng.uniform(0.01, 2.0)))
pts = np.array(pts)
inp = "\n".join(f"{x:.17e} {y:.17e} {t:.17e}" for x, y, t in pts)
out = subprocess.run([exe], input=inp, capture_output=True, text=True).stdout
C = np.array([[float(v) for v in line.split()] for line in out.strip().split("\n")])
P = np.array(io.segment_ramp_stress_open(pts[:, 0], pts[:, 1], pts[:, 2], a, 6000.0, 3464.0, 2670.0)).T
rel = np.abs(C - P).max()/np.abs(P).max()
print(f"test_fd_open.py: 200 receivers, max |C - prototype| relative {rel:.2e}: {'PASS' if rel < 1e-8 else 'FAIL'}")
