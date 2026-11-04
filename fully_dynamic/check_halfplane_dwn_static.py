"""Static check of the DWN half-plane kernel (halfplane_dwn.py, fd_plan.md M10): the transfer
function at a small frequency against Crouch and Starfield's half-plane solution (interact's
tdd_coeff through fs_static2d.py; TDD_LIB names the library) for horizontal and dipping segments.
Prints the relative difference per case and PASS if all are below 1e-5."""
import sys, os
import numpy as np
os.environ.setdefault('TDD_LIB', '/tmp/tdd/libtdd.so')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import halfplane_dwn as hd
import fs_static2d as fs
al, be, rho = 6000., 3464., 2670.
a = 250.
fs.MU = rho*be**2
fs.NU = (al**2 - 2*be**2)/(2*(al**2 - be**2))
worst = 0.0
for (xs, zs, xr, zr, ang) in ((0., 5000., 800., 4400., 0.), (0., 2000., 1500., 800., 0.), (0., 1000., 300., 200., 0.),
                              (0., 1500., -700., 900., 30.), (0., 800., 400., 300., -50.)):
    s = (np.cos(np.radians(ang)), np.sin(np.radians(ang)))          # tangent in (x, z), z down
    w = 2*np.pi*1e-4 + 1j*1e-5
    S = hd.spectrum(xr, zr, xs, zs, s, a, w, al, be, rho, L=4e6, M=40000, part="total").real
    st, op = fs.tdd(xr, -zr, xs, -zs, a, s[0], -s[1], 1)            # tdd: y up, D = u(-) - u(+)
    ref = -np.array([st[0], -st[1], st[2]])
    rel = np.abs(S - ref).max()/np.abs(ref).max()
    worst = max(worst, rel)
    print(f"source depth {zs:g} m, dip {ang:g}: DWN {S} Crouch-Starfield {ref} relative {rel:.1e}")
print(f"check_halfplane_dwn_static: largest relative difference {worst:.1e}: {'PASS' if worst < 1e-5 else 'FAIL'}")
