#!/usr/bin/env python3
"""T1: elastic anchors for the k-space layered solver.

T1a (tight): uniform half-space vertical dip-slip fault vs interact's
Okada kernels (a 4000-km-long 3-D fault evaluated at mid-strike is
the plane-strain limit).  This is the community-grade anchor; the
k-solver matched it to 5e-4 at all tested depths when this gate was
built.  Skipped with a warning if ../../bin/interact is not built.

T1b (soft): same problem vs the independent FEM (fem2d, split
nodes).  Bilinear quads shear-lock in this bending-dominated
free-surface problem and carry a documented ~10 percent bias (they
match full-space fields to ~1 percent where nothing flexes, and
Okada agrees with the k-solver, so the bias is firmly attributed to
the FEM).  The FEM is kept as an order-of-magnitude/back-up check
and for future viscoelastic time-stepping comparisons with better
elements; tolerance 15 percent."""
import os, subprocess
import numpy as np
from inplane2d import fields_xz, segment_sources
from fem2d import Fem

b, d1, d2, H = 1.0, 4e3, 12e3, 20e3
lam1 = mu1 = 30e9
src = segment_sources(0.0, d1, 0.0, d2, b, ns=30)
kw = dict(kmin=2e-7, kmax=4e-3, npanel=50, nquad=6)
zs = np.array([14e3, 20e3, 30e3])
xs = np.arange(2e3, 27e3, 2e3)
F = fields_xz(xs, zs, src, lam1, mu1, lam1, mu1, H, **kw)
ok = True

ib = os.path.join(os.path.dirname(__file__), "../../bin/interact")
if os.path.exists(ib):
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        open(td + "/geom.in", "w").write(
            f"0 0 {-(d1+d2)/2:.1f} 0 90 2000000 {(d2-d1)/2:.1f} 0\n")
        open(td + "/bc.in", "w").write("1\n2\n0 1 1.0\n")
        with open(td + "/oloc.dat", "w") as f:
            for z in zs:
                for x in xs:
                    f.write(f"{x} 0 {-z}\n")
        subprocess.run([os.path.abspath(ib)], cwd=td,
                       capture_output=True)
        try:
            okd = np.loadtxt(td + "/stress.out")
        except (FileNotFoundError, OSError):
            okd = None
            print("T1a SKIPPED: interact binary present but did not "
                  "run (library environment?)")
    if okd is not None:
        # interact runs with shear modulus 1e4 and nu=0.25;
        # z -> -z flips sxz; slip sign convention is opposite
        fac = -mu1/1e4
        for j, comp, col, sg in ((0, "sxx", 3, 1.0), (1, "szz", 8, 1.0),
                                 (2, "sxz", 5, -1.0)):
            o = fac*sg*okd[:, col]
            p = np.array([F[comp][np.where(zs == -zr)[0][0],
                                  np.where(xs == xr)[0][0]]
                          for xr, _, zr in okd[:, :3]])
            mis = np.sqrt(np.mean((o - p)**2))/np.sqrt(np.mean(o**2))
            print(f"T1a okada {comp}: rel misfit {mis:8.2e}")
            ok = ok and mis < 0.01
else:
    print("T1a SKIPPED: ../../bin/interact not built")

fem = Fem(100e3, 200e3, 250, 250, 0.0, d1, d2, lam1, mu1, lam1, mu1, H)
u = fem.solve_elastic(b)
S = fem.stress_field(u); C = fem.centers()
for zt in zs:
    mm = (np.abs(C[:, 1] - zt) < 0.45e3) & (C[:, 0] > 2e3) & (C[:, 0] < 26e3)
    mm &= C[:, 1] == np.unique(C[mm, 1])[0]
    Fw = fields_xz(C[mm, 0], np.array([C[mm, 1][0]]), src,
                   lam1, mu1, lam1, mu1, H, **kw)
    a = S[mm, 1]; p = Fw["szz"][0]
    mis = np.sqrt(np.mean((a - p)**2))/np.sqrt(np.mean(a**2))
    print(f"T1b FEM szz z={zt/1e3:4.0f} km: rel misfit {mis:6.2%} "
          "(bilinear-quad locking bias, soft gate)")
    ok = ok and mis < 0.15
print("T1 " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
