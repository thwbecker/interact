#!/usr/bin/env python3
"""demo_lineage: tie the plane-strain plate-over-Maxwell(+gravity)
machinery to the classical results it must reproduce:
 P1 coseismic surface displacement of a 30-deg thrust vs interact's
    Okada kernels (the community form of the Freund & Barnett 1976 /
    Rani & Singh 1992 elastic dip-slip solutions);
 P2 uniform-Maxwell stress relaxation vs the exact correspondence-
    principle factor (Singh & Rosenman 1974; Nur & Mavko 1974 class);
 P3 relaxed limit with gravity vs thin-plate flexure-isostasy
    scaling (Turcotte & Schubert): the k-domain relaxed surface
    response follows 1/(D k^4 + drho g) below kH ~ 1;
 P4 the postseismic surface-uplift transient connecting P1 to the
    relaxed state on the Maxwell timescale (the regime computed by
    Rundle 1978, 1982 and Thatcher & Rundle 1984).
writes inplane_lineage.png"""
import os, subprocess, tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inplane2d import (fields_xz, ve_fields_xz, segment_sources,
                       jumps, solve_k, eval_k)

yr = 3.15576e7
mu1 = lam1 = mu2 = lam2 = 30e9
K2 = lam2 + 2*mu2/3
rho1, rho2, g = 2800.0, 3300.0, 9.81
tau = 10.0*yr
H = 20e3
# 30-degree thrust, top at 2 km, down-dip length 15 km
dip = np.deg2rad(30.0)
x1, z1 = 0.0, 2e3
x2, z2 = x1 + 15e3*np.cos(dip), z1 + 15e3*np.sin(dip)
b = -1.0        # thrust sense (hanging wall up)
src = segment_sources(x1, z1, x2, z2, b, ns=16)
kw = dict(kmin=2e-7, kmax=5e-3, npanel=32, nquad=5)
xs = np.linspace(-40e3, 40e3, 81)
zsurf = np.array([1.0])

fig, axs = plt.subplots(2, 2, figsize=(12.5, 9))

# ---- P1: coseismic vs Okada --------------------------------------------------
F = fields_xz(xs, zsurf, src, lam1, mu1, lam2, mu2, H, **kw)
ib = "../../bin/interact"
okx = oku = None
if os.path.exists(ib):
    W = 7.5e3          # half down-dip width
    zc = 0.5*(z1 + z2)
    xc = 0.5*(x1 + x2)
    with tempfile.TemporaryDirectory() as td:
        open(td + "/geom.in", "w").write(
            f"{xc:.1f} 0 {-zc:.1f} 0 30 2000000 {W:.1f} 0\n")
        open(td + "/bc.in", "w").write("1\n2\n0 1 -1.0\n")
        with open(td + "/oloc.dat", "w") as f:
            for x in xs:
                f.write(f"{x:.1f} 0 -1.0\n")
        subprocess.run([os.path.abspath(ib)], cwd=td,
                       capture_output=True)
        try:
            dd = np.loadtxt(td + "/disp.out")
            okx, oku = dd[:, 0], dd[:, 3:6]
        except OSError:
            pass
ax = axs[0, 0]
ax.plot(xs/1e3, -F["uz"][0], "b-", label="k-solver uz (up)")
ax.plot(xs/1e3, F["ux"][0], "r-", label="k-solver ux")
if oku is not None:
    # interact: z up, slip sign convention opposite (T1a); uz_up = -u
    ax.plot(okx/1e3, -oku[:, 2], "bo", ms=3.5, mfc="none",
            label="interact/Okada uz")
    ax.plot(okx/1e3, -oku[:, 0], "ro", ms=3.5, mfc="none",
            label="interact/Okada ux")
ax.set_title("P1  coseismic 30-deg thrust, surface displacement\n"
             "(elastic anchor: Okada = Freund-Barnett/Rani-Singh)",
             fontsize=10)
ax.set_xlabel("x [km]"); ax.set_ylabel("u [m]")
ax.legend(fontsize=8); ax.axhline(0, color="0.8", lw=0.5)

# ---- P2: uniform-Maxwell relaxation vs exact correspondence ------------------
import sympy as sp
s_, t_ = sp.symbols("s t", positive=True)
mus = mu1*s_/(s_ + 1.0/tau)
lams = (lam1 + 2*mu1/3) - 2*mus/3
D0 = mu1*(lam1 + mu1)/(lam1 + 2*mu1)
fex = sp.inverse_laplace_transform(
    sp.nsimplify(sp.simplify(mus*(lams + mus)/((lams + 2*mus)*s_)/D0)),
    s_, t_)
fl = sp.lambdify(t_, fex, "numpy")
dsrc = segment_sources(0.0, 98e3, 0.0, 102e3, 1.0, ns=8)
pt = (np.array([10e3]), np.array([88e3]))
kwd = dict(kmin=1e-7, kmax=2e-3, npanel=30, nquad=5)
e0 = fields_xz(pt[0], pt[1], dsrc, lam1, mu1, lam1, mu1, 150e3, **kwd)
ts = np.array([0.1, 0.5, 1.5, 3.0])
vv = [float(np.real(ve_fields_xz(pt[0], pt[1], dsrc, lam1, mu1, lam1,
      mu1, 150e3, tau, tf*tau, tau1=tau, M=10, **kwd)["sxz"][0, 0]))
      for tf in ts]
ax = axs[0, 1]
tt = np.linspace(0.01, 4, 200)
ax.plot(tt, fl(tt*tau), "k-", label="exact correspondence factor")
ax.plot(ts, np.array(vv)/float(e0["sxz"][0, 0]), "rs", ms=6,
        mfc="none", label="k-solver + Talbot")
ax.set_title("P2  uniform Maxwell: stress relaxation of a held\n"
             "dislocation (Singh-Rosenman 1974 correspondence class)",
             fontsize=10)
ax.set_xlabel("t / tau_M"); ax.set_ylabel("sigma(t)/sigma(0)")
ax.legend(fontsize=8)

# ---- P3: relaxed limit vs flexure-isostasy, source-free test -----------------
# the ratio of relaxed responses at two gravities cancels the source
# spectrum exactly; thin-plate flexure with FULL isostatic restoring
# predicts r(k) = (D k^4 + rho2 g2)/(D k^4 + rho2 g1): at k -> 0 the
# whole displaced column counts (surface rho1 g plus interface
# (rho2-rho1) g add up to rho2 g), not just the interface contrast.
ks = np.geomspace(3e-7, 3e-4, 25)
mur = 1e-10*mu2; lr = K2 - 2*mur/3
g2 = 0.5*g
uz1, uz2 = [], []
for k in ks:
    J = jumps(k, lam1, mu1, 1.0, (0, 1), (1, 0), 0.0)
    c1 = solve_k(k, lam1, mu1, lr, mur, H, 10e3, J, rho1, rho2, g)
    c2 = solve_k(k, lam1, mu1, lr, mur, H, 10e3, J, rho1, rho2, g2)
    uz1.append(abs(eval_k(0.0, k, lam1, mu1, lr, mur, H, 10e3, c1)[1]))
    uz2.append(abs(eval_k(0.0, k, lam1, mu1, lr, mur, H, 10e3, c2)[1]))
r = np.array(uz2)/np.array(uz1)
nu = lam1/(2*(lam1 + mu1))
E = 2*mu1*(1 + nu)
Dfl = E*H**3/(12*(1 - nu**2))
ax = axs[1, 0]
ax.semilogx(ks, r, "b-", lw=2, label="k-solver: uz(g/2)/uz(g), relaxed")
ax.semilogx(ks, (Dfl*ks**4 + rho2*g)/(Dfl*ks**4 + rho2*g2), "k--",
            label="flexure, restoring rho2 g (full isostasy)")
ax.semilogx(ks, (Dfl*ks**4 + (rho2-rho1)*g)/(Dfl*ks**4 + (rho2-rho1)*g2),
            "r:", label="flexure, restoring (rho2-rho1) g")
kf = (rho2*g/Dfl)**0.25
ax.axvline(kf, color="0.6", lw=0.8)
ax.text(kf*1.15, 1.35, "flexural k", fontsize=8, rotation=90)
ax.set_title("P3  relaxed limit with gravity vs thin-plate\n"
             "flexure-isostasy (Turcotte-Schubert)", fontsize=10)
ax.set_xlabel("k [1/m]"); ax.set_ylabel("response ratio")
ax.legend(fontsize=8)

# ---- P4: postseismic surface uplift transient --------------------------------
ax = axs[1, 1]
xp = np.array([-15e3, 5e3, 25e3])
u0 = fields_xz(xp, zsurf, src, lam1, mu1, lam2, mu2, H,
               rho1=rho1, rho2=rho2, g=g, **kw)
tfs = np.array([0.1, 0.5, 2.0, 8.0, 30.0])
U = [-u0["uz"][0]]
for tf in tfs:
    v = ve_fields_xz(xp, zsurf, src, lam1, mu1, lam2, mu2,
                     H, tau, tf*tau, M=10, rho1=rho1, rho2=rho2,
                     g=g, **kw)
    U.append(-np.real(v["uz"][0]))
U = np.array(U)
for i, x in enumerate(xp):
    ax.semilogx(np.concatenate([[0.03], tfs]), U[:, i], "o-", ms=4,
                label=f"x = {x/1e3:.0f} km")
ax.set_title("P4  postseismic uplift above the thrust\n"
             "(the Rundle 1978/Thatcher-Rundle 1984 regime)",
             fontsize=10)
ax.set_xlabel("t / tau_M"); ax.set_ylabel("uplift [m]")
ax.legend(fontsize=8)
fig.suptitle("inplane_ve_proto vs the classical results "
             "(elastic plate H = 20 km over Maxwell half-space, "
             "30-deg thrust 2-9.5 km depth, b = 1 m, with gravity)",
             fontsize=11)
fig.tight_layout()
fig.savefig("inplane_lineage.png", dpi=150, bbox_inches="tight")
print("wrote inplane_lineage.png")
