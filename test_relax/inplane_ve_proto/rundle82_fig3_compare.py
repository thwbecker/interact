#!/usr/bin/env python3
"""rundle82_fig3_compare: the GRAVITY case, Rundle (1982) Figure 3:
same 30-deg surface-breaking thrust (W = H = 30 km, D = 0) but in an
elastic-GRAVITATIONAL layer over a viscoelastic-gravitational half
space: rhoL = 3.3, rhoH = 3.8 g/cc (Rundle 1982), muL = lamL = 30 GPa
(Rundle 1981), lambda = const.  Rundle's Fig 2 -> Fig 3 shows gravity
cutting the 45 tau_a basin from -75 to -45 and narrowing the profile;
this must emerge from the advected-boundary buoyancy rows."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inplane2d import fields_xz, ve_fields_xz, segment_sources

H = 30e3
mu = lam = 30e9
rho1, rho2, g = 3300.0, 3800.0, 9.8
tau = 10.0*3.15576e7
b = -1.0
dip = np.deg2rad(30.0)
src = segment_sources(0.0, 1.0, H*np.cos(dip), 0.5*H, b, ns=40)
kw = dict(kmin=1e-7, kmax=9.0/H, npanel=40, nquad=6,
          rho1=rho1, rho2=rho2, g=g)
ys = np.linspace(-90e3, 160e3, 101)
zs = np.array([1.0])
el = fields_xz(ys, zs, src, lam, mu, lam, mu, H, **kw)
out = {}
for lab, tfac in (("5ta", 10.0), ("45ta", 90.0)):
    ve = ve_fields_xz(ys, zs, src, lam, mu, lam, mu, H, tau,
                      tfac*tau, M=14, lam_mode="lambda", **kw)
    out[lab] = -100.0*(np.real(ve["uz"][0]) - el["uz"][0])

dig5 = np.array([[-90,5],[-60,5],[-40,1],[-30,0],[-20,-3],[0,-13],
                 [10,-20],[20,-25],[30,-24],[40,-17],[50,-10],
                 [60,-5],[70,0],[90,4],[120,5],[160,3]], float)
dig45 = np.array([[-90,10],[-70,9],[-50,3],[-43,0],[-30,-9],
                  [-20,-17],[-10,-27],[0,-37],[10,-43],[20,-45],
                  [30,-42],[40,-34],[50,-24],[60,-13],[70,-5],
                  [80,0],[100,8],[130,11],[160,9]], float)

fig, ax = plt.subplots(figsize=(9.5, 6))
ax.plot(ys/1e3, out["5ta"], "b-", label="this work (2-D), 5 tau_a")
ax.plot(dig5[:, 0], dig5[:, 1], "bo", mfc="none",
        label="Rundle 1982 Fig. 3 (digitized)")
ax.plot(ys/1e3, out["45ta"], "r-", label="this work (2-D), 45 tau_a")
ax.plot(dig45[:, 0], dig45[:, 1], "rs", mfc="none",
        label="Rundle 1982 Fig. 3 (digitized)")
ax.axhline(0, color="0.7", lw=0.6)
ax.set_xlabel("y [km]"); ax.set_ylabel("change in (uz/U) x 100")
ax.set_title("WITH GRAVITY: postseismic change, 30-deg thrust, "
             "H = 30 km\n(Rundle 1982 Fig. 3; rhoL 3.3, rhoH 3.8, "
             "mu = lam = 30 GPa)", fontsize=11)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("rundle82_fig3_compare.png", dpi=150, bbox_inches="tight")
for lab, dig in (("5ta", dig5), ("45ta", dig45)):
    p = np.interp(dig[:, 0], ys/1e3, out[lab])
    mx = np.max(np.abs(dig[:, 1]))
    mc = np.abs(dig[:, 0]) <= 75.0
    rmsc = np.sqrt(np.mean((p[mc] - dig[mc, 1])**2))
    print(f"{lab}: core rms {rmsc:5.2f} ({rmsc/mx:5.1%} of peak); "
          f"model min {out[lab].min():6.1f} vs Rundle "
          f"{dig[:, 1].min():6.1f}")
print("wrote rundle82_fig3_compare.png")
