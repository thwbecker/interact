#!/usr/bin/env python3
"""rundle82_compare: quantitative comparison with Rundle (JGR 1982),
Figure 2: 30-deg dipping thrust, surface-breaking (D = 0), down-dip
width W = H, in an elastic layer over a Maxwell half-space, NO
gravity, lambda = const.  Rundle plots the CHANGE in surface uz
(x100, per unit slip) after 5 tau_a and 45 tau_a, tau_a = 2 eta/G =
2 tau_M.  The change field is controlled by the substrate and is
smooth on scale H, so both epochs and the elastic reference share a
k-grid truncated at kH ~ 9 with cancelling truncation errors.

Rundle's fault is 3-D with half-length L = 10H/3 (2L = 6.7H); we
compute the 2-D (infinitely long) limit and quote the finite-length
elastic correction at mid-strike from interact's Okada kernels for
context.  Digitized points read from Figure 2 (estimated +-2 units).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inplane2d import fields_xz, ve_fields_xz, segment_sources

H = 30e3
mu = lam = 30e9
tau = 10.0*3.15576e7         # arbitrary; times are in units of tau
b = -1.0                     # thrust sense (hanging wall up, +x dip)
dip = np.deg2rad(30.0)
src = segment_sources(0.0, 1.0, H*np.cos(dip), 0.5*H, b, ns=40)
kw = dict(kmin=1e-7/1.0, kmax=9.0/H, npanel=40, nquad=6)
ys = np.linspace(-3*H, 5*H, 81)
zs = np.array([1.0])

el = fields_xz(ys, zs, src, lam, mu, lam, mu, H, **kw)
out = {}
for lab, tfac in (("5ta", 10.0), ("45ta", 90.0)):
    ve = ve_fields_xz(ys, zs, src, lam, mu, lam, mu, H, tau,
                      tfac*tau, M=14, lam_mode="lambda", **kw)
    # uz positive DOWN in the solver; Rundle plots up, per unit slip
    out[lab] = -100.0*(np.real(ve["uz"][0]) - el["uz"][0])

# digitized from Rundle (1982) Fig. 2 (y in units of H, d_uz x 100)
dig5 = np.array([[-3.0, 5], [-2.0, 4], [-1.5, 2], [-0.85, 0],
                 [-0.5, -6], [0.0, -15], [0.5, -24], [0.65, -26],
                 [1.0, -23], [1.5, -12], [2.0, -4], [2.2, 0],
                 [3.0, 5], [4.0, 6], [5.0, 5]], dtype=float)
dig45 = np.array([[-3.0, 5], [-2.5, 4], [-2.0, 0], [-1.5, -13],
                  [-1.0, -30], [-0.5, -57], [0.0, -70], [0.5, -75],
                  [1.0, -70], [1.5, -52], [2.0, -30], [2.5, -13],
                  [3.0, -2], [3.5, 5], [4.0, 9], [4.5, 11],
                  [5.0, 13]], dtype=float)

fig, ax = plt.subplots(figsize=(9.5, 6))
ax.plot(ys/H, out["5ta"], "b-", label="this work (2-D), t = 5 tau_a")
ax.plot(dig5[:, 0], dig5[:, 1], "bo", mfc="none",
        label="Rundle 1982 Fig. 2, 5 tau_a (digitized)")
ax.plot(ys/H, out["45ta"], "r-", label="this work (2-D), t = 45 tau_a")
ax.plot(dig45[:, 0], dig45[:, 1], "rs", mfc="none",
        label="Rundle 1982 Fig. 2, 45 tau_a (digitized)")
ax.axhline(0, color="0.7", lw=0.6)
ax.set_xlabel("y / H"); ax.set_ylabel("change in (uz/U) x 100")
ax.set_title("Postseismic surface-displacement change: 30-deg "
             "surface-breaking thrust,\nelastic layer over Maxwell "
             "half-space, no gravity (Rundle 1982 Fig. 2 setup)",
             fontsize=11)
ax.legend(fontsize=9)
ax.text(0.02, 0.03, "Rundle fault is 3-D (2L = 6.7H): Okada elastic\n"
        "finite-length factors are 0.8-1.1 over the basin and\n"
        "collapse beyond ~2.5H; expect the 2-D curves to\n"
        "overshoot mildly at 45 tau_a and far-field lobes to differ",
        transform=ax.transAxes, fontsize=8, va="bottom")
fig.tight_layout()
fig.savefig("rundle82_fig2_compare.png", dpi=150, bbox_inches="tight")

# statistics: full range, and the |y| <= 2.5H core where the 2-D
# (infinite-length) limit is meaningful.  Rundle's fault has
# 2L = 6.7H; interact-Okada elastic factors uz(6.7H)/uz(inf) are
# 0.8-1.1 over the basin and collapse (even changing sign) beyond
# ~2.5H, so far-field lobes are intrinsically 3-D.
for lab, dig in (("5ta", dig5), ("45ta", dig45)):
    p = np.interp(dig[:, 0], ys/H, out[lab])
    mx = np.max(np.abs(dig[:, 1]))
    rms = np.sqrt(np.mean((p - dig[:, 1])**2))
    mc = np.abs(dig[:, 0]) <= 2.5
    rmsc = np.sqrt(np.mean((p[mc] - dig[mc, 1])**2))
    print(f"{lab}: rms {rms:5.2f} ({rms/mx:5.1%} of peak) full; "
          f"{rmsc:5.2f} ({rmsc/mx:5.1%}) core |y|<=2.5H; model min "
          f"{out[lab].min():6.1f} vs Rundle {dig[:, 1].min():6.1f}")
print("wrote rundle82_fig2_compare.png")
