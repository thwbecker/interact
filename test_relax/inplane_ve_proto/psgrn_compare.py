#!/usr/bin/env python3
"""psgrn_compare: cross-validation against Wang's PSGRN/PSCMP (2020),
the modern community standard for layered viscoelastic-gravitational
dislocations.  Case: elastic layer H = 30 km (mu = lam = 30 GPa,
rho 3300) over a Maxwell half-space (mu 30 GPa unrelaxed, rho 3800,
eta = 9.468e17 Pa s -> tau_M = 1.0 yr), 30-deg thrust, W = 30 km
down-dip, surface-breaking, slip 1 m; PSCMP fault 600 km long,
profile at mid-strike (2-D limit).  Gravity on AND off (PSGRN has a
gravity factor switch; the prototype toggles its buoyancy rows).
Compares postseismic CHANGE profiles at 10 and 90 tau_M.
Usage: psgrn_compare.py <psgrn_test_dir> [stage]
  stage compute: run one missing prototype config and cache to npz
  stage plot:    assemble comparison figure + stats from caches"""
import sys, os
import numpy as np
from inplane2d import fields_xz, ve_fields_xz, segment_sources

td = sys.argv[1]
stage = sys.argv[2] if len(sys.argv) > 2 else "compute"
H = 30e3
mu = lam = 30e9
rho1, rho2 = 3300.0, 3800.0
yr = 3.15576e7
tau = 9.468e17/30e9            # Maxwell time [s]
tau_yr = tau/yr                # = 1.0003 yr
dip = np.deg2rad(30.0)
b = -1.0
KM = 111.195
lat = np.linspace(90.0/KM, -160.0/KM, 101)
xs = -lat*KM*1e3               # our x (down-dip horizontal, m)
src = segment_sources(0.0, 1.0, H*np.cos(dip), 0.5*H, b, ns=40)
kw = dict(kmin=1e-7, kmax=9.0/H, npanel=40, nquad=6)
zs = np.array([1.0])

if stage == "compute":
    for gg, tag in ((9.81, "g1"), (0.0, "g0")):
        f = f"psgrn_cache_{tag}.npz"
        if os.path.exists(f):
            continue
        kwg = dict(kw, rho1=rho1, rho2=rho2, g=gg)
        el = fields_xz(xs, zs, src, lam, mu, lam, mu, H, **kwg)
        d = {"el": el["uz"][0]}
        for lab, tf in (("t10", 10.0), ("t90", 90.0)):
            ve = ve_fields_xz(xs, zs, src, lam, mu, lam, mu, H,
                              tau*yr, tf*tau*yr, M=14, **kwg)
            d[lab] = np.real(ve["uz"][0])
        np.savez(f, xs=xs, **d)
        print(f"cached {f}")
        break                   # one config per invocation (time cap)
    else:
        print("all configs cached")
    raise SystemExit(0)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4), sharey=True)
stats = []
for ax, tag, ttl in ((axs[0], "g1", "WITH gravity"),
                     (axs[1], "g0", "gravity OFF")):
    z = np.load(f"psgrn_cache_{tag}.npz")
    for lab, tf, col in (("t10", 10.0, "b"), ("t90", 90.0, "r")):
        # early times from the short-window GF run (fine dt), late
        # times from the long-window run (FFT time sampling)
        sub = tag + ("s" if tf <= 15.0 else "")
        U = np.loadtxt(os.path.join(td, f"out_{sub}", "U_down.dat"),
                       skiprows=1)
        t_d = U[:, 0]/365.25
        # pscmp: rows times, cols receivers; uplift = -U_down; change
        pj = np.array([np.interp(tf*tau_yr, t_d, U[:, 1 + i])
                       for i in range(101)])
        p0 = U[0, 1:102]
        dps = -(pj - p0)*100.0            # uplift change x100, b=1
        mine = -(z[lab] - z["el"])*100.0
        ax.plot(xs/1e3, mine, col + "-",
                label=f"prototype, {tf:.0f} tau_M")
        ax.plot(xs/1e3, dps, col + "o", ms=3, mfc="none",
                label=f"PSGRN/PSCMP, {tf:.0f} tau_M")
        sc = np.max(np.abs(dps))
        rms = np.sqrt(np.mean((mine - dps)**2))
        stats.append((tag, lab, rms, sc))
    ax.set_title(ttl, fontsize=11)
    ax.set_xlabel("x (down-dip horizontal) [km]")
    ax.axhline(0, color="0.7", lw=0.6)
    ax.legend(fontsize=8)
axs[0].set_ylabel("uplift change x 100 [cm/m of slip]")
fig.suptitle("Prototype vs PSGRN/PSCMP: postseismic uplift change, "
             "30-deg thrust, H = 30 km plate over Maxwell half-space",
             fontsize=12)
fig.tight_layout()
fig.savefig("psgrn_compare.png", dpi=150, bbox_inches="tight")
for tag, lab, rms, sc in stats:
    print(f"{tag} {lab}: rms {rms:6.3f} = {rms/sc:6.2%} of peak")
print("wrote psgrn_compare.png")
