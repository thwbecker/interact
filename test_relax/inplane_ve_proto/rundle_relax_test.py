#!/usr/bin/env python3
"""
rundle_relax_test.py: viscoelastic THRUST-FAULT RELAXATION test,
with and without gravity, against Rundle (JGR 1982).

Configuration (his Figures 2 and 3): a 30-degree dipping thrust,
surface-breaking (D = 0), down-dip width W = H, in an elastic
(-gravitational) layer of thickness H over a viscoelastic
(-gravitational) Maxwell half-space; lambda = const (his figure
convention); times in units of tau_a = 2 eta/G = 2 tau_M.  His
purely viscoelastic case is nondimensional (only geometry and
tau_a enter); the gravitational case additionally uses
rho_L = 3.3, rho_H = 3.8 g/cc and mu_L = lam_L = 30 GPa (the latter
from Rundle 1981, J. Phys. Earth 29, 173).

What is tested, quantitatively (misfit versus points digitized from
his figures, estimated reading error +-2 units):

  no gravity  (his Fig. 2): postseismic surface uplift CHANGE at
              5 tau_a and 45 tau_a
  gravity     (his Fig. 3): the same, and the basin-depth
              SUPPRESSION by buoyancy (his -75 -> -45 units), which
              is the actual gravity signal

Beyond his two epochs, the profiles are computed on a dense time
ladder to show the relaxation evolution and, at fixed stations, as
time series (his Fig. 5 style), which is what the earthquake-cycle
kernels ultimately sample.

usage: rundle_relax_test.py [nt] [nx] [outfig]
       nt     epochs per case on a log ladder 0.5-90 tau_a (default 8)
       nx     profile samples (default 121)
       outfig default rundle_relax_test.png
exit code 0 if both cases match the digitized data within tolerance
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inplane2d import fields_xz, ve_fields_xz, segment_sources

nt = int(sys.argv[1]) if len(sys.argv) > 1 else 8
nx = int(sys.argv[2]) if len(sys.argv) > 2 else 121
outfig = sys.argv[3] if len(sys.argv) > 3 else "rundle_relax_test.png"

H = 30e3
mu = lam = 30e9
tau = 10.0*3.15576e7          # tau_M; his tau_a = 2 tau_M
dip = np.deg2rad(30.0)
b = -1.0                      # thrust sense
src = segment_sources(0.0, 1.0, H*np.cos(dip), 0.5*H, b, ns=16)
kw = dict(kmin=1e-7, kmax=9.0/H, npanel=28, nquad=5)
ys = np.linspace(-3*H, 5.4*H, nx)
zs = np.array([1.0])
# station probes live ON the profile grid (one solve set, not two)
probe_y = np.array([-1.0, 0.5, 1.5, 3.0])*H
ys = np.unique(np.concatenate([ys, probe_y]))
ip = np.searchsorted(ys, probe_y)
ta = np.geomspace(0.5, 90.0, nt)          # epochs in tau_a
# his two digitized epochs are always included
ta = np.unique(np.concatenate([ta, [5.0, 45.0]]))

# digitized from Rundle (1982), y in units of H, change x 100
DIG = {
 ("nog", 5.0): np.array([[-3,5],[-2,4],[-1.5,2],[-0.85,0],[-0.5,-6],
   [0,-15],[0.5,-24],[0.65,-26],[1,-23],[1.5,-12],[2,-4],[2.2,0],
   [3,5],[4,6],[5,5]], float),
 ("nog", 45.0): np.array([[-3,5],[-2.5,4],[-2,0],[-1.5,-13],[-1,-30],
   [-0.5,-57],[0,-70],[0.5,-75],[1,-70],[1.5,-52],[2,-30],[2.5,-13],
   [3,-2],[3.5,5],[4,9],[4.5,11],[5,13]], float),
 ("grav", 5.0): np.array([[-3,5],[-2,5],[-1.33,1],[-1,0],[-0.67,-3],
   [0,-13],[0.33,-20],[0.67,-25],[1,-24],[1.33,-17],[1.67,-10],
   [2,-5],[2.33,0],[3,4],[4,5],[5.33,3]], float),
 ("grav", 45.0): np.array([[-3,10],[-2.33,9],[-1.67,3],[-1.43,0],
   [-1,-9],[-0.67,-17],[-0.33,-27],[0,-37],[0.33,-43],[0.67,-45],
   [1,-42],[1.33,-34],[1.67,-24],[2,-13],[2.33,-5],[2.67,0],
   [3.33,8],[4.33,11],[5.33,9]], float)}
CASES = (("nog", "no gravity (Fig. 2)", 0.0, 0.0, 0.0),
         ("grav", "with gravity (Fig. 3)", 3300.0, 3800.0, 9.8))

fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
ok = True
summary = []
for ic, (tag, ttl, r1, r2, g) in enumerate(CASES):
    kwg = dict(kw, rho1=r1, rho2=r2, g=g)
    el = fields_xz(ys, zs, src, lam, mu, lam, mu, H, **kwg)
    prof, pts = {}, {}
    for t in ta:
        ve = ve_fields_xz(ys, zs, src, lam, mu, lam, mu, H, tau,
                          t*2.0*tau, M=12, lam_mode="lambda", **kwg)
        prof[t] = -100.0*(np.real(ve["uz"][0]) - el["uz"][0])
        pts[t] = prof[t][ip]
    ax = axs[0, ic]
    cm = plt.get_cmap("viridis")
    for i, t in enumerate(ta):
        lab = None
        kwl = dict(lw=1.0, color=cm(i/max(len(ta)-1, 1)))
        if t in (5.0, 45.0):
            kwl.update(lw=2.0, color=("tab:blue" if t == 5 else "tab:red"))
            lab = f"{t:.0f} tau_a"
        ax.plot(ys/H, prof[t], label=lab, **kwl)
    for t, col, mk in ((5.0, "tab:blue", "o"), (45.0, "tab:red", "s")):
        d = DIG[(tag, t)]
        ax.plot(d[:, 0], d[:, 1], mk, color=col, mfc="none", ms=5,
                label=f"Rundle 1982, {t:.0f} tau_a")
        p = np.interp(d[:, 0], ys/H, prof[t])
        m = np.abs(d[:, 0]) <= 2.5
        rms = np.sqrt(np.mean((p[m] - d[m, 1])**2))
        pk = np.max(np.abs(d[:, 1]))
        summary.append((tag, t, rms, rms/pk, prof[t].min(), d[:, 1].min()))
        if rms/pk > 0.15:
            ok = False
    ax.axhline(0, color="0.8", lw=0.6)
    ax.set_title(f"30-deg thrust, {ttl}", fontsize=10)
    ax.set_xlabel("y / H"); ax.set_ylabel("uplift change x 100")
    ax.legend(fontsize=7)
    ax = axs[1, ic]
    for k, xp in enumerate(probe_y):
        ax.semilogx(ta, [pts[t][k] for t in ta], "o-", ms=3,
                    label=f"y = {xp/H:+.1f} H")
    ax.axhline(0, color="0.8", lw=0.6)
    ax.set_xlabel("t / tau_a"); ax.set_ylabel("uplift change x 100")
    ax.set_title(f"station histories, {ttl}", fontsize=10)
    ax.legend(fontsize=7)

print("case  epoch   core rms   of peak   model min   Rundle min")
for tag, t, rms, rel, mn, dm in summary:
    print(f"{tag:5s} {t:5.0f}  {rms:9.2f}  {rel:7.1%}  {mn:10.1f}  {dm:10.1f}")
sup = dict((tag, mn) for tag, t, r, rl, mn, dm in summary if t == 45.0)
print(f"gravity basin suppression: {sup['nog']:.1f} -> {sup['grav']:.1f} "
      f"units (Rundle: -75 -> -45)")
if not (0.45 < sup["grav"]/sup["nog"] < 0.75):
    ok = False
    print("  SUPPRESSION RATIO OUT OF RANGE")
fig.suptitle("Viscoelastic thrust-fault relaxation vs Rundle (JGR 1982): "
             "elastic layer over Maxwell half-space", fontsize=12)
fig.tight_layout()
fig.savefig(outfig, dpi=150, bbox_inches="tight")
print(f"wrote {outfig}")
print("rundle_relax_test: " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
