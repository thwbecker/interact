#!/usr/bin/env python3
"""
plot_xsection.py: reconstruct and plot the OFF-FAULT stress field
sigma_xy(x, z, t) (and sigma_yz) in x-depth cross section for a
run_xsection output directory: an antiplane rate-and-state fault
(0 <= z <= D, x = 0) in an elastic plate (0 <= z <= H, modulus G) over
a Maxwell half-space of the same elastic modulus (equal rigidities;
tM_yr = 0 means uniform elastic half-space).

METHOD.  The stress anywhere is the hereditary response to the fault's
slip-deficit history q_j(t) = v_j - vpl, read from the -field_slip
frames.  Per image-reflection order n the time dependence is carried by
Erlang memory chains built ONCE over the frames (S_0 = deficit D_j,
dS_m/dt = b (S_{m-1} - S_m); b = 1/(2 tM) for equal rigidities):

  PLATE receivers (z < H): sigma_1 = K_0 D + sum_{n>=1} K_n S_n, where
  K_0 is the elastic half-space dipole field (direct + free-surface
  mirror) and K_n the same field of the image pair at depths 2nH +- d
  (each with its own surface mirror); S_n carries the validated Erlang
  CDF weight A_n(bt) = L^-1[Gamma(s)^n / s], Gamma(s) = b/(s+b).

  SUBSTRATE receivers (z > H): the transmitted field.  With T(s) =
  1 + Gamma(s) and the Maxwell modulus mu_2(s) = G s/(s+2b), the
  per-order stress weight is W_n(s) = (1/s) T Gamma^n mu_2(s)/G =
  b^n/(s+b)^{n+1}, i.e. the POISSON mass p_n(bt) = (bt)^n e^-bt/n! =
  A_n - A_{n+1}, obtained for free from the same chains (P_n = S_n -
  S_{n+1}); the spatial kernel is the plain elastic dipole field of
  the image pair at depths +-d - 2nH (rays reaching the substrate
  after n internal reflections, with and without one surface bounce;
  the n = 0 kernel is identically the direct + mirror field, so at
  t = 0, where p_n = delta_n0, the whole medium is the uniform elastic
  half-space, and as t -> infinity substrate stresses relax to zero).

BUILT-IN CONSISTENCY CHECKS (run every time; abort on failure):
(1) the on-fault sigma_xy of the elastic kernel matches the
independent closed-form fault kernel; (2) the free surface is traction
free; (3) the interface traction sigma_yz is continuous across z = H
at EVERY snapshot time, which jointly tests the image positions and
the time weights of both media (sigma_xy may jump once the substrate
relaxes; sigma_yz must not).

WHAT IS PLOTTED.  Snapshots through the LAST full recurrence of the
run: the change of stress relative to the immediately post-event state
(this is "how stress evolves over the cycle" and is independent of the
reconstruction's virgin start), at phases post, 0.25, 0.50, 0.75, and
pre(next event).  A separate figure shows the coseismic change of the
cycle-closing event.  Fields are symmetric (sigma_xy) / antisymmetric
(sigma_yz) in x; only x >= 0 is shown.

usage: plot_xsection.py rundir [nx nz xfac zfac sat]
       grid defaults 90 x 90 over (0, xfac*H] x (0, zfac*H],
       xfac = 1.5, zfac = 1.75 (extra substrate depth in view); sat is
       the saturation quantile of the color scale (default 0.99: the
       top 1 percent of |values|, i.e. the near-fault extremes, are
       clipped so that the far-field differences read better;
       colorbars show arrows where saturated)
"""
import sys, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib import cm

NLEV = 21          # target discrete color subdivisions (nice-number levels)
from matplotlib.ticker import MaxNLocator
def discrete(vmax):
    """symmetric discrete levels with round values, about NLEV bins"""
    lev = MaxNLocator(nbins=NLEV, symmetric=True).tick_values(-vmax, vmax)
    return lev, BoundaryNorm(lev, len(lev)-1), plt.get_cmap("RdBu_r", len(lev)-1)

rundir = sys.argv[1]
nx  = int(sys.argv[2]) if len(sys.argv) > 2 else 90
nz  = int(sys.argv[3]) if len(sys.argv) > 3 else 90
xfac = float(sys.argv[4]) if len(sys.argv) > 4 else 1.5
zfac = float(sys.argv[5]) if len(sys.argv) > 5 else 1.75
sat  = float(sys.argv[6]) if len(sys.argv) > 6 else 0.99
os.chdir(rundir)

par = dict((l.split()[0], l.split()[1]) for l in open("xsect_params.txt") if l.split())
H   = float(par["H_km"])*1e3
G   = float(par["G_Pa"])
vpl = float(par["vpl_ms"])
tM  = float(par["tM_yr"])
Df  = float(par["D_km"])
label = par["label"]
yr  = 3.15576e7
elastic = tM <= 0.0
b   = 0.0 if elastic else 1.0/(2.0*tM*yr)
pref = -G/(2.0*np.pi)   # global dipole orientation, pinned by self-test 1

# fault cells
geo = np.loadtxt("rsf_geom.g000.dat", comments="#")
zc_f = -geo[:, 4]; hl = geo[:, 8]
ncell = zc_f.size
ze1, ze2 = zc_f - hl, zc_f + hl

# frames
T = np.loadtxt("rsf_vel.times", comments="#")
t_fr = T[:, 3]
fl = sorted(glob.glob("tmp_rsf/rsf_slip.g000.*.bin"))
nfr = min(len(fl), t_fr.size)
slip = np.empty((nfr, ncell))
for i in range(nfr):
    slip[i] = np.fromfile(fl[i], dtype="<f4").reshape(-1, 3)[:, 2]
t_fr = t_fr[:nfr]
D_fr = slip - vpl*t_fr[:, None]

# image order count
if elastic:
    N = 0
else:
    tcyc = 0.2*(t_fr[-1] - t_fr[0])          # ~ one recurrence
    N = int(min(48, max(8, np.ceil(3.0*b*tcyc) + 4)))

# Erlang chains (one exponential step per frame interval, averaged forcing)
S = np.zeros((N+2, nfr, ncell))
S[0] = D_fr
for i in range(1, nfr):
    if N == 0:
        break
    e = np.exp(-b*(t_fr[i] - t_fr[i-1]))
    for m in range(1, N+2):
        fbar = 0.5*(S[m-1, i-1] + S[m-1, i])
        S[m, i] = S[m, i-1]*e + fbar*(1.0 - e)

# ---------------------------------------------------------------------------
def line_stress(x, z, zc, sgn):
    R2 = x*x + (z - zc)**2
    return -sgn*(z - zc)/R2, sgn*x/R2

def cell_stress(x, z, z1, z2, mirror=True):
    """dipole = line(z2, +) + line(z1, -); free-surface mirror pair at
    (-z2, -), (-z1, +)."""
    s1 = line_stress(x, z, z2, +1.0)
    s2 = line_stress(x, z, z1, -1.0)
    sxy, syz = s1[0] + s2[0], s1[1] + s2[1]
    if mirror:
        s3 = line_stress(x, z, -z2, -1.0)
        s4 = line_stress(x, z, -z1, +1.0)
        sxy += s3[0] + s4[0]; syz += s3[1] + s4[1]
    return sxy, syz

# grid
xg = np.linspace(0.015*H, xfac*H, nx)
zg = np.linspace(0.01*H, zfac*H, nz)
XX, ZZ = np.meshgrid(xg, zg)
mpl_ = zg < H; msb_ = ~mpl_
XXp, ZZp = XX[mpl_], ZZ[mpl_]
XXs, ZZs = XX[msb_], ZZ[msb_]

def accum(xx, zz, weights_by_order, sub):
    """stress (sxy, syz) at points (xx, zz): sum over orders n of
    weight-vector[n] (ncell) applied to the order-n kernels"""
    sxy = np.zeros(xx.shape); syz = np.zeros(xx.shape)
    for n, w in weights_by_order:
        if not np.any(w):
            continue
        for j in range(ncell):
            if w[j] == 0.0:
                continue
            if not sub:
                if n == 0:
                    a = cell_stress(xx, zz, ze1[j], ze2[j])
                else:
                    o = 2.0*n*H
                    a1 = cell_stress(xx, zz, o + ze1[j], o + ze2[j])
                    b1 = cell_stress(xx, zz, o - ze2[j], o - ze1[j])
                    a = (a1[0] + b1[0], a1[1] + b1[1])
            else:
                o = -2.0*n*H
                a1 = cell_stress(xx, zz, o + ze1[j], o + ze2[j], mirror=False)
                b1 = cell_stress(xx, zz, o - ze2[j], o - ze1[j], mirror=False)
                a = (a1[0] + b1[0], a1[1] + b1[1])
            sxy += pref*w[j]*a[0]; syz += pref*w[j]*a[1]
    return sxy, syz

def field_at(ifr):
    sxy = np.zeros((nz, nx)); syz = np.zeros((nz, nx))
    wpl = [(0, S[0, ifr])] + [(n, S[n, ifr]) for n in range(1, N+1)]
    a, c = accum(XXp, ZZp, wpl, sub=False)
    sxy[mpl_] = a; syz[mpl_] = c
    if elastic:
        a, c = accum(XXs, ZZs, [(0, S[0, ifr])], sub=True)
    else:
        wsb = [(n, S[n, ifr] - S[n+1, ifr]) for n in range(0, N+1)]
        a, c = accum(XXs, ZZs, wsb, sub=True)
    sxy[msb_] = a; syz[msb_] = c
    return sxy, syz

# ---------------------------------------------------------------------------
# self-test 1: on-fault sigma_xy vs the closed-form fault kernel
jt = ncell//3
zt = np.array([[zc_f[2*ncell//3]]]); xt = np.array([[1.0]])
sxy_t, _ = cell_stress(xt, zt, ze1[jt], ze2[jt])
sxy_t = (pref/(G/(2.0*np.pi)))*np.asarray(sxy_t)   # apply the global sign
def F0(z, zp): return 1.0/(z - zp) - 1.0/(z + zp)
ref = F0(zt[0, 0], ze2[jt]) - F0(zt[0, 0], ze1[jt])
if abs(sxy_t[0, 0] - ref)/abs(ref) > 1e-3:
    sys.exit(f"plot_xsection: SELF-TEST 1 FAILED "
             f"(on-fault kernel: {sxy_t[0,0]:.6e} vs {ref:.6e})")

# self-test 2: free surface traction-free
xs = np.linspace(0.1*H, H, 20); zs = np.full(20, 1e-3)
a = cell_stress(xs, zs, ze1[ncell//2], ze2[ncell//2])
if np.max(np.abs(a[1])) > 1e-6*np.max(np.abs(a[0])):
    sys.exit("plot_xsection: SELF-TEST 2 FAILED (free surface traction)")

# ---------------------------------------------------------------------------
# snapshots: last full cycle
ev = [l.split() for l in open("rsf_events.dat") if not l.startswith("#")]
on  = np.array([float(r[0]) for r in ev if int(r[2]) == 1])
off = np.array([float(r[0]) for r in ev if int(r[2]) == -1])
# cluster the threshold crossings into events (gaps > 1 yr)
keep = np.concatenate([[True], np.diff(on) > 1.0*yr])
on = on[keep]
# per clustered event, the arrest is the LAST -1 crossing before the
# next onset
offc = []
for i, t0 in enumerate(on):
    t1 = on[i+1] if i+1 < on.size else np.inf
    sel = off[(off > t0) & (off < t1)]
    offc.append(sel[-1] if sel.size else t0)
off = np.array(offc)
if on.size < 2:
    sys.exit("plot_xsection: fewer than 2 clustered events; run longer")
if on.size == 2:
    print("plot_xsection: WARNING: only 2 clustered events; the single "
          "available cycle is used (consider running longer)")
t_a = off[-2]          # arrest of the penultimate event
t_b = on[-1]           # onset of the final event
t_c = off[-1]          # arrest of the final event
phases = [("post", t_a + 1e-4*yr), ("0.25", t_a + 0.25*(t_b - t_a)),
          ("0.50", t_a + 0.50*(t_b - t_a)), ("0.75", t_a + 0.75*(t_b - t_a)),
          ("pre", t_b - 1e-4*yr)]
idx = [min(max(np.searchsorted(t_fr, tt), 0), nfr-1) for _, tt in phases]
i_end = min(np.searchsorted(t_fr, t_c) + 1, nfr-1)

fields = [field_at(i) for i in idx]
f_end  = field_at(i_end)

# save the snapshot fields for plot_xsection_overview.py
np.savez_compressed(
    "xsect_fields.npz",
    xg=xg, zg=zg, H=H, D=Df, tM=tM, label=label,
    t_snap=np.array([t_fr[i] for i in idx]),
    phases=np.array([p for p, _ in phases]),
    sxy=np.array([f[0] for f in fields]),
    syz=np.array([f[1] for f in fields]),
    sxy_end=f_end[0], syz_end=f_end[1], t_end=t_fr[i_end])

# self-test 3: interface traction continuity at every snapshot
iH1 = np.argmin(np.abs(zg - 0.985*H)); iH2 = np.argmin(np.abs(zg - 1.015*H))
for (sxy, syz), (ph, _) in zip(fields, phases):
    num = np.max(np.abs(syz[iH1] - syz[iH2]))
    den = max(np.max(np.abs(syz)), 1.0)
    if num/den > 0.08:
        sys.exit(f"plot_xsection: SELF-TEST 3 FAILED (interface traction "
                 f"jump {num/den:.2%} at phase {ph})")

# ---------------------------------------------------------------------------
ref_xy, ref_yz = fields[0]
def panelplot(comp, cbl, outname, title):
    """panel 1: the ACTUAL post-event field (the reconstructed stress
    anomaly of the deficit history, own color scale); panels 2-5: the
    change since that state through the cycle (shared scale)"""
    fig, axs = plt.subplots(1, 5, figsize=(17.5, 3.9), sharey=True)
    diffs = [fields[k][comp] - fields[0][comp] for k in range(1, 5)]
    vmax = max(np.quantile(np.abs(np.array(diffs)), sat)/1e6, 1e-12)
    va = max(np.quantile(np.abs(fields[0][comp]), sat)/1e6, 1e-12)
    lev0, nrm0, cmp0 = discrete(va)
    im0 = axs[0].pcolormesh(xg/1e3, zg/1e3, fields[0][comp]/1e6,
                            cmap=cmp0, norm=nrm0, shading="auto")
    axs[0].set_title(f"post-event state, t = {t_fr[idx[0]]/yr:.1f} yr",
                     fontsize=9)
    lev, nrm, cmpd = discrete(vmax)
    im = None
    for k in range(1, 5):
        ax = axs[k]
        im = ax.pcolormesh(xg/1e3, zg/1e3, diffs[k-1]/1e6, cmap=cmpd,
                           norm=nrm, shading="auto")
        ax.set_title(f"change, phase {phases[k][0]}, "
                     f"t = {t_fr[idx[k]]/yr:.1f} yr", fontsize=9)
    for ax in axs:
        ax.axhline(H/1e3, color="k", lw=0.8, ls="--")
        ax.plot([0, 0], [0, Df], color="k", lw=2.5)
        ax.set_xlabel("x [km]")
    axs[0].set_ylabel("depth [km]")
    axs[0].set_ylim(zfac*H/1e3, 0)
    fig.colorbar(im0, ax=axs[0], label=cbl + " (state)", shrink=0.8,
                 pad=0.03, extend="both")
    fig.colorbar(im, ax=axs[1:], label="d " + cbl, shrink=0.8, pad=0.01,
                 extend="both")
    fig.suptitle(f"{title}  [{label}]   dashed: plate base H = "
                 f"{H/1e3:.0f} km; thick: fault (D = {Df:.0f} km)",
                 fontsize=10)
    fig.savefig(outname, dpi=150, bbox_inches="tight")
    plt.close(fig)

panelplot(0, "sigma_xy [MPa]", "xsect_sxy.png",
          "stress evolution, sigma_xy")
panelplot(1, "sigma_yz [MPa]", "xsect_syz.png",
          "stress evolution, sigma_yz")

co = f_end[0] - fields[-1][0]
vmax = max(np.quantile(np.abs(co), sat)/1e6, 1e-12)
lev, nrm, cmpd = discrete(vmax)
fig, ax = plt.subplots(figsize=(5.8, 4.8))
im = ax.pcolormesh(xg/1e3, zg/1e3, co/1e6, cmap=cmpd, norm=nrm,
                   shading="auto")
ax.axhline(H/1e3, color="k", lw=0.8, ls="--")
ax.plot([0, 0], [0, Df], color="k", lw=2.5)
ax.set_ylim(zfac*H/1e3, 0)
ax.set_xlabel("x [km]"); ax.set_ylabel("depth [km]")
fig.colorbar(im, ax=ax, label="coseismic d sigma_xy [MPa]", extend="both")
ax.set_title(f"coseismic stress change  [{label}]", fontsize=10)
fig.savefig("xsect_coseis.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"plot_xsection: {label}: N = {N} orders, cycle "
      f"[{t_a/yr:.1f}, {t_b/yr:.1f}] yr; wrote xsect_sxy/syz/coseis.png; "
      "self-tests passed")
