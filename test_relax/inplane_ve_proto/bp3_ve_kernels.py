#!/usr/bin/env python3
"""
bp3_ve_kernels.py: build a -ve_prony_file (rsf_solve -ve_mode 3) for
a 2-D PLANE-STRAIN dipping fault in an ELASTIC PLATE over a MAXWELL
half-space (+ optional gravity), from the validated kernel machinery
in inplane2d.py (PSGRN-cross-checked at the 1-3 percent level; the
same file contract is intended for PSGRN-derived 3-D kernels later).

Method: for each source segment j (unit slip), the RELAXATION part of
the shear-traction kernel row, dK_ij(t) = K_ij(t) - K_ij(0+), is
computed as the DIFFERENCE of the viscoelastic and elastic solutions
on one common truncated k-grid (the singular direct parts cancel, so
dK is regular everywhere including the diagonal; kmax ~ 8/H
suffices).  The relaxed limit dK(inf) comes from an elastic solve
with the substrate shear modulus removed.  Per-pair amplitudes C_p
are least-squares fitted on a FIXED shared ladder tau_p with the
relaxed limit as an exact constraint row and one held-out sample as
the fit gate:

    K(t) = C_inf + sum_p C_p e^{-t/tau_p},   dK(t) = sum_p C_p (e^{-t/tau_p} - 1)

The file's K0 block carries the generator's elastic kernel at
|i-j| > 2 only (near-diagonal entries are written as 0 and skipped
by rsf_solve's consistency gate; the singular elastic self terms
live in the assembled Is = implicit C_inf).

usage: bp3_ve_kernels.py geom_file H_m tM_yr out_file [np] [sign] [g]
       geom_file  bp3_geom.in from gen_bp3.py (x y 0 strike 90 l 0 g)
       H_m        plate thickness [m] (must exceed fault bottom)
       tM_yr      substrate Maxwell time [yr]
       np         ladder size (default 5, max 6)
       sign       +1/-1 global kernel sign (default -1: interact assembles
                  Is in the stress-DROP convention)
                  convention (default +1; fix once via the runtime
                  K0 gate)
       g          gravity flag 0/1 (default 1: rho 2800/3300, 9.81)
"""
import sys
import numpy as np
from inplane2d import fields_xz, ve_fields_xz, segment_sources

geomf = sys.argv[1]
H = float(sys.argv[2])
tM = float(sys.argv[3])*3.15576e7
outf = sys.argv[4]
npr = int(sys.argv[5]) if len(sys.argv) > 5 else 5
sgn = float(sys.argv[6]) if len(sys.argv) > 6 else -1.0
gflag = int(sys.argv[7]) if len(sys.argv) > 7 else 1

mu = lam = 32.04e9           # BP3 elastic parameters
K2 = lam + 2*mu/3
rho1, rho2, g = 2800.0, 3300.0, (9.81 if gflag else 0.0)

# geometry: centers, tangents, normals in the proto frame (z down)
G = np.loadtxt(geomf)
xc, yc = G[:, 0], G[:, 1]
zc = -yc
alpha = np.deg2rad(90.0 - G[:, 3])
tx, tz = np.cos(alpha), -np.sin(alpha)     # tangent, z-down frame
nx, nz = -tz, tx                           # normal (90 deg CCW from t)
hl = G[:, 5]
n = G.shape[0]
if zc.max() + hl.max() >= H:
    sys.exit("bp3_ve_kernels: fault reaches below the plate; increase H")

# shared ladder and sample times
taus = np.geomspace(tM/2.0, 20.0*tM, npr)
ts = np.geomspace(tM/10.0, 60.0*tM, 9)
ihold = 4                                   # held-out sample
kw = dict(kmin=1e-7, kmax=8.0/H, npanel=24, nquad=4)
mur = 1e-10*mu

def shear_rows(F):
    """resolve the field dict onto receiver shear tractions"""
    sxx = F["sxx"]; sxz = F["sxz"]; szz = F["szz"]
    # traction t_i = sigma_ij n_j; shear = t . that
    tx_ = sxx*nxr + sxz*nzr
    tz_ = sxz*nxr + szz*nzr
    return tx_*txr + tz_*tzr

# receivers: evaluate at the segment centers; fields_xz wants separate
# x list and z list -> evaluate per-receiver-depth batches would be
# n^2; instead evaluate all receivers per source using matched (x,z)
# pairs: loop depths is expensive, so use the diagonal trick: one
# call per source with zs = all receiver depths and pick matched x.
# Cost is nz*nx evals internally; for the demonstration sizes
# (n <= 100) this stays manageable.
zs = zc.copy()
xs = xc.copy()
txr, tzr, nxr, nzr = tx, tz, nx, nz

import os
cachef = outf + ".cache.npz"
if os.path.exists(cachef):
    cz = np.load(cachef)
    K0, DK, DKinf, done = (cz["K0"], cz["DK"], cz["DKinf"],
                           cz["done"].astype(bool))
    print(f"bp3_ve_kernels: resuming, {done.sum()}/{n} sources done")
else:
    K0 = np.zeros((n, n))
    DK = np.zeros((len(ts), n, n))
    DKinf = np.zeros((n, n))
    done = np.zeros(n, dtype=bool)
import time as _time
t0_ = _time.time()
for j in range(n):
    if done[j]:
        continue
    x1, z1 = xc[j] - hl[j]*tx[j], zc[j] - hl[j]*tz[j]
    x2, z2 = xc[j] + hl[j]*tx[j], zc[j] + hl[j]*tz[j]
    src = segment_sources(x1, z1, x2, z2, 1.0, ns=2)
    def rows(F):
        # matched (x_i, z_i): fields_xz returns (nz, nx); take diag
        S = {c: np.diag(F[c]) for c in ("sxx", "sxz", "szz")}
        tx_ = S["sxx"]*nxr + S["sxz"]*nzr
        tz_ = S["sxz"]*nxr + S["szz"]*nzr
        return tx_*txr + tz_*tzr
    el = fields_xz(xs, zs, src, lam, mu, lam, mu, H,
                   rho1=rho1, rho2=rho2, g=g, **kw)
    r0 = rows(el)
    # K0 (gate block) needs the high-k content the smooth-dK grid
    # truncates: separate elastic pass with kmax set by the minimum
    # off-diagonal receiver separation (|i-j| > 2 entries only)
    elhi = fields_xz(xs, zs, src, lam, mu, lam, mu, H,
                     rho1=rho1, rho2=rho2, g=g,
                     kmin=1e-7, kmax=25.0/(3.0*2.0*np.min(hl)),
                     npanel=56, nquad=6)
    K0[:, j] = sgn*rows(elhi)
    elr = fields_xz(xs, zs, src, lam, mu, K2 - 2*mur/3, mur, H,
                    rho1=rho1, rho2=rho2, g=g, **kw)
    DKinf[:, j] = sgn*(rows(elr) - r0)
    for k, t in enumerate(ts):
        ve = ve_fields_xz(xs, zs, src, lam, mu, lam, mu, H, tM, t,
                          M=10, rho1=rho1, rho2=rho2, g=g, **kw)
        vr = {c: np.real(ve[c]) for c in ("sxx", "sxz", "szz")}
        DK[k, :, j] = sgn*(rows(vr) - r0)
    done[j] = True
    np.savez(cachef, K0=K0, DK=DK, DKinf=DKinf, done=done)
    print(f"  source {j+1}/{n} done  ({_time.time()-t0_:.0f} s)",
          flush=True)

# per-pair LSQ on the fixed ladder, relaxed limit as constraint row
E = np.exp(-np.outer(ts, 1.0/taus)) - 1.0        # (nt, np)
fit_rows = [k for k in range(len(ts)) if k != ihold]
A = np.vstack([E[fit_rows], -np.ones((1, npr))])  # last row: t = inf
C = np.zeros((npr, n, n))
worst = 0.0
scl = np.max(np.abs(DKinf)) + 1e-30
for i in range(n):
    for j in range(n):
        b = np.concatenate([DK[fit_rows, i, j], [DKinf[i, j]]])
        w = np.ones(len(b)); w[-1] = 3.0          # weight the limit
        cij, *_ = np.linalg.lstsq(A*w[:, None], b*w, rcond=None)
        C[:, i, j] = cij
        pred = E[ihold] @ cij
        r = abs(pred - DK[ihold, i, j])/scl
        if r > worst:
            worst = r
print(f"bp3_ve_kernels: held-out worst residual {worst:.3e} "
      f"(of max |dK_inf| = {scl:.3e} Pa/m)")

with open(outf, "w") as f:
    f.write("# -ve_prony_file for rsf_solve -ve_mode 3\n")
    f.write(f"# generator: inplane_ve_proto plate-over-Maxwell"
            f" H={H/1e3:g} km tM={tM/3.15576e7:g} yr g={g:g}"
            f" np={npr} heldout={worst:.3e}\n")
    f.write("# K(t) = C_inf + sum C_p exp(-t/tau_p); C_inf implicit\n")
    f.write(f"{n} {npr}\n")
    f.write(" ".join(f"{t:.8e}" for t in taus) + "\n")
    K0g = K0.copy()
    for i in range(n):
        for j in range(max(0, i-2), min(n, i+3)):
            K0g[i, j] = 0.0                       # gate skips these
    for M_ in [K0g] + [C[p] for p in range(npr)]:
        for i in range(n):
            f.write(" ".join(f"{v:.8e}" for v in M_[i]) + "\n")
print(f"bp3_ve_kernels: wrote {outf}")
