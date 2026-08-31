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

NORMAL TRACTIONS: with normal=1 (default) the same machinery is
applied a second time to the FAULT-NORMAL traction n.sigma.n, giving
a second family of amplitude matrices on the SAME tau ladder (the
relaxation spectrum belongs to the medium, not to the traction
component).  The file then declares two families in its header and
carries an In gate block plus the normal amplitudes; rsf_solve
activates them exactly when the elastic normal path is on
(-calc_sigma_dot), and refuses that flag with a shear-only file so
that a relaxing-shear / frozen-normal medium can never be run by
accident.  The extra cost is one projection per already-computed
field, i.e. essentially free.

usage: bp3_ve_kernels.py geom_file H_m tM_yr out_file [np] [sign] [g] [normal]
       geom_file  bp3_geom.in from gen_bp3.py (x y 0 strike 90 l 0 g)
       H_m        plate thickness [m] (must exceed fault bottom)
       tM_yr      substrate Maxwell time [yr]
       np         ladder size (default 5, max 6)
       sign       +1/-1 global kernel sign (default -1: interact assembles
                  Is in the stress-DROP convention)
                  convention (default +1; fix once via the runtime
                  K0 gate)
       g          gravity flag 0/1 (default 1: rho 2800/3300, 9.81)
       normal     0/1 (default 1): also fit the normal-traction family
       ipart npart parallel column split (see below); default -1 1
                  worker:  ... [normal] <i> <n>   for i = 0..n-1
                  merge:   ... [normal] -1 <n>    once all parts exist
                  wrapper: bp3_ve_kernels_par does both
"""
import sys
import numpy as np
from inplane2d import (fields_xz, ve_fields_xz, fields_pairs,
                       ve_fields_pairs, segment_sources)

geomf = sys.argv[1]
H = float(sys.argv[2])
tM = float(sys.argv[3])*3.15576e7
outf = sys.argv[4]
npr = int(sys.argv[5]) if len(sys.argv) > 5 else 5
sgn = float(sys.argv[6]) if len(sys.argv) > 6 else -1.0
gflag = int(sys.argv[7]) if len(sys.argv) > 7 else 1
donorm = int(sys.argv[8]) if len(sys.argv) > 8 else 1
# PARALLEL over source columns (trivially parallel: columns are
# independent).  ipart >= 0 with npart > 1 computes only the columns
# j with j % npart == ipart and stores them in <out>.part<i>.npz
# WITHOUT fitting or writing the kernel file; a final call with
# ipart < 0 merges all parts, fits, and writes.  npart == 1 is the
# ordinary serial path (unchanged, single <out>.cache.npz).
ipart = int(sys.argv[9]) if len(sys.argv) > 9 else -1
npart = int(sys.argv[10]) if len(sys.argv) > 10 else 1
if (npart > 1) and (ipart >= npart):
    sys.exit("bp3_ve_kernels: ipart must be < npart")
# interact assembles the NORMAL interaction matrix In compression
# positive (In is scaled by -1 internally), while the projection
# n.sigma.n below is tension positive, so the normal family carries
# the opposite global sign to the shear family.  Both are pinned by
# the runtime gates (K0 vs Is, N0 vs In), which is how this was
# found: a pure convention flip shows up as a gate deviation of
# exactly 2.0
sgnn = -sgn

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
if npart > 1 and ipart >= 0:
    cachef = f"{outf}.part{ipart}.npz"
else:
    cachef = outf + ".cache.npz"
if os.path.exists(cachef):
    cz = np.load(cachef)
    K0, DK, DKinf, done = (cz["K0"], cz["DK"], cz["DKinf"],
                           cz["done"].astype(bool))
    N0, DN, DNinf = cz["N0"], cz["DN"], cz["DNinf"]
    print(f"bp3_ve_kernels: resuming, {done.sum()}/{n} sources done")
else:
    K0 = np.zeros((n, n))
    DK = np.zeros((len(ts), n, n))
    DKinf = np.zeros((n, n))
    N0 = np.zeros((n, n))
    DN = np.zeros((len(ts), n, n))
    DNinf = np.zeros((n, n))
    done = np.zeros(n, dtype=bool)
import time as _time
t0_ = _time.time()
cols = range(n) if (npart <= 1 or ipart < 0) else \
    range(ipart, n, npart)
for j in cols:
    if done[j]:
        continue
    x1, z1 = xc[j] - hl[j]*tx[j], zc[j] - hl[j]*tz[j]
    x2, z2 = xc[j] + hl[j]*tx[j], zc[j] + hl[j]*tz[j]
    src = segment_sources(x1, z1, x2, z2, 1.0, ns=2)
    def rows(F, comp="shear"):
        """matched (x_i, z_i) traction resolved on the receiver:
        'shear' = t.sigma.n projected on the tangent, 'normal' =
        n.sigma.n (positive in tension here; the sign convention is
        fixed once by the runtime gate against In)"""
        S = {c: F[c] for c in ("sxx", "sxz", "szz")}
        tx_ = S["sxx"]*nxr + S["sxz"]*nzr
        tz_ = S["sxz"]*nxr + S["szz"]*nzr
        if comp == "normal":
            return tx_*nxr + tz_*nzr
        return tx_*txr + tz_*tzr
    el = fields_pairs(xs, zs, src, lam, mu, lam, mu, H,
                      rho1=rho1, rho2=rho2, g=g, **kw)
    r0 = rows(el)
    # K0 (gate block) needs the high-k content the smooth-dK grid
    # truncates: separate elastic pass with kmax set by the minimum
    # off-diagonal receiver separation (|i-j| > 2 entries only)
    elhi = fields_pairs(xs, zs, src, lam, mu, lam, mu, H,
                        rho1=rho1, rho2=rho2, g=g,
                        kmin=1e-7, kmax=25.0/(3.0*2.0*np.min(hl)),
                        npanel=56, nquad=6)
    K0[:, j] = sgn*rows(elhi)
    n0 = rows(el, "normal")
    N0[:, j] = sgnn*rows(elhi, "normal")
    elr = fields_pairs(xs, zs, src, lam, mu, K2 - 2*mur/3, mur, H,
                       rho1=rho1, rho2=rho2, g=g, **kw)
    DKinf[:, j] = sgn*(rows(elr) - r0)
    DNinf[:, j] = sgnn*(rows(elr, "normal") - n0)
    for k, t in enumerate(ts):
        ve = ve_fields_pairs(xs, zs, src, lam, mu, lam, mu, H, tM, t,
                             M=10, rho1=rho1, rho2=rho2, g=g, **kw)
        vr = {c: np.real(ve[c]) for c in ("sxx", "sxz", "szz")}
        DK[k, :, j] = sgn*(rows(vr) - r0)
        DN[k, :, j] = sgnn*(rows(vr, "normal") - n0)
    done[j] = True
    np.savez(cachef, K0=K0, DK=DK, DKinf=DKinf, N0=N0, DN=DN,
             DNinf=DNinf, done=done)
    print(f"  source {j+1}/{n} done  ({_time.time()-t0_:.0f} s)",
          flush=True)

if npart > 1 and ipart >= 0:
    nd = int(done[ipart::npart].sum())
    print(f"bp3_ve_kernels: part {ipart}/{npart} done "
          f"({nd} of {len(range(ipart, n, npart))} columns) -> {cachef}")
    raise SystemExit(0)
if npart > 1:
    # MERGE: columns are disjoint between parts, so summing the
    # arrays and OR-ing the done masks reconstructs the whole kernel
    K0 = np.zeros((n, n)); N0 = np.zeros((n, n))
    DK = np.zeros((len(ts), n, n)); DN = np.zeros((len(ts), n, n))
    DKinf = np.zeros((n, n)); DNinf = np.zeros((n, n))
    done = np.zeros(n, dtype=bool)
    for ip in range(npart):
        f = f"{outf}.part{ip}.npz"
        if not os.path.exists(f):
            sys.exit(f"bp3_ve_kernels: missing {f}; run that worker first")
        z = np.load(f)
        K0 += z["K0"]; N0 += z["N0"]
        DK += z["DK"]; DN += z["DN"]
        DKinf += z["DKinf"]; DNinf += z["DNinf"]
        done |= z["done"].astype(bool)
    if not done.all():
        sys.exit(f"bp3_ve_kernels: {int((~done).sum())} columns still "
                 "missing across the parts")
    print(f"bp3_ve_kernels: merged {npart} parts, all {n} columns present")

# per-pair LSQ on the fixed ladder, relaxed limit as constraint row
E = np.exp(-np.outer(ts, 1.0/taus)) - 1.0        # (nt, np)
fit_rows = [k for k in range(len(ts)) if k != ihold]
A = np.vstack([E[fit_rows], -np.ones((1, npr))])  # last row: t = inf

def fit_family(D, Dinf, label):
    """per-pair amplitudes C_p from sampled dK(t) and dK(inf).
    ALL pairs are solved at once through the (small, well
    conditioned) weighted normal equations A^T W A c = A^T W b,
    which is the same least-squares solution as a per-pair lstsq
    but O(n^2) numpy work instead of n^2 python calls; returns
    (C, worst held-out residual relative to the relaxation
    scale)"""
    scl = np.max(np.abs(Dinf)) + 1e-30
    w = np.ones(A.shape[0]); w[-1] = 3.0          # weight the limit
    Aw = A*w[:, None]
    G = Aw.T @ Aw                                  # (np, np)
    # right-hand sides for every pair: (nfit+1, n, n)
    B = np.concatenate([D[fit_rows], Dinf[None]], axis=0)*w[:, None, None]
    rhs = np.einsum("kp,kij->pij", Aw, B)          # (np, n, n)
    C = np.linalg.solve(G, rhs.reshape(npr, -1)).reshape(npr, n, n)
    worst = float(np.max(np.abs(np.einsum("p,pij->ij", E[ihold], C)
                                - D[ihold])))/scl
    print(f"bp3_ve_kernels: {label}: held-out worst residual "
          f"{worst:.3e} (of max |d{label[0].upper()}_inf| = "
          f"{scl:.3e} Pa/m)")
    return C, worst

C, worst = fit_family(DK, DKinf, "shear")
if donorm:
    Cn, worstn = fit_family(DN, DNinf, "normal")

def gate_block(M):
    """near-diagonal entries are set to zero: the runtime gate skips
    them (the singular elastic self terms live in the assembled
    operator = implicit C_inf)"""
    G = M.copy()
    for i in range(n):
        for j in range(max(0, i-2), min(n, i+3)):
            G[i, j] = 0.0
    return G

nfam = 2 if donorm else 1
with open(outf, "w") as f:
    f.write("# -ve_prony_file for rsf_solve -ve_mode 3\n")
    f.write(f"# generator: inplane_ve_proto plate-over-Maxwell"
            f" H={H/1e3:g} km tM={tM/3.15576e7:g} yr g={g:g}"
            f" np={npr} nfam={nfam} heldout={worst:.3e}"
            + (f"/{worstn:.3e}" if donorm else "") + "\n")
    f.write("# K(t) = C_inf + sum C_p exp(-t/tau_p); C_inf implicit\n")
    f.write("# blocks: K0 gate, C_1..C_np"
            + (", N0 gate, Cn_1..Cn_np\n" if donorm else "\n"))
    f.write(f"{n} {npr} {nfam}\n")
    f.write(" ".join(f"{t:.8e}" for t in taus) + "\n")
    blocks = [gate_block(K0)] + [C[p] for p in range(npr)]
    if donorm:
        blocks += [gate_block(N0)] + [Cn[p] for p in range(npr)]
    for M_ in blocks:
        for i in range(n):
            f.write(" ".join(f"{v:.8e}" for v in M_[i]) + "\n")
print(f"bp3_ve_kernels: wrote {outf} ({nfam} traction "
      f"{'families' if nfam > 1 else 'family'})")
