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

MAXWELL-TIME SCALING: the substrate enters the wavenumber-domain
system only through mu2(s) = mu2 s/(s + 1/tau) and lam2(s) =
K2 - 2 mu2(s)/3, and the buoyancy rows carry no s at all, so with
s = sigma/tau the system depends on s only through s*tau.  The step
response then satisfies

    dK(t; tau) = dKtilde(t/tau)

exactly.  The sampling times and the tau ladder are accordingly
fixed multiples of tM here, so in units of tM they are the same set
for every Maxwell time: the sampled cache carries no Maxwell time,
and a file for any tM follows by scaling the ladder, with the
amplitude blocks unchanged.  Two consequences:

  * cache and part files are named after a SIGNATURE of what the
    samples actually depend on (geometry, H, gravity, sign,
    sampling), not after the output file, so a sweep over Maxwell
    times writing into one directory samples once and the remaining
    files cost only the fit;
  * tM_yr may be a comma-separated LIST, and then out_file must
    contain %TM%, which is replaced by each Maxwell time as written
    on the command line.

The scaling is checked against direct generation by
run_tau_scaling_test, and on the sampled arrays by t6_tau_scaling.py;
agreement is at roundoff level rather than at fit-noise level,
because the sample set in units of tM is the same one.  A file
written for the Maxwell time that was sampled is bit-identical to
what a single-Maxwell-time run produced before this scaling path
existed; a file written for another Maxwell time differs from a
direct generation of it only through the rounding of the scaled
ladder, at the 1e-16 level.  It holds for a SINGLE Maxwell time with an
elastic plate.  A second relaxation time (Maxwell plate, Burgers or
bi-viscous substrate, depth-dependent viscosity) leaves only the
ratios free and would need one sampling pass per case.

usage: bp3_ve_kernels.py geom_file H_m tM_yr out_file [np] [sign] [g] [normal]
       geom_file  bp3_geom.in from gen_bp3.py (x y 0 strike 90 l 0 g)
       H_m        plate thickness [m] (must exceed fault bottom)
       tM_yr      substrate Maxwell time [yr], or a comma-separated
                  list, e.g. 15,45,200 (then out_file needs %TM%)
       out_file   output file; %TM% is replaced by the Maxwell time
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
import os
import glob
import hashlib
import numpy as np
from inplane2d import (fields_xz, ve_fields_xz, fields_pairs,
                       ve_fields_pairs, segment_sources)

YR = 3.15576e7
# bump when anything that changes the SAMPLED arrays changes (grids,
# sample times, projections, sign folding): old caches are then not
# picked up by accident
CACHE_VERSION = "v2-dimensionless"

geomf = sys.argv[1]
H = float(sys.argv[2])
tmarg = sys.argv[3]
outf = sys.argv[4]
npr = int(sys.argv[5]) if len(sys.argv) > 5 else 5
sgn = float(sys.argv[6]) if len(sys.argv) > 6 else -1.0
gflag = int(sys.argv[7]) if len(sys.argv) > 7 else 1
donorm = int(sys.argv[8]) if len(sys.argv) > 8 else 1
# PARALLEL over source columns (trivially parallel: columns are
# independent).  ipart >= 0 with npart > 1 computes only the columns
# j with j % npart == ipart and stores them in a .part<i>.npz
# WITHOUT fitting or writing the kernel file; a final call with
# ipart < 0 merges all parts, fits, and writes.  npart == 1 is the
# ordinary serial path (single .cache.npz).
ipart = int(sys.argv[9]) if len(sys.argv) > 9 else -1
npart = int(sys.argv[10]) if len(sys.argv) > 10 else 1
if (npart > 1) and (ipart >= npart):
    sys.exit("bp3_ve_kernels: ipart must be < npart")

# Maxwell times: one sampling pass, one output file per requested tM
tmtok = [t for t in tmarg.replace(",", " ").split() if t]
if not tmtok:
    sys.exit("bp3_ve_kernels: need at least one Maxwell time")
try:
    tmyr = [float(t) for t in tmtok]
except ValueError:
    sys.exit(f"bp3_ve_kernels: cannot read '{tmarg}' as Maxwell time(s)")
if min(tmyr) <= 0.0:
    sys.exit("bp3_ve_kernels: Maxwell times must be positive")
if len(tmtok) > 1 and "%TM%" not in outf:
    sys.exit("bp3_ve_kernels: a list of Maxwell times needs %TM% in "
             "out_file, e.g. bp3_tm%TM%_g1n1.prony")
outs = [outf.replace("%TM%", t) for t in tmtok]
if len(set(outs)) != len(outs):
    sys.exit("bp3_ve_kernels: the Maxwell times give duplicate file names")
tM = tmyr[0]*YR                 # reference: what the sampling uses

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

# shared ladder and sample times, at the REFERENCE Maxwell time.
# Both are fixed multiples of tM, so in units of tM they are the same
# set for every Maxwell time; that is what makes the sampled arrays
# and the fitted amplitudes Maxwell-time independent (see the
# docstring).  They are formed here in absolute time, exactly as a
# single-Maxwell-time run has always formed them, so that a file
# written for the reference Maxwell time is unchanged to the last
# digit; the other files carry the same ladder scaled by tM/tM_ref.
TS_SPAN = (0.1, 60.0, 9)                    # t/tM: first, last, count
taus = np.geomspace(tM/2.0, 20.0*tM, npr)
ts = np.geomspace(TS_SPAN[0]*tM, TS_SPAN[1]*tM, TS_SPAN[2])
ihold = 4                                   # held-out sample
nsrc = 2                                    # subsegments per source
kw = dict(kmin=1e-7, kmax=8.0/H, npanel=24, nquad=4)
kwhi = dict(kmin=1e-7, kmax=25.0/(3.0*2.0*np.min(hl)),
            npanel=56, nquad=6)
murel = 1e-10
mur = murel*mu

# receivers: the segment centers, as matched (x, z) pairs
zs = zc.copy()
xs = xc.copy()
txr, tzr, nxr, nzr = tx, tz, nx, nz

def rows(F, comp="shear"):
    """matched (x_i, z_i) traction resolved on the receiver:
    'shear' = t.sigma.n projected on the tangent, 'normal' =
    n.sigma.n (positive in tension here; the sign convention is
    fixed once by the runtime gate against In)"""
    tx_ = F["sxx"]*nxr + F["sxz"]*nzr
    tz_ = F["sxz"]*nxr + F["szz"]*nzr
    if comp == "normal":
        return tx_*nxr + tz_*nzr
    return tx_*txr + tz_*tzr

# ----------------------------------------------------------------------
# cache naming: a signature of everything the SAMPLED arrays depend
# on.  The Maxwell time is deliberately NOT part of it (the samples
# are taken at fixed t/tau), so Maxwell times written into the same
# directory share one cache.
with open(geomf, "rb") as _f:
    _geom_bytes = _f.read()
# note: the sample times enter as their t/tM SPAN, not as values,
# so that rounding in tM*(t/tM) cannot make two Maxwell times miss
# each other's cache
_sig_src = repr((CACHE_VERSION, _geom_bytes, H, sgn, gflag, nsrc,
                 TS_SPAN, tuple(sorted(kw.items())),
                 tuple(sorted(kwhi.items())), murel, rho1, rho2,
                 mu, lam))
sig = hashlib.blake2b(_sig_src.encode("utf-8", "surrogateescape"),
                      digest_size=6).hexdigest()
cdir = os.path.dirname(os.path.abspath(outs[0]))
cachebase = os.path.join(cdir, f"ve_kernels_{sig}")
fullf = cachebase + ".cache.npz"
partf = cachebase + ".part%d.npz"

CKEYS = ("K0", "DK", "DKinf", "N0", "DN", "DNinf", "done")
# the Maxwell time the samples were taken at.  It does not enter the
# fit (only t/tau does), and is kept for provenance and for the
# message when a later call asks for a different one.  Absent from
# caches written before the scaling path existed.
tm_sampled = None

def load_cache(path):
    """read a cache/part file, returning None if it is absent or does
    not match this geometry and sampling"""
    if not os.path.exists(path):
        return None
    try:
        z = np.load(path)
        d = {k: z[k] for k in CKEYS}
        d["tm_sampled"] = (float(z["tm_sampled"])
                           if "tm_sampled" in z.files else None)
    except Exception as e:
        print(f"bp3_ve_kernels: ignoring {path} ({e})")
        return None
    if (d["K0"].shape != (n, n) or d["DK"].shape != (len(ts), n, n)
            or d["done"].shape != (n,)):
        print(f"bp3_ve_kernels: ignoring {path} (shape mismatch)")
        return None
    d["done"] = d["done"].astype(bool)
    return d

def save_cache(path, d):
    """write through a temporary file: the cache is rewritten after
    every column, and a run killed during the write would otherwise
    leave a truncated archive (hours of sampling)"""
    tmp = path + ".tmp.npz"
    out = {k: d[k] for k in CKEYS}
    out["tm_sampled"] = (tmyr[0] if d.get("tm_sampled") is None
                         else d["tm_sampled"])
    np.savez(tmp, **out)
    os.replace(tmp, path)

def new_cache():
    return dict(K0=np.zeros((n, n)), DK=np.zeros((len(ts), n, n)),
                DKinf=np.zeros((n, n)), N0=np.zeros((n, n)),
                DN=np.zeros((len(ts), n, n)),
                DNinf=np.zeros((n, n)),
                done=np.zeros(n, dtype=bool), tm_sampled=None)

def merge_into(dst, src):
    """columns are disjoint between parts, so summing the arrays and
    OR-ing the done masks reconstructs the whole kernel"""
    for k in CKEYS[:-1]:
        dst[k] = dst[k] + src[k]
    dst["done"] = dst["done"] | src["done"]
    if dst.get("tm_sampled") is None:
        dst["tm_sampled"] = src.get("tm_sampled")
    return dst

def collect_parts(paths):
    """assemble a full cache from part files, or None if any of them
    is missing, unreadable or leaves a column uncovered"""
    if not paths:
        return None
    out = new_cache()
    for f in paths:
        d = load_cache(f)
        if d is None:
            return None
        out = merge_into(out, d)
    return out if out["done"].all() else None

def full_cache():
    """the shared serial cache, only if it covers every column"""
    d = load_cache(fullf)
    return d if (d is not None and d["done"].all()) else None

# the pre-scaling version of this script named the cache after the
# OUTPUT file.  Adopt such a cache rather than resampling: the
# sampling itself is unchanged, only the name is.
def adopt_legacy(path, target):
    if os.path.exists(target):
        return False
    d = load_cache(path)
    if d is None:
        return False
    save_cache(target, d)
    print(f"bp3_ve_kernels: adopted {os.path.basename(path)} as "
          f"{os.path.basename(target)} ({int(d['done'].sum())}/{n} "
          "columns; delete it to resample from scratch)")
    return True

worker = (npart > 1 and ipart >= 0)
if worker:
    adopt_legacy(f"{outs[0]}.part{ipart}.npz", partf % ipart)
    cachef = partf % ipart
    cache = load_cache(cachef)
    # a complete serial cache covers this worker's columns as well
    if (cache is None or not cache["done"].all()) and \
            full_cache() is not None:
        print(f"bp3_ve_kernels: part {ipart}/{npart}: complete "
              f"{os.path.basename(fullf)} present, nothing to sample")
        raise SystemExit(0)
else:
    adopt_legacy(outs[0] + ".cache.npz", fullf)
    cachef = fullf
    cache = load_cache(cachef)
    if npart > 1:
        # merge path: the parts if they are all there, otherwise a
        # complete serial cache (the two paths are interchangeable)
        cache = collect_parts([partf % ip for ip in range(npart)]) \
            or full_cache()
        if cache is None:
            missing = [ip for ip in range(npart)
                       if load_cache(partf % ip) is None]
            sys.exit("bp3_ve_kernels: cannot assemble the kernel: "
                     f"missing or unreadable parts {missing} and no "
                     f"complete {os.path.basename(fullf)}")
        print(f"bp3_ve_kernels: assembled all {n} columns")
    elif cache is None or not cache["done"].all():
        # serial call in a directory where a PARALLEL run left parts
        cache = collect_parts(sorted(glob.glob(cachebase + ".part*.npz"))) \
            or cache

if cache is None:
    cache = new_cache()
else:
    nd = int(cache["done"].sum())
    if nd and nd < n:
        print(f"bp3_ve_kernels: resuming, {nd}/{n} sources done")
tm_sampled = cache.get("tm_sampled") or tmyr[0]
if tm_sampled != tmyr[0]:
    print(f"bp3_ve_kernels: cache was sampled at tM = {tm_sampled:g} "
          f"yr; reused for tM = {tmyr[0]:g} yr by scaling the ladder "
          "(dK depends on t/tau)")
K0, DK, DKinf = cache["K0"], cache["DK"], cache["DKinf"]
N0, DN, DNinf = cache["N0"], cache["DN"], cache["DNinf"]
done = cache["done"]

# ----------------------------------------------------------------------
# sampling pass
import time as _time
t0_ = _time.time()
cols = range(n) if not worker else range(ipart, n, npart)
for j in cols:
    if done[j]:
        continue
    x1, z1 = xc[j] - hl[j]*tx[j], zc[j] - hl[j]*tz[j]
    x2, z2 = xc[j] + hl[j]*tx[j], zc[j] + hl[j]*tz[j]
    src = segment_sources(x1, z1, x2, z2, 1.0, ns=nsrc)
    el = fields_pairs(xs, zs, src, lam, mu, lam, mu, H,
                      rho1=rho1, rho2=rho2, g=g, **kw)
    r0 = rows(el)
    # K0 (gate block) needs the high-k content the smooth-dK grid
    # truncates: separate elastic pass with kmax set by the minimum
    # off-diagonal receiver separation (|i-j| > 2 entries only)
    elhi = fields_pairs(xs, zs, src, lam, mu, lam, mu, H,
                        rho1=rho1, rho2=rho2, g=g, **kwhi)
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
    save_cache(cachef, dict(K0=K0, DK=DK, DKinf=DKinf, N0=N0, DN=DN,
                            DNinf=DNinf, done=done))
    print(f"  source {j+1}/{n} done  ({_time.time()-t0_:.0f} s)",
          flush=True)

if worker:
    nd = int(done[ipart::npart].sum())
    print(f"bp3_ve_kernels: part {ipart}/{npart} done "
          f"({nd} of {len(range(ipart, n, npart))} columns) -> {cachef}")
    raise SystemExit(0)
if not done.all():
    sys.exit(f"bp3_ve_kernels: {int((~done).sum())} columns missing")

# ----------------------------------------------------------------------
# per-pair LSQ on the fixed ladder, relaxed limit as constraint row.
# Done ONCE for all requested Maxwell times: t and tau enter only
# through t/tau, which is the same set for each of them, so the
# amplitudes are Maxwell-time independent (see the docstring).
# only t/tau enters: if the samples come from a cache taken at
# another Maxwell time, these ts and taus are that time's own values
# up to the rounding of the two geomspace calls (1e-16 relative)
E = np.exp(-np.outer(ts, 1.0/taus)) - 1.0        # (nt, np)
fit_rows = [k for k in range(len(ts)) if k != ihold]
A = np.vstack([E[fit_rows], -np.ones((1, npr))])    # last row: t = inf

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
    Gb = M.copy()
    for i in range(n):
        for j in range(max(0, i-2), min(n, i+3)):
            Gb[i, j] = 0.0
    return Gb

nfam = 2 if donorm else 1
blocks = [gate_block(K0)] + [C[p] for p in range(npr)]
if donorm:
    blocks += [gate_block(N0)] + [Cn[p] for p in range(npr)]
for tok, tm_yr, out in zip(tmtok, tmyr, outs):
    fac = tm_yr/tmyr[0]                 # only the ladder scales
    taus_out = fac*taus
    with open(out, "w") as f:
        f.write("# -ve_prony_file for rsf_solve -ve_mode 3\n")
        f.write(f"# generator: inplane_ve_proto plate-over-Maxwell"
                f" H={H/1e3:g} km tM={tm_yr:g} yr g={g:g}"
                f" np={npr} nfam={nfam} heldout={worst:.3e}"
                + (f"/{worstn:.3e}" if donorm else "") + "\n")
        if fac != 1.0 or tm_sampled != tmyr[0]:
            f.write(f"# sampled at tM={tm_sampled:g} yr; ladder"
                    f" scaled by {tm_yr/tm_sampled:.8g} (dK is a"
                    " function of t/tau, see run_tau_scaling_test)\n")
        f.write("# K(t) = C_inf + sum C_p exp(-t/tau_p); C_inf implicit\n")
        f.write("# blocks: K0 gate, C_1..C_np"
                + (", N0 gate, Cn_1..Cn_np\n" if donorm else "\n"))
        f.write(f"{n} {npr} {nfam}\n")
        f.write(" ".join(f"{t:.8e}" for t in taus_out) + "\n")
        for M_ in blocks:
            for i in range(n):
                f.write(" ".join(f"{v:.8e}" for v in M_[i]) + "\n")
    print(f"bp3_ve_kernels: wrote {out} (tM = {tok} yr, {nfam} traction "
          f"{'families' if nfam > 1 else 'family'})")
