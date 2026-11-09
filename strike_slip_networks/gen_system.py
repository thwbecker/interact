#!/usr/bin/env python3
"""
gen_system.py: write rsf_solve input files for a quasi 2-D ("2.5-D")
strike-slip fault system.

Each fault is a strip of ncell rectangular Okada patches along strike,
one patch deep (from the free surface to the locking depth W), so slip
varies only along strike while stress transfer is computed in the 3-D
half space with the finite depth extent.  All faults are vertical and
slip in mode II (strike slip, -rsf_slip_mode 0); each fault has its
own strike (--strikes).  Loading is a uniform far-field horizontal
stressing-rate tensor, resolved onto each fault as shear (along
strike, positive = left-lateral drive, interact's slip convention) and
normal (compression positive) rates, written to -rsf_stress_rate_file
and used with -vpl 0 (no backslip).  Normal stress evolves
(-calc_sigma_dot) from the far-field normal rate and from the slip on
the other faults through the normal-stress interaction matrix, with
the limiter -limit_sigma; --sigma_evolve 0 switches both off, and
--sigma_interact 0 keeps the far-field rate but drops the change from
slip on the faults.

Far-field tensor (compression positive, horizontal, vertical rate 0):
principal rates --sr_s1 >= --sr_s3 with the S1 axis at azimuth --sr_az
(degrees clockwise from north = +y).  With the defaults (sr_s1 = 0,
sr_s3 = 0) a pure shear field s1 = tau_ref, s3 = -tau_ref is used,
where tau_ref = --taudot or G v_eq/(pi W), and sr_az = 135 (S1 axis
NW-SE), so a fault of strike 0 receives tau_dot = +tau_ref
(left-lateral) and no normal rate, as in the first version of this
generator; --rot (default 90) rotates the whole system including
sr_az, so the default layout has the faults along x (strike 90),
driven left-lateral by S1 at azimuth 225.  In general, with
alpha = phi - sr_az,
  tau_dot   = (s1 - s3)/2 sin(2 alpha)                 (> 0: left-lateral)
  sigma_dot = (s1 + s3)/2 - (s1 - s3)/2 cos(2 alpha)   (compression positive)
(a fault striking 45 degrees clockwise from S1 is driven left-lateral,
45 degrees counterclockwise right-lateral; e.g. the Garlock and San
Andreas faults with S1 roughly north-south).  Faults oblique to the
optimal orientation accumulate normal stress linearly in time, which
has no steady state (the limiter keeps sigma within --sigma_min,
--sigma_max).  --ff_normal 0 drops the far-field normal rate but keeps
the interaction-induced normal stress changes.

Files written to --out directory:
  geom.in    x y z strike dip hlength hwidth group   (m, deg; group = fault index)
  rsf.dat    a b per cell
  ic.in      tau[Pa] v[m/s] per cell (steady state at v_init with
             optional log-normal perturbation, see --v_pert)
  srate.in   tau_dot sigma_dot [Pa/s] per cell (resolved far-field rates)
  groups.in  fault index per cell (same as last geom column)
  params.json  all parameters plus derived scales (Lb, h*, N)
  run.opts   rsf_solve command line fragment matching these files

Parameter notes
  --taudot gives the stressing rate directly in Pa/s.  If it is 0 (the
  default), tau_dot = G v_eq/(pi W) is used, the stressing rate at the
  surface above an infinitely long strike-slip fault creeping at v_eq
  below depth W.  The long-term slip rate of a finite, stress-driven
  fault is an outcome of the model, not v_eq.

  --vs_margin puts velocity-strengthening cells (a = --a_vs) over that
  fraction of each fault's length at both ends.  0 = whole fault VW.
  --tip_creep L_t grades the outer L_t of each fault end from a_vs at
  the inner edge to --tip_a at the end (creeping tips that taper the
  slip toward the end); it overrides the profile there.
  --ab_noise adds along-strike Gaussian noise to a (std given as a
  fraction of a), seeded with --seed, for heterogeneous models.

Fault layout for nfault > 1 (en echelon): faults are numbered along
+y; fault k is shifted by k*--step_x in x (perpendicular to strike)
and starts at y = k*(L - overlap).  --overlap > 0 overlapping tips,
< 0 a gap.  step_x should not be 0 when overlap > 0 (coplanar
overlapping patches).  --stepover right (default) puts each next
fault to the right when looking along strike, --stepover left to the
left (the sign of step_x is flipped).  For left-lateral slip a right
step is restraining and a left step releasing, for right-lateral slip
the reverse; the type of each step is printed and written to
params.json.  The layout is built for strike 0 and then rotated by
--rot.

--layout chain instead places the faults tip to tip: fault k starts
--gap beyond the end of fault k-1 (along the strike of k-1, or with
--gap_dir bisect along the bisector of the strikes of k-1 and k, which
keeps point-symmetric chains symmetric) and runs
along its own strike, so --strikes gives the shape of a segmented
trace; --lengths sets one length per fault (chain only).  --step_x,
--overlap, --stepover and --pivot are ignored for the chain.

Faults with a strike other than 0 in --strikes are rotated about their
centre (--pivot center, default) or about their end nearest to the
middle of the system (--pivot inner), which keeps the stepover
geometry at the inner tips unchanged; with --pivot inner the middle
fault of an odd number of faults is rotated about its centre.
"""
import argparse
import json
import math
import os

import numpy as np

SEC_PER_YEAR = 365.25 * 86400.0


def steady_state_tau(v, a, b, sigma, f0, v0):
    """shear stress at steady state for slip rate v with the regularized law
    v = 2 v0 exp(-psi/a) sinh(tau/(sigma a)), psi_ss = f0 - b ln(v/v0)"""
    psi = f0 - b * np.log(v / v0)
    return sigma * a * np.arcsinh(v / (2.0 * v0) * np.exp(psi / a))


def a_profile(ncell, dx, args, rng):
    """direct-effect parameter a along strike of one fault: base value,
    velocity-strengthening margins, n_asp asperities separated by VS
    barriers, and optional (correlated) multiplicative noise"""
    frac = (np.arange(ncell) + 0.5) / ncell
    a = np.full(ncell, args.a)
    vs = (frac < args.vs_margin) | (frac > 1.0 - args.vs_margin)
    if args.n_asp > 1:
        inner = 1.0 - 2.0 * args.vs_margin
        nb = args.n_asp - 1
        wb = args.barrier_frac * inner / nb
        # asperity lengths: equal, or random with relative spread asp_rand
        wa = np.full(args.n_asp, (inner - nb * wb) / args.n_asp)
        if args.asp_rand > 0:
            wa = wa * np.exp(args.asp_rand * rng.standard_normal(args.n_asp))
            wa = wa / wa.sum() * (inner - nb * wb)
        pos = args.vs_margin
        for k in range(args.n_asp):
            pos += wa[k]
            if k < nb:
                vs |= (frac >= pos) & (frac < pos + wb)
                pos += wb
    a[vs] = args.a_vs
    if args.tip_creep > 0:
        # creeping tips: over --tip_creep from each end, a increases
        # linearly from a_vs at the inner edge to tip_a at the end
        dist = np.minimum(frac, 1.0 - frac) * ncell * dx
        tip = dist < args.tip_creep
        w = 1.0 - dist[tip] / args.tip_creep
        a[tip] = args.a_vs + (args.tip_a - args.a_vs) * w
    if args.ab_noise > 0:
        noise = rng.standard_normal(ncell)
        if args.ab_corr > 0:
            # Gaussian smoothing with correlation length ab_corr, then
            # restore unit standard deviation
            nk = int(math.ceil(3.0 * args.ab_corr / dx))
            kx = np.arange(-nk, nk + 1) * dx
            ker = np.exp(-0.5 * (kx / args.ab_corr) ** 2)
            ker /= ker.sum()
            noise = np.convolve(noise, ker, mode="same")
            sd = noise.std()
            if sd > 0:
                noise /= sd
        a = a * (1.0 + args.ab_noise * noise)
        a = np.maximum(a, 1e-4)
    return a


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default=".", help="output directory")
    p.add_argument("--nfault", type=int, default=1)
    p.add_argument("--L", type=float, default=60e3, help="fault length [m]")
    p.add_argument("--W", type=float, default=15e3, help="locking depth = patch depth extent [m]")
    p.add_argument("--dx", type=float, default=0.0,
                   help="along-strike cell size [m]; 0: Lb/--cells_per_lb, rounded to divide L")
    p.add_argument("--cells_per_lb", type=float, default=3.0)
    p.add_argument("--step_x", type=float, default=5e3, help="en echelon step perpendicular to strike [m]")
    p.add_argument("--overlap", type=float, default=0.0, help="along-strike overlap (>0) or gap (<0) [m]")
    p.add_argument("--layout", choices=("echelon", "chain"), default="echelon",
                   help="echelon (default): parallel offsets by --step_x, --overlap; chain: each fault "
                        "starts --gap beyond the end of the previous one, along its own strike")
    p.add_argument("--gap", type=float, default=2e3, help="tip-to-tip gap for --layout chain [m]")
    p.add_argument("--gap_dir", choices=("previous", "bisect"), default="previous",
                   help="direction of the chain gaps: along the previous fault (default) or along the "
                        "bisector of the two faults (point-symmetric chains stay symmetric)")
    p.add_argument("--lengths", default="",
                   help="comma separated length per fault [m] (--layout chain only); default all --L")
    p.add_argument("--stepover", choices=("right", "left"), default="right",
                   help="each next fault to the right (default) or left looking along strike")
    p.add_argument("--pivot", choices=("center", "inner"), default="center",
                   help="rotation point for faults with a non-zero entry in --strikes")
    p.add_argument("--strikes", default="",
                   help="comma separated strike per fault [deg clockwise from north]; default all 0")
    p.add_argument("--sr_s1", type=float, default=0.0, help="far-field S1 stressing rate [Pa/s], compression positive")
    p.add_argument("--sr_s3", type=float, default=0.0, help="far-field S3 stressing rate [Pa/s]; both 0: pure shear tau_ref")
    p.add_argument("--sr_az", type=float, default=135.0, help="azimuth of the S1 axis [deg clockwise from north]")
    p.add_argument("--ff_normal", type=int, default=1, help="1: include the far-field normal stressing rate")
    p.add_argument("--rot", type=float, default=90.0,
                   help="rotate the whole system (geometry, strikes, S1 azimuth) clockwise by this [deg]; "
                        "default 90 lays the faults out left to right (strike 90 for --strikes 0); "
                        "a pure rotation must not change the results")
    p.add_argument("--sigma_evolve", type=int, default=1,
                   help="1: -calc_sigma_dot -limit_sigma in run.opts (normal stress evolves)")
    p.add_argument("--sigma_interact", type=int, default=1,
                   help="with --sigma_evolve 1: 1 normal stress also changes from slip on the faults, "
                        "0 only from the far-field rate (rsf_solve -sigma_dot_interact 0)")
    p.add_argument("--sigma_min", type=float, default=5e6, help="normal stress limiter floor [Pa]")
    p.add_argument("--sigma_max", type=float, default=200e6, help="normal stress limiter ceiling [Pa]")
    p.add_argument("--G", type=float, default=32.04e9)
    p.add_argument("--sigma", type=float, default=50e6)
    p.add_argument("--a", type=float, default=0.010)
    p.add_argument("--b", type=float, default=0.015)
    p.add_argument("--a_vs", type=float, default=0.020, help="a in velocity-strengthening margins")
    p.add_argument("--vs_margin", type=float, default=0.0, help="fraction of L per end that is VS")
    p.add_argument("--tip_creep", type=float, default=0.0,
                   help="length [m] of graded velocity-strengthening tips at each fault end, a from "
                        "--a_vs to --tip_a; 0: off")
    p.add_argument("--tip_a", type=float, default=0.04, help="a at the fault ends with --tip_creep")
    p.add_argument("--ab_noise", type=float, default=0.0, help="std of relative noise on a")
    p.add_argument("--ab_corr", type=float, default=0.0, help="correlation length of that noise [m], 0: white")
    p.add_argument("--n_asp", type=int, default=1, help="number of VW asperities per fault")
    p.add_argument("--barrier_frac", type=float, default=0.2,
                   help="total fraction of the inner fault length taken by VS barriers between asperities")
    p.add_argument("--asp_rand", type=float, default=0.0, help="relative spread of asperity lengths (log-normal)")
    p.add_argument("--dc", type=float, default=0.02)
    p.add_argument("--b2", type=float, default=0.0,
                   help="b of a second state variable (rsf_solve -rsf_two_state); 0: one state variable")
    p.add_argument("--dc2", type=float, default=0.0, help="dc of the second state variable [m]")
    p.add_argument("--state_law", type=int, default=0,
                   help="rsf_solve -state_law (1 aging, 2 slip, ...); 0: rsf_solve default (aging)")
    p.add_argument("--f0", type=float, default=0.6)
    p.add_argument("--v0", type=float, default=1e-6)
    p.add_argument("--v_eq", type=float, default=1e-9, help="equivalent deep slip rate for tau_dot [m/s]")
    p.add_argument("--taudot", type=float, default=0.0, help="far-field shear stressing rate [Pa/s], 0: from v_eq")
    p.add_argument("--v_init", type=float, default=0.0, help="initial slip rate [m/s], 0: v_eq")
    p.add_argument("--load", choices=("farfield", "backslip", "both"), default="farfield",
                   help="farfield (default): far-field stressing rate; backslip: uniform plate-rate "
                        "backslip v_eq on all cells (rsf_solve -vpl), no far-field rate, so that "
                        "creeping velocity-strengthening ends load the locked patch from its edges; "
                        "both: the two together")
    p.add_argument("--ic_kick", type=float, default=0.0,
                   help="add this shear stress [Pa] to the initial stress of one cell (twin runs)")
    p.add_argument("--ic_kick_cell", type=int, default=-1,
                   help="cell for --ic_kick; -1: middle cell of fault 0")
    p.add_argument("--v_pert", type=float, default=0.1, help="std of ln(v_init) perturbation per cell")
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    G = args.G
    nu = 0.25
    Lb = G * args.dc / (args.b * args.sigma)
    # Rubin and Ampuero (2005) aging-law nucleation length, in-plane (mode II) modulus G/(1-nu)
    if args.b > args.a:
        hstar = 2.0 / math.pi * G / (1.0 - nu) * args.dc * args.b / ((args.b - args.a) ** 2 * args.sigma)
    if args.b2 > 0:
        # two state variables: velocity weakening if a < b + b2; the
        # one-state h* formula is not meaningful there, an estimate from
        # the second variable's weakening (dc2, b + b2) is reported instead
        bt = args.b + args.b2
        hstar = (2.0 / math.pi * G / (1.0 - nu) * args.dc2 * bt / ((bt - args.a) ** 2 * args.sigma)
                 if bt > args.a else float("inf"))
    else:
        hstar = float("inf")
    if args.lengths:
        lengths = [float(v) for v in args.lengths.split(",")]
        if len(lengths) != args.nfault:
            raise SystemExit("--lengths needs %i values" % args.nfault)
        if args.layout != "chain":
            raise SystemExit("--lengths needs --layout chain")
    else:
        lengths = [args.L] * args.nfault
    ncell_list = []
    dx_list = []
    for Lk in lengths:
        if args.dx > 0:
            nk = max(1, int(round(Lk / args.dx)))
        else:
            nk = max(1, int(math.ceil(Lk / (Lb / args.cells_per_lb))))
        ncell_list.append(nk)
        dx_list.append(Lk / nk)
    ncell = ncell_list[0]
    dx = dx_list[0]

    taudot = args.taudot
    if taudot <= 0.0:
        taudot = G * args.v_eq / (math.pi * args.W)
    v_init = args.v_init if args.v_init > 0 else args.v_eq
    # far-field stressing-rate tensor: pure shear of magnitude tau_ref unless given
    if args.sr_s1 == 0.0 and args.sr_s3 == 0.0:
        s1, s3 = taudot, -taudot
    else:
        s1, s3 = args.sr_s1, args.sr_s3
    if args.strikes:
        strikes = [float(v) for v in args.strikes.split(",")]
        if len(strikes) != args.nfault:
            raise SystemExit("--strikes needs %i values" % args.nfault)
    else:
        strikes = [0.0] * args.nfault
    strikes_rel = list(strikes)
    strikes = [s + args.rot for s in strikes]
    step_x = args.step_x if args.stepover == "right" else -args.step_x
    kmid = 0.5 * (args.nfault - 1)
    # chain layout: fault centres before the rotation (strike 0 along +y)
    chain_centres = []
    if args.layout == "chain":
        px, py = 0.0, 0.0
        for k in range(args.nfault):
            dl = math.radians(strikes_rel[k])
            ux, uy = math.sin(dl), math.cos(dl)
            if k > 0:
                dlp = math.radians(strikes_rel[k - 1])
                if args.gap_dir == "bisect":
                    gx = math.sin(dlp) + ux
                    gy = math.cos(dlp) + uy
                    gn = math.hypot(gx, gy)
                    px += args.gap * gx / gn
                    py += args.gap * gy / gn
                else:
                    px += args.gap * math.sin(dlp)
                    py += args.gap * math.cos(dlp)
            chain_centres.append((px + 0.5 * lengths[k] * ux, py + 0.5 * lengths[k] * uy))
            px += lengths[k] * ux
            py += lengths[k] * uy
    args.sr_az += args.rot
    rot = math.radians(args.rot)
    fault_rates = []

    geom, ab, ic, srate, groups = [], [], [], [], []
    for k in range(args.nfault):
        # fault centre on the en echelon layout (built along +y, then the
        # whole layout rotated clockwise by --rot), each fault at its own strike
        xc0 = k * step_x
        yc1 = k * (args.L - args.overlap) + 0.5 * args.L
        if args.layout == "chain":
            xc0, yc1 = chain_centres[k]
        elif args.pivot == "inner" and k != kmid:
            # keep the end nearest to the middle of the system fixed
            side = 1.0 if k < kmid else -1.0
            dl = math.radians(strikes_rel[k])
            xc0 = xc0 - side * 0.5 * args.L * math.sin(dl)
            yc1 = yc1 + side * 0.5 * args.L - side * 0.5 * args.L * math.cos(dl)
        xc = xc0 * math.cos(rot) + yc1 * math.sin(rot)
        yc0 = -xc0 * math.sin(rot) + yc1 * math.cos(rot)
        phi = math.radians(strikes[k])
        sx, sy = math.sin(phi), math.cos(phi)          # along-strike unit vector
        alpha = phi - math.radians(args.sr_az)
        tau_k = 0.5 * (s1 - s3) * math.sin(2.0 * alpha)
        sig_k = 0.5 * (s1 + s3) - 0.5 * (s1 - s3) * math.cos(2.0 * alpha)
        if not args.ff_normal:
            sig_k = 0.0
        fault_rates.append((strikes[k], tau_k, sig_k))
        nk = ncell_list[k]
        dxk = dx_list[k]
        aprof = a_profile(nk, dxk, args, rng)
        for i in range(nk):
            d = (i + 0.5 - 0.5 * nk) * dxk
            geom.append((xc + d * sx, yc0 + d * sy, -0.5 * args.W, strikes[k], 90.0,
                         0.5 * dxk, 0.5 * args.W, k))
            a = float(aprof[i])
            ab.append((a, args.b))
            vi = v_init * math.exp(args.v_pert * rng.standard_normal()) if args.v_pert > 0 else v_init
            # steady state with the total state dependence b + b2
            tau = float(steady_state_tau(vi, a, args.b + args.b2, args.sigma, args.f0, args.v0))
            # signs of the initial stress and slip rate follow the drive, so
            # the fault starts at steady state in its driven sense (the
            # regularized law is odd in tau and v; rsf_solve inverts
            # psi = a ln(2 v0/v sinh(tau/(sigma a))), which needs equal signs)
            sgn = 1.0 if tau_k >= 0.0 else -1.0
            ic.append((sgn * tau, sgn * vi))
            srate.append((tau_k, sig_k))
            groups.append(k)

    # stepover type between consecutive faults from the sense of slip
    step_types = []
    for k in range(args.nfault - 1):
        if args.layout == "chain":
            step_types.append("n/a (chain)")
            continue
        s0 = fault_rates[k][1]
        s1_ = fault_rates[k + 1][1]
        if s0 * s1_ <= 0.0:
            step_types.append("mixed sense")
        elif (s0 > 0.0) == (args.stepover == "right"):
            step_types.append("restraining")
        else:
            step_types.append("releasing")

    if args.ic_kick != 0.0:
        kc = args.ic_kick_cell if args.ic_kick_cell >= 0 else ncell_list[0] // 2
        tau_k, v_k = ic[kc]
        ic[kc] = (tau_k + args.ic_kick, v_k)
    with open(os.path.join(args.out, "geom.in"), "w") as f:
        for g in geom:
            f.write("%.6e %.6e %.6e %.1f %.1f %.6e %.6e %i\n" % g)
    with open(os.path.join(args.out, "rsf.dat"), "w") as f:
        for a, b in ab:
            f.write("%.8e %.8e\n" % (a, b))
    with open(os.path.join(args.out, "ic.in"), "w") as f:
        for tau, v in ic:
            f.write("%.10e %.8e\n" % (tau, v))
    with open(os.path.join(args.out, "srate.in"), "w") as f:
        for td, sd in srate:
            if args.load == "backslip":
                td, sd = 0.0, 0.0
            f.write("%.8e %.8e\n" % (td, sd))
    with open(os.path.join(args.out, "groups.in"), "w") as f:
        for k in groups:
            f.write("%i\n" % k)

    N = len(geom)
    # closeness check between faults: cells of different faults closer
    # than two cell sizes give large, poorly resolved interaction terms
    if args.nfault > 1:
        gx = np.array([g[0] for g in geom])
        gy = np.array([g[1] for g in geom])
        gg = np.array([g[7] for g in geom])
        dmin = float("inf")
        for k in range(args.nfault):
            sk = gg == k
            for m in range(k + 1, args.nfault):
                sm = gg == m
                dd = np.hypot(gx[sk][:, None] - gx[sm][None, :], gy[sk][:, None] - gy[sm][None, :])
                dmin = min(dmin, float(dd.min()))
        if dmin < 2.0 * max(dx_list):
            print("WARNING: closest cells of different faults are %.0f m apart (cell size %.0f m)" %
                  (dmin, max(dx_list)))
    a_arr = np.array([x[0] for x in ab])
    n_vw = int(np.sum(a_arr < args.b + args.b2))
    par = vars(args).copy()
    par.update(dict(ncell_per_fault=ncell, dx=dx, N=N, N_vw=n_vw, Lb=Lb, hstar=hstar,
                    lengths=lengths, ncell_list=ncell_list, dx_list=dx_list,
                    L_over_hstar=args.L / hstar if hstar < float("inf") else 0.0,
                    taudot=taudot, taudot_kPa_per_yr=taudot * SEC_PER_YEAR / 1e3,
                    v_init=v_init, sr_s1=s1, sr_s3=s3, strikes=strikes, step_types=step_types,
                    fault_rates=[dict(strike=s, tau_dot=t, sigma_dot=n,
                                      tau_dot_kPa_per_yr=t * SEC_PER_YEAR / 1e3,
                                      sigma_dot_kPa_per_yr=n * SEC_PER_YEAR / 1e3)
                                 for s, t, n in fault_rates]))
    with open(os.path.join(args.out, "params.json"), "w") as f:
        json.dump(par, f, indent=1)
    # option fragment: everything that must match the files
    opts = ("-geom_file geom.in -rsf_file rsf.dat -rsf_ic_file ic.in "
            "-rsf_stress_rate_file srate.in -vpl %s -vel_init %.6e -rsf_slip_mode 0 "
            "-shear_modulus %.6e -sigma_init %.6e -dc %.6e -f0 %g -v0 %.3e"
            % ("0" if args.load == "farfield" else "%.6e" % args.v_eq, v_init, G, args.sigma,
               args.dc, args.f0, args.v0))
    if args.b2 > 0:
        if args.dc2 <= 0:
            raise SystemExit("--b2 needs --dc2 > 0")
        opts += " -rsf_two_state 1 -rsf_b2 %.8e -rsf_dc2 %.8e" % (args.b2, args.dc2)
    if args.state_law > 0:
        opts += " -state_law %i" % args.state_law
    if args.sigma_evolve:
        opts += (" -calc_sigma_dot 1 -limit_sigma 1 -min_sigma %.6e -max_sigma %.6e"
                 % (args.sigma_min, args.sigma_max))
        if not args.sigma_interact:
            opts += " -sigma_dot_interact 0"
    with open(os.path.join(args.out, "run.opts"), "w") as f:
        f.write(opts + "\n")
    print("N %i cells (%i per fault, dx %.0f m, %i VW), Lb %.0f m (%.1f cells), h* %.1f km, L/h* %.1f, "
          "tau_ref %.3e Pa/s = %.1f kPa/yr, S1 %.2e S3 %.2e Pa/s at az %.0f, sigma evolve %i" %
          (N, ncell, dx, n_vw, Lb, Lb / dx, hstar / 1e3, par["L_over_hstar"], taudot,
           par["taudot_kPa_per_yr"], s1, s3, args.sr_az, args.sigma_evolve))
    for k, st in enumerate(step_types):
        if args.layout == "chain":
            break
        print("  step %i-%i: %s %s step" % (k, k + 1, args.stepover, st))
    for s, t, n in fault_rates:
        print("  strike %6.1f: tau_dot %8.2f (%s) sigma_dot %8.2f kPa/yr" %
              (s, t * SEC_PER_YEAR / 1e3, "left-lateral" if t >= 0 else "right-lateral",
               n * SEC_PER_YEAR / 1e3))


if __name__ == "__main__":
    main()
