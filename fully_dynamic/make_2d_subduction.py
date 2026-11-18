#!/usr/bin/env python3
"""Idealized 2-D (in-plane, plane strain) subduction-zone sections for rsf_solve: a megathrust
whose dip increases with depth, optionally with splay faults branching from it towards the
surface, written in the input format of the Cascadia sections in test_twod/cascadia (geom_*.in,
geom_*.in.idx, rsf_ab.in, rsf_dc.in, rsf_sigma.in, rsf_vpl.in, rsf_ic.in), so that
run_2d_cascadia.sh (and coarsen_2d.py, check_2d_resolution.py) can use them directly.

Frame: x landward, y up, the free surface y = 0; the megathrust starts at the trench, x = 0, at
depth z_top (a buried top: the exact half-plane kernels allow a top at the surface too, but the
reflected tables are cheapest and most accurate for buried tops, fd_plan.md M10).

Megathrust: the dip increases from dip0 at the top to dip1 at depth zkink, linearly with depth,
and stays at dip1 below, down to depth zmax; elements of about h m (resample with coarsen_2d.py,
or use h directly).
Friction (rate and state, aging law): b constant; a - b = ab_vs (velocity strengthening) above
z_vw0 and below z_vw1, ab_vw (weakening) between, with linear tapers of width taper on both
sides; Dc constant. Effective normal stress sigma = min(sigma_max, sigma_grad z) (at least
sigma_min). Loading: vpl on the megathrust, vpl_splay on the splays (0: loaded only by the
megathrust). Initial state: steady state at vinit, tau0 = a sigma asinh(vinit/(2 v0) exp(psi/a)),
psi = f0 + b ln(v0/vinit), as in the Cascadia inputs.
Splays: --splay x_top,z_top,dip (m, m, degrees, dipping landward like the megathrust): a planar
fault from (x_top, -z_top) down to within gap (default one element) of the megathrust; its
friction is velocity strengthening (ab_splay) unless --splay_vw is given.

Usage example (megathrust and two splays, elements of 100 m):
  make_2d_subduction.py --out test_twod/ideal/IM01 --name IM01 --h 100 \
      --splay 15000,500,30 --splay 25000,500,35
Then: fully_dynamic/run_2d_cascadia.sh "q o c" 1000 32 mpirun 700 test_twod/ideal/IM01 2d_ideal_IM01_1km
"""
import argparse
import os
import numpy as np


def megathrust_line(a):
    """points (x, z) along the megathrust from z_top to zmax with the dip profile, fine steps"""
    pts = [(0.0, a.z_top)]
    ds = min(a.h/10.0, 10.0)
    x, z = 0.0, a.z_top
    while z < a.zmax:
        f = min(max((z - a.z_top)/max(a.zkink - a.z_top, 1e-9), 0.0), 1.0)
        dip = np.radians(a.dip0 + f*(a.dip1 - a.dip0))
        x += ds*np.cos(dip)
        z += ds*np.sin(dip)
        pts.append((x, z))
    return np.array(pts)


def resample(pts, h):
    """elements of nearly equal length about h along a polyline (x, z): centres, unit tangents
    (x, y up), half-lengths"""
    seg = np.diff(pts, axis=0)
    s = np.concatenate([[0.0], np.cumsum(np.hypot(seg[:, 0], seg[:, 1]))])
    n = max(1, int(round(s[-1]/h)))
    sb = np.linspace(0.0, s[-1], n + 1)
    xb = np.interp(sb, s, pts[:, 0])
    zb = np.interp(sb, s, pts[:, 1])
    cx = 0.5*(xb[1:] + xb[:-1])
    cz = 0.5*(zb[1:] + zb[:-1])
    dx = xb[1:] - xb[:-1]
    dz = zb[1:] - zb[:-1]
    ln = np.hypot(dx, dz)
    tx = dx/ln
    ty = -dz/ln                                  # y up
    return cx, -cz, tx, ty, 0.5*ln


def ab_profile(z, a):
    """a - b of the megathrust with depth: VS above z_vw0 and below z_vw1, VW between, linear tapers"""
    w = max(a.taper, 1e-9)
    up = np.clip((z - (a.z_vw0 - w))/w, 0.0, 1.0)        # 0 above the upper taper, 1 below it
    dn = np.clip(((a.z_vw1 + w) - z)/w, 0.0, 1.0)        # 1 above the lower taper, 0 below it
    f = np.minimum(up, dn)                               # 1 in the VW zone
    return a.ab_vs + f*(a.ab_vw - a.ab_vs)


def splay_line(spec, mt, a):
    """planar splay from (x_top, z_top) down at its dip to within gap of the megathrust polyline"""
    xt, zt, dip = [float(v) for v in spec.split(",")]
    d = np.radians(dip)
    gap = a.gap if a.gap > 0 else a.h
    step = min(a.h/10.0, 10.0)
    pts = [(xt, zt)]
    x, z = xt, zt
    for k in range(200000):
        x += step*np.cos(d)
        z += step*np.sin(d)
        # depth of the megathrust below this x
        zm = np.interp(x, mt[:, 0], mt[:, 1], left=np.nan, right=np.nan)
        if np.isnan(zm):
            raise SystemExit(f"make_2d_subduction: splay {spec} leaves the megathrust's x range")
        dist = np.min(np.hypot(mt[:, 0] - x, mt[:, 1] - z))
        if z >= zm or dist <= gap:
            break
        pts.append((x, z))
    if len(pts) < 2:
        raise SystemExit(f"make_2d_subduction: splay {spec} starts below or within {gap} m of the megathrust")
    return np.array(pts)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, help="output directory (a section like test_twod/cascadia/PD13_nosigma)")
    p.add_argument("--name", default="IM", help="geometry name: geom_<name>.in")
    p.add_argument("--h", type=float, default=100.0, help="element length [m]")
    p.add_argument("--z_top", type=float, default=1000.0, help="depth of the megathrust top at the trench [m]")
    p.add_argument("--zmax", type=float, default=60e3, help="depth of the megathrust bottom [m]")
    p.add_argument("--dip0", type=float, default=5.0, help="dip at the top [deg]")
    p.add_argument("--dip1", type=float, default=25.0, help="dip below zkink [deg]")
    p.add_argument("--zkink", type=float, default=30e3, help="depth where the dip reaches dip1 [m]")
    p.add_argument("--b", type=float, default=0.015)
    p.add_argument("--ab_vs", type=float, default=0.003, help="a - b, velocity strengthening")
    p.add_argument("--ab_vw", type=float, default=-0.0035, help="a - b, velocity weakening")
    p.add_argument("--z_vw0", type=float, default=4e3, help="top of the VW zone [m]")
    p.add_argument("--z_vw1", type=float, default=40e3, help="bottom of the VW zone [m]")
    p.add_argument("--taper", type=float, default=3e3, help="width of the a - b tapers [m]")
    p.add_argument("--dc", type=float, default=0.05, help="Dc [m]")
    p.add_argument("--sigma_grad", type=float, default=18e3, help="effective normal stress gradient [Pa/m]")
    p.add_argument("--sigma_max", type=float, default=50e6, help="[Pa]")
    p.add_argument("--sigma_min", type=float, default=5e6, help="[Pa]")
    p.add_argument("--vpl", type=float, default=1.2e-9, help="megathrust loading rate [m/s]")
    p.add_argument("--vpl_splay", type=float, default=0.0, help="splay loading rate [m/s]")
    p.add_argument("--vinit", type=float, default=1e-9, help="initial slip rate [m/s]")
    p.add_argument("--f0", type=float, default=0.6)
    p.add_argument("--v0", type=float, default=1e-6)
    p.add_argument("--splay", action="append", default=[], help="x_top,z_top,dip of a splay [m, m, deg]")
    p.add_argument("--ab_splay", type=float, default=0.003, help="a - b of the splays (VS)")
    p.add_argument("--splay_vw", action="store_true", help="splays with the megathrust's a - b profile")
    p.add_argument("--gap", type=float, default=0.0, help="closest approach of a splay to the megathrust [m] (default h)")
    a = p.parse_args()

    mt = megathrust_line(a)
    groups = [("megathrust", mt)]
    for k, spec in enumerate(a.splay):
        groups.append((f"splay{k + 1}", splay_line(spec, mt, a)))
    rows = []
    idx = []
    ab = []
    sig = []
    vpl = []
    names = []
    for gi, (name, pts) in enumerate(groups):
        cx, cy, tx, ty, hl = resample(pts, a.h)
        z = -cy
        strike = (90.0 - np.degrees(np.arctan2(ty, tx))) % 360.0
        if gi == 0 or a.splay_vw:
            amb = ab_profile(z, a)
        else:
            amb = np.full(len(z), a.ab_splay)
        s = np.clip(a.sigma_grad*z, a.sigma_min, a.sigma_max)
        for i in range(len(cx)):
            rows.append((cx[i], cy[i], strike[i], hl[i], gi))
            idx.append((len(rows) - 1, gi, name))
            ab.append((a.b + amb[i], a.b))
            sig.append(s[i])
            vpl.append(a.vpl if gi == 0 else a.vpl_splay)
        names.append((name, len(cx), z.min(), z.max()))
    os.makedirs(a.out, exist_ok=True)
    gfile = os.path.join(a.out, f"geom_{a.name}.in")
    with open(gfile, "w") as f:
        for (x, y, st, hl, gi) in rows:
            f.write(f"{x:18.10e} {y:18.10e} {0.0:18.10e} {st:11.6f} {90.0:10.6f} {hl:17.10e} {0.0:18.10e} {gi:10d}\n")
    with open(gfile + ".idx", "w") as f:
        for (i, gi, name) in idx:
            f.write(f"{i} {gi} {name}\n")
    ab = np.array(ab)
    sig = np.array(sig)
    psi = a.f0 + ab[:, 1]*np.log(a.v0/a.vinit)
    tau0 = ab[:, 0]*sig*np.arcsinh(a.vinit/(2.0*a.v0)*np.exp(psi/ab[:, 0]))
    np.savetxt(os.path.join(a.out, "rsf_ab.in"), ab, fmt="%.6f")
    np.savetxt(os.path.join(a.out, "rsf_dc.in"), np.full(len(rows), a.dc), fmt="%.6e")
    np.savetxt(os.path.join(a.out, "rsf_sigma.in"), sig, fmt="%.6e")
    np.savetxt(os.path.join(a.out, "rsf_vpl.in"), np.array(vpl), fmt="%.6e")
    np.savetxt(os.path.join(a.out, "rsf_ic.in"), np.stack([tau0, np.full(len(rows), a.vinit)], axis=1), fmt="%.6e")
    with open(os.path.join(a.out, "make_2d_subduction.args"), "w") as f:
        f.write(" ".join(f"--{k} {v}" for k, v in vars(a).items()) + "\n")
    print(f"make_2d_subduction: {gfile}: {len(rows)} elements of about {a.h:g} m: " +
          ", ".join(f"{n} {m} (depth {z0/1e3:.2f} to {z1/1e3:.2f} km)" for n, m, z0, z1 in names))


if __name__ == "__main__":
    main()
