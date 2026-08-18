#!/usr/bin/env python3
"""
Convert the simplified PD06/PD07/PD13 fault curves written by
make_faults_cascadia.py into interact 2D half-plane geometry input, one
constant-displacement element per row,

    x  y  0  strike  90  half_length  0  0  group

with x horizontal, y = -depth (negative down) and strike = 90 + local
dip in degrees, the same convention makefault uses in test_twod/run_test.

The group ID in the last column is 0 for the plate interface and N for
splay N, splays numbered seaward to landward.  A companion <out>.idx
file carries element_index group_id fault_name for building the per
fault rate and state property files.

Usage: faults2geom_cascadia.py [PD06|PD07|PD13]

Parameters are plain assignments below; no environment variables.
"""
import sys
import numpy as np

profile = sys.argv[1] if len(sys.argv) > 1 else "PD06"
infile = "faults_%s.d" % profile
outgeom = "geom_%s.in" % profile
outplot = "geom_%s_overview.png" % profile
make_plot = 1     # write outplot (needs matplotlib); 0: off
dx = 0.025
#dx = 0.1          # target element length [km].  Set by the SPLAYS, not by
                  # the interface or the nucleation length: the shortest
                  # splay is 2.07 km (PD13 splay4), so 0.1 km is what it
                  # takes to give every splay at least twenty elements.
                  # The generator prints the per-splay element count, so
                  # check it after changing anything upstream that alters
                  # splay length.  For reference L_b = G D_c / (b sigma)
                  # is 2.0 km at D_c = 0.05 m and sigma 50 MPa, so the
                  # interface is heavily over-resolved at this dx; 0.25
                  # would be ample for the interface alone
min_splay_elements = 20
                  # warn if any splay ends up with fewer elements than
                  # this.  A splay carrying only a handful of elements
                  # cannot resolve slip on itself whatever its friction
min_depth = 0.0   # clamp shallow tips [km]; make_faults_cascadia.py
                  # already trims them, so this should stay a no-op
scale = 1000.0    # km -> m.  interact's kernels scale as G/length, so
                  # geometry in km with SI material properties weakens
                  # the elastic coupling by 1000 and stretches every
                  # recurrence interval by the same factor.  Keep this
                  # at 1000.0 and keep coord_unit_km = 0.001 in the rsf
                  # input generator
monotone_x = 0    # the PD06 and PD07 splays dip seaward while the
                  # interface dips landward, so x is not monotone in a
                  # common sense across the model; each curve is already
                  # single valued along its own arclength, so leave off

faults = {}
name = None
for ln in open(infile):
    ln = ln.strip()
    if not ln or ln.startswith('#'):
        continue
    if ln.startswith('>'):
        name = ln[1:].strip()
        faults[name] = []
        continue
    x, d = map(float, ln.split()[:2])
    faults[name].append((x, d))

rows = []
idx = []
for name, pts in faults.items():
    p = np.array(pts)
    p[:, 1] = np.maximum(p[:, 1], min_depth)
    seg = np.hypot(*np.diff(p, axis=0).T)
    p = p[np.hstack([True, seg > 1e-9])]
    if monotone_x:
        sgn = 1.0 if p[-1, 0] >= p[0, 0] else -1.0
        kept = [p[0]]
        for q in p[1:]:
            if sgn*(q[0] - kept[-1][0]) > 1e-6:
                kept.append(q)
        p = np.array(kept)
    s = np.hstack([0.0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
    n = max(int(round(s[-1]/dx)), 1)
    si = np.linspace(0.0, s[-1], n + 1)
    xi = np.interp(si, s, p[:, 0])
    di = np.interp(si, s, p[:, 1])
    gid = 0 if name == "megathrust" else (int(name[5:]) if name.startswith("splay") else 99)
    for k in range(n):
        x0, x1, d0, d1 = xi[k], xi[k+1], di[k], di[k+1]
        if x1 < x0:                       # orient +x
            x0, x1, d0, d1 = x1, x0, d1, d0
        rows.append((0.5*(x0+x1)*scale, -0.5*(d0+d1)*scale,
                     90.0 + np.degrees(np.arctan2(d1-d0, x1-x0)),
                     0.5*np.hypot(x1-x0, d1-d0)*scale, gid))
        idx.append((name, gid))

with open(outgeom, 'w') as f:
    for cx, cy, strike, hl, gid in rows:
        f.write(f" {cx:.10e} {cy:.10e}  0.0000000000e+00 {strike:10.6f}  90.000000 "
                f"{hl:.10e}  0.0000000000e+00 {gid:10d}\n")
with open(outgeom + ".idx", 'w') as f:
    for i, (name, gid) in enumerate(idx):
        f.write(f"{i} {gid} {name}\n")
counts = {}
for name, gid in idx:
    counts[name] = counts.get(name, 0) + 1
print(f"wrote {outgeom} ({len(rows)} elements, dx = {dx} km, scale = {scale}) and {outgeom}.idx")
for name, c in counts.items():
    print(f"  {name}: {c} elements")
sp = {k: v for k, v in counts.items() if k != "megathrust"}
if sp:
    worst = min(sp, key=sp.get)
    print(f"  shortest splay is {worst} with {sp[worst]} elements", end="")
    print("" if sp[worst] >= min_splay_elements else
          f"  <-- below the {min_splay_elements} wanted; reduce dx to about "
          f"{dx*sp[worst]/min_splay_elements:.3f} km")

# ---------------------------------------------------------------- plot
# the elements as rsf_solve will see them, plus dip against depth, which
# is where shallow-tip damage and fitting artefacts show up first
if make_plot:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping " + outplot)
    else:
        R = np.array([(r[0], r[1], r[2], r[3]) for r in rows])/np.array([scale, scale, 1.0, scale])
        names = np.array([q[0] for q in idx])
        cols = plt.cm.tab10(np.arange(10))
        fig = plt.figure(figsize=(12, 7.5))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.2])
        axc = fig.add_subplot(gs[0, :])
        axz = fig.add_subplot(gs[1, 0])
        axd = fig.add_subplot(gs[1, 1])
        for k, nm in enumerate(counts):
            s = names == nm
            cx, cy, st, hl = R[s, 0], R[s, 1], R[s, 2], R[s, 3]
            th = np.radians(st - 90.0)
            xs = np.vstack([cx - hl*np.cos(th), cx + hl*np.cos(th)])
            ys = np.vstack([cy + hl*np.sin(th), cy - hl*np.sin(th)])
            for a in (axc, axz):
                a.plot(xs, ys, "-", color=cols[k % 10], lw=1.6)
            axc.plot([], [], "-", color=cols[k % 10], lw=2, label=f"{nm} ({int(s.sum())})")
            axd.plot(st - 90.0, -cy, ".", ms=2.5, color=cols[k % 10])
        for a in (axc, axz):
            a.axhline(0.0, color="tab:blue", lw=1.0)
            a.set_aspect("equal")
            a.set_xlabel("x [km]")
            a.set_ylabel("y [km]")
        axc.legend(fontsize=7, loc="lower left", ncol=6)
        axc.set_title(f"{outgeom}: {len(rows)} elements, dx = {dx} km")
        sp = names != "megathrust"
        if sp.any():
            axz.set_xlim(R[sp, 0].min() - 4, R[sp, 0].max() + 4)
            axz.set_ylim(R[sp, 1].min() - 4, 2)
            axz.set_title("splay detail, free surface in blue")
        axd.invert_yaxis()
        axd.set_xlabel("local dip [deg]")
        axd.set_ylabel("depth [km]")
        axd.set_title("element dip against depth")
        axd.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(outplot, dpi=130)
        print("wrote " + outplot)
