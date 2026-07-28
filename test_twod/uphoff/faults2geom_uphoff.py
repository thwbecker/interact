#!/usr/bin/env python3
"""
Convert labeled distance-depth fault polylines (faults_zelst.d format:
"> name" headers followed by "distance_km depth_km" lines, depth
positive down) into interact 2D half-plane geometry input as used in
test_twod: one constant-displacement element per row,
    x  y  0  strike  90  half_length  0  0
with x horizontal [km -> model units], y = -depth (negative down), and
strike = 90 + local dip (degrees, clockwise from horizontal, measured
along the local element direction with positive x orientation), the
same convention makefault uses in test_twod/run_test.

The geom.in group ID (last column) encodes the fault: megathrust = 0,
splay N = N (splays numbered by trench distance of their megathrust
junction). Also writes <out>.idx with one line per element:
element_index group_id fault_name, for convenience when building per
fault rsf property files.

Parameters are plain assignments below; no environment variables.
"""
import sys
import numpy as np

infile = "faults_uphoff.d"
outgeom = "geom_uphoff.in"
outplot = "faults_uphoff_overview.png"
make_plot = 1     # write outplot (needs matplotlib); 0: off
dx = 0.25         # target element length [km]; their production res_f 0.25
smooth_km = 0.0   # exact spline samples, no digitization noise
                  # digitization wiggle from the skeleton trace (0: off)
end_trim_km = 0.0
                  # (skeleton endpoints carry bulb artifacts)
min_depth = 0.0   # clamp shallow fault tips to at least this depth [km]
                  # (digitization can put a tip a fraction of a km above
                  # the surface; the half-plane needs y < 0)
scale = 1000.0    # km -> m.  interact's kernels scale as G/length, so the
                  # geometry must be in the same length unit as -shear_modulus
                  # [Pa], vpl [m/s] and D_c [m]; km geometry with SI inputs
                  # weakens the elastic coupling by 1000x and stretches every
                  # recurrence time by the same factor
monotone_x = 1    # enforce monotonically increasing x along each fault
                  # (true for this geometry; removes digitization
                  # double-backs at junctions)
r_junc = 0.0      # splays intentionally do NOT touch the main fault (7 km gap along the 40 deg chord); never snap
                  # splay root, the skeleton trace weaves; the megathrust
                  # is bridged smoothly through the zone and each splay is
                  # extended along its last clean tangent to an exact
                  # intersection with the repaired megathrust

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

def fine_resample(p, ds=0.1):
    s0 = np.hstack([0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
    if s0[-1] < 4*ds or len(p) < 3:
        return p
    sf = np.linspace(0, s0[-1], max(int(s0[-1]/ds), 4))
    return np.column_stack([np.interp(sf, s0, p[:,0]), np.interp(sf, s0, p[:,1])])

# junction-zone repair (see r_junc above)
if r_junc > 0 and "megathrust" in faults and len(faults) > 1:
    mega = np.array(faults["megathrust"])
    roots = {nm: np.array(pts[-1]) for nm, pts in faults.items() if nm != "megathrust"}
    mega = fine_resample(mega)
    keep = np.ones(len(mega), bool)
    for rt in roots.values():
        keep &= np.hypot(mega[:,0]-rt[0], mega[:,1]-rt[1]) > r_junc*1.2
    mega = mega[keep]                      # bridge gaps by interpolation
    faults["megathrust"] = [tuple(q) for q in mega]
    for nm in roots:
        p = fine_resample(np.array(faults[nm]))
        rt = roots[nm]
        sa = np.hstack([0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
        cut = sa[-1] - r_junc
        p = p[sa < cut] if cut > 2.0 else p[:max(3, len(p)//2)]
        # tangent from the last 2 km
        sb = np.hstack([0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
        tail = p[sb > sb[-1]-2.0]
        t = tail[-1] - tail[0]
        t = t/np.hypot(*t)
        # intersect the ray p[-1] + a*t with the megathrust polyline
        best = None
        for k in range(len(mega)-1):
            a, bpt = mega[k], mega[k+1]
            M = np.array([[t[0], a[0]-bpt[0]], [t[1], a[1]-bpt[1]]])
            if abs(np.linalg.det(M)) < 1e-12:
                continue
            al, be = np.linalg.solve(M, a - p[-1])
            if al > 0 and -0.05 <= be <= 1.05:
                if best is None or al < best[0]:
                    best = (al, p[-1] + al*t)
        if best is not None and best[0] < 3.0*r_junc:
            p = np.vstack([p, best[1]])
        faults[nm] = [tuple(q) for q in p]

rows = []
idx = []
for name, pts in faults.items():
    p = np.array(pts)
    p[:,1] = np.maximum(p[:,1], min_depth)
    seg = np.hypot(*np.diff(p, axis=0).T)
    keep = np.hstack([True, seg > 1e-9])
    p = p[keep]
    if monotone_x:
        sgn = 1.0 if p[-1,0] >= p[0,0] else -1.0
        kept = [p[0]]
        for q in p[1:]:
            if sgn*(q[0] - kept[-1][0]) > 1e-6:
                kept.append(q)
        p = np.array(kept)
    if end_trim_km > 0 and len(p) > 4:
        s0 = np.hstack([0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
        if s0[-1] > 4*end_trim_km:
            sf = np.linspace(end_trim_km, s0[-1]-end_trim_km, max(int(s0[-1]/0.1),4))
            p = np.column_stack([np.interp(sf, s0, p[:,0]), np.interp(sf, s0, p[:,1])])
    if smooth_km > 0 and len(p) > 4:
        # resample finely, then boxcar along arclength
        s0 = np.hstack([0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
        fine = max(int(s0[-1]/0.1), 4)
        sf = np.linspace(0, s0[-1], fine)
        xf = np.interp(sf, s0, p[:,0]); df = np.interp(sf, s0, p[:,1])
        w = max(int(smooth_km/0.1), 1)
        ker = np.ones(w)/w
        xs = np.convolve(xf, ker, mode='same'); ds = np.convolve(df, ker, mode='same')
        # boxcar edge bias: keep the exact original endpoints
        xs[0], ds[0] = xf[0], df[0]; xs[-1], ds[-1] = xf[-1], df[-1]
        m = w//2
        if len(xs) > 2*m+2:
            xs[1:m+1] = np.linspace(xs[0], xs[m+1], m+2)[1:-1]
            ds[1:m+1] = np.linspace(ds[0], ds[m+1], m+2)[1:-1]
            xs[-m-1:-1] = np.linspace(xs[-m-2], xs[-1], m+2)[1:-1]
            ds[-m-1:-1] = np.linspace(ds[-m-2], ds[-1], m+2)[1:-1]
        p = np.column_stack([xs, ds])
    s = np.hstack([0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
    L = s[-1]
    n = max(int(round(L/dx)), 1)
    si = np.linspace(0, L, n+1)
    xi = np.interp(si, s, p[:,0])
    di = np.interp(si, s, p[:,1])
    # group ID: megathrust = 0, splayN = N
    if name == "megathrust":
        gid = 0
    elif name.startswith("splay"):
        gid = int(name[5:])
    else:
        gid = 99
    for k in range(n):
        x0, x1 = xi[k], xi[k+1]
        d0, d1 = di[k], di[k+1]
        if x1 < x0:               # orient +x
            x0, x1, d0, d1 = x1, x0, d1, d0
        cx = 0.5*(x0+x1)*scale
        cy = -0.5*(d0+d1)*scale
        hl = 0.5*np.hypot(x1-x0, d1-d0)*scale
        dip_loc = np.degrees(np.arctan2(d1-d0, x1-x0))
        strike = 90.0 + dip_loc
        rows.append((cx, cy, strike, hl, gid))
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
print(f"wrote {outgeom} ({len(rows)} elements, dx = {dx} km) and {outgeom}.idx")
for name, c in counts.items():
    print(f"  {name}: {c} elements")

# ---------------------------------------------------------------- plot
# quick look at what was written: the elements as they will be handed to
# rsf_solve, and their local dip against depth.  The dip panel is the one
# that catches shallow-tip damage, since a clamped rather than trimmed
# polyline shows up there as a few elements at the wrong dip near the
# surface while the cross section still looks fine at this scale
if make_plot:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping " + outplot)
    else:
        R = np.array([(r[0], r[1], r[2], r[3], r[4]) for r in rows])
        names = np.array([q[0] for q in idx])
        cols = plt.cm.tab10(np.arange(10))
        fig = plt.figure(figsize=(12, 7.5))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25])
        axc = fig.add_subplot(gs[0, :])
        axd = fig.add_subplot(gs[1, 0])
        axz = fig.add_subplot(gs[1, 1])
        for k, name in enumerate(counts):
            s = names == name
            cx, cy, st, hl = R[s, 0], R[s, 1], R[s, 2], R[s, 3]
            th = np.radians(st - 90.0)
            xseg = np.vstack([cx - hl*np.cos(th), cx + hl*np.cos(th)])
            yseg = np.vstack([cy + hl*np.sin(th), cy - hl*np.sin(th)])
            for q in (axc, axz):
                q.plot(xseg, yseg, "-", color=cols[k], lw=1.6)
            axc.plot([], [], "-", color=cols[k], lw=2,
                     label=f"{name} ({int(s.sum())})")
            axd.plot(st - 90.0, -cy, ".", ms=2.5, color=cols[k])
        for q in (axc, axz):
            q.axhline(0.0, color="k", lw=0.8)
            q.set_aspect("equal")
            q.set_xlabel("x [km]")
            q.set_ylabel("y [km]")
        axc.set_title(f"{outgeom}: {len(rows)} elements, dx = {dx} km")
        axc.legend(fontsize=7, loc="lower left", ncol=5)
        # zoom on the shallowest splay tip and the gap it leaves to the
        # main fault, which is the part of the geometry that is easiest
        # to break and hardest to see at full scale
        s1 = names == "splay1"
        if s1.any():
            th1 = np.radians(R[s1, 2] - 90.0)
            ex = R[s1, 0] + R[s1, 3]*np.cos(th1)
            ey = R[s1, 1] - R[s1, 3]*np.sin(th1)
            m = int(np.argmax(ex))
            axz.set_xlim(ex[m] - 6.0, ex[m] + 6.0)
            axz.set_ylim(ey[m] - 5.0, ey[m] + 5.0)
            # perpendicular distance from the splay tip to the planar main
            # fault, their f2 = 7 km along the 40 deg chord
            gap = abs(ey[m] + ex[m]*np.tan(np.radians(20.0)))*np.cos(np.radians(20.0))
            axz.set_title(f"splay1 deep tip, perpendicular gap to main fault {gap:.3f} km")
        axd.invert_yaxis()
        axd.set_xlabel("local dip [deg]")
        axd.set_ylabel("depth [km]")
        axd.set_title("element dip against depth")
        axd.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(outplot, dpi=140)
        print(f"wrote {outplot}")
