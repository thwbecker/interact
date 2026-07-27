#!/usr/bin/env python3
"""
Write faults_uphoff.d for the Uphoff, May & Gabriel (2023 GJI, ggac467)
section 8.1 splay scenario, from the exact curves of their splays.geo.

NumPy only; gmsh is not needed.  gmsh's OpenCASCADE "Spline" through
three points is a quadratic Bezier that interpolates all three, with the
middle point placed at the chord-length parameter
t* = |P1P2| / (|P1P2| + |P2P3|).  That form was read off the OCC BREP
record for these curves (degree 2, three poles, knots of multiplicity 3
at 0 and at |P1P2|+|P2P3|) and then compared against gmsh's own sampling
of all four splays, agreeing to about 6e-11 m, i.e. to roundoff.  If a
future gmsh or OpenCASCADE release changes what "Spline" means for three
points, this closed form would no longer follow the .geo, so the
--check option below re-runs the comparison when gmsh is importable.

This replaces the earlier hand-built polyline, which had two defects:

  1. the deep splay tips fell short of the exact curve end by 0.18 to
     0.65 km, growing with splay length, and
  2. the shallow ends started at depth 0, so the min_depth clamp in
     faults2geom_uphoff.py flattened the first few nodes instead of
     trimming them, leaving a horizontal element under the free surface
     on splay1 and dips of 10 to 45 degrees where 20 to 55 are correct.

Here the shallow end is trimmed along the true curve, so
faults2geom_uphoff.py can be run with min_depth = 0.0 and its clamp
becomes a no-op.

Parameters are plain assignments below; no environment variables.
"""
import sys
import numpy as np

outfile     = "faults_uphoff.d"
dip         = 20.0      # main normal fault dip [deg], their "dip"
splay_dip1  = 50.0      # their splay_dip1
splay_dip2  = 40.0      # their splay_dip2
f1          = 0.3       # their Macro Splay f1
f2          = 7.0       # their Macro Splay f2 [km], the intentional gap
offsets     = [30.0, 50.0, 70.0, 90.0]   # splay surface offsets [km]
H           = 60.0      # main fault depth [km], their H
min_depth   = 0.2       # trim the shallow tip to this depth [km]; the
                        # half-plane kernels need y < 0 and an element
                        # centre a few hundred m below the free surface
ds          = 0.25      # polyline sampling along the curve [km]; keep at
                        # or below the target element length dx of
                        # faults2geom_uphoff.py so no curve detail is lost
main_extend = 0.0       # extra down-dip length of the main fault beyond
                        # H [km along dip], to emulate their fault
                        # extension towards 1000 km depth (0.0: off).
                        # Their .geo continues the main fault to the
                        # domain corner; that segment creeps and carries
                        # part of the loading


def occ_spline3(P1, P2, P3, n=20001):
    """quadratic Bezier through P1, P2, P3 as OCC builds it (see docstring),
    sampled densely in its own parameter"""
    d1 = float(np.hypot(*(P2 - P1)))
    d2 = float(np.hypot(*(P3 - P2)))
    t = d1/(d1 + d2)
    C = (P2 - (1.0 - t)**2*P1 - t**2*P3)/(2.0*(1.0 - t)*t)
    tt = np.linspace(0.0, 1.0, n)[:, None]
    B = (1.0 - tt)**2*P1 + 2.0*(1.0 - tt)*tt*C + tt**2*P3
    return B[:, 0], B[:, 1]


def trim_shallow(x, d, zmin):
    """cut a polyline (x, depth) at its first crossing of depth = zmin,
    interpolating the crossing, and drop everything shallower"""
    if d[0] >= zmin:
        return x, d
    k = int(np.argmax(d >= zmin))
    if k == 0:
        return x, d
    w = (zmin - d[k-1])/(d[k] - d[k-1])
    return (np.hstack([x[k-1] + w*(x[k] - x[k-1]), x[k:]]),
            np.hstack([zmin, d[k:]]))


def resample(x, d, step):
    s = np.hstack([0.0, np.cumsum(np.hypot(np.diff(x), np.diff(d)))])
    n = max(int(round(s[-1]/step)), 2)
    si = np.linspace(0.0, s[-1], n + 1)
    return np.interp(si, s, x), np.interp(si, s, d)


def control_points(off):
    """their Macro Splay, with det/l2 written out; for dip 20 and
    splay_dip2 40 this gives l2 = off exactly"""
    det = (np.sin(np.radians(dip))*np.cos(np.radians(splay_dip2))
           - np.cos(np.radians(dip))*np.sin(np.radians(splay_dip2)))
    l2 = -np.sin(np.radians(dip))*off/det
    P1 = np.array([off, 0.0])
    P2 = np.array([off + f1*l2*np.cos(np.radians(splay_dip1)),
                   f1*l2*np.sin(np.radians(splay_dip1))])
    P3 = np.array([off + (l2 - f2)*np.cos(np.radians(splay_dip2)),
                   (l2 - f2)*np.sin(np.radians(splay_dip2))])
    return P1, P2, P3


if "--check" in sys.argv:
    try:
        import gmsh
    except ImportError:
        sys.exit("--check needs the gmsh module (pip install gmsh); the "
                 "generator itself does not")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("check")
    worst = 0.0
    for off in offsets:
        P1, P2, P3 = control_points(off)
        a = gmsh.model.occ.addPoint(P1[0], -P1[1], 0.0)
        b = gmsh.model.occ.addPoint(P2[0], -P2[1], 0.0)
        c = gmsh.model.occ.addPoint(P3[0], -P3[1], 0.0)
        cu = gmsh.model.occ.addSpline([a, b, c])
        gmsh.model.occ.synchronize()
        t0, t1 = gmsh.model.getParametrizationBounds(1, cu)
        uu = np.linspace(t0[0], t1[0], 5001)
        xyz = np.array(gmsh.model.getValue(1, cu, uu)).reshape(-1, 3)
        bx, bd = occ_spline3(P1, P2, P3, n=5001)
        e = float(np.max(np.hypot(xyz[:, 0] - bx, -xyz[:, 1] - bd)))
        worst = max(worst, e)
        print("offset %5.1f km: max |closed form - gmsh| = %.3e km" % (off, e))
    gmsh.finalize()
    print("worst %.3e km; the closed form is %s"
          % (worst, "consistent with this gmsh" if worst < 1e-9
             else "NOT consistent with this gmsh, do not use"))
    sys.exit(0 if worst < 1e-9 else 1)

faults = {}
Lmain = H/np.sin(np.radians(dip)) + main_extend      # main fault is planar
smain = np.linspace(0.0, Lmain, max(int(round(Lmain/ds)), 2) + 1)
faults["megathrust"] = (smain*np.cos(np.radians(dip)),
                        smain*np.sin(np.radians(dip)))
for j, off in enumerate(offsets):
    P1, P2, P3 = control_points(off)
    bx, bd = occ_spline3(P1, P2, P3)
    faults["splay%d" % (j + 1)] = resample(bx, bd, ds)

with open(outfile, "w") as f:
    f.write("# Uphoff, May & Gabriel (2023 GJI, ggac467) splay fault SEAS scenario\n")
    f.write("# geometry from the exact curves of their splays.geo (make_faults_uphoff.py):\n")
    f.write("# main normal fault planar, dip %g deg, to %g km depth" % (dip, H))
    f.write(", plus %g km along dip.\n" % main_extend if main_extend > 0 else ".\n")
    f.write("# splays: quadratic Bezier through (off,0), %g along the %g deg direction,\n"
            % (f1, splay_dip1))
    f.write("# and %g km short of the main fault along the %g deg chord; offsets %s km.\n"
            % (f2, splay_dip2, " ".join("%g" % o for o in offsets)))
    f.write("# shallow tips trimmed at %g km depth along the true curve, sampled at %g km.\n"
            % (min_depth, ds))
    f.write("# columns: x_km depth_km (depth positive down)\n")
    for name in ["megathrust"] + ["splay%d" % (j + 1) for j in range(len(offsets))]:
        x, d = faults[name]
        x, d = trim_shallow(np.asarray(x), np.asarray(d), min_depth)
        f.write("> %s\n" % name)
        for xx, dd in zip(x, d):
            f.write("%10.5f %10.5f\n" % (xx, dd))
        dips = np.degrees(np.arctan2(np.diff(d), np.diff(x)))
        print("%-11s n=%5d  depth %6.3f-%7.3f km  dip %5.2f-%5.2f deg  L=%8.3f km"
              % (name, len(x), d.min(), d.max(), dips.min(), dips.max(),
                 np.sum(np.hypot(np.diff(x), np.diff(d)))))
print("wrote %s ; now run faults2geom_uphoff.py with min_depth = 0.0" % outfile)
