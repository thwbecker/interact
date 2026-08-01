#!/usr/bin/env python3
"""
Build simplified fault geometries for the PD06, PD07 and PD13 margin
cross sections and write them in the faults_*.d format read by
faults2geom_cascadia.py.

Input is cascadia_picks.d, the colour-separated digitization of the
three figures (seafloor, plate interface, splay faults, and a second
interpreted horizon that is carried for reference only and not used).

Three things happen here, and each is an approximation worth stating:

1. FLATTENING.  interact is a half plane with its free surface on
   y = 0, so the sloping seafloor has to become a flat line.  Since sea
   water carries no shear stress, the seafloor and not the sea surface
   is the mechanical free surface, so y = 0 is placed on the seafloor.
   The section is then rotated (not sheared, not simply shifted) by the
   best-fitting seafloor slope.  A rotation is a rigid motion: every
   fault dip and every fault length is preserved exactly, and the dips
   that come out are dips relative to the free surface, which is the
   angle the elastic problem actually depends on.  What is lost is the
   residual seafloor relief about that straight line, printed below as
   an RMS.  The straight line is fitted at the horizontal positions of
   the fault picks rather than uniformly along the profile, so the
   error is minimized where the faults actually are.

2. SMOOTHING.  The interface is replaced by a degree deg_mega
   polynomial in x and each splay by a degree deg_splay polynomial in
   depth.  The RMS misfit to the picks is printed per curve; raise the
   degree if it is large compared with the digitization scatter (a few
   tens of metres), lower it if the fitted dip at the downdip end looks
   unreasonable, since that dip is what the extension below is built
   on.

3. DOWNDIP EXTENSION.  The interface is continued from the last pick to
   ext_depth_km with a dip that ramps smoothly from the fitted end dip
   to dip_deep_deg.  This part is NOT constrained by the sections: on
   these profiles the picks stop at 8 to 13 km depth, so most of the
   extended fault is an assumption.  The printed summary reports how
   much of the final interface length is extrapolated.

Parameters are plain assignments below; no environment variables.
"""
import sys
import numpy as np

picks = "cascadia_picks.d"
profiles = ["PD06", "PD07", "PD13"]
outfmt = "faults_%s.d"
outplot = "faults_cascadia_overview.png"
make_plot = 1        # write outplot (needs matplotlib); 0: off

deg_mega = 3         # polynomial degree for interface depth against x
deg_splay = 2        # polynomial degree for splay x against depth
ds_out = 0.05        # output sampling along each curve [km]; keep at or
                     # below the element length dx of faults2geom_cascadia.py

ext_depth_km = 70.0  # continue the interface down to this depth [km].
                     # The a-b profiles put the downdip end of the
                     # velocity-weakening band at 55 km (zab_TzT.dat) or
                     # 52 km (the figure), so a 40 km fault would be
                     # weakening all the way to its tip and would
                     # nucleate on the mesh edge rather than at a
                     # frictional transition.  Set 40.0 only if you want
                     # the downdip limit to be the mesh
dip_deep_deg = 25.0  # dip the extension asymptotes to [deg]
ext_zmid_km = 20.0   # depth at which the ramp is half way [km]
ext_zwidth_km = 8.0  # ramp width in depth [km]
dip_end_override = None
                     # None: start the extension at the dip the
                     # polynomial ends on.  A number [deg] forces it,
                     # which is what the warning below is about
updip_km = 0.0       # continue the interface updip beyond the first pick
                     # by this arclength [km] along its end tangent
                     # (0.0: off; the deformation front is off-section
                     # on these profiles)

min_depth_km = 0.3   # trim splay tops to at least this depth below the
                     # flattened free surface, so no element centre ends
                     # up at or above y = 0
junction_gap_km = 0.5
                     # minimum perpendicular distance to leave between a
                     # splay's deep tip and the interface.  A splay that
                     # is closer than this, or that crosses, is pulled
                     # back along itself.  Do not set this to zero:
                     # elements that touch make the stiffness matrix
                     # near singular and the rate and state solve runs
                     # away.  It should be a few element lengths
extend_to_junction = 0
                     # 0: leave splays that already clear the interface
                     # at their drawn length, which is what the sections
                     # show.  1: also extend those splays down their end
                     # tangent until the gap equals junction_gap_km.
                     # Where a splay and the interface converge slowly
                     # (PD13) that extension is many km of invented
                     # fault, so check the printed adjustments
max_adjust_km = 6.0  # refuse to move a splay tip further than this

# --------------------------------------------------------------- input
data = {}
key = None
for ln in open(picks):
    ln = ln.strip()
    if not ln or ln.startswith("#"):
        continue
    if ln.startswith(">"):
        key = tuple(ln[1:].split())
        data[key] = []
        continue
    data[key].append(tuple(map(float, ln.split()[:2])))
for k in data:
    data[k] = np.array(data[k])


def rotate(p, ang, x0, z0):
    """rigid rotation of (x, depth) points by -ang about (x0, z0)"""
    c, s = np.cos(ang), np.sin(ang)
    dx, dz = p[:, 0] - x0, p[:, 1] - z0
    return np.column_stack([dx*c + dz*s, -dx*s + dz*c])


def dip_ramp(z, z0, dip0, dip1, zmid, zw):
    """dip that equals dip0 at z0 and tends to dip1 with depth, C1 in z"""
    t0 = np.tanh((z0 - zmid)/zw)
    f = (np.tanh((z - zmid)/zw) - t0)/(1.0 - t0)
    return dip0 + (dip1 - dip0)*np.clip(f, 0.0, 1.0)


def resample(p, ds):
    s = np.hstack([0.0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
    n = max(int(round(s[-1]/ds)), 2)
    si = np.linspace(0.0, s[-1], n + 1)
    return np.column_stack([np.interp(si, s, p[:, 0]), np.interp(si, s, p[:, 1])])


def perp_gap(pt, curve):
    """signed perpendicular distance from pt to a polyline, positive above"""
    d = np.hypot(curve[:, 0] - pt[0], curve[:, 1] - pt[1])
    k = int(np.argmin(d))
    return d[k]*(1.0 if curve[k, 1] > pt[1] else -1.0)


out = {}
for prof in profiles:
    surf = data[(prof, "surface")]
    mega = data[(prof, "interface")]
    splays = [data[k] for k in sorted(data) if k[0] == prof and k[1].startswith("splay")]
    splays.sort(key=lambda p: p[0, 0])

    # -- 1. flattening: fit the seafloor where the faults are, then rotate
    xf = np.hstack([mega[:, 0]] + [s[:, 0] for s in splays])
    xf = xf[(xf >= surf[0, 0]) & (xf <= surf[-1, 0])]
    zf = np.interp(xf, surf[:, 0], surf[:, 1])
    m, b = np.polyfit(xf, zf, 1)
    ang = np.arctan(m)
    x0, z0 = 0.0, b
    rms_flat = float(np.sqrt(np.mean((zf - (m*xf + b))**2)))
    surf_r = rotate(surf, ang, x0, z0)
    mega_r = rotate(mega, ang, x0, z0)
    splays_r = [rotate(s, ang, x0, z0) for s in splays]

    # -- 2. smoothing: polynomial interface, evaluated over the pick span
    cm = np.polyfit(mega_r[:, 0], mega_r[:, 1], deg_mega)
    rms_mega = float(np.sqrt(np.mean((mega_r[:, 1] - np.polyval(cm, mega_r[:, 0]))**2)))
    xs = np.linspace(mega_r[0, 0], mega_r[-1, 0], 400)
    fit = np.column_stack([xs, np.polyval(cm, xs)])
    dip_end = float(np.degrees(np.arctan(np.polyval(np.polyder(cm), xs[-1]))))
    if dip_end_override is not None:
        dip_end = float(dip_end_override)
    # a robust comparison value: straight line through the last 20 km
    tail = mega_r[mega_r[:, 0] > mega_r[-1, 0] - 20.0]
    dip_rob = float(np.degrees(np.arctan(np.polyfit(tail[:, 0], tail[:, 1], 1)[0])))
    if abs(dip_end - dip_rob) > 3.0:
        print("  %s WARNING: the degree %d fit ends at %.1f deg but a straight line "
              "through the last 20 km gives %.1f deg; the extension is built on the "
              "former, so consider a lower degree" % (prof, deg_mega, dip_end, dip_rob))

    # -- 3. downdip extension, and optionally updip
    seg = []
    x, z, step = fit[-1, 0], fit[-1, 1], 0.05
    while z < ext_depth_km:
        dd = np.radians(dip_ramp(z, fit[-1, 1], dip_end, dip_deep_deg,
                                 ext_zmid_km, ext_zwidth_km))
        x += step*np.cos(dd)
        z += step*np.sin(dd)
        seg.append((x, z))
    ext = np.array(seg)
    if updip_km > 0.0:
        t = fit[1] - fit[0]
        t = t/np.hypot(*t)
        up = fit[0] - np.outer(np.arange(step, updip_km + step, step), t)[::-1]
        fit = np.vstack([up, fit])
    curve = np.vstack([fit, ext])
    L_tot = float(np.sum(np.hypot(*np.diff(curve, axis=0).T)))
    L_ext = float(np.sum(np.hypot(*np.diff(ext, axis=0).T)))

    # -- splays: polynomial in depth, top trimmed, deep tip gapped
    sp_out = []
    rms_sp = []
    adjust = []
    for s in splays_r:
        # fit in the splay's own chord frame: u along the tip-to-tip
        # chord, v perpendicular to it.  Fitting x against depth instead
        # puts the vertex of the polynomial inside the depth range
        # whenever a splay steepens upward, and the curve then doubles
        # back in x and reports dips of 50 to 70 degrees where the picks
        # show 15 to 40
        e = s[-1] - s[0]
        e = e/np.hypot(*e)
        nrm = np.array([-e[1], e[0]])
        u = (s - s[0]) @ e
        v = (s - s[0]) @ nrm
        cs = np.polyfit(u, v, deg_splay)
        rms_sp.append(float(np.sqrt(np.mean((v - np.polyval(cs, u))**2))))
        uu = np.linspace(u[0], u[-1], 200)
        q = s[0] + np.outer(uu, e) + np.outer(np.polyval(cs, uu), nrm)
        keep = q[:, 1] >= min_depth_km        # trim the top, do not clamp it
        if keep.sum() > 4:
            q = q[keep]
        # walk the deep tip along its own end tangent, one small step at
        # a time, until the perpendicular gap to the interface is right.
        # Stepping rather than solving keeps this correct whatever the
        # angle between the splay and the interface
        t = q[-1] - q[-5]
        t = t/np.hypot(*t)
        step_j, moved = 0.05, 0.0
        g0 = perp_gap(q[-1], curve)
        if g0 < junction_gap_km:                # too close, or crossing
            sl = np.hstack([0.0, np.cumsum(np.hypot(*np.diff(q, axis=0).T))])
            while (len(q) > 4 and moved < max_adjust_km
                   and perp_gap(q[-1], curve) < junction_gap_km):
                sl = np.hstack([0.0, np.cumsum(np.hypot(*np.diff(q, axis=0).T))])
                q = q[sl < sl[-1] - step_j]
                moved -= step_j
        elif extend_to_junction:
            while moved < max_adjust_km and perp_gap(q[-1], curve) > junction_gap_km:
                q = np.vstack([q, q[-1] + step_j*t])
                moved += step_j
        adjust.append(moved)
        sp_out.append(resample(q, ds_out))
    out[prof] = dict(surf=surf_r, picks=mega_r, curve=resample(curve, ds_out),
                     splays=sp_out, sp_picks=splays_r, ang=ang, rms_flat=rms_flat)

    with open(outfmt % prof, "w") as f:
        f.write("# %s: simplified fault geometry for interact (make_faults_cascadia.py)\n" % prof)
        f.write("# section rotated by %.3f deg so the best-fitting seafloor is flat;\n"
                "# y = 0 is that seafloor, residual relief %.3f km RMS.\n" % (np.degrees(ang), rms_flat))
        f.write("# interface: degree %d polynomial (%.3f km RMS to the picks) to %.1f km\n"
                "# distance, then extended to %.1f km depth with the dip ramping from\n"
                "# %.1f to %.1f deg.  The extension is an assumption, not a pick.\n"
                % (deg_mega, rms_mega, fit[-1, 0], ext_depth_km, dip_end, dip_deep_deg))
        f.write("# splays: degree %d polynomials, numbered seaward to landward, deep tips\n"
                "# held %.2f km clear of the interface.\n" % (deg_splay, junction_gap_km))
        f.write("# columns: distance_km depth_km (depth positive down)\n")
        f.write("> megathrust\n")
        for xx, zz in out[prof]["curve"]:
            f.write("%10.4f %10.4f\n" % (xx, zz))
        for k, q in enumerate(sp_out):
            f.write("> splay%d\n" % (k + 1))
            for xx, zz in q:
                f.write("%10.4f %10.4f\n" % (xx, zz))

    gaps = [perp_gap(q[-1], out[prof]["curve"]) for q in sp_out]
    print("%s  rotated %+.2f deg, seafloor residual %.3f km RMS" % (prof, np.degrees(ang), rms_flat))
    print("   interface: %d picks, degree %d fit %.3f km RMS, end dip %.1f deg" % (len(mega_r), deg_mega, rms_mega, dip_end))
    print("   interface: %.1f km long, %.1f km (%.0f%%) of it extrapolated below the deepest pick at %.1f km"
          % (L_tot, L_ext, 100*L_ext/L_tot, fit[-1, 1]))
    tops = [q[0, 1] for q in sp_out]
    print("   %d splays, fit %.3f to %.3f km RMS, lengths %.1f to %.1f km, tops at %.2f to %.2f km depth"
          % (len(sp_out), min(rms_sp), max(rms_sp),
             min(np.sum(np.hypot(*np.diff(q, axis=0).T)) for q in sp_out),
             max(np.sum(np.hypot(*np.diff(q, axis=0).T)) for q in sp_out),
             min(tops), max(tops)))
    print("   junction gaps %.2f to %.2f km after tip adjustments of %+.2f to %+.2f km"
          % (min(gaps), max(gaps), min(adjust), max(adjust)))
    print("   wrote " + outfmt % prof)

# ---------------------------------------------------------------- plot
if make_plot:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping " + outplot)
    else:
        fig, ax = plt.subplots(len(profiles), 2, figsize=(14, 3.4*len(profiles)),
                               gridspec_kw=dict(width_ratios=[1.7, 1.0]))
        ax = np.atleast_2d(ax)
        for i, prof in enumerate(profiles):
            o = out[prof]
            for j in (0, 1):
                a = ax[i, j]
                a.axhline(0.0, color="tab:blue", lw=1.5)
                a.plot(o["surf"][:, 0], o["surf"][:, 1], color="tab:blue", lw=0.8, alpha=0.6)
                a.plot(o["picks"][:, 0], o["picks"][:, 1], "k.", ms=1.5, alpha=0.5)
                for q in o["sp_picks"]:
                    a.plot(q[:, 0], q[:, 1], ".", color="0.6", ms=1.5)
                a.plot(o["curve"][:, 0], o["curve"][:, 1], "-", color="k", lw=1.8)
                for k, q in enumerate(o["splays"]):
                    a.plot(q[:, 0], q[:, 1], "-", color=plt.cm.tab10((k + 1) % 10), lw=1.6)
                a.invert_yaxis()
                a.set_aspect("equal")
                a.set_ylabel("depth below seafloor [km]")
                a.grid(alpha=0.25)
            ax[i, 0].set_title("%s: rotated %+.2f deg, seafloor residual %.2f km RMS"
                               % (prof, np.degrees(o["ang"]), o["rms_flat"]))
            xr = o["curve"][:, 0]
            ax[i, 1].set_xlim(min(q[:, 0].min() for q in o["splays"]) - 5,
                              max(q[:, 0].max() for q in o["splays"]) + 5)
            ax[i, 1].set_ylim(max(q[:, 1].max() for q in o["splays"]) + 3, -1)
            ax[i, 1].set_title("splay detail (dots: picks, lines: fits)")
            ax[i, 0].set_xlabel("distance [km]")
            ax[i, 1].set_xlabel("distance [km]")
        fig.tight_layout()
        fig.savefig(outplot, dpi=130)
        print("wrote " + outplot)
