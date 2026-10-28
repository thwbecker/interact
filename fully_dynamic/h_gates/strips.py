"""H5 prototype: surface layouts for faults of different strikes.

A small dense time-domain integral solver (first order, piecewise-constant rates, square
elements, tada_quad.rect_ramp_stress) for faults with prescribed slip rate and a virtual
free surface of square cells (3 slip components each) whose rates follow from zero
traction at every step. Surface layouts are compared by the traction histories they
produce on the faults:
  ref    one uniform grid aligned with x, y over the union of the bands around the traces
  strip  for each fault a grid aligned with its own strike, cells kept where that fault's
         trace is the nearest (a Voronoi-type partition), so no surface is represented twice
  over   the per-fault strips without the partition (overlapping: surface counted twice)
  none   no surface (whole space)
usage: python3 strips.py [angle of fault 2's strike, degrees] [trace distance at the near end, m]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import tada_quad as tq
cL, cT, rho = 6000.0, 3464.0, 2670.0
dxf, hs = 2000.0, 500.0                      # fault cells 2 km, surface cells 500 m (cT dt/hs = 2, stable)
dt = 0.5*dxf/cT
ang2 = float(sys.argv[1]) if len(sys.argv) > 1 else 30.0
dnear = float(sys.argv[2]) if len(sys.argv) > 2 else 3000.0
band = 2500.0; nsteps = 20; trise = 5     # prescribed slip rate 1 m/s on both faults for 5 steps
def frame(strike_deg, dip_deg):
    s = np.radians(strike_deg); d = np.radians(dip_deg)
    e1 = np.array([np.sin(s), np.cos(s), 0.0])                     # strike (azimuth from y)
    h = np.array([np.cos(s), -np.sin(s), 0.0])                     # horizontal, perpendicular, dip direction
    e2 = -np.cos(d)*h + np.sin(d)*np.array([0, 0, 1.0])            # up dip
    e3 = np.cross(e1, e2)
    return np.array([e1, e2, e3])
def fault(origin, strike, dip, nx, nz):
    E = frame(strike, dip); cells = []
    for i in range(nx):
        for j in range(nz):
            cells.append(origin + (i + 0.5)*dxf*E[0] - (j + 0.5)*dxf*E[1])
    return dict(c=np.array(cells), E=E, a=0.5*dxf, trace=(origin, origin + nx*dxf*E[0]))
# fault 1 along y from (0, 0); fault 2 turned by ang2, its trace starting dnear off fault 1's start
F1 = fault(np.array([0.0, 0.0, 0.0]), 0.0, 60.0, 4, 2)
F2 = fault(np.array([-dnear, 0.0, 0.0]), ang2, 60.0, 4, 2)
faults = [F1, F2]
if os.environ.get("SINGLE"):
    faults = [F2]           # one fault only (the turned one): trace-aligned strip against the x-y grid
def seg_dist(p, a, b):
    ab = b - a; t = np.clip(np.dot(p - a, ab)/np.dot(ab, ab), 0, 1); return np.linalg.norm(p - (a + t*ab))
def surface_grid(E, center, half_len, half_wid):
    # cells of size hs aligned with frame E (e1 along strike), z = 0, normal +z
    cells = []
    nl, nw = int(round(2*half_len/hs)), int(round(2*half_wid/hs))
    for i in range(nl):
        for j in range(nw):
            p = center + (-half_len + (i + 0.5)*hs)*E[0] + (-half_wid + (j + 0.5)*hs)*np.array([E[0][1], -E[0][0], 0.0])
            p[2] = 0.0; cells.append(p)
    return np.array(cells)
def near_band(p):
    return min(seg_dist(p[:2], f["trace"][0][:2], f["trace"][1][:2]) for f in faults) <= band
def layout(kind):
    S = []
    Es = np.eye(3)                               # surface cells: frame x, y, z (normal +z)
    if kind == "none": return []
    if kind == "ref":
        lo = np.min([f["trace"][k][:2] for f in faults for k in (0, 1)], axis=0) - band
        hi = np.max([f["trace"][k][:2] for f in faults for k in (0, 1)], axis=0) + band
        xs = np.arange(lo[0], hi[0], hs) + 0.5*hs; ys = np.arange(lo[1], hi[1], hs) + 0.5*hs
        for x in xs:
            for y in ys:
                p = np.array([x, y, 0.0])
                if near_band(p): S.append(dict(c=p, E=Es))
        return S
    for k, f in enumerate(faults):
        mid = 0.5*(f["trace"][0] + f["trace"][1]); L = np.linalg.norm(f["trace"][1] - f["trace"][0])
        e1 = f["E"][0].copy(); e1[2] = 0; e1 /= np.linalg.norm(e1)
        Ek = np.array([e1, [-e1[1], e1[0], 0.0], [0, 0, 1.0]])
        for p in surface_grid(np.array([e1, [-e1[1], e1[0], 0], [0, 0, 1]]), mid, 0.5*L + band, band):
            if seg_dist(p[:2], f["trace"][0][:2], f["trace"][1][:2]) > band: continue
            if kind == "strip":
                dk = seg_dist(p[:2], f["trace"][0][:2], f["trace"][1][:2])
                if any(seg_dist(p[:2], g["trace"][0][:2], g["trace"][1][:2]) < dk for g in faults if g is not f): continue
            S.append(dict(c=p, E=Ek))
    return S
def tensor_from(src_c, E, a, comp, X, t):
    x = (X - src_c) @ E.T
    small = np.abs(x) < 1e-9*a; x = np.where(small, 0.0, x)
    if t <= 0: return np.zeros((len(X), 3, 3))
    S = tq.rect_ramp_stress(x[:, 0], x[:, 1], x[:, 2], np.full(len(X), t), a, a, cL, cT, rho, comp=comp)
    s = np.zeros((len(X), 3, 3))
    s[:, 0, 0], s[:, 1, 1], s[:, 2, 2] = S["s11"], S["s22"], S["s33"]
    s[:, 0, 1] = s[:, 1, 0] = S["s12"]; s[:, 0, 2] = s[:, 2, 0] = S["s31"]; s[:, 1, 2] = s[:, 2, 1] = S["s32"]
    return np.einsum("ai,nab,bj->nij", E, s, E)
def run(kind):
    surf = layout(kind)
    # sources: fault cells (dip slip, comp 2) and surface cells (3 comps: along x, y of the cell frame, opening)
    src = []
    for f in faults:
        for c in f["c"]: src.append((c, f["E"], f["a"], 2, None))
    for s in surf:
        for comp in (1, 2, 3): src.append((s["c"], s["E"], 0.5*hs, comp, None))
    nfs = sum(len(f["c"]) for f in faults); nsrc = len(src)
    # receivers: fault cells (shear along up dip, normal), surface cells (3 traction components on +z)
    rcvX = np.array([c for f in faults for c in f["c"]] + [s["c"] for s in surf]) if surf else np.array([c for f in faults for c in f["c"]])
    rcv_rows = []
    for f in faults:
        for c in f["c"]: rcv_rows.append(("f", f["E"]))
    for s in surf: rcv_rows.append(("s", None))
    nrow = 2*nfs + 3*len(surf)
    Ef = np.array([f["E"] for f in faults for c in f["c"]])           # receiver frames of the fault cells
    def rows_of(sig):
        out = np.empty(nrow)
        tf = np.einsum("nab,nb->na", sig[:nfs], Ef[:, 2])                # traction on the fault plane
        out[0:2*nfs:2] = np.einsum("na,na->n", tf, Ef[:, 1])             # shear along up dip
        out[1:2*nfs:2] = np.einsum("na,na->n", tf, Ef[:, 2])             # normal (tension positive)
        out[2*nfs:] = sig[nfs:, :, 2].reshape(-1)                        # surface: sigma_xz, sigma_yz, sigma_zz
        return out
    dS = np.zeros((nsteps + 1, nrow, nsrc))
    for j, (c, E, a, comp, _) in enumerate(src):
        if comp == 3: E_use, comp_use = None, None
        prev = np.zeros(nrow)
        for l in range(1, nsteps + 1):
            if comp == 3:     # opening of a horizontal surface cell: slip normal to the cell -> use comp 3 of tada if available
                sig = tensor_from(c, E, a, 3, rcvX, l*dt)
            else:
                sig = tensor_from(c, E, a, comp, rcvX, l*dt)
            cur = rows_of(sig); dS[l, :, j] = cur - prev; prev = cur
    # time stepping: fault rates prescribed, surface rates from zero surface traction
    nf_rows = 2*nfs
    V = np.zeros((nsteps + 1, nsrc)); trac = np.zeros((nsteps + 1, nf_rows))
    ns_src = nsrc - nfs
    if ns_src:
        A = dS[1][nf_rows:, nfs:]
        Ainv = np.linalg.inv(A)
    for n in range(1, nsteps + 1):
        V[n, :nfs] = 1.0 if n <= trise else 0.0
        known = np.zeros(nrow)
        for l in range(1, n + 1):
            known += dS[l] @ V[n - l + 1] if l > 1 else dS[1][:, :nfs] @ V[n, :nfs]
        if ns_src:
            V[n, nfs:] = -Ainv @ known[nf_rows:]
            known[:nf_rows] += dS[1][:nf_rows, nfs:] @ V[n, nfs:]
        trac[n] = known[:nf_rows]
    return trac, len(surf)
if __name__ == "__main__":
    res = {}
    for kind in ("ref", "strip", "over", "none"):
        res[kind] = run(kind)
        print(f"{kind:6s}: {res[kind][1]} surface cells", flush=True)
    ref = res["ref"][0]; sc = np.abs(ref).max()
    for kind in ("strip", "over", "none"):
        d = res[kind][0] - ref
        print(f"fault 2 strike turned {ang2:g} deg, traces {dnear/1e3:g} km apart at one end: {kind:6s} max |traction - ref| {np.abs(d).max()/sc:.3f} of the largest, shear {np.abs(d[:, 0::2]).max()/sc:.3f}, normal {np.abs(d[:, 1::2]).max()/sc:.3f}")
