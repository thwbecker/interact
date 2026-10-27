"""Gate 3: low-rank compressibility of a fault -> surface block in Ando's domains.
Source: 16 x 8 squares of 500 m of a thrust dipping 30 degrees (dip slip), top at 1 km depth.
Receivers: 32 x 32 surface cells of 250 m (8 x 8 km) either right above the fault's top
(near) or 10 km away (far). Traction on the surface plane (s13, s23, s33 global) per unit
slip rate of each source cell. Quantities: static slope (Domain S), cubic coefficients
between the fronts (Domain I), wavefront strips (Domain F, values per step after the strip
start); singular-value ranks for relative Frobenius tolerances, and the storage against
the present dense window."""
import sys, numpy as np
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), ".."))
import tada_quad as tq
cL, cT, rho, dxf, dxs = 6000.0, 3464.0, 2670.0, 500.0, 250.0
a = 0.5*dxf; dt = 0.5*dxf/cT; dip = np.radians(30.0); ztop = -1000.0
e1 = np.array([1.0, 0, 0]); e2 = np.array([0, np.cos(dip), np.sin(dip)]); e3 = np.cross(e1, e2)   # e2 up dip
E = np.array([e1, e2, e3])
src = np.array([[ (i+0.5)*dxf - 4000.0, -(j+0.5)*dxf*np.cos(dip), ztop - (j+0.5)*dxf*np.sin(dip)] for i in range(16) for j in range(8)])
def receivers(x0, y0):
    return np.array([[x0 + (i+0.5)*dxs, y0 + (j+0.5)*dxs, 0.0] for i in range(32) for j in range(32)])
def tensor(X, t):
    # global stress tensor at points X (n x 3) relative to a source cell at the origin of frame E, comp 2 (dip slip)
    x = X @ E.T
    S = tq.rect_ramp_stress(x[:,0], x[:,1], x[:,2], t, a, a, cL, cT, rho, comp=2)
    s = np.zeros((X.shape[0], 3, 3))
    s[:,0,0], s[:,1,1], s[:,2,2] = S["s11"], S["s22"], S["s33"]
    s[:,0,1] = s[:,1,0] = S["s12"]; s[:,0,2] = s[:,2,0] = S["s31"]; s[:,1,2] = s[:,2,1] = S["s32"]
    return np.einsum("ai,nab,bj->nij", E, s, E)        # E^T s E
def block(rcv):
    nr, ns = len(rcv), len(src)
    out = {}
    nmax = 0
    for j, c in enumerate(src):
        X = rcv - c; r = np.linalg.norm(X, axis=1)
        rmin = np.maximum(r - np.sqrt(2)*a, 0.0); rmax = r + np.sqrt(2)*a
        lP0 = np.floor(rmin/cL/dt).astype(int); lP1 = np.ceil(rmax/cL/dt).astype(int) + 1
        lS0 = np.floor(rmin/cT/dt).astype(int); lS1 = np.ceil(rmax/cT/dt).astype(int) + 1
        nmax = max(nmax, lS1.max())
        T1, T2 = 60.0, 120.0
        Ks = (tensor(X, np.full(nr, T2))[:, :, 2] - tensor(X, np.full(nr, T1))[:, :, 2])/(T2 - T1)  # traction on z plane: s_i3
        for comp in range(3):
            out.setdefault(("S", comp), np.zeros((nr, ns)))[:, j] = Ks[:, comp]
        # Domain F: 5 steps from each strip start (P and S)
        for name, l0 in (("FP", lP0), ("FS", lS0)):
            for k in range(5):
                t1 = (l0 + k + 1)*dt; t0 = (l0 + k)*dt
                W = tensor(X, t1)[:, :, 2] - tensor(X, np.maximum(t0, 1e-12))[:, :, 2]*(t0 > 0)[:, None] - Ks*dt
                for comp in range(3):
                    out.setdefault((name + str(k), comp), np.zeros((nr, ns)))[:, j] = W[:, comp]
        # Domain I: cubic fit on [lP1, lS0] (needs >= 5 points; else zero, the strips then cover it)
        cub = np.zeros((nr, 3, 4))
        ok = (lS0 - lP1) >= 5
        if ok.any():
            idx = np.where(ok)[0]
            for q in idx:
                tt = np.arange(lP1[q], lS0[q] + 1)*dt
                v = tensor(np.repeat(X[q:q+1], tt.size, axis=0), tt)[:, :, 2]
                for comp in range(3): cub[q, comp] = np.polyfit(tt, v[:, comp], 3)
        for comp in range(3):
            for k in range(4):
                out.setdefault(("I" + str(k), comp), np.zeros((nr, ns)))[:, j] = cub[:, comp, k]
    return out, nmax
def rank(M, tol):
    s = np.linalg.svd(M, compute_uv=False)
    if s[0] == 0: return 0
    e = np.sqrt(np.cumsum(s[::-1]**2))[::-1]/np.sqrt((s**2).sum())     # relative Frobenius error of rank-k truncation
    k = int(np.argmax(e <= tol)) if np.any(e <= tol) else len(s)
    return k
for label, (x0, y0) in (("near (over the top of the fault)", (-4000.0, -4000.0)), ("far (10 km off)", (-4000.0, 10000.0))):
    rcv = receivers(x0, y0)
    out, nmax = block(rcv)
    nr, ns = len(rcv), len(src)
    print(f"{label}: block {nr} x {ns}, present dense window about {nmax} lags")
    tot = {1e-3: 0, 1e-5: 0}
    for key in sorted(out):
        M = out[key]
        r3, r5 = rank(M, 1e-3), rank(M, 1e-5)
        tot[1e-3] += r3*(nr + ns); tot[1e-5] += r5*(nr + ns)
        if key[1] == 2: print(f"   {key[0]:4s} (normal comp): rank {r3:3d} (1e-3), {r5:3d} (1e-5) of {min(nr, ns)}")
    nq = len(out)
    dense_fdpm = nq*nr*ns; dense_now = nmax*3*nr*ns
    print(f"   storage: present dense window {dense_now:.2e}, domain-partitioned dense {dense_fdpm:.2e} ({dense_now/dense_fdpm:.1f}x less), compressed {tot[1e-3]:.2e} at 1e-3 ({dense_now/tot[1e-3]:.0f}x less), {tot[1e-5]:.2e} at 1e-5 ({dense_now/tot[1e-5]:.0f}x less)")
