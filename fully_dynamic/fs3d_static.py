"""Static check of the free-surface options for a vertical strike-slip fault
reaching the surface (3-D), with interact's Okada routine (dc3dmod.F, built
as a shared library) for all building blocks:
  exact      : Okada half space
  full       : Okada full space, no surface
  mirror     : full space, fault plus its image above the surface
  mirror+open: mirror plus a layer of opening elements on the surface plane,
               solved for zero normal traction there (the shear tractions
               vanish by symmetry)
Measure: relative Frobenius error of the fault's shear traction operator
(strike slip to strike traction) against the exact half space.
Build first:  gfortran -cpp -shared -fPIC -O2 -DUSE_DOUBLE_PRECISION \\
   -Isrc/includes -o /tmp/okd/libdc3d.so src/green/dc3dmod.F"""
import ctypes, os, numpy as np
MU, NU = 3.2e10, 0.25
LAM = 2*MU*NU/(1 - 2*NU)
_lib = ctypes.CDLL(os.environ.get("DC3D_LIB", "/tmp/okd/libdc3d.so"))
D = ctypes.c_double

def okada_stress(x, y, z, depth, dip, al, aw, disl, full):
    """stress tensor at (x, y, z) for a rectangle with reference point (0, 0, -depth),
    length al = (al1, al2) along x, width aw = (aw1, aw2) up dip, dislocation (strike, dip, tensile)"""
    alpha = (LAM + MU)/(LAM + 2*MU)
    ins = [D(alpha), D(x), D(y), D(z), D(depth), D(dip), D(al[0]), D(al[1]), D(aw[0]), D(aw[1]), D(disl[0]), D(disl[1]), D(disl[2])]
    outs = [D() for _ in range(12)]; iret = ctypes.c_int(); ifs = ctypes.c_int(1 if full else 0)
    _lib.dc3dm_(*[ctypes.byref(v) for v in ins], *[ctypes.byref(v) for v in outs], ctypes.byref(iret), ctypes.byref(ifs))
    u = [o.value for o in outs]
    # derivative order: uxx uyx uzx uxy uyy uzy uxz uyz uzz (u_i,j = d u_i / d x_j)
    G = np.array([[u[3], u[6], u[9]], [u[4], u[7], u[10]], [u[5], u[8], u[11]]])   # G[i][j] = du_i/dx_j
    e = 0.5*(G + G.T)
    return LAM*np.trace(e)*np.eye(3) + 2*MU*e, iret.value

def run(L=8e3, W=4e3, h=1e3, ext=8e3, hy_fine=250.0, w_fine=2e3, hy_out=1e3, margin=4e3, hx_s=None, verbose=True):
    hx_s = h if hx_s is None else hx_s
    nx = int(L/h); nz = int(W/h)
    xf = -L/2 + (np.arange(nx) + 0.5)*h; zf = -(np.arange(nz) + 0.5)*h
    F = [(x, z) for x in xf for z in zf]; N = len(F)
    def fault_src(i, X, full, mirror=False):
        xc, zc = F[i]
        if mirror: zc = -zc
        # vertical rectangle centred at (xc, 0, zc): reference at its bottom edge, width h up dip
        return okada_stress(X[0] - xc, X[1], X[2], -(zc - h/2), 90.0, (-h/2, h/2), (0.0, h), (1.0, 0.0, 0.0), full)[0]
    def tau_f(S): return S[0, 1]                 # strike traction on the plane y = 0
    Kex = np.zeros((N, N)); Kfu = np.zeros((N, N)); Kim = np.zeros((N, N))
    for j in range(N):
        for i, (x, z) in enumerate(F):
            X = (x, 0.0, z)
            Kex[i, j] = tau_f(fault_src(j, X, False)); Kfu[i, j] = tau_f(fault_src(j, X, True))
            Kim[i, j] = tau_f(fault_src(j, X, True, mirror=True))
    # surface opening layer: columns as the fault plus margins, rows across strike (fine strip
    # of width 2 w_fine around the trace, y = 0 on a row boundary, coarse rows out to +-ext)
    xs = -L/2 - margin + (np.arange(int(round((L + 2*margin)/hx_s))) + 0.5)*hx_s
    yb = np.concatenate([np.arange(-ext, -w_fine, hy_out), np.arange(-w_fine, w_fine, hy_fine), np.arange(w_fine, ext + 1e-6, hy_out)])
    yc = 0.5*(yb[:-1] + yb[1:]); dy = np.diff(yb)
    S = [(x, y, d) for x in xs for (y, d) in zip(yc, dy)]; M = len(S)
    def surf_src(k, X):
        xc, yc_, d = S[k]
        # horizontal rectangle (dip 0) at z = 0, x in xc +- h/2, y in yc +- d/2, tensile dislocation
        return okada_stress(X[0] - xc, X[1] - (yc_ - d/2), X[2], 0.0, 0.0, (-hx_s/2, hx_s/2), (0.0, d), (0.0, 0.0, 1.0), True)[0]
    A = np.zeros((M, M)); B = np.zeros((M, N)); C = np.zeros((N, M))
    for k in range(M):
        for m, (x, y, d) in enumerate(S):
            A[m, k] = surf_src(k, (x, y, 0.0))[2, 2]
        for i, (x, z) in enumerate(F):
            C[i, k] = tau_f(surf_src(k, (x, 0.0, z)))
    for j in range(N):
        for m, (x, y, d) in enumerate(S):
            B[m, j] = fault_src(j, (x, y, 0.0), True)[2, 2] + fault_src(j, (x, y, 0.0), True, mirror=True)[2, 2]
    Dop = np.linalg.solve(A, -B)
    Kmo = Kfu + Kim + C @ Dop
    nrm = np.linalg.norm(Kex)
    errs = dict(full=np.linalg.norm(Kfu - Kex)/nrm, mirror=np.linalg.norm(Kfu + Kim - Kex)/nrm, mirror_open=np.linalg.norm(Kmo - Kex)/nrm)
    if verbose:
        print(f"fault {L/1e3:g} x {W/1e3:g} km, h {h:g} m ({N}); surface {M} opening elements ({len(xs)} x {len(yc)}, dx {hx_s:g} m, fine dy {hy_fine:g} m, +-{ext/1e3:g} km): "
              + ", ".join(f"{k} {v:.2e}" for k, v in errs.items()))
    return errs

if __name__ == "__main__":
    run()
