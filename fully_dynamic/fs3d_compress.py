"""How compressible are the fault-to-layer coupling tables of the mirror +
opening-layer free surface (fd_plan.md, stage A)? Builds W[l] for strike slip
on a vertical fault (plus its image) to the normal traction s33 on the
surface layer, Toeplitz along strike (offsets on the layer lattice), FFT
along strike, and reports the numerical rank of each (fault row) x (layer
row) block per wavenumber and lag at a relative tolerance."""
import sys, numpy as np, time
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tada_quad as tq
AL, BE, RHO = 6000.0, 3464.0, 2670.0

def build(Lx=20e3, W=10e3, h=1e3, hs=250.0, fine=2e3, ext=20e3, h_out=1e3, margin=10e3, L=None):
    nz = int(W/h); zc = -(np.arange(nz) + 0.5)*h
    # layer rows across strike: fine near the trace, graded (h_out) beyond; y = 0 on a boundary
    yb = np.unique(np.concatenate([np.arange(-ext, -fine, h_out), np.arange(-fine, fine + 1e-6, hs), np.arange(fine, ext + 1e-6, h_out)]))
    yc = 0.5*(yb[:-1] + yb[1:]); ny = len(yc)
    nxl = int((Lx + 2*margin)/hs); offs = (np.arange(-(nxl - 1), nxl))*hs     # receiver - source along strike
    dt = 0.5*h/BE
    rmax = np.hypot(Lx + 2*margin, np.hypot(ext, W)); Lw = int(np.ceil(rmax/BE/dt)) + 2 if L is None else L
    times = np.arange(Lw + 1)*dt; Tl = 3*rmax/BE; tl = np.array([Tl, 2*Tl])
    tt = np.concatenate([times, tl])
    Wt = np.zeros((Lw, len(offs), nz, ny))
    t0 = time.time()
    X, T = np.meshgrid(offs, tt, indexing="ij")
    for iz in range(nz):
        for iy in range(ny):
            s = 0.0
            for z_src in (zc[iz], -zc[iz]):                        # fault element and its image
                # receiver (x, y, 0) in the source frame e1 = x, e2 = up, e3 = -y; s33(global) = s22(local)
                r = tq.rect_ramp_stress(X, np.full_like(X, -z_src), -yc[iy], T, h/2, h/2, AL, BE, RHO, comp=1)
                s = s + r["s22"]
            ks = (s[:, -1] - s[:, -2])/(tl[1] - tl[0])
            K = np.diff(s[:, :Lw + 1], axis=1)
            Wt[:, :, iz, iy] = (K - ks[:, None]*dt).T
    return Wt, dict(nz=nz, ny=ny, noff=len(offs), L=Lw, dt=dt, build_s=time.time() - t0)

def ranks(Wt, tol=1e-4):
    L, noff, nz, ny = Wt.shape
    M = int(2**np.ceil(np.log2(noff)))
    F = np.fft.rfft(Wt, n=M, axis=1)                                # (L, M/2+1, nz, ny)
    gmax = np.abs(F).max()
    full = 0; comp = 0; rk = []
    for l in range(L):
        for k in range(F.shape[1]):
            B = F[l, k]; full += nz*ny
            sv = np.linalg.svd(B, compute_uv=False)
            r = int(np.sum(sv > tol*gmax))                        # absolute tolerance relative to the global maximum
            rk.append(r); comp += min(r*(nz + ny), nz*ny)
    rk = np.array(rk)
    return dict(ratio=comp/full, rank_mean=rk.mean(), rank_max=rk.max(), frac_zero=np.mean(rk == 0))

if __name__ == "__main__":
    Wt, info = build()
    print(info)
    for tol in (1e-3, 1e-4, 1e-5, 1e-6):
        r = ranks(Wt, tol)
        print(f"tol {tol:g}: storage ratio {r['ratio']:.3f}, mean rank {r['rank_mean']:.1f} (max {r['rank_max']}), blocks below tolerance {r['frac_zero']:.2f}")
