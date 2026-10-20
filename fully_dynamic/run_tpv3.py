"""SCEC TPV3 (Day et al. 2005): 30 x 15 km fault in a whole space, slip weakening.
usage: python3 run_tpv3.py dx_m [tend]"""
import sys, time, numpy as np
import fd3d
dx = float(sys.argv[1]) if len(sys.argv) > 1 else 500.0
tend = float(sys.argv[2]) if len(sys.argv) > 2 else 12.0
cL, cT, rho = 6000.0, 3464.0, 2670.0
nx, nz = int(round(30e3/dx)), int(round(15e3/dx))
dt = 0.5*dx/cT
t0 = time.time()
tab = fd3d.Tables3D(nx, nz, dx, dt, cL, cT, rho, cdtype=(np.complex64 if dx < 200 else complex))
print(f"tables {time.time() - t0:.1f} s")
xc = -15e3 + (np.arange(nx) + 0.5)*dx; zc = -7.5e3 + (np.arange(nz) + 0.5)*dx
X, Z = np.meshgrid(xc, zc, indexing="ij")
tau0 = np.zeros((2, nx, nz)); tau0[0] = 70e6
nuc = (np.abs(X) < 1.5e3) & (np.abs(Z) < 1.5e3); tau0[0][nuc] = 81.6e6
f = fd3d.Fault3D(tab, tau0, np.full((nx, nz), 120e6), 0.677, 0.525, 0.40)
def stn(x, z):
    i = np.argsort(np.abs(xc - x))[:2] if abs(((x + 15e3)/dx) % 1) < 1e-9 else [np.argmin(np.abs(xc - x))]
    j = np.argsort(np.abs(zc - z))[:2] if abs(((z + 7.5e3)/dx) % 1) < 1e-9 else [np.argmin(np.abs(zc - z))]
    return [(a, b) for a in i for b in j]
st = {"PI": stn(7.5e3, 0.0), "PA": stn(0.0, 6.0e3)}
rec = {k: [] for k in st}
rt = np.full((nx, nz), np.nan)
nsteps = int(round(tend/dt)); t0 = time.time()
for n in range(nsteps):
    its = f.step()
    vm = np.sqrt(f.V[0]**2 + f.V[1]**2)
    newly = np.isnan(rt) & (vm > 1e-3); rt[newly] = f.time
    for k, cells in st.items():
        rec[k].append([f.time, np.mean([f.delta[0][c] for c in cells]), np.mean([f.V[0][c] for c in cells]),
                       np.mean([f.tau[0][c] for c in cells])])
    if n % 20 == 0:
        print(f"t {f.time:6.3f}  max V {vm.max():.3f}  ruptured {np.sum(~np.isnan(rt))*dx*dx/1e6:.0f} km2  its {its}  wall {time.time()-t0:.0f} s", flush=True)
out = {k: np.array(v) for k, v in rec.items()}
for k, a in out.items():
    on = a[a[:, 2] > 1e-3, 0]
    print(f"{k}: rupture time {on[0] if len(on) else np.nan:.3f} s, peak slip rate {a[:,2].max():.3f} m/s at {a[np.argmax(a[:,2]),0]:.3f} s, final slip {a[-1,1]:.3f} m")
np.savez(f"tpv3_dx{int(dx)}.npz", rt=rt, xc=xc, zc=zc, delta=f.delta, **out)
