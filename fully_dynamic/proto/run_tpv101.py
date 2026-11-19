"""SCEC TPV101: rate and state (aging law), whole space. usage: run_tpv101.py dx [tend]
Computational fault: x in [-21, 21] km, y in [-6, 21] km (3 km of fully
velocity-strengthening fault beyond the transition layer on all sides);
beyond it the fault is taken as locked."""
import sys, time, numpy as np
import fd3d, fd3d_rs
dx = float(sys.argv[1]) if len(sys.argv) > 1 else 500.0
tend = float(sys.argv[2]) if len(sys.argv) > 2 else 12.0
cL, cT, rho = 6000.0, 3464.0, 2670.0
xlo, xhi, ylo, yhi = -21e3, 21e3, -6e3, 21e3
nx, nz = int(round((xhi - xlo)/dx)), int(round((yhi - ylo)/dx))
dt = 0.5*dx/cT
t0 = time.time()
tab = fd3d.Tables3D(nx, nz, dx, dt, cL, cT, rho, cdtype=np.complex64)
print(f"tables {time.time() - t0:.1f} s", flush=True)
xc = xlo + (np.arange(nx) + 0.5)*dx; yc = ylo + (np.arange(nz) + 0.5)*dx
X, Y = np.meshgrid(xc, yc, indexing="ij")
W, w, da0 = 15e3, 3e3, 0.008
def B(x, Wh):
    """smooth boxcar of half-width Wh and transition width w (spec eq. 5)"""
    x = np.abs(x); out = np.zeros_like(x)
    out[x <= Wh] = 1.0
    m = (x > Wh) & (x < Wh + w); xm = x[m]
    out[m] = 0.5*(1 + np.tanh(w/(xm - Wh - w) + w/(xm - Wh)))
    return out
a = 0.008 + da0*(1 - B(X, W)*B(Y - 7.5e3, W/2))     # spec eq. 4: B(x; W, w) B(y - y0; W/2, w)
b, Lc, f0, V0 = 0.012, 0.02, 0.6, 1e-6
Vini, tauini, sig = 1e-12, 75e6, 120e6
theta0 = Lc/V0*np.exp((a*np.log(2*np.sinh(tauini/(a*sig))) - f0 - a*np.log(Vini/V0))/b)
R, T, dtau0 = 3e3, 1.0, 25e6
r2 = X**2 + (Y - 7.5e3)**2
Fr = np.where(r2 < R*R, np.exp(r2/np.minimum(r2 - R*R, -1e-300)), 0.0)
def load(t):
    G = np.exp((t - T)**2/(t*(t - 2*T))) if 0 < t < T else (1.0 if t >= T else 0.0)
    out = np.zeros((2, nx, nz)); out[0] = dtau0*Fr*G; return out
tau0 = np.zeros((2, nx, nz)); tau0[0] = tauini
f = fd3d_rs.FaultRS3D(tab, tau0, np.full((nx, nz), sig), a, b, Lc, f0, V0, theta0, Vini, load)
import os
f.fixed_rake = bool(os.environ.get("FIXED_RAKE"))
if f.fixed_rake: print("fixed rake: slip along strike only", flush=True)
stations = {"st000dp030": (0, 3e3), "st000dp075": (0, 7.5e3), "st000dp120": (0, 12e3), "st090dp075": (9e3, 7.5e3),
            "st120dp030": (12e3, 3e3), "st120dp120": (12e3, 12e3)}
def cells(x, y):
    fi = (x - xlo)/dx - 0.5; fj = (y - ylo)/dx - 0.5
    i0, j0 = int(np.floor(fi)), int(np.floor(fj)); wi, wj = fi - i0, fj - j0
    return [(i0, j0, (1 - wi)*(1 - wj)), (i0 + 1, j0, wi*(1 - wj)), (i0, j0 + 1, (1 - wi)*wj), (i0 + 1, j0 + 1, wi*wj)]
cw = {k: cells(*v) for k, v in stations.items()}
rec = {k: [] for k in stations}
rt = np.full((nx, nz), np.nan)
nsteps = int(round(tend/dt)); t0 = time.time()
for n in range(nsteps):
    its = f.step()
    vm = np.sqrt(f.V[0]**2 + f.V[1]**2)
    newly = np.isnan(rt) & (vm > 1e-3); rt[newly] = f.time
    for k, cc in cw.items():
        g = lambda A: sum(wt*A[i, j] for i, j, wt in cc)
        rec[k].append([f.time, g(f.delta[0]), g(f.V[0]), g(f.tau[0]), g(np.log10(f.theta))])
    if n % 25 == 0:
        print(f"t {f.time:6.3f}  max V {vm.max():.3e}  ruptured {np.sum(~np.isnan(rt))*dx*dx/1e6:.0f} km2  its {its}  wall {time.time()-t0:.0f} s", flush=True)
out = {k: np.array(v) for k, v in rec.items()}
for k, A in out.items():
    on = A[A[:, 2] > 1e-3, 0]
    print(f"{k}: rupture time {on[0] if len(on) else np.nan:.3f} s, peak slip rate {A[:,2].max():.3f} m/s, final slip {A[-1,1]:.3f} m")
np.savez(f"tpv101_dx{int(dx)}{'_fr' if f.fixed_rake else ''}.npz", rt=rt, xc=xc, yc=yc, delta=f.delta, **out)
