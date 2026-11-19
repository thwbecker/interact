"""TPV205-2D (SCEC, 2-D slice of TPV5) with the space-time BIEM prototype."""
import numpy as np, time, sys
import shells as sh, stepper as st

def run(dx=50.0, tend=12.0, tol=1e-2, verbose=True, coupled=True):
    AL, BE, RHO = 6000.0, 3464.0, 2670.0
    L = 30e3; n = int(L / dx)
    f = sh.Fault([-L / 2, 0.0], [1.0, 0.0], L, n)
    dt = 0.5 * dx / BE
    nsteps = int(round(tend / dt))
    t0 = time.time()
    import os
    sys_ = sh.System([f], dt, AL, BE, RHO, tol=tol, lmax=nsteps + 2, cap_ok=True, et=float(os.environ.get('FD_ET', '1.0')))
    if verbose: print(f"dx={dx} n={n} dt={dt*1e3:.2f} ms steps={nsteps} L={sys_.L} tables {time.time()-t0:.1f}s mem {sys_.nbytes/1e6:.0f} MB")
    x = f.centers[:, 0]
    tau0 = np.full(n, 70e6)
    tau0[np.abs(x) < 1.5e3] = 81.6e6
    tau0[(x > -9e3) & (x < -6e3)] = 78e6
    tau0[(x > 6e3) & (x < 9e3)] = 62e6
    sigma0 = 120e6
    fr = st.SlipWeakening(0.677, 0.525, 0.40)
    S = st.Stepper(sys_, fr, tau0=tau0, sigma0=sigma0, normal_coupling=False)
    if os.environ.get('FD_NL'): S.enable_nl()
    S.coupled = coupled
    stations = {"st-045": -4.5e3, "st-075": -7.5e3, "st120": 12.0e3}
    idx = {k: int(np.argmin(np.abs(x - v))) for k, v in stations.items()}
    rec = {k: [] for k in stations}
    t0 = time.time()
    for k in range(nsteps):
        S.step()
        for name, i in idx.items():
            rec[name].append((S.t, S.delta[i], S.V[i], S.tau[i]))
        if verbose and k % 400 == 0:
            print(f"  t={S.t:5.2f} s  max V={S.V.max():.3f} m/s  ruptured {np.sum(S.delta>0.01)*dx/1e3:.1f} km  ({time.time()-t0:.0f}s)")
    out = {k: np.array(v) for k, v in rec.items()}
    out["x"] = x; out["delta"] = S.delta.copy(); out["dx"] = dx
    return out

if __name__ == "__main__":
    dx = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
    coupled = (sys.argv[2] != "pc") if len(sys.argv) > 2 else True
    out = run(dx, coupled=coupled)
    import os
    np.savez(f"tpv205_2d_dx{int(dx)}{'' if coupled else '_pc'}_et{os.environ.get('FD_ET', '1.0')}{'_nl' if os.environ.get('FD_NL') else ''}.npz", **out)
