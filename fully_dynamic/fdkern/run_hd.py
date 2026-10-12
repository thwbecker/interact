"""Harris and Day (1993) stepover experiments with the space-time BIEM prototype.
Right-lateral convention (positive slip); rupture nucleated at the centre of
fault 1 and considered as it reaches the +x end. Fault 2 on the -y side is a
right step (dilational for right-lateral slip), on the +y side a left step
(compressional). Parameters from their Table 1 (cases A-D)."""
import numpy as np, time, sys
import shells as sh, stepper as st

CASES = {"A": (0.75, 0.30), "B": (0.645, 0.51), "C": (1.10, 0.30), "D": (0.75, 0.51)}

def run(case="D", W=1000.0, overlap=5e3, dx=100.0, tend=12.0, rc=2.5e3, vnuc=0.5, t0=0.2,
        tol=1e-2, verbose=True):
    AL, BE, RHO = 6000.0, 3464.0, 2700.0
    mus, mud = CASES[case]
    L = 28e3; n = int(L / dx)
    f1 = sh.Fault([-L / 2, 0.0], [1.0, 0.0], L, n, "f1")
    f2 = sh.Fault([L / 2 - overlap, W], [1.0, 0.0], L, n, "f2")
    dt = 0.5 * dx / BE
    nsteps = int(round(tend / dt))
    tb = time.time()
    sys_ = sh.System([f1, f2], dt, AL, BE, RHO, tol=tol, lmax=nsteps + 2, cap_ok=True)
    x = np.concatenate([f1.centers[:, 0], f2.centers[:, 0]])
    on2 = np.zeros(sys_.N, bool); on2[n:] = True
    # forced nucleation on fault 1 within rc of the centre, front speed vnuc*beta
    ft = np.full(sys_.N, np.inf)
    inner = (~on2) & (np.abs(x) < rc)
    ft[inner] = np.abs(x[inner]) / (vnuc * BE)
    fr = st.SlipWeakening(mus, mud, 0.10, force_time=ft, force_t0=t0)
    S = st.Stepper(sys_, fr, tau0=20e6, sigma0=33.3e6, normal_coupling=True)
    if verbose:
        print(f"case {case} W={W/1e3:+.1f} km overlap={overlap/1e3:.1f} km dx={dx} N={sys_.N} steps={nsteps} L={sys_.L} tables {time.time()-tb:.1f}s mem {sys_.nbytes/1e6:.0f} MB")
    rec = []
    t_end1 = None; t_trig2 = None; x_trig2 = None
    tb = time.time()
    for k in range(nsteps):
        S.step()
        d1 = S.delta[:n]; d2 = S.delta[n:]
        if t_end1 is None and d1[-1] > 1e-2:
            t_end1 = S.t
        if t_trig2 is None and (S.V[n:] > 1e-3).any():
            t_trig2 = S.t; x_trig2 = f2.centers[np.argmax(S.V[n:] > 1e-3), 0]
        if k % 100 == 0:
            rec.append((S.t, (d1 > 1e-2).sum() * dx / 1e3, (d2 > 1e-2).sum() * dx / 1e3, S.V.max(),
                        (S.sigma[n:]).min() / 1e6, (S.sigma[n:]).max() / 1e6))
            if verbose:
                print(f"  t={S.t:5.2f} ruptured f1 {rec[-1][1]:5.1f} km, f2 {rec[-1][2]:5.1f} km, maxV {rec[-1][3]:.2f}, sigma_eff on f2 [{rec[-1][4]:.1f},{rec[-1][5]:.1f}] MPa ({time.time()-tb:.0f}s)")
    len2 = (S.delta[n:] > 1e-2).sum() * dx / 1e3
    out = dict(case=case, W=W, overlap=overlap, dx=dx, t_end1=t_end1, t_trig2=t_trig2, x_trig2=x_trig2,
               len2=len2, final_slip1=S.delta[:n].copy(), final_slip2=S.delta[n:].copy(), x1=f1.centers[:, 0], x2=f2.centers[:, 0],
               rec=np.array(rec))
    return out

if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "D"
    W = float(sys.argv[2]) * 1e3 if len(sys.argv) > 2 else 1000.0
    dx = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0
    out = run(case=case, W=W, dx=dx)
    print(f"RESULT case {case} W={W/1e3:+.1f} km: fault1 end reached at {out['t_end1']} s; fault 2 first slip at {out['t_trig2']} s, x={None if out['x_trig2'] is None else out['x_trig2']/1e3} km; fault 2 ruptured length {out['len2']:.1f} km")
    np.savez(f"hd_{case}_W{int(W)}_dx{int(dx)}.npz", **{k: v for k, v in out.items() if v is not None})
