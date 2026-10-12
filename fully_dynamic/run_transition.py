"""Supershear transition length versus S (Andrews 1976b; Dunham 2007)."""
import numpy as np, sys, time
import shells as sh, stepper as st

AL, BE, RHO = 6000.0, 3464.0, 2700.0
MU = RHO * BE**2; NU = (AL**2 - 2 * BE**2) / (2 * (AL**2 - BE**2)); MUP = MU / (1 - NU)

def run(S, dx=100.0, dtau=10e6, d0=0.2, sigma=100e6, rc_over_Lc=1.5, vnuc=0.5, t0=0.3, half=40e3, verbose=True):
    tau_f = 0.4 * sigma
    tau_u = tau_f + dtau
    tau0 = tau_f + dtau / (1 + S)
    Lc = MUP * dtau * d0 / (np.pi * (tau0 - tau_f) ** 2)
    L = 2 * half; n = int(L / dx)
    f = sh.Fault([-half, 0.0], [1.0, 0.0], L, n)
    dt = 0.5 * dx / BE
    tend = (half - 2e3) / BE + 1.0
    nsteps = int(tend / dt) + 1
    sys_ = sh.System([f], dt, AL, BE, RHO, tol=1e-2, lmax=nsteps + 2, cap_ok=True)
    x = f.centers[:, 0]
    rc = rc_over_Lc * Lc
    ft = np.full(n, np.inf); m = np.abs(x) < rc; ft[m] = np.abs(x[m]) / (vnuc * BE)
    fr = st.SlipWeakening(tau_u / sigma, tau_f / sigma, d0, force_time=ft, force_t0=t0)
    S_ = st.Stepper(sys_, fr, tau0=tau0, sigma0=sigma, normal_coupling=False)
    crit = {"V>1e-3": lambda V, d: V > 1e-3, "slip>0.1d0": lambda V, d: d > 0.1 * d0,
            "slip>0.5d0": lambda V, d: d > 0.5 * d0, "slip>d0": lambda V, d: d > d0}
    trup = {k: np.full(n, np.inf) for k in crit}
    t1 = time.time()
    while S_.t < tend:
        V = S_.step()
        for k, fn in crit.items():
            new = fn(V, S_.delta) & ~np.isfinite(trup[k]); trup[k][new] = S_.t
        if np.isfinite(trup["slip>d0"][-1]):
            break
    sel = x > rc + 2 * dx; xs = x[sel]
    res = {}
    for k in crit:
        ts = trup[k][sel]
        v = np.diff(xs) / np.diff(ts)
        i = [j for j in range(len(v) - 2) if (v[j:j + 3] > BE).all()]
        res[k] = xs[i[0]] if i else np.nan
    if verbose:
        # Andrews (1976b) uses L_c(Andrews) = 2 L_c as defined here (eq. in his section 2 with lambda = mu)
        print(f"S={S}: Lc={Lc/1e3:.2f} km (Andrews' L_c = {2*Lc/1e3:.2f} km), rc={rc/1e3:.2f} km; transition by criterion: " +
              ", ".join(f"{k}: {res[k]/Lc:.1f} Lc = {res[k]/(2*Lc):.1f} Lc_A" for k in crit) + f"  ({time.time()-t1:.0f}s)")
    return dict(S=S, Lc=Lc, Ltr=res, x=xs, trup={k: trup[k][sel] for k in crit})

if __name__ == "__main__":
    S = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    half = float(sys.argv[2]) * 1e3 if len(sys.argv) > 2 else 40e3
    rc = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5
    out = run(S, half=half, rc_over_Lc=rc)
    np.savez(f"transition_S{S:g}_rc{rc:g}.npz", **out)
