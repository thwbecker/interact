"""Dunham (2007) Fig. 5: supershear transition length under Dugdale friction,
rupture nucleated from a static Dugdale crack in equilibrium plus a 20 percent
stress perturbation inside it. Lengths in L_fric = mu d0 / (tau_p - tau_r)."""
import numpy as np, sys, time
import shells as sh, stepper as st

AL, BE, RHO = 6000.0, 3464.0, 2700.0
MU = RHO * BE**2; NU = (AL**2 - 2 * BE**2) / (2 * (AL**2 - BE**2)); MUP = MU / (1 - NU)

def run(S, dx=200.0, dtau=10e6, d0=1.0, sigma=100e6, half=70e3, a_over_crit=1.0, pert=0.2, verbose=True):
    tau_f = 0.4 * sigma; tau_u = tau_f + dtau; tau0 = tau_f + dtau / (1 + S)
    Lfric = MU * d0 / dtau
    # Dugdale (BCS) critical crack: slip at the inner edge c equals d0,
    #   d0 = (4 dtau c / (pi mu')) ln sec(theta),  theta = pi/(2(1+S)),  a = c / cos(theta)
    theta = np.pi / (2 * (1 + S))
    ccrit = np.pi * MUP * d0 / (4 * dtau * np.log(1 / np.cos(theta)))
    acrit = ccrit / np.cos(theta)
    a0 = a_over_crit * acrit
    L = 2 * half; n = int(L / dx)
    f = sh.Fault([-half, 0.0], [1.0, 0.0], L, n)
    dt = 0.5 * dx / BE
    tend = (half - 3e3) / BE + 2.0
    nsteps = int(tend / dt) + 1
    sys_ = sh.System([f], dt, AL, BE, RHO, tol=1e-2, lmax=nsteps + 2, cap_ok=True)
    x = f.centers[:, 0]
    fr = st.SlipWeakening(tau_u / sigma, tau_f / sigma, d0, dugdale=True)
    # 1. static Dugdale crack (Bilby-Cottrell-Swinden geometry): slip zone |x| < a,
    #    tau = tau_f for |x| < c (slip > d0) and tau = tau_u in the cohesive zones
    #    c < |x| < a, with c = a cos(pi/(2(1+S))) for finite stress at a.
    a = a0; c = a * np.cos(np.pi / (2 * (1 + S)))
    slipset = np.abs(x) < a
    tau_target = np.where(np.abs(x) < c, tau_f, tau_u)
    d_eq = np.zeros(n)
    idx = np.where(slipset)[0]
    d_eq[idx] = np.linalg.solve(sys_.Ks[np.ix_(idx, idx)], tau_target[idx] - tau0)
    tau_eq = tau0 + sys_.Ks @ d_eq
    if verbose:
        ic = np.argmin(np.abs(np.abs(x) - c))
        print(f"S={S}: Lfric={Lfric/1e3:.2f} km, a_crit={acrit/1e3:.2f} km ({acrit/Lfric:.2f} Lfric), a={a/1e3:.2f} km, c={c/1e3:.2f} km; "
              f"static crack: min slip {d_eq[idx].min()/d0:.3f} d0, slip at c {d_eq[ic]/d0:.3f} d0, max tau outside/tau_u {tau_eq[~slipset].max()/tau_u:.3f}")
    # 2. perturb: extra stress inside the crack, then fully dynamic
    sys_.reset_history()
    S_ = st.Stepper(sys_, fr, tau0=tau0 + np.where(np.abs(x) < c, pert * (tau0 - tau_f), 0.0), sigma0=sigma, normal_coupling=False)
    S_.delta = d_eq.copy()
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
    sel = x > a0 + 3 * dx; xs = x[sel]
    res = {}
    for k in crit:
        ts = trup[k][sel]
        v = np.diff(xs) / np.diff(ts)
        i = [j for j in range(len(v) - 2) if (v[j:j + 3] > BE).all()]
        res[k] = xs[i[0]] if i else np.nan
    if verbose:
        print("   transition (3 consecutive pairs faster than beta), by criterion: " +
              ", ".join(f"{k}: {res[k]/Lfric:.1f} Lfric" for k in crit) + f"  ({time.time()-t1:.0f}s)")
    return dict(S=S, Lfric=Lfric, Ltr=res, x=xs, trup={k: trup[k][sel] for k in crit})

if __name__ == "__main__":
    S = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    run(S, dx=float(sys.argv[2]) if len(sys.argv) > 2 else 200.0, half=float(sys.argv[3]) * 1e3 if len(sys.argv) > 3 else 70e3)
