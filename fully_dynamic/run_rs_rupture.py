"""Short in-plane rate-and-state rupture (aging law) for time-scheme convergence:
20 km fault, uniform velocity weakening, steady state at 1e-9 m/s with a 2 km
overstressed patch; stations at x = 4 and 7 km. usage: run_rs_rupture.py dx [tend]
FD_NL=1 selects the second-order (Noda et al. 2020) time marching."""
import sys, os, numpy as np
import shells as sh, stepper as st
AL, BE, RHO = 6000.0, 3464.0, 2670.0
dx = float(sys.argv[1]); tend = float(sys.argv[2]) if len(sys.argv) > 2 else 7.0
L = 20e3; n = int(round(L/dx)); dt = 0.5*dx/BE; nsteps = int(round(tend/dt))
f = sh.Fault([-L/2, 0.0], [1.0, 0.0], L, n)
sys_ = sh.System([f], dt, AL, BE, RHO, tol=1e-2, lmax=nsteps + 2, cap_ok=True)
x = f.centers[:, 0]
a, b, dc, V0, f0, sig, Vi = 0.008, 0.012, 0.01, 1e-6, 0.6, 50e6, 1e-9
eta = RHO*BE/2
psi0 = f0 + b*np.log(V0/Vi)                                  # steady state at Vi
tau_ss = sig*a*np.arcsinh(Vi/(2*V0)*np.exp(psi0/a)) + eta*Vi
tau0 = np.full(n, tau_ss); tau0[np.abs(x) < 1e3] += 8e6
fr = st.RateState(np.full(n, a), np.full(n, b), dc, V0, f0)
S = st.Stepper(sys_, fr, tau0=tau0, sigma0=sig, normal_coupling=False)
S.state = np.full(n, psi0); S.V = np.full(n, Vi)
if os.environ.get("FD_NL"): S.enable_nl()
ist = {"x4": np.argmin(np.abs(x - 4e3)), "x7": np.argmin(np.abs(x - 7e3))}
rec = {k: [] for k in ist}
for m in range(nsteps):
    S.step()
    for k, i in ist.items(): rec[k].append([S.t, S.delta[i], S.V[i]])
out = {k: np.array(v) for k, v in rec.items()}
tag = "nl" if os.environ.get("FD_NL") else "cm"
np.savez(f"rs_rupture_dx{dx:g}_{tag}.npz", **out)
for k, A in out.items():
    print(f"{tag} dx {dx:g} {k}: V>0.1 at {A[np.argmax(A[:,2]>0.1),0]:.4f} s, peak {A[:,2].max():.4f} m/s, slip {A[-1,1]:.4f} m")
