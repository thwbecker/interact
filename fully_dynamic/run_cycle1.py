"""Single planar in-plane fault with a velocity-weakening patch: QD/FD cycle test."""
import numpy as np, time, sys
import shells as sh, stepper as st, cycle as cy

def setup(dx=100.0, Lf=30e3, Lvw=16e3, V_sw=1e-3, V_back=1e-4, tol=1e-2, lmax=None, dynamic=True):
    AL, BE, RHO = 6000.0, 3464.0, 2700.0
    n = int(Lf / dx)
    f = sh.Fault([-Lf / 2, 0.0], [1.0, 0.0], Lf, n)
    dt = 0.5 * dx / BE
    sys_ = sh.System([f], dt, AL, BE, RHO, tol=tol, lmax=lmax, cap_ok=lmax is not None)
    x = f.centers[:, 0]
    a = np.where(np.abs(x) < Lvw / 2, 0.010, 0.020); b = np.full(n, 0.015)
    dc = 0.016; V0 = 1e-6; f0 = 0.6; Vpl = 1e-9; sigma0 = 50e6
    fr = st.RateState(a, b, dc, V0, f0)
    psi_ss = f0 + b * np.log(V0 / Vpl)
    # start slightly above steady state in the VW centre to break symmetry in time
    Vinit = Vpl * np.where(np.abs(x) < 2e3, 1.5, 1.0)
    tau0 = fr.tau_of_V(Vinit, psi_ss, sigma0)
    model = cy.CycleModel(sys_, fr, tau0, sigma0, Vpl, normal_coupling=True, V_sw=V_sw, V_back=V_back, dynamic=dynamic, verbose=False)
    y0 = np.concatenate([np.zeros(n), psi_ss])
    return model, y0, x

if __name__ == "__main__":
    V_sw = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3
    nev = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    dynamic = (sys.argv[3] != "qd") if len(sys.argv) > 3 else True
    model, y0, x = setup(V_sw=V_sw, dynamic=dynamic)
    print(f"N={model.N} L={model.sys.L} dt={model.sys.dt*1e3:.2f} ms")
    t0 = time.time()
    y, t = model.run(y0, t_end=400 * 3.15576e7, max_events=nev)
    print(f"done in {time.time()-t0:.0f}s; events: {len(model.events)}")
    for e in model.events:
        print(f"  event at {e['t_start']/86400:.3f} d, duration {e['t_end']-e['t_start']:.1f} s ({e['nsteps']} steps), Vmax {e['Vmax']:.2f} m/s, mean slip {e['slip'].mean():.3f} m, max slip {e['slip'].max():.3f} m, ruptured {np.sum(e['slip']>1e-2)*100/1e3:.1f} km")
    np.savez(f"cycle1_Vsw{V_sw:g}{'' if dynamic else '_qd'}.npz", x=x, events_slip=np.array([e["slip"] for e in model.events]),
             t_start=np.array([e["t_start"] for e in model.events]), log=np.array([(l[0], l[1]) for l in model.log]))
