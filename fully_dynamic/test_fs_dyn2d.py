"""Dynamic tests of the virtual free surface in 2-D (densesys.py).
1. Antiplane, vertical fault reaching the surface, kinematic slip (a rupture
   front running up from the bottom at 0.8 cS): fault shear traction with
   a virtual surface against the exact image solution (fault plus mirrored
   fault in a full space), both schemes."""
import sys, numpy as np, time
import densesys as ds
AL, BE, RHO = 6000.0, 3464.0, 2670.0

def vertical_fault(Df, hf, top=0.0):
    nf = int(round(Df/hf)); y = -top - (np.arange(nf) + 0.5)*hf
    return np.column_stack([np.zeros(nf), y]), np.tile([0.0, -1.0], (nf, 1)), np.full(nf, hf/2)

def surface(Ls, hs, x0=0.0):
    ns = int(round(2*Ls/hs)); x = x0 - Ls + (np.arange(ns) + 0.5)*hs
    return np.column_stack([x, np.zeros(ns)]), np.tile([1.0, 0.0], (ns, 1)), np.full(ns, hs/2)

def antiplane_case(Df=4e3, hf=250.0, Ls=8e3, hs=125.0, T=4.0, nl=True, verbose=True):
    dt = 0.5*min(hf, hs)/BE; nsteps = int(round(T/dt)); L = nsteps + 1
    Xf, sf, af = vertical_fault(Df, hf); nf = len(Xf)
    Xs, ss, as_ = surface(Ls, hs)
    vr = 0.8*BE
    def vfun_factory(yel):
        return lambda t, f: (t > (Df + yel[f])/vr).astype(float)      # front from the bottom (y = -Df) upward
    # virtual surface
    t0 = time.time()
    dv = ds.DenseSys(np.vstack([Xf, Xs]), np.vstack([sf, ss]), np.concatenate([af, as_]), np.r_[np.zeros(nf), np.ones(len(Xs))], dt, AL, BE, RHO, L, mode="ap")
    yel = np.r_[Xf[:, 1], Xs[:, 1]]
    mv = ds.Marcher(dv, nl=nl, vf_fun=vfun_factory(yel))
    # image reference: fault plus mirrored fault, both prescribed with the same slip rate
    Xi = Xf*np.array([1.0, -1.0])
    di = ds.DenseSys(np.vstack([Xf, Xi]), np.vstack([sf, sf]), np.concatenate([af, af]), np.zeros(2*nf), dt, AL, BE, RHO, L, mode="ap")
    yim = np.r_[Xf[:, 1], Xf[:, 1]]
    mi = ds.Marcher(di, nl=nl, vf_fun=vfun_factory(yim))
    tb = time.time() - t0
    ev = []; ref = []
    for n in range(nsteps):
        mv.step(); mi.step()
        ev.append(mv.tau[:nf].copy()); ref.append(mi.tau[:nf].copy())
    ev = np.array(ev); ref = np.array(ref)
    err = np.abs(ev - ref).max()/np.abs(ref).max()
    if verbose:
        print(f"antiplane {'NL' if nl else 'first order'}: fault {Df/1e3:g} km / {hf:g} m, surface +-{Ls/1e3:g} km / {hs:g} m, {nsteps} steps "
              f"(tables {tb:.0f} s): max |tau diff| / max |tau| = {err:.2e}; surface slip rate max {np.abs(mv.V[nf:]).max():.3f}")
    return err


def inplane_kinematic(dip=30.0, Lf=3e3, hf=250.0, top=500.0, Ls=6e3, hs=125.0, T=5.0, trise=1.0, nl=True, verbose=True):
    """In-plane, buried dipping fault with prescribed uniform slip rate 1 for
    0 < t < trise: (i) before the surface reflection arrives the fault traction
    equals the full-space one; (ii) at the end the traction equals the exact
    half-plane static traction of the final slip (tdd_coeff); (iii) no
    instability on the surface elements (no filter)."""
    import fs_static2d as fs
    dt = 0.5*min(hf, hs)/BE; nsteps = int(round(T/dt)); L = nsteps + 1
    ang = -np.radians(dip); sfv = np.array([np.cos(ang), np.sin(ang)])
    nf = int(round(Lf/hf)); Xf = np.array([[0.0, -top]]) + ((np.arange(nf) + 0.5)*hf)[:, None]*sfv[None, :]
    sf = np.tile(sfv, (nf, 1)); af = np.full(nf, hf/2)
    Xs, ss, as_ = surface(Ls, hs, x0=Lf*np.cos(ang)/2)
    t0 = time.time()
    dv = ds.DenseSys(np.vstack([Xf, Xs]), np.vstack([sf, ss]), np.concatenate([af, as_]), np.r_[np.zeros(nf), np.ones(len(Xs))], dt, AL, BE, RHO, L, mode="ip")
    df = ds.DenseSys(Xf, sf, af, np.zeros(nf), dt, AL, BE, RHO, L, mode="ip")
    tb = time.time() - t0
    vf = lambda t, f: np.full(len(f), 1.0 if t <= trise + 1e-12 else 0.0)
    mv = ds.Marcher(dv, nl=nl, vf_fun=vf); mf = ds.Marcher(df, nl=nl, vf_fun=vf)
    # first arrival of the surface reflection at the fault: distance fault -> surface -> fault, P wave
    dmin = np.min(-Xf[:, 1]); t_refl = 2*dmin/AL
    early = 0.0; scale = 0.0; vsmax = []
    for n in range(nsteps):
        mv.step(); mf.step()
        if mv.t < t_refl - 2*dt: early = max(early, np.abs(mv.tau[:nf] - mf.tau[:nf]).max())
        scale = max(scale, np.abs(mf.tau[:nf]).max())
        vsmax.append(np.abs(mv.V[nf:]).max())
    # static exact: half-plane traction of the final slip on the fault
    MU = RHO*BE**2; fs.MU = MU; fs.NU = (AL**2 - 2*BE**2)/(2*(AL**2 - BE**2))
    Th = np.zeros(nf); Tf = np.zeros(nf)
    for i in range(nf):
        for j in range(nf):
            st, so = fs.tdd(Xf[i, 0], Xf[i, 1], Xf[j, 0], Xf[j, 1], hf/2, sfv[0], sfv[1], 1)
            Th[i] += fs.traction(st, sfv)[0]*mv.delta[j]
    g = df.Ks[:nf, :nf] @ mf.delta[:nf]
    sgn = np.sign(np.sum(g*Th)) if np.sum(np.abs(Th)) > 0 else 1.0
    stat = np.abs(mv.tau[:nf] - sgn*Th).max()/np.abs(Th).max()
    statf = np.abs(mf.tau[:nf] - sgn*Th).max()/np.abs(Th).max()
    vs = np.array(vsmax)
    if verbose:
        print(f"in-plane {'NL' if nl else 'first order'}: dip {dip:g}, fault {Lf/1e3:g} km / {hf:g} m, top {top:g} m, surface +-{Ls/1e3:g} km / {hs:g} m, "
              f"{nsteps} steps, tables {tb:.0f} s")
        print(f"   before the reflection (t < {t_refl:.2f} s): max |tau - tau_fullspace| / max|tau| = {early/scale:.1e}")
        print(f"   at t = {T:g} s: |tau - exact half-plane static| / max = {stat:.2e} (full space: {statf:.2e})")
        print(f"   surface slip-rate max over the last 20 % of the run {vs[int(0.8*len(vs)):].max():.2e} (peak {vs.max():.2e})")
    return early/scale, stat

def inplane_rupture(dip=30.0, Lf=6e3, hf=200.0, top=400.0, Ls=8e3, hs=200.0, T=5.0, nl=True, surface_on=True, verbose=True, tag="", sense=-1.0, stagger=True, nucl=1.2e3, normal_coupling=True):
    """Spontaneous rupture (slip weakening, normal-stress coupling) on a dipping
    thrust whose top is 'top' below the surface; nucleation by overstress at the
    bottom 1.2 km. Returns final slip along the fault and diagnostics."""
    # the step follows the FAULT elements (cS dt/dx = 0.5 there); finer surface elements then
    # see cS dt/h_s = 0.5 hf/hs (1 for hs = hf/2, 2 for hf/4), inside the NL stability range
    # (0.4 to 2, stab_random.py); a step set by the surface elements would put the fault at
    # 0.5 hs/hf, below 0.4, where NL is unstable
    dt = 0.5*hf/BE; nsteps = int(round(T/dt)); L = nsteps + 1
    ang = -np.radians(dip); sfv = np.array([np.cos(ang), np.sin(ang)])
    nf = int(round(Lf/hf)); Xf = np.array([[0.0, -top]]) + ((np.arange(nf) + 0.5)*hf)[:, None]*sfv[None, :]
    sf = np.tile(sfv, (nf, 1)); af = np.full(nf, hf/2)
    if surface_on:
        Xs, ss, as_ = surface(Ls, hs, x0=Lf*np.cos(ang)/2)
        if stagger and top == 0.0:
            # a fault reaching the surface: put the fault trace (x = 0) on a boundary between two
            # surface elements (a collocation point on the fault tip would see its singular field)
            k = np.argmin(np.abs(Xs[:, 0])); Xs[:, 0] -= Xs[k, 0] - hs/2
        X = np.vstack([Xf, Xs]); S = np.vstack([sf, ss]); A = np.concatenate([af, as_]); K = np.r_[np.zeros(nf), np.ones(len(Xs))]
    else:
        # a fault alone has no surface unknowns: add one tiny far-away surface element (zero traction, negligible)
        X = np.vstack([Xf, [[1e9, 0.0]]]); S = np.vstack([sf, [[1.0, 0.0]]]); A = np.r_[af, 1.0]; K = np.r_[np.zeros(nf), 1.0]
    t0 = time.time()
    d = ds.DenseSys(X, S, A, K, dt, AL, BE, RHO, L, mode="ip")
    tb = time.time() - t0
    tau0 = np.zeros(d.Nc); s_along = (np.arange(nf) + 0.5)*hf
    # sense -1: thrust (the tangent points down dip, so negative slip moves the hanging wall up)
    tau0[:nf] = sense*np.where(s_along > Lf - nucl, 81.6e6, 70e6)            # nucleation at the bottom
    sigma0 = np.full(nf, 120e6)
    m = ds.Marcher(d, nl=nl, tau0=tau0, sigma0=sigma0, sw=(0.677, 0.525, 0.40), normal_coupling=normal_coupling)
    vmax_s = []; vf_hist = []
    for n in range(nsteps):
        m.step()
        vmax_s.append(np.abs(m.V[d.surf]).max() if len(d.surf) else 0.0)
        vf_hist.append(m.V[:nf].copy())
    vf_hist = np.array(vf_hist); tt = (np.arange(nsteps) + 1)*dt
    av = np.abs(vf_hist)
    rt = np.array([tt[np.argmax(av[:, i] > 0.1)] if (av[:, i] > 0.1).any() else np.nan for i in range(nf)])
    if verbose:
        vs = np.array(vmax_s)
        print(f"rupture {tag}{'thrust ' if sense < 0 else 'normal '}{'NL' if nl else 'first order'} {'with' if surface_on else 'without'} surface: dip {dip:g}, fault {Lf/1e3:g} km / {hf:g} m, top {top:g} m, "
              f"surface +-{Ls/1e3:g} km / {hs:g} m, {nsteps} steps, tables {tb:.0f} s")
        print(f"   |slip| top / middle / bottom: {abs(m.delta[0]):.3f} / {abs(m.delta[nf//2]):.3f} / {abs(m.delta[nf-1]):.3f} m; arrival at the top {rt[0]:.3f} s; "
              f"peak fault slip rate {np.abs(vf_hist).max():.2f} m/s; surface rate max (last 20 %) {vs[int(0.8*len(vs)):].max():.2e}")
    return np.abs(m.delta[:nf]).copy(), rt


if __name__ == "__main__":
    # quick checks (about 2 minutes); the longer runs recorded in fd_plan.md use the same functions
    e = antiplane_case(nl=True, T=3.0)
    assert e < 5e-3, e
    early, stat = inplane_kinematic(T=3.0, Ls=4e3, hs=250.0)
    assert early < 1e-12 and stat < 1e-2, (early, stat)
    full = inplane_rupture(Lf=3e3, Ls=4e3, T=3.0, surface_on=False)
    thr = inplane_rupture(Lf=3e3, Ls=4e3, T=3.0)
    nrm = inplane_rupture(Lf=3e3, Ls=4e3, T=3.0, sense=1.0)
    assert thr[0][0] > nrm[0][0] > 0.0 and thr[0][0] > full[0][0]      # free surface favours thrusts near the top
    print("test_fs_dyn2d: PASS")
