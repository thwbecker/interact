import numpy as np, pytest, time
import inplane as ip
import shells as sh

AL, BE, RHO = 6000.0, 3464.0, 2700.0
MU = RHO * BE**2
NU = (AL**2 - 2 * BE**2) / (2 * (AL**2 - BE**2))


def small_system(n1=24, n2=16, W=600.0, gap=-200.0, dx=50.0, tol=1e-2):
    f1 = sh.Fault([0.0, 0.0], [1.0, 0.0], n1 * dx, n1, "f1")
    f2 = sh.Fault([n1 * dx + gap, W], [1.0, 0.0], n2 * dx, n2, "f2")
    dt = 0.5 * dx / BE
    return sh.System([f1, f2], dt, AL, BE, RHO, tol=tol), dt


def test_tables_match_direct_kernels():
    sys_, dt = small_system()
    f1, f2 = sys_.faults
    for (a_, b_), p in sys_.pairs.items():
        fs, fr = sys_.faults[a_], sys_.faults[b_]
        for (i, j) in [(0, 0), (3, 7), (fr.n - 1, 0), (2, fs.n - 1)]:
            lags = np.arange(1, min(p.L, 40) + 1)
            ks, kn = ip.lag_kernel(fs.centers[j], fs.s, fs.a, fr.centers[i], fr.s, lags, dt, AL, BE, RHO)
            Ss, Sn = ip.static_traction(fs.centers[j], fs.s, fs.a, fr.centers[i], fr.s, MU, NU)
            for l in lags:
                Ws, Wn = p.W_dense(l)
                assert abs(Ws[i, j] - (ks[l - 1] - Ss * dt)) < 1e-9 * (abs(Ss * dt) + abs(ks[l - 1]))
                assert abs(Wn[i, j] - (kn[l - 1] - Sn * dt)) < 1e-9 * (abs(Ss * dt) + abs(kn[l - 1]) + 1e-30)
            ksb, knb = p.static_blocks()
            assert abs(ksb[i, j] - Ss) < 1e-12 * abs(Ss) + 1e-20
            assert abs(knb[i, j] - Sn) < 1e-12 * (abs(Ss) + abs(Sn))


def test_history_convolution_matches_dense():
    sys_, dt = small_system()
    rng = np.random.default_rng(0)
    nsteps = sys_.L + 30
    Vs = rng.standard_normal((nsteps, sys_.N)) * 1e-3
    # dense reference tables
    Wd = {}
    for key, p in sys_.pairs.items():
        Wd[key] = [p.W_dense(l) for l in range(1, p.L + 1)]
    for n in range(nsteps):
        sys_.push(Vs[n])
        if n in (0, 5, sys_.L // 2, nsteps - 1):
            ts, tn = sys_.history_stress(lag0=1)
            rs = np.zeros(sys_.N); rn = np.zeros(sys_.N)
            for (a_, b_), p in sys_.pairs.items():
                for l in range(1, p.L + 1):
                    m = n + 1 - l
                    if m < 0:
                        break
                    Ws, Wn = Wd[(a_, b_)][l - 1]
                    v = Vs[m][sys_.off[a_]:sys_.off[a_ + 1]]
                    rs[sys_.off[b_]:sys_.off[b_ + 1]] += Ws @ v
                    rn[sys_.off[b_]:sys_.off[b_ + 1]] += Wn @ v
            assert np.allclose(ts, rs, rtol=0, atol=1e-10 * np.abs(rs).max())
            assert np.allclose(tn, rn, rtol=0, atol=1e-10 * (np.abs(rn).max() + np.abs(rs).max()))
            # lag0 = 2 excludes the newest slice
            ts2, tn2 = sys_.history_stress(lag0=2)
            w1s = sys_.K1s - sys_.Ks * dt; w1n = sys_.K1n - sys_.Kn * dt
            assert np.allclose(ts2 + w1s @ Vs[n], ts, rtol=0, atol=1e-10 * np.abs(rs).max())
            assert np.allclose(tn2 + w1n @ Vs[n], tn, rtol=0, atol=1e-10 * (np.abs(rn).max() + np.abs(rs).max()))


def test_window_and_first_lag_properties():
    sys_, dt = small_system()
    # first-lag self term close to -eta (P from the ends arrives within dt for beta = 1/2)
    d = np.diag(sys_.K1s)
    assert np.all(d < 0) and np.all(np.abs(d + sys_.eta) < 0.5 * sys_.eta)
    # first-lag coupling only within the P distance of one step
    f1 = sys_.faults[0]
    row = sys_.K1s[5, :f1.n]
    far = np.abs(np.arange(f1.n) - 5) * f1.dx - f1.dx > AL * dt
    assert np.all(np.abs(row[far]) < 1e-12 * np.abs(row).max())   # FFT round-trip noise only
    # the static Kn block vanishes on the same line, and is nonzero across the step
    assert np.abs(sys_.Kn[:f1.n, :f1.n]).max() < 1e-9 * np.abs(sys_.Ks).max()
    assert np.abs(sys_.Kn[f1.n:, :f1.n]).max() > 1e-3 * np.abs(sys_.Ks).max()
    # window grows when the tolerance tightens
    s2, _ = small_system(tol=1e-3)
    assert s2.L > sys_.L


def kostrov_slip_avg(x, dx, t, v, C, nsub=41):
    """element-averaged self-similar slip C sqrt(v^2 t^2 - x^2)"""
    xs = x[:, None] + np.linspace(-0.5, 0.5, nsub)[None, :] * dx
    r = v * v * t * t - xs * xs
    return np.where(r > 0, C * np.sqrt(np.maximum(r, 0)), 0.0).mean(axis=1)


@pytest.mark.parametrize("vfrac", [0.5, 0.8])
def test_kostrov_self_similar_crack_kinematic(vfrac):
    """Impose the self-similar slip delta = C sqrt(v^2 t^2 - x^2) (Kostrov /
    Burridge-Willis form) on a single fault and check the computed stress:
    the mean shear stress change over the central part of the crack is
    constant in time (self-similarity), the normal stress change on the plane
    vanishes, and the implied dynamic factor C(v)/C(0) lies below 1 and
    decreases with v (energy radiated, as for the stress intensity factor
    universal function). The piecewise-constant slip-rate representation of
    the square-root tip makes the stress profile ragged (Tada and Madariaga
    2001 remark on this), so uniformity in x is not asserted here; it is the
    stepper test with prescribed stress that checks the profile."""
    dx = 100.0; n = 120
    f = sh.Fault([-n * dx / 2, 0.0], [1.0, 0.0], n * dx, n, "f")
    dt = 0.5 * dx / BE
    sys_ = sh.System([f], dt, AL, BE, RHO, tol=1e-2)
    v = vfrac * BE
    x = f.centers[:, 0]
    C = 1e-3
    nsteps = int(0.8 * (n * dx / 2) / v / dt)
    delta_prev = np.zeros(n)
    rows = []
    for k in range(nsteps):
        t1 = (k + 1) * dt
        delta = kostrov_slip_avg(x, dx, t1, v, C)
        sys_.push((delta - delta_prev) / dt)
        delta_prev = delta
        if k % 10 == 9 and v * t1 > 12 * dx:
            ts, tn = sys_.history_stress(lag0=1)
            tau = sys_.Ks @ delta + ts
            mid = np.abs(x) < 0.25 * v * t1
            rows.append((t1, tau[mid].mean(), np.abs(tn[mid]).max()))
    rows = np.array(rows)
    means = rows[:, 1]
    assert np.std(means) < 0.02 * abs(means.mean())        # self-similar: constant in time (ragged tip)
    assert np.all(rows[:, 2] < 1e-9 * abs(means.mean()))   # no normal stress change on the plane
    assert means.mean() < 0
    Cstat = 2 * (1 - NU) / MU
    factor = (C / Cstat) / abs(means.mean())
    assert 0.5 < factor < 1.0
    np.save(f"/tmp/kostrov_factor_{vfrac}.npy", np.array([abs(means.mean()), factor]))


def test_kostrov_factor_decreases_with_speed():
    import os
    f5 = "/tmp/kostrov_factor_0.5.npy"; f8 = "/tmp/kostrov_factor_0.8.npy"
    if not (os.path.exists(f5) and os.path.exists(f8)):
        pytest.skip("run the kinematic tests first")
    assert np.load(f8)[1] < np.load(f5)[1]
