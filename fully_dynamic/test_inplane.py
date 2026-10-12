import numpy as np, pytest
from scipy.integrate import quad
from inplane_kernel import sigma_point_ramp
import inplane as ip

AL, BE, RHO = 6000.0, 3464.0, 2700.0
MU = RHO * BE**2
NU = (AL**2 - 2 * BE**2) / (2 * (AL**2 - BE**2))
ETA = MU / (2 * BE)


def kelvin_dc(x, y):
    """static double couple (point slip patch) = d/dx of Hirth-Lothe field"""
    D = MU / (2 * np.pi * (1 - NU)); r6 = (x * x + y * y) ** 3
    return (2 * x * y * (3 * x * x - y * y) / r6 * D,
            (-x**4 + 6 * x * x * y * y - y**4) / r6 * D,
            2 * x * y * (-x * x + 3 * y * y) / r6 * D)


@pytest.mark.parametrize("xy", [(1000.0, 300.0), (-200.0, 900.0), (500.0, -500.0), (2000.0, 10.0), (0.0, 50.0)])
def test_point_static_limit(xy):
    x, y = xy
    # convergence to the static limit goes like (R/(beta t))^-2 while the
    # roundoff of the generated expressions grows with beta t / R; choose
    # beta t / R = 1e3 (truncation 1e-6, roundoff well below that)
    t = 1e3 * np.hypot(x, y) / BE
    s = sigma_point_ramp(x, y, t, AL, BE, RHO)
    k = kelvin_dc(x, y)
    scale = max(abs(v) for v in k)
    for si, ki in zip(s, k):
        assert abs(si / t - ki) < 1e-5 * scale


def test_point_causality_and_regions():
    x, y = 300.0, 400.0
    R = 500.0
    assert all(v == 0 for v in sigma_point_ramp(x, y, R / AL * 0.999, AL, BE, RHO))
    sp = sigma_point_ramp(x, y, 0.5 * (R / AL + R / BE), AL, BE, RHO)
    assert any(v != 0 for v in sp)


@pytest.mark.parametrize("y,t", [(100.0, 0.05), (1000.0, 0.2), (1000.0, 1.0)])
def test_plane_wave_identity(y, t):
    """Integral of the point field over the whole wave cone equals the plane
    S-wave field of an infinite fault with unit slip rate: (0, mu/(2 beta), 0)
    behind the S front, zero ahead of it (sign: generated field is minus the
    physical one)."""
    L = AL * t * 1.0001
    bps = sorted(set(sum([[-np.sqrt(c * c * t * t - y * y), np.sqrt(c * c * t * t - y * y)]
                          for c in (AL, BE) if c * t > abs(y)], [])))
    vals = []
    for comp in range(3):
        f = lambda xi: sigma_point_ramp(xi, y, t, AL, BE, RHO)[comp].item()
        vals.append(quad(f, -L, L, points=bps, limit=500)[0])
    ref = (0.0, ETA if BE * t > abs(y) else 0.0, 0.0)
    for v, r in zip(vals, ref):
        assert abs(v - r) < 1e-4 * ETA   # adaptive quad accuracy near the 1/sqrt fronts


@pytest.mark.parametrize("method", ["closed", "quad"])
@pytest.mark.parametrize("xr,yr", [(0.0, 400.0), (900.0, 400.0), (-2500.0, 150.0), (300.0, 0.0), (0.0, 0.0), (120.0, 0.0), (-45.0, 0.0)])
def test_segment_static_limit(xr, yr, method):
    a = 50.0
    t = 1e3 * max(np.hypot(abs(xr) + a, yr), a) / BE
    s = ip.segment_ramp_stress(xr, yr, t, a, AL, BE, RHO, method=method)
    k = ip.static_segment_stress(xr, yr, a, MU, NU)
    scale = max(abs(float(v)) for v in k)
    for si, ki in zip(s, k):
        assert abs(float(si) / t - float(ki)) < 1e-5 * scale, (si / t, ki)


@pytest.mark.parametrize("method", ["closed", "quad"])
def test_segment_causality(method):
    a = 50.0
    xr, yr = 800.0, 300.0
    R = np.hypot(xr - a, yr)  # nearest source point
    s = ip.segment_ramp_stress(xr, yr, R / AL * 0.999, a, AL, BE, RHO, method=method)
    assert all(float(v) == 0 for v in s)
    s = ip.segment_ramp_stress(xr, yr, R / AL * 1.01, a, AL, BE, RHO, method=method)
    assert any(float(v) != 0 for v in s)


def test_p_radiates_from_segment_ends_only():
    """Uniform slip radiates P waves only from the slip-gradient
    discontinuities at the segment ends; above the interior nothing arrives
    before the P wave from the nearer end. Exact in the closed form; the
    quadrature engine cancels the interior contributions only numerically."""
    a = 50.0
    xr, yr = 0.0, 30.0
    # after the P wave from the point directly below (yr/AL = 0.0050 s), before
    # both the plane S wave (yr/BE = 0.0087 s) and the P wave from the ends
    # (hypot(a, yr)/AL = 0.0097 s)
    t = 0.007
    assert yr / AL < t < min(yr / BE, np.hypot(a, yr) / AL)
    s = ip.segment_ramp_stress(xr, yr, t, a, AL, BE, RHO, method="closed")
    assert all(float(v) == 0 for v in s)
    q = ip.segment_ramp_stress(xr, yr, t, a, AL, BE, RHO, method="quad")
    assert all(abs(float(v)) < 1e-6 * ETA for v in q)


@pytest.mark.parametrize("method", ["closed", "quad"])
def test_self_element_early_time_is_radiation_damping(method):
    a = 50.0
    for x in (0.0, 20.0, -40.0):
        t = (a - abs(x)) / AL * 0.9
        s = ip.segment_ramp_stress(x, 0.0, t, a, AL, BE, RHO, method=method)
        assert abs(float(s[1]) + ETA) < 1e-9 * ETA
        assert abs(float(s[0])) < 1e-9 * ETA and abs(float(s[2])) < 1e-9 * ETA
    # after the P wave from the ends arrives the value departs from -eta
    s = ip.segment_ramp_stress(0.0, 0.0, 3 * a / AL, a, AL, BE, RHO, method=method)
    assert abs(float(s[1]) + ETA) > 1e-3 * ETA


def test_closed_form_against_independent_derivation():
    """External validation of the in-plane kernel: the Tada-Madariaga /
    Ando closed form (production engine) against the quadrature of the point
    field derived from the 3-D Green's function in gen_inplane.py."""
    a = 50.0
    rng = np.random.default_rng(1)
    pts = [(rng.uniform(-3000, 3000), rng.uniform(-3000, 3000)) for _ in range(120)]
    pts = [(x, (y if abs(y) > 1 else 1.0)) for x, y in pts]
    pts += [(0.0, 0.0), (20.0, 0.0), (-45.0, 0.0), (120.0, 0.0), (-300.0, 0.0), (1000.0, 0.0), (0.0, 5.0), (0.0, 30.0)]
    errs = []; errs_front = []
    for x, y in pts:
        Rn = np.hypot(max(abs(x) - a, 0), y); Rn = Rn if Rn > 0 else a
        for frac in [0.5, 0.9, 0.999, 1.001, 1.01, 1.05, 1.3, 1.7, 2.0, 5.0, 20.0, 200.0]:
            for c in (AL, BE):
                t = frac * Rn / c
                o = ip.segment_ramp_stress(x, y, t, a, AL, BE, RHO, method="quad")
                r = ip.segment_ramp_stress(x, y, t, a, AL, BE, RHO, method="closed")
                scale = max(abs(float(v)) for v in r)
                if scale < 1e-9 * ETA:
                    assert all(abs(float(v)) < 1e-6 * ETA for v in o)
                    continue
                e = max(abs(float(u) - float(v)) for u, v in zip(o, r)) / scale
                (errs_front if abs(frac - 1) < 0.02 else errs).append(e)
    errs = np.array(errs); errs_front = np.array(errs_front)
    assert np.median(errs) < 1e-10
    assert errs.max() < 1e-5
    assert errs_front.max() < 5e-2     # tiny absolute values just after an arrival


@pytest.mark.parametrize("xr,yr,t", [(300.0, 200.0, 0.08), (300.0, 200.0, 0.0605), (1500.0, 40.0, 0.3), (200.0, 0.0, 0.1)])
def test_quadrature_against_adaptive(xr, yr, t):
    a = 50.0
    s = ip.segment_ramp_stress(xr, yr, t, a, AL, BE, RHO, ngauss=12)
    bps = []
    for c in (AL, BE):
        if c * t > abs(yr):
            h = np.sqrt(c * c * t * t - yr * yr)
            bps += [xr - h, xr + h]
    bps = [b for b in bps if -a < b < a]
    ref = []
    for comp in range(3):
        f = lambda xp: sigma_point_ramp(xr - xp, yr, t, AL, BE, RHO)[comp].item()
        ref.append(-quad(f, -a, a, points=bps or None, limit=400, epsabs=0, epsrel=1e-10)[0])
    scale = max(abs(v) for v in ref)
    for si, ri in zip(s, ref):
        assert abs(float(si) - ri) < 1e-6 * scale, (si, ri)   # limited by roundoff near fronts


def test_lag_kernel_telescopes_and_projection_invariance():
    a = 25.0; dt = 0.5 * 2 * a / BE
    Xs = np.array([0.0, 0.0]); ss = np.array([1.0, 0.0])
    Xr = np.array([600.0, 250.0]); sr = np.array([np.cos(0.3), np.sin(0.3)])
    lags = np.arange(1, 400)
    ks, kn = ip.lag_kernel(Xs, ss, a, Xr, sr, lags, dt, AL, BE, RHO)
    xr, yr = ip.receiver_in_source_frame(Xs, ss, Xr)
    s = ip.segment_ramp_stress(xr, yr, lags[-1] * dt, a, AL, BE, RHO)
    # both engines give the same lag kernels
    ks2q, kn2q = ip.lag_kernel(Xs, ss, a, Xr, sr, lags, dt, AL, BE, RHO, method="quad")
    assert np.allclose(ks, ks2q, rtol=0, atol=1e-5 * np.abs(ks).max())   # quadrature front error
    assert np.allclose(kn, kn2q, rtol=0, atol=1e-5 * np.abs(kn).max())
    ts, tn = ip.project(s[0], s[1], s[2], ss, sr)
    assert abs(ks.sum() - ts) < 1e-9 * abs(ts) + 1e-12
    assert abs(kn.sum() - tn) < 1e-9 * abs(tn) + 1e-12
    # rotate the whole configuration: tractions must not change
    th = 0.7; Q = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    ks2, kn2 = ip.lag_kernel(Q @ Xs, Q @ ss, a, Q @ Xr, Q @ sr, lags, dt, AL, BE, RHO)
    assert np.allclose(ks, ks2, rtol=1e-10, atol=1e-9 * np.abs(ks).max())
    assert np.allclose(kn, kn2, rtol=1e-10, atol=1e-9 * np.abs(kn).max())
    # static traction consistency of the projection
    st_s, st_n = ip.static_traction(Xs, ss, a, Xr, sr, MU, NU)
    st_s2, st_n2 = ip.static_traction(Q @ Xs, Q @ ss, a, Q @ Xr, Q @ sr, MU, NU)
    assert abs(st_s - st_s2) < 1e-9 * abs(st_s) and abs(st_n - st_n2) < 1e-9 * abs(st_n)


def test_lag_kernel_static_sum():
    """sum_l K^(l) over a long window approaches K_s * dt (per unit slip rate step)."""
    a = 25.0; dt = 0.5 * 2 * a / BE
    Xs = np.array([0.0, 0.0]); ss = np.array([1.0, 0.0])
    Xr = np.array([300.0, 150.0]); sr = np.array([1.0, 0.0])
    T = 400.0
    L = int(T / dt)
    xr, yr = ip.receiver_in_source_frame(Xs, ss, Xr)
    s = ip.segment_ramp_stress(xr, yr, L * dt, a, AL, BE, RHO)
    ts, tn = ip.project(s[0], s[1], s[2], ss, sr)
    st_s, st_n = ip.static_traction(Xs, ss, a, Xr, sr, MU, NU)
    assert abs(ts / (L * dt) - st_s) < 1e-5 * abs(st_s)
    assert abs(tn / (L * dt) - st_n) < 1e-5 * abs(st_n)


def test_reciprocity_parallel_equal_elements():
    a = 25.0; dt = 0.5 * 2 * a / BE
    s = np.array([1.0, 0.0])
    X1 = np.array([0.0, 0.0]); X2 = np.array([700.0, 300.0])
    lags = np.arange(1, 300)
    k12 = ip.lag_kernel(X1, s, a, X2, s, lags, dt, AL, BE, RHO)
    k21 = ip.lag_kernel(X2, s, a, X1, s, lags, dt, AL, BE, RHO)
    for u, v in zip(k12, k21):
        assert np.allclose(u, v, rtol=1e-9, atol=1e-9 * np.abs(u).max())


def test_time_step_consistency():
    """Response to a smooth slip-rate history must converge as dt is refined."""
    a = 25.0
    Xs = np.array([0.0, 0.0]); ss = np.array([1.0, 0.0])
    Xr = np.array([400.0, 200.0]); sr = np.array([1.0, 0.0])
    tc, w = 0.25, 0.04
    V = lambda t: np.exp(-0.5 * ((t - tc) / w) ** 2)
    tobs = 0.6
    res = []
    for n in (1, 2, 4):
        dt = 0.5 * 2 * a / BE / n
        L = int(tobs / dt)
        lags = np.arange(1, L + 1)
        ks, kn = ip.lag_kernel(Xs, ss, a, Xr, sr, lags, dt, AL, BE, RHO)
        # V^m is the slip rate during [m dt, (m+1) dt); midpoint sampling
        Vm = V((L - lags + 0.5) * dt)
        res.append((np.sum(ks * Vm), np.sum(kn * Vm)))
    e1 = abs(res[0][0] - res[2][0]); e2 = abs(res[1][0] - res[2][0])
    assert e2 < 0.6 * e1
    assert abs(res[2][0]) > 0


def test_lag_kernel_table_matches_lag_kernel():
    a = 25.0; dt = 0.5 * 2 * a / BE
    Xs = np.array([0.0, 0.0]); ss = np.array([1.0, 0.0])
    Xr = np.array([[300.0, 150.0], [-120.0, 40.0], [700.0, 0.0], [0.0, 0.0]])
    sr = np.array([1.0, 0.0])
    ks, kn = ip.lag_kernel_table(Xs, ss, a, Xr, sr, 60, dt, AL, BE, RHO)
    lags = np.arange(1, 61)
    for m in range(Xr.shape[0]):
        k1, k2 = ip.lag_kernel(Xs, ss, a, Xr[m], sr, lags, dt, AL, BE, RHO)
        assert np.allclose(ks[m], k1, rtol=1e-12, atol=1e-10 * np.abs(k1).max())
        assert np.allclose(kn[m], k2, rtol=1e-12, atol=1e-10 * (np.abs(k2).max() + 1e-30))
    # self element: first lag close to -eta (P from the ends arrives within dt here)
    assert ks[3, 0] < 0 and abs(ks[3, 0] + ETA) < 0.5 * ETA
