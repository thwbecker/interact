import numpy as np, pytest, warnings
from scipy.integrate import quad
import antiplane as ap
import ref_ando2007 as ref
import inplane as ip

AL, BE, RHO = 6000.0, 3464.0, 2700.0
MU = RHO * BE**2
ETA = MU / (2 * BE)


@pytest.mark.parametrize("y,t", [(100.0, 0.05), (1000.0, 1.0), (1000.0, 0.2)])
def test_plane_sh_wave_identity(y, t):
    L = BE * t * 1.0001
    pts = [-np.sqrt((BE * t) ** 2 - y * y), np.sqrt((BE * t) ** 2 - y * y)] if BE * t > y else None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = [quad(lambda xi: ap.sigma_point_ramp_ap(xi, y, t, AL, BE, RHO)[c].item(), -L, L, points=pts, limit=500)[0] for c in (0, 1)]
    ref_v = (0.0, ETA if BE * t > y else 0.0)
    assert abs(v[0] - ref_v[0]) < 1e-5 * ETA and abs(v[1] - ref_v[1]) < 1e-5 * ETA


@pytest.mark.parametrize("method", ["closed", "quad"])
@pytest.mark.parametrize("xr,yr", [(300.0, 200.0), (0.0, 0.0), (120.0, 0.0), (-45.0, 0.0), (0.0, 30.0), (-900.0, -400.0)])
def test_antiplane_segment_static_limit(xr, yr, method):
    a = 50.0
    t = 1e3 * max(np.hypot(abs(xr) + a, yr), a) / BE
    o = ap.segment_ramp_stress_ap(xr, yr, t, a, AL, BE, RHO, method=method)
    st = ap.static_segment_stress_ap(xr, yr, a, MU)
    scale = max(abs(float(v)) for v in st)
    for oi, si in zip(o, st):
        assert abs(float(oi) / t - float(si)) < 1e-5 * scale


def test_reference_transcription_static_limit():
    a = 50.0
    for xr, yr in [(300.0, 200.0), (0.0, 0.0), (120.0, 0.0)]:
        t = 1e3 * max(np.hypot(abs(xr) + a, yr), a) / BE
        r = ref.antiplane_segment_ramp_stress_ref(xr, yr, t, a, BE, RHO)
        st = ap.static_segment_stress_ap(xr, yr, a, MU)
        scale = max(abs(float(v)) for v in st)
        for ri, si in zip(r, st):
            assert abs(float(ri) / t - float(si)) < 1e-5 * scale


def test_against_ando2007_antiplane():
    """Entry-by-entry comparison of the derivation pipeline (double time
    integral, double couple, segment quadrature, complement identity; the
    "quad" engine) with the closed-form Ando et al. (2007) eqs (30)-(31)
    (the production engine)."""
    a = 50.0
    rng = np.random.default_rng(2)
    pts = [(rng.uniform(-3000, 3000), rng.uniform(-3000, 3000)) for _ in range(120)]
    pts = [(x, (y if abs(y) > 1 else 1.0)) for x, y in pts]
    pts += [(0.0, 0.0), (20.0, 0.0), (-45.0, 0.0), (120.0, 0.0), (-300.0, 0.0), (1000.0, 0.0), (0.0, 5.0), (0.0, 30.0)]
    errs = []; errs_front = []
    for x, y in pts:
        Rn = np.hypot(max(abs(x) - a, 0), y); Rn = Rn if Rn > 0 else a
        for frac in [0.5, 0.9, 0.999, 1.001, 1.01, 1.05, 1.3, 1.7, 2.0, 5.0, 20.0, 200.0]:
            t = frac * Rn / BE
            o = ap.segment_ramp_stress_ap(x, y, t, a, AL, BE, RHO, method="quad")
            r = ap.segment_ramp_stress_ap(x, y, t, a, AL, BE, RHO, method="closed")
            scale = max(abs(float(v)) for v in r)
            if scale < 1e-9 * ETA:
                assert all(abs(float(v)) < 1e-9 * ETA for v in o)
                continue
            e = max(abs(float(u) - float(v)) for u, v in zip(o, r)) / scale
            (errs_front if abs(frac - 1) < 0.02 else errs).append(e)
    errs = np.array(errs); errs_front = np.array(errs_front)
    assert np.median(errs) < 1e-10
    assert errs.max() < 1e-6
    assert errs_front.max() < 1e-3     # limited by front cancellation in double precision


def test_antiplane_lag_kernel_self_first_lag():
    """With beta*dt < a (slow S wave relative to the element) the first-lag
    self kernel is exactly -eta (no end effects within one step)."""
    a = 50.0
    dt = 0.5 * a / BE
    s1 = ap.segment_ramp_stress_ap(0.0, 0.0, dt, a, AL, BE, RHO)
    assert abs(float(s1[1]) + ETA) < 1e-9 * ETA and abs(float(s1[0])) < 1e-9 * ETA
