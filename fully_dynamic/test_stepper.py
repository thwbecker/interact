import numpy as np, pytest, os
import shells as sh, stepper as st

AL, BE, RHO = 6000.0, 3464.0, 2700.0
MU = RHO * BE**2
NU = (AL**2 - 2 * BE**2) / (2 * (AL**2 - BE**2))

from scipy.special import ellipe, ellipk


def kostrov_factor_III(vfrac):
    """Exact dynamic factor of the symmetrically expanding mode III crack
    (Kostrov 1964; Broberg 1999 eq. 6.9.145): slip = (2 tau/mu) sqrt(a^2-x^2) / E(k),
    k^2 = 1 - (v/beta)^2, E the complete elliptic integral of the second kind."""
    return 1.0 / ellipe(1 - vfrac**2)


def kostrov_factor_II(vfrac, alpha=AL, beta=BE):
    """Exact dynamic factor of the symmetrically expanding mode II crack
    (Broberg 1999 eqs 6.9.85-6.9.87): C(v)/C(0) = 2 (1-k^2)(k^2-b^2)/g2(b),
    b = v/alpha, k = beta/alpha."""
    k = beta / alpha; b = vfrac * k
    E1 = ellipe(1 - b * b); K1 = ellipk(1 - b * b)
    E2 = ellipe(1 - b * b / (k * k)); K2 = ellipk(1 - b * b / (k * k))
    g2 = ((8 * k * k * (k * k - b * b) / (b * b)) * E1 - 4 * k * k * (k * k - b * b) * K1
          - (8 * k**4 - 8 * k * k * b * b + b**4) / (b * b) * E2 + (4 * k * k - 3 * b * b) * K2)
    return 2 * (1 - k * k) * (k * k - b * b) / g2


# the independent spectral solver (spectral_ap.py) gave 0.8258 and 0.7042 for
# vfrac = 0.5 and 0.8, against the exact 0.8257 and 0.7052


def run_kostrov(mode, vfrac, dx=100.0, n=120, tol=1e-2):
    f = sh.Fault([-n * dx / 2, 0.0], [1.0, 0.0], n * dx, n)
    dt = 0.5 * dx / BE
    sys_ = sh.System([f], dt, AL, BE, RHO, tol=tol, mode=mode)
    v = vfrac * BE; x = f.centers[:, 0]
    S = st.Stepper(sys_, st.PrescribedTip(x, v, 1e6), tau0=0.0, sigma0=1e8)
    for _ in range(int(0.8 * (n * dx / 2) / v / dt)):
        S.step()
    return S, x, v


@pytest.mark.parametrize("vfrac", [0.5, 0.8])
def test_self_similar_crack_antiplane_against_spectral(vfrac):
    S, x, v = run_kostrov("antiplane", vfrac)
    t = S.t
    mid = np.abs(x) < 0.25 * v * t
    C = np.mean(S.delta[mid] / np.sqrt(v * v * t * t - x[mid] ** 2)) / (2 / MU * 1e6)
    assert abs(C - kostrov_factor_III(vfrac)) < 3e-3
    # elliptical profile over the inner 90 percent of the crack
    ins = np.abs(x) < 0.9 * v * t
    ref = C * (2 / MU * 1e6) * np.sqrt(v * v * t * t - x[ins] ** 2)
    assert np.abs(S.delta[ins] - ref).max() < 1e-2 * ref.max()
    # stress on the free part equals the prescribed drop; locked part does not slip
    assert np.abs(S.tau[np.abs(x) < v * t] + 1e6).max() < 1e-6
    assert np.all(S.V[np.abs(x) >= v * t] == 0)


@pytest.mark.parametrize("vfrac", [0.5, 0.8])
def test_self_similar_crack_inplane_against_broberg(vfrac):
    """In-plane: stepper against the exact Broberg (1999) mode II factor;
    the slip profile must be elliptical."""
    fac = kostrov_factor_II(vfrac)
    S, x, v = run_kostrov("inplane", vfrac)
    t = S.t
    Cstat = 2 * (1 - NU) / MU
    ref = fac * Cstat * 1e6 * np.sqrt(np.maximum(v * v * t * t - x**2, 0))
    mid = np.abs(x) < 0.25 * v * t; ins = np.abs(x) < 0.9 * v * t
    assert np.abs(S.delta[mid] - ref[mid]).max() < 6e-3 * ref[mid].max()   # dx = 100 m; 2e-3 at 25 m
    assert np.abs(S.delta[ins] - ref[ins]).max() < 1.5e-2 * ref[ins].max()
    assert np.abs(S.tau[np.abs(x) < v * t] + 1e6).max() < 1e-6
    # no normal-stress change on a planar fault
    assert np.abs(S.sigma - 1e8).max() < 1e-6 * 1e8


def test_qd_mode_static_limit():
    """QD mode with a slowly expanding crack approaches the static elliptical
    slip (2(1-nu) dsigma/mu) sqrt(a^2 - x^2)."""
    dx = 100.0; n = 120
    f = sh.Fault([-n * dx / 2, 0.0], [1.0, 0.0], n * dx, n)
    dt = 0.5 * dx / BE
    sys_ = sh.System([f], dt, AL, BE, RHO, tol=1e-2)
    v = 0.01 * BE; x = f.centers[:, 0]
    S = st.Stepper(sys_, st.PrescribedTip(x, v, 1e6), tau0=0.0, sigma0=1e8, mode="qd")
    for _ in range(int(0.6 * (n * dx / 2) / v / dt)):
        S.step()
    t = S.t; a = v * t
    mid = np.abs(x) < 0.5 * a
    ref = 2 * (1 - NU) / MU * 1e6 * np.sqrt(a * a - x[mid] ** 2)
    assert np.abs(S.delta[mid] - ref).max() < 2e-2 * ref.max()


def test_kinematic_factors_match_exact():
    """The dynamic factors obtained kinematically in test_shells (imposed
    elliptical slip, measured stress drop) against the exact values."""
    for vfrac in (0.5, 0.8):
        fac_file = f"/tmp/kostrov_factor_{vfrac}.npy"
        if not os.path.exists(fac_file):
            pytest.skip("run test_shells first")
        fac = np.load(fac_file)[1]
        assert abs(fac - kostrov_factor_II(vfrac)) < 4e-3


def test_tpv205_2d_coarse():
    """TPV205-2D at 100 m against Dunham's MDSBI 6.25 m reference: final
    slips within 2 percent, rupture arrivals late by a first-order
    discretisation delay (0.27, 0.23, 0.45 s at 100 m; halved at 50 m; the
    Richardson-extrapolated values match the reference to 1 ms)."""
    import run_tpv205
    out = run_tpv205.run(dx=100.0, verbose=False)
    ref_dir = os.path.join(os.path.dirname(__file__), "ref")
    delays = {"st-045": 0.271, "st-075": 0.229, "st120": 0.455}
    for name in ("st-045", "st-075", "st120"):
        r = np.load(os.path.join(ref_dir, f"tpv205_2d_dunham_{name}.npy"))
        o = out[name]
        ti = o[np.argmax(o[:, 2] > 1e-3), 0]; tr = r[np.argmax(r[:, 2] > 1e-3), 0]
        assert abs((ti - tr) - delays[name]) < 0.03
        assert abs(o[-1, 1] - r[-1, 1]) < 0.02 * r[-1, 1]
        assert abs(o[-1, 3] - r[-1, 3] * 1e6) < 1e6
