import numpy as np, pytest, os
import shells as sh, stepper as st

AL, BE, RHO = 6000.0, 3464.0, 2700.0
MU = RHO * BE**2
NU = (AL**2 - 2 * BE**2) / (2 * (AL**2 - BE**2))

# reference dynamic factors C(v)/C(0) of the self-similar antiplane crack from the
# independent spectral solver (spectral_ap.py), converged to about 1e-3
SPECTRAL_AP = {0.5: 0.8258, 0.8: 0.7042}


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
    assert abs(C - SPECTRAL_AP[vfrac]) < 3e-3
    # elliptical profile over the inner 90 percent of the crack
    ins = np.abs(x) < 0.9 * v * t
    ref = C * (2 / MU * 1e6) * np.sqrt(v * v * t * t - x[ins] ** 2)
    assert np.abs(S.delta[ins] - ref).max() < 1e-2 * ref.max()
    # stress on the free part equals the prescribed drop; locked part does not slip
    assert np.abs(S.tau[np.abs(x) < v * t] + 1e6).max() < 1e-6
    assert np.all(S.V[np.abs(x) >= v * t] == 0)


@pytest.mark.parametrize("vfrac", [0.5, 0.8])
def test_self_similar_crack_inplane_consistent_with_kinematic(vfrac):
    """In-plane: the stepper's dynamic factor must agree with the one obtained
    kinematically (test_shells), and the slip profile must be elliptical.
    The two routes share the kernels, so this is a consistency check; the
    time-marching convention itself is validated by the antiplane test."""
    fac_file = f"/tmp/kostrov_factor_{vfrac}.npy"
    if not os.path.exists(fac_file):
        pytest.skip("run test_shells first")
    fac = np.load(fac_file)[1]
    S, x, v = run_kostrov("inplane", vfrac)
    t = S.t
    Cstat = 2 * (1 - NU) / MU
    ref = fac * Cstat * 1e6 * np.sqrt(np.maximum(v * v * t * t - x**2, 0))
    mid = np.abs(x) < 0.25 * v * t; ins = np.abs(x) < 0.9 * v * t
    assert np.abs(S.delta[mid] - ref[mid]).max() < 5e-3 * ref[mid].max()
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
