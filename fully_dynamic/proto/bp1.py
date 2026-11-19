"""
SEAS BP1 (antiplane, vertical strike-slip fault in a half-space, rate-and-state
with the aging law) with the space-time BIEM prototype, QD and FD.

Half-space by the image method: the fault is doubled to z in [-Wf, Wf] in a
full space and the mirror half carries the same slip as the real half. The
proxy below exposes the symmetrised operators (N real degrees of freedom) to
stepper.Stepper and cycle.CycleModel.

Depth z is the along-fault coordinate (positive down); stations at 2.5, 7.5,
12.5 km match the uploaded reference series.
"""
import numpy as np
import shells as sh, stepper as st, cycle as cy

YR = 3.15576e7
P = dict(rho=2670.0, beta=3464.0, sigma=50e6, a0=0.010, amax=0.025, b=0.015, Dc=0.008,
         Vp=1e-9, Vinit=1e-9, V0=1e-6, f0=0.6, H=15e3, h=3e3, Wf=40e3)


class HalfSpaceAntiplane:
    """Symmetrised view of a doubled antiplane fault."""
    def __init__(self, full, N):
        self.full = full; self.N = N
        self.dt = full.dt; self.L = full.L; self.eta = full.eta; self.mu = full.mu
        self.beta = full.beta; self.nu = full.nu
        R = slice(N, 2 * N)
        mir = np.arange(N - 1, -1, -1)            # full index of the mirror of real j is N-1-j
        self.Ks = full.Ks[R, R] + full.Ks[R, :][:, mir]
        self.Kn = np.zeros((N, N))
        self.K1s = full.K1s[R, R] + full.K1s[R, :][:, mir]
        self.K1n = np.zeros((N, N))
        self._R, self._mir = R, mir
        self._Kk = {1: (self.K1s, self.K1n)}
        self.nbytes = full.nbytes

    def _mirror(self, V):
        return np.concatenate([V[::-1], V])

    def push(self, V):
        self.full.push(self._mirror(V))

    def reset_history(self):
        self.full.reset_history()

    def history_stress(self, lag0=1):
        ts, tn = self.full.history_stress(lag0=lag0)
        return ts[self._R], np.zeros(self.N)

    def cumulative_kernel(self, k):
        if k not in self._Kk:
            Ks_k, _ = self.full.cumulative_kernel(k)
            R, mir = self._R, self._mir
            self._Kk[k] = (Ks_k[R, R] + Ks_k[R, :][:, mir], np.zeros((self.N, self.N)))
        return self._Kk[k]


def build(dz=50.0, tol=1e-2, V_sw=1e-3, V_back=1e-4, dynamic=True, lmax=None):
    Wf = P["Wf"]; N = int(Wf / dz)
    f = sh.Fault([0.0, -Wf], [0.0, 1.0], 2 * Wf, 2 * N, "bp1")     # z from -Wf to +Wf, real half is z > 0
    dt = 0.5 * dz / P["beta"]
    full = sh.System([f], dt, 2 * P["beta"], P["beta"], P["rho"], tol=tol, mode="antiplane",
                     lmax=lmax, cap_ok=lmax is not None)
    sys_ = HalfSpaceAntiplane(full, N)
    z = f.centers[N:, 1]
    a = np.where(z < P["H"], P["a0"], np.where(z < P["H"] + P["h"], P["a0"] + (P["amax"] - P["a0"]) * (z - P["H"]) / P["h"], P["amax"]))
    b = np.full(N, P["b"])
    fr = st.RateState(a, b, P["Dc"], P["V0"], P["f0"])
    eta = sys_.eta
    # initial shear stress: steady state with Vinit at the depth of Wf (a = amax), eq. (14)
    Vi = P["Vinit"]; sig = P["sigma"]
    tau0 = sig * P["amax"] * np.arcsinh(Vi / (2 * P["V0"]) * np.exp((P["f0"] + P["b"] * np.log(P["V0"] / Vi)) / P["amax"])) + eta * Vi
    # initial state consistent with tau0 and Vinit for the local a, eq. (15) in psi form
    psi0 = a * np.log(2 * P["V0"] / Vi * np.sinh((tau0 - eta * Vi) / (a * sig)))
    model = cy.CycleModel(sys_, fr, np.full(N, tau0), sig, P["Vp"], normal_coupling=False,
                          V_sw=V_sw, V_back=V_back, dynamic=dynamic, verbose=False)
    y0 = np.concatenate([np.zeros(N), psi0])
    return model, y0, z, fr


def station_series(model, z, depths=(2.5e3, 7.5e3, 12.5e3)):
    """Interpolate recorded snapshots (FD) and log entries to station depths."""
    idx = {d: int(np.argmin(np.abs(z - d))) for d in depths}
    return idx
