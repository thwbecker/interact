"""
WP3: fixed-step dynamic time marching on a shells.System.

At each step n (slip rate V^n during [t_n, t_{n+1}), collocation at t_{n+1}):

    tau(t_{n+1}) = H + K1 V^n,
    H = tau0 + taudot t_{n+1} + Ks delta^n + sum_{l>=2} W[l] V^{n+1-l},

    sigma(t_{n+1}) = sigma0 - (Kn delta^n + sum_{l>=2} Wn[l] V^{n+1-l} + K1n V^n)
                   (effective normal stress, compression positive; the K1n V^n
                    term is taken from the predictor),

and the friction law closes the system for V^n. K1 = W[1] + Ks dt couples
each element to the few neighbours reached by the P wave within one step;
the stepper treats the diagonal implicitly and the off-diagonal part with a
predictor, except for the linear prescribed-tip model which solves the
coupled system exactly.

Quasi-dynamic mode (mode="qd"): the history term is dropped and K1 is
replaced by -eta on the diagonal (radiation damping), which is the same
model as rsf_solve's QD formulation on a fixed step.
"""
import numpy as np
import shells as sh


class PrescribedTip:
    """Frictionless crack with uniform stress drop dsigma expanding at speed v
    from x = 0 (self-similar benchmark); elements outside the tip are locked."""
    def __init__(self, x, v, dsigma, t0=0.0):
        self.x = np.asarray(x, float); self.v = v; self.dsigma = dsigma; self.t0 = t0

    def free_mask(self, t):
        return np.abs(self.x) < self.v * (t - self.t0)

    def solve(self, H, K1, sigma, t, state, Vguess):
        free = self.free_mask(t)
        V = np.zeros_like(H)
        if free.any():
            idx = np.where(free)[0]
            A = K1[np.ix_(idx, idx)]
            rhs = -self.dsigma - H[idx]           # tau(t_{n+1}) - tau0 = -dsigma  (H already contains tau0-free part)
            V[idx] = np.linalg.solve(A, rhs)
        return V, state


class RateState:
    """Regularized rate-and-state friction with the aging law, psi form
    (as in rsf_solve):
        V = 2 V0 exp(-psi/a) sinh(tau/(a sigma))
        dpsi/dt = (b/dc) (V0 exp((f0 - psi)/b) - |V|)
    State array: psi."""
    def __init__(self, a, b, dc, V0, f0, newton_tol=1e-10, max_iter=50):
        self.a = np.asarray(a, float); self.b = np.asarray(b, float); self.dc = np.asarray(dc, float)
        self.V0 = V0; self.f0 = f0; self.tol = newton_tol; self.max_iter = max_iter

    def tau_of_V(self, V, psi, sigma):
        return self.a * sigma * np.arcsinh(V / (2 * self.V0) * np.exp(psi / self.a))

    def dtau_dV(self, V, psi, sigma):
        z = V / (2 * self.V0) * np.exp(psi / self.a)
        return self.a * sigma * np.exp(psi / self.a) / (2 * self.V0) / np.sqrt(1 + z * z)

    def psidot(self, V, psi):
        return self.b / self.dc * (self.V0 * np.exp((self.f0 - psi) / self.b) - np.abs(V))

    def _solve_local(self, Hloc, k11, sigma, psi, Vg):
        """Newton on g(V) = Hloc + k11 V - tau_f(V) = 0, elementwise; k11 < 0 and
        tau_f increasing in V, so g is monotone decreasing: unique root."""
        V = np.maximum(Vg, 1e-30)
        for _ in range(self.max_iter):
            g = Hloc + k11 * V - self.tau_of_V(V, psi, sigma)
            dg = k11 - self.dtau_dV(V, psi, sigma)
            dV = -g / dg
            # damp steps in log space to keep V positive and robust
            Vn = V + dV
            bad = Vn <= 0.1 * V
            Vn = np.where(bad, 0.1 * V, Vn)
            Vn = np.where(Vn > 10 * V, 10 * V, Vn)
            conv = np.max(np.abs(Vn - V) / np.maximum(V, 1e-30)) < self.tol
            V = Vn
            if conv:
                break
        return V

    def solve(self, H, K1, sigma_fun, t, state, Vguess, dt):
        """Predictor-corrector. sigma_fun(V) returns the effective normal stress
        for a trial current slip rate. state = psi at t_n. Returns V^n, psi(t_{n+1})."""
        psi = state
        k11 = np.diag(K1)
        Koff = K1 - np.diag(k11)
        # predictor
        psd0 = self.psidot(Vguess, psi)
        psi_s = psi + dt * psd0
        sig_s = sigma_fun(Vguess)
        Vs = self._solve_local(H + Koff @ Vguess, k11, sig_s, psi_s, Vguess)
        # corrector
        psi_ss = psi + 0.5 * dt * (psd0 + self.psidot(Vs, psi_s))
        sig_ss = sigma_fun(Vs)
        Vss = self._solve_local(H + Koff @ Vs, k11, sig_ss, psi_ss, Vs)
        return Vss, psi_ss


class Stepper:
    def __init__(self, system, friction, tau0, sigma0, taudot=0.0, mode="fd", normal_coupling=True):
        self.sys = system; self.fr = friction
        self.N = system.N
        self.tau0 = np.broadcast_to(np.asarray(tau0, float), (self.N,)).copy()
        self.sigma0 = np.broadcast_to(np.asarray(sigma0, float), (self.N,)).copy()
        self.taudot = np.broadcast_to(np.asarray(taudot, float), (self.N,)).copy()
        self.mode = mode
        self.normal_coupling = normal_coupling
        self.delta = np.zeros(self.N)
        self.V = np.zeros(self.N)
        self.t = 0.0; self.n = 0
        self.state = None
        if mode == "fd":
            self.K1 = system.K1s
            self.K1n = system.K1n
        else:
            self.K1 = -system.eta * np.eye(self.N)
            self.K1n = np.zeros((self.N, self.N))

    def known_parts(self):
        dt = self.sys.dt
        t1 = self.t + dt
        H = self.tau0 + self.taudot * t1 + self.sys.Ks @ self.delta
        Hn = self.sys.Kn @ self.delta
        if self.mode == "fd":
            hs, hn = self.sys.history_stress(lag0=2)
            H = H + hs; Hn = Hn + hn
        return H, Hn, t1

    def sigma_fun_factory(self, Hn):
        if not self.normal_coupling:
            return lambda V: self.sigma0
        return lambda V: self.sigma0 - (Hn + self.K1n @ V)

    def step(self):
        dt = self.sys.dt
        H, Hn, t1 = self.known_parts()
        if isinstance(self.fr, PrescribedTip):
            # H here includes tau0; the prescribed model wants tau(t_{n+1}) = tau0 - dsigma
            V, self.state = self.fr.solve(H - self.tau0, self.K1, None, t1, self.state, self.V)
        else:
            V, self.state = self.fr.solve(H, self.K1, self.sigma_fun_factory(Hn), t1, self.state, self.V, dt)
        self.tau = H + self.K1 @ V
        self.sigma = self.sigma0 - (Hn + self.K1n @ V)
        if self.mode == "fd":
            self.sys.push(V)
        self.V = V
        self.delta = self.delta + V * dt
        self.t = t1; self.n += 1
        return V
