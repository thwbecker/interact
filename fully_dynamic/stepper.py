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
    def __init__(self, x, v, dsigma, t0=0.0, allowed=None):
        self.x = np.asarray(x, float); self.v = v; self.dsigma = dsigma; self.t0 = t0
        self.allowed = allowed            # optional mask of elements that may ever slip

    def free_mask(self, t):
        m = np.abs(self.x) < self.v * (t - self.t0)
        return m if self.allowed is None else (m & self.allowed)

    def solve(self, H, K1, sigma, t, state, Vguess):
        free = self.free_mask(t)
        V = np.zeros_like(H)
        if free.any():
            idx = np.where(free)[0]
            A = K1[np.ix_(idx, idx)]
            rhs = -self.dsigma - H[idx]           # tau(t_{n+1}) - tau0 = -dsigma  (H already contains tau0-free part)
            V[idx] = np.linalg.solve(A, rhs)
        return V, state


def banded_from_dense(A):
    """Extract the band of a dense matrix in scipy solve_banded (ab) layout."""
    from scipy.linalg import solve_banded  # noqa
    N = A.shape[0]
    nz = np.nonzero(A)
    bw = int(np.max(np.abs(nz[0] - nz[1]))) if nz[0].size else 0
    ab = np.zeros((2 * bw + 1, N))
    for d in range(-bw, bw + 1):
        diag = np.diagonal(A, offset=d)
        if d >= 0:
            ab[bw - d, d:] = diag
        else:
            ab[bw - d, :N + d] = diag
    return ab, bw


def banded_solve(ab, bw, rhs):
    from scipy.linalg import solve_banded
    return solve_banded((bw, bw), ab, rhs)


class SlipWeakening:
    """Linear slip-weakening friction (TPV5 type), no cohesion, no back slip:
        strength = sigma * (mu_s - (mu_s - mu_d) min(delta/d0, 1)),
    with the strength evaluated at the end-of-step slip delta^n + V dt.
    Elements outside [xmin, xmax] (or where mu_s is huge) stay locked.
    Predictor-corrector: diagonal of K1 implicit, neighbours explicit."""
    def __init__(self, mu_s, mu_d, d0, locked=None, force_time=None, force_t0=0.0, dugdale=False):
        """dugdale: strength stays at mu_s sigma until slip reaches d0, then drops
        to mu_d sigma (Dugdale law; Dunham 2007). Otherwise linear weakening.
        force_time: per-element time at which a forced (nucleation) weakening
        begins; the strength then decreases linearly from mu_s to mu_d over
        force_t0 seconds regardless of slip (SCEC-style forced rupture)."""
        self.mu_s = np.asarray(mu_s, float); self.mu_d = np.asarray(mu_d, float); self.d0 = float(d0)
        self.locked = locked
        self.force_time = None if force_time is None else np.asarray(force_time, float)
        self.force_t0 = float(force_t0)
        self.dugdale = dugdale

    def _mu_forced(self, t):
        if self.force_time is None:
            return None
        f = np.clip((t - self.force_time) / max(self.force_t0, 1e-12), 0.0, 1.0)
        return self.mu_s - (self.mu_s - self.mu_d) * f

    def _solve_local(self, Hloc, k11, sigma, delta, dt):
        """Per element: find V >= 0 with Hloc + k11 V = strength(delta + V dt),
        or V = 0 if Hloc <= strength(delta) (k11 < 0)."""
        mus, mud, d0 = self.mu_s, self.mu_d, self.d0
        str0 = sigma * (mus - (mus - mud) * np.minimum(delta / d0, 1.0))
        slope = np.where(delta < d0, sigma * (mus - mud) * dt / d0, 0.0)   # strength decrease per unit V
        V = np.where(Hloc > str0, (Hloc - str0) / (-k11 + slope), 0.0)
        # if the step crosses d0, strength saturates at mu_d sigma
        cross = (delta < d0) & (delta + V * dt > d0) & (V > 0)
        Vc = (Hloc - sigma * mud) / (-k11)
        V = np.where(cross, np.maximum(Vc, 0.0), V)
        if self.locked is not None:
            V = np.where(self.locked, 0.0, V)
        return V

    def solve(self, H, K1, sigma_fun, t, state, Vguess, dt, delta, coupled=True, K1band=None):
        if not coupled:
            k11 = np.diag(K1); Koff = K1 - np.diag(k11)
            sig = sigma_fun(Vguess)
            Vs = self._solve_local(H + Koff @ Vguess, k11, sig, delta, dt)
            sig = sigma_fun(Vs)
            Vss = self._solve_local(H + Koff @ Vs, k11, sig, delta, dt)
            return Vss, state
        return self._solve_coupled(H, K1, sigma_fun, Vguess, dt, delta, K1band, t=t), state

    def _solve_coupled(self, H, K1, sigma_fun, Vguess, dt, delta, K1band, max_iter=12, t=None):
        import scipy.sparse as sps
        if not sps.issparse(K1):
            K1 = sps.csr_matrix(K1)
        """Exact solution of the coupled (banded) system with an active-set
        loop: rows of locked elements are replaced by V_i = 0. Within a
        step the strength is linear in V (slope) until delta reaches d0,
        after which it is mu_d sigma; regime changes are iterated."""
        ab0, bw = K1band
        N = H.size
        mus, mud, d0 = self.mu_s, self.mu_d, self.d0
        sig = sigma_fun(Vguess)
        if self.dugdale:
            mu_slip = np.where(delta < d0, mus, mud)
            slope0 = np.zeros_like(delta)
        else:
            mu_slip = mus - (mus - mud) * np.minimum(delta / d0, 1.0)
            slope0 = np.where(delta < d0, sig * (mus - mud) * dt / d0, 0.0)
        muf = self._mu_forced(t) if t is not None else None
        if muf is not None:
            forced = muf < mu_slip
            mu_slip = np.where(forced, muf, mu_slip)
            slope0 = np.where(forced, 0.0, slope0)
        str0 = sig * mu_slip
        # initial active set: trial stress with V = guess exceeds strength, or currently slipping
        active = (H + K1 @ Vguess > str0) | (Vguess > 0)
        saturate = np.zeros(N, bool)
        if self.locked is not None:
            active &= ~self.locked
        V = np.zeros(N)
        for it in range(max_iter):
            strength0 = np.where(saturate, sig * mud, str0)
            slope = np.where(saturate, 0.0, slope0)
            ab = ab0.copy()
            # diagonal row index in ab layout is bw
            ab[bw, :] = ab0[bw, :] + np.where(active, slope, 0.0)
            rhs = np.where(active, strength0 - H, 0.0)
            # inactive rows -> identity
            for d in range(-bw, bw + 1):
                if d == 0:
                    continue
                r = bw - d
                # entry (i, i+d) stored at ab[r, i+d]
                i = np.arange(max(0, -d), min(N, N - d))
                mask = ~active[i]
                ab[r, i[mask] + d] = 0.0
            ab[bw, :] = np.where(active, ab[bw, :], 1.0)
            Vn = banded_solve(ab, bw, rhs)
            Vn = np.where(active, Vn, 0.0)
            # checks: negative slip rates -> deactivate; locked elements overstressed -> activate
            tau = H + K1 @ Vn
            neg = active & (Vn < 0)
            over = (~active) & (tau > strength0)
            if self.locked is not None:
                over &= ~self.locked
            cross = active & (~saturate) & (delta < d0) & (delta + Vn * dt > d0)
            uncross = active & saturate & (delta < d0) & (delta + Vn * dt < d0)
            changed = neg | over | cross | uncross
            V = Vn
            if not changed.any():
                break
            active = (active & ~neg) | over
            saturate = (saturate | cross) & ~uncross
        return np.maximum(V, 0.0)


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
        with np.errstate(over="ignore", invalid="ignore"):
            return self.a * sigma * np.arcsinh(V / (2 * self.V0) * np.exp(psi / self.a))

    def dtau_dV(self, V, psi, sigma):
        with np.errstate(over="ignore", invalid="ignore"):
            z = V / (2 * self.V0) * np.exp(psi / self.a)
            return self.a * sigma * np.exp(psi / self.a) / (2 * self.V0) / np.sqrt(1 + z * z)

    def psidot(self, V, psi):
        with np.errstate(over="ignore", invalid="ignore"):
            return self.b / self.dc * (self.V0 * np.exp((self.f0 - psi) / self.b) - np.abs(V))

    def _solve_local(self, Hloc, k11, sigma, psi, Vg):
        """Elementwise solution of Hloc + k11 V = tau_f(V) for V > 0 (k11 < 0).
        g(u) = Hloc + k11 e^u - tau_f(e^u) is strictly decreasing in u = ln V,
        so a bracketed Newton converges from any start. Elements with
        Hloc <= 0 (no forward driving stress) get V = 0 for the regularized law's
        positive branch; back slip is not modelled here."""
        Hloc = np.asarray(Hloc, float); k11 = np.broadcast_to(k11, Hloc.shape)
        sigma = np.broadcast_to(sigma, Hloc.shape); psi = np.broadcast_to(psi, Hloc.shape)
        pos = Hloc > 0
        u = np.log(np.maximum(np.where(pos, Vg, 1e-30), 1e-300))
        lo = np.full(Hloc.shape, -700.0); hi = np.full(Hloc.shape, 50.0)
        def g(uu):
            V = np.exp(uu)
            return Hloc + k11 * V - self.tau_of_V(V, psi, sigma)
        for _ in range(100):
            gu = g(u)
            # update bracket (g decreasing: g > 0 means root is above u)
            lo = np.where(gu > 0, np.maximum(lo, u), lo)
            hi = np.where(gu <= 0, np.minimum(hi, u), hi)
            V = np.exp(u)
            dg = (k11 - self.dtau_dV(V, psi, sigma)) * V
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                un = u - gu / dg
            bad = ~np.isfinite(un) | (un <= lo) | (un >= hi)
            un = np.where(bad, 0.5 * (lo + hi), un)
            conv = np.all(np.abs(un - u)[pos] < self.tol) if pos.any() else True
            u = un
            if conv:
                break
        return np.where(pos, np.exp(u), 0.0)

    def solve(self, H, K1, sigma_fun, t, state, Vguess, dt, coupled=True, K1band=None):
        """Predictor-corrector in the state variable; the slip rate is solved
        from the coupled first-lag system (banded Newton in log V) when
        coupled=True, otherwise elementwise with the neighbour terms lagged.
        sigma_fun(V) returns the effective normal stress for a trial current
        slip rate. state = psi at t_n. Returns V^n, psi(t_{n+1})."""
        psi = state
        if coupled:
            solve_V = lambda Hh, sig, ps, Vg: self._solve_banded(Hh, K1, K1band, sig, ps, Vg)
            Koff = None
        else:
            k11 = np.diag(K1); Koff = K1 - np.diag(k11)
            solve_V = lambda Hh, sig, ps, Vg: self._solve_local(Hh + Koff @ Vg, k11, sig, ps, Vg)
        psi_s = self.advance_psi(psi, Vguess, Vguess, dt)
        Vs = solve_V(H, sigma_fun(Vguess), psi_s, Vguess)
        psi_ss = self.advance_psi(psi, Vguess, Vs, dt)
        Vss = solve_V(H, sigma_fun(Vs), psi_ss, Vs)
        return Vss, psi_ss

    def advance_psi(self, psi, V_a, V_b, dt):
        """Exact solution of the aging law over one step for a constant slip rate
        V = mean(|V_a|, |V_b|):  theta(t+dt) = dc/V + (theta - dc/V) exp(-V dt/dc),
        with theta = (dc/V0) exp((psi - f0)/b). Unconditionally stable, which the
        explicit Heun step is not once V dt / dc or the healing rate becomes large
        (Lapusta et al. 2000)."""
        V = np.maximum(0.5 * (np.abs(V_a) + np.abs(V_b)), 1e-30)
        with np.errstate(over="ignore", invalid="ignore"):
            theta = self.dc / self.V0 * np.exp((psi - self.f0) / self.b)
            theta_n = self.dc / V + (theta - self.dc / V) * np.exp(-V * dt / self.dc)
            return self.f0 + self.b * np.log(self.V0 * theta_n / self.dc)

    def _solve_banded(self, H, K1, K1band, sigma, psi, Vg):
        """Newton in u = ln V on g(u) = H + K1 e^u - tau_f(e^u) = 0 with the
        banded Jacobian (K1 - diag(dtau_f/dV)) diag(e^u)."""
        import scipy.sparse as sps
        ab0, bw = K1band
        K1s = K1 if sps.issparse(K1) else sps.csr_matrix(K1)
        u = np.clip(np.log(np.maximum(Vg, 1e-30)), -70.0, 10.0)
        converged = False
        for _ in range(self.max_iter):
            V = np.exp(u)
            g = H + K1s @ V - self.tau_of_V(V, psi, sigma)
            dgdV = -self.dtau_dV(V, psi, sigma)          # diagonal part beyond K1
            ab = ab0 * V[None, :]                          # column scaling by e^u
            # band layout: entry (i, j) at ab[bw + i - j, j]; diagonal row bw
            ab[bw, :] += dgdV * V
            if not np.all(np.isfinite(ab)) or not np.all(np.isfinite(g)):
                break
            du = -banded_solve(ab, bw, g)
            if not np.all(np.isfinite(du)):
                break
            du = np.clip(du, -2.0, 2.0)
            u = np.clip(u + du, -70.0, 10.0)
            if np.max(np.abs(du)) < self.tol:
                converged = True
                break
        if converged:
            return np.exp(u)
        # fallback: elementwise bracketed solve with the neighbour terms lagged (fixed point);
        # each inner solve is exact and monotone, the outer iteration converges because the
        # neighbour coupling is a fraction of the diagonal
        k11 = K1s.diagonal().copy()
        Koff = K1s - sps.diags(k11)
        V = np.exp(np.clip(u, -70.0, 10.0))
        for _ in range(60):
            Vn = self._solve_local(H + Koff @ V, k11, sigma, psi, V)
            if np.max(np.abs(Vn - V) / np.maximum(V, 1e-30)) < self.tol:
                V = Vn
                break
            V = Vn
        return V


class FlashRateState(RateState):
    """Strong rate weakening (flash heating) in the form of Noda, Dunham and
    Rice (2009): same regularized tau(V, psi), slip-law evolution towards a
    steady state that weakens strongly above V_w,
        f_ss(V)  = f_w + (f_LV(V) - f_w) / (1 + (V/V_w)^8)^(1/8),
        f_LV(V)  = f0 - (b - a) ln(V/V0),
        psi_ss(V) = a ln[(2 V0 / V) sinh(f_ss(V)/a)],
        dpsi/dt  = -(V/L)(psi - psi_ss(V)).
    No additional state variable; L plays the role of dc."""
    def __init__(self, a, b, dc, V0, f0, Vw=0.1, fw=0.2, newton_tol=1e-10, max_iter=50):
        super().__init__(a, b, dc, V0, f0, newton_tol, max_iter)
        self.Vw = float(Vw); self.fw = float(fw)

    def f_ss(self, V):
        V = np.maximum(np.abs(V), 1e-30)
        fLV = self.f0 - (self.b - self.a) * np.log(V / self.V0)
        return self.fw + (fLV - self.fw) / (1 + (V / self.Vw) ** 8) ** 0.125

    def psi_ss(self, V):
        V = np.maximum(np.abs(V), 1e-30)
        with np.errstate(over="ignore", invalid="ignore"):
            arg = self.f_ss(V) / self.a
            # ln(2 sinh(x)) = x + ln(1 - exp(-2x)), stable for large x
            ln2sinh = arg + np.log1p(-np.exp(-2 * arg))
            return self.a * (np.log(self.V0 / V) + ln2sinh)

    def psidot(self, V, psi):
        V = np.maximum(np.abs(V), 1e-30)
        return -(V / self.dc) * (psi - self.psi_ss(V))

    def advance_psi(self, psi, V_a, V_b, dt):
        """Exact relaxation towards psi_ss at the mean slip rate: the slip law
        is linear in psi with rate V/L, so this is unconditionally stable."""
        Vm = 0.5 * (np.abs(V_a) + np.abs(V_b))
        Vm = np.maximum(Vm, 1e-30)
        pss = self.psi_ss(Vm)
        return pss + (psi - pss) * np.exp(-Vm * dt / self.dc)


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
        self.Hn_dot = np.zeros(self.N)      # rate of the normal-stress change (e.g. backslip), physical sign
        if mode == "fd":
            thr = 1e-13 * np.abs(system.K1s).max()
            self.K1 = np.where(np.abs(system.K1s) > thr, system.K1s, 0.0)
            self.K1n = np.where(np.abs(system.K1n) > thr, system.K1n, 0.0)
        else:
            self.K1 = -system.eta * np.eye(self.N)
            self.K1n = np.zeros((self.N, self.N))
        self.K1band = banded_from_dense(self.K1)
        import scipy.sparse as sps
        self.K1sp = sps.csr_matrix(self.K1)           # banded: use sparse products in the solvers
        self.coupled = True
        self._ops = {1: (self.K1, self.K1n, self.K1band, self.K1sp)}
        self.k = 1

    def set_substeps(self, k):
        """Use a step of k dt_min (FD mode). Operators for each k are cached."""
        import scipy.sparse as sps
        k = int(k)
        if k not in self._ops:
            if self.mode == "fd":
                Ks_k, Kn_k = self.sys.cumulative_kernel(k)
            else:
                Ks_k, Kn_k = -self.sys.eta * np.eye(self.N), np.zeros((self.N, self.N))
            # tiny entries from FFT round trips would widen the band: prune
            thr = 1e-13 * np.abs(Ks_k).max()
            Ks_k = np.where(np.abs(Ks_k) > thr, Ks_k, 0.0)
            Kn_k = np.where(np.abs(Kn_k) > thr, Kn_k, 0.0)
            self._ops[k] = (Ks_k, Kn_k, banded_from_dense(Ks_k), sps.csr_matrix(Ks_k))
        self.K1, self.K1n, self.K1band, self.K1sp = self._ops[k]
        self.k = k

    def known_parts(self):
        dt = self.sys.dt * self.k
        t1 = self.t + dt
        H = self.tau0 + self.taudot * t1 + self.sys.Ks @ self.delta
        Hn = self.sys.Kn @ self.delta + self.Hn_dot * t1
        if self.mode == "fd":
            hs, hn = self.sys.history_stress(lag0=self.k + 1)
            H = H + hs; Hn = Hn + hn
        return H, Hn, t1

    def sigma_fun_factory(self, Hn):
        if not self.normal_coupling:
            return lambda V: self.sigma0
        return lambda V: self.sigma0 - (Hn + self.K1n @ V)

    def step(self):
        dt = self.sys.dt * self.k
        H, Hn, t1 = self.known_parts()
        if isinstance(self.fr, PrescribedTip):
            # H here includes tau0; the prescribed model wants tau(t_{n+1}) = tau0 - dsigma
            V, self.state = self.fr.solve(H - self.tau0, self.K1, None, t1, self.state, self.V)
        elif isinstance(self.fr, SlipWeakening):
            V, self.state = self.fr.solve(H, self.K1sp if self.coupled else self.K1, self.sigma_fun_factory(Hn), t1, self.state, self.V, dt, self.delta,
                                          coupled=self.coupled, K1band=self.K1band)
        else:
            V, self.state = self.fr.solve(H, self.K1sp if self.coupled else self.K1, self.sigma_fun_factory(Hn), t1, self.state, self.V, dt,
                                          coupled=self.coupled, K1band=self.K1band)
        self.tau = H + self.K1sp @ V
        self.sigma = self.sigma0 - (Hn + self.K1n @ V)
        if self.mode == "fd":
            for _ in range(self.k):
                self.sys.push(V)
        self.V = V
        self.delta = self.delta + V * dt
        self.t = t1; self.n += 1
        return V
