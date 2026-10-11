"""
WP4: sequences of earthquakes with a quasi-dynamic (QD) adaptive integrator
for the slow phases and the fixed-step fully dynamic (FD) stepper for events,
on the same static operators, friction law and loading.

Loading (backslip): tau(t) = tau0 + Ks (delta - V_pl t) [+ history],
                    sigma(t) = sigma0 - Kn (delta - V_pl t) [+ history].

QD phase: y = [delta, psi]; V from tau = tau_f(V, psi, sigma) with
tau = tau0 + Ks(delta - V_pl t) - eta V (elementwise Newton), integrated with
scipy's DOP853 and an event on max V = V_sw.

Hand-off QD -> FD: the dense output of the QD solution is sampled on the
fine grid at t_sw - l dt, l = L..1, slip rates recomputed from the QD
closure and pushed into the history buffers (the correction they produce is
of order V_slow, negligible, but the buffer is then consistent).

Hand-off FD -> QD: after max V has stayed below V_back for at least the
longest window L dt, so that the history term has flushed, the FD state
(delta, psi) seeds the QD integrator.

Note: rsf_solve integrates tau as a state variable; here delta is the state
and tau is algebraic, which is equivalent for these operators and avoids the
impulsive terms the hereditary kernel would put into dtau/dt.
"""
import numpy as np
from scipy.integrate import solve_ivp
import shells as sh
import stepper as st


class CycleModel:
    def __init__(self, system, friction, tau0, sigma0, Vpl, normal_coupling=True,
                 V_sw=1e-3, V_back=1e-4, rtol=1e-8, atol_delta=1e-9, atol_psi=1e-8, verbose=True,
                 xi=0.25, kmax=64, dynamic=True):
        self.xi, self.kmax = xi, kmax
        self.dynamic = dynamic
        self.record = 0            # if > 0: store snapshots every `record` steps and every step above 1e-2 m/s
        self.sys = system; self.fr = friction
        self.N = system.N
        self.tau0 = np.broadcast_to(np.asarray(tau0, float), (self.N,)).copy()
        self.sigma0 = np.broadcast_to(np.asarray(sigma0, float), (self.N,)).copy()
        self.Vpl = float(Vpl)
        self.normal_coupling = normal_coupling
        self.V_sw, self.V_back = V_sw, V_back
        self.rtol, self.atol = rtol, np.concatenate([np.full(self.N, atol_delta), np.full(self.N, atol_psi)])
        self.eta = system.eta
        self.verbose = verbose
        self.log = []          # (t, max V, phase)
        self.events = []       # dicts with t_start, t_end, max V, moment-like slip

    # ---- QD closure ----
    def qd_rates(self, t, y, Vguess=None):
        d, psi = y[:self.N], y[self.N:]
        back = d - self.Vpl * t
        tau = self.tau0 + self.sys.Ks @ back
        sig = self.sigma0 - (self.sys.Kn @ back) if self.normal_coupling else self.sigma0
        if Vguess is None:
            Vguess = np.full(self.N, self.Vpl)
        V = self.fr._solve_local(tau, -self.eta * np.ones(self.N), sig, psi, Vguess)
        return V, self.fr.psidot(V, psi), tau - self.eta * V, sig

    def _rhs(self, t, y):
        V, psd, _, _ = self.qd_rates(t, y, self._Vlast)
        self._Vlast = V
        return np.concatenate([V, psd])

    def run_qd(self, y0, t0, t_max):
        self._Vlast = np.full(self.N, self.Vpl)
        def ev(t, y):
            V, _, _, _ = self.qd_rates(t, y, self._Vlast)
            return np.log(V.max() / self.V_sw)
        ev.terminal = True; ev.direction = 1
        sol = solve_ivp(self._rhs, (t0, t_max), y0, method="DOP853", rtol=self.rtol, atol=self.atol,
                        events=ev, dense_output=True, max_step=t_max - t0)
        return sol

    # ---- hand-offs ----
    def seed_history_from_qd(self, sol, t_sw):
        dt = self.sys.dt
        self.sys.reset_history()
        L = self.sys.L
        for l in range(L, 0, -1):
            t = t_sw - l * dt
            if t < sol.t[0]:
                V = np.zeros(self.N)
            else:
                y = sol.sol(t)
                V, _, _, _ = self.qd_rates(t, y)
            self.sys.push(V)

    def make_fd_stepper(self, y, t, V):
        d, psi = y[:self.N], y[self.N:]
        S = st.Stepper(self.sys, self.fr, tau0=self.tau0, sigma0=self.sigma0,
                       taudot=-(self.sys.Ks @ np.full(self.N, self.Vpl)), mode="fd" if self.dynamic else "qd",
                       normal_coupling=self.normal_coupling)
        # backslip also enters the normal stress: sigma = sigma0 - Kn(delta - Vpl t)
        S.Hn_dot = -(self.sys.Kn @ np.full(self.N, self.Vpl))
        S.delta = d.copy(); S.state = psi.copy(); S.V = V.copy(); S.t = t
        return S

    def run(self, y0, t_end, max_events=None):
        t = 0.0; y = np.asarray(y0, float).copy()
        while t < t_end:
            sol = self.run_qd(y, t, t_end)
            tq = sol.t[-1]; y = sol.y[:, -1]
            Vq, _, _, _ = self.qd_rates(tq, y)
            self.log.append((tq, Vq.max(), "qd"))
            if sol.status != 1:          # reached t_end without an event
                t = tq
                break
            if self.verbose:
                print(f"  QD -> FD at t = {tq/3.15576e7:10.4f} yr (max V = {Vq.max():.2e})")
            self.seed_history_from_qd(sol, tq)
            S = self.make_fd_stepper(y, tq, Vq)
            ev = dict(t_start=tq, Vmax=0.0, slip0=S.delta.copy(), snaps=[])
            below = 0.0; n = 0
            dcmin = np.min(self.fr.dc)
            Vmax = Vq.max()
            while True:
                # Lapusta-style step choice: dt <= xi dc / Vmax, as a power of two times dt_min
                k = int(min(self.kmax, max(1, 2 ** np.floor(np.log2(max(1.0, self.xi * dcmin / (Vmax * self.sys.dt)))))))
                S.set_substeps(k)
                V = S.step()
                Vmax = V.max()
                n += 1
                if self.record and (n % self.record == 0 or Vmax > 1e-2):
                    ev["snaps"].append((S.t, V.copy(), S.tau.copy()))
                ev["Vmax"] = max(ev["Vmax"], V.max())
                if V.max() < self.V_back:
                    below += self.sys.dt * k
                else:
                    below = 0.0
                if n % 200 == 0:
                    self.log.append((S.t, V.max(), "fd"))
                    if self.verbose:
                        print(f"    fd step {n}: t-t_sw = {S.t - tq:8.2f} s, k = {k:2d}, max V = {Vmax:.3e}")
                    if self.verbose:
                        print(f"    fd step {n}: t - t_sw = {S.t - tq:7.2f} s, max V = {V.max():.3e}, below for {below:.1f} s", flush=True)
                if below >= self.sys.L * self.sys.dt:
                    break
            ev["t_end"] = S.t; ev["slip"] = S.delta - ev["slip0"]; ev["nsteps"] = n
            self.events.append(ev)
            if self.verbose:
                print(f"  FD -> QD at t = {S.t/3.15576e7:10.4f} yr after {n} steps; max V {ev['Vmax']:.2f} m/s, mean slip {ev['slip'].mean():.3f} m")
            y = np.concatenate([S.delta, S.state]); t = S.t
            if max_events is not None and len(self.events) >= max_events:
                break
        return y, t
