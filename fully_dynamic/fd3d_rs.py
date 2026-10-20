"""Rate-and-state (aging law, regularized) friction for the 3-D prototype
(fd3d.py): vector slip rate parallel to the shear traction, exact aging-law
state update over a step at the mean slip rate, predictor-corrector in the
state, coupled first-lag solve by Jacobi iteration with an exact local
solve (Newton in ln v with bracketing) along the traction direction."""
import numpy as np
import fd3d

class FaultRS3D(fd3d.Fault3D):
    def __init__(self, tab, tau0, sigma0, a, b, Lc, f0, V0, theta0, Vini, load=None):
        super().__init__(tab, tau0, sigma0, 0.0, 0.0, 1.0)
        self.a, self.b, self.Lc, self.f0, self.V0 = a, b, Lc, f0, V0
        self.theta = theta0.copy()
        self.V = np.zeros((2,) + a.shape); self.V[0] = Vini
        self.load = load                  # load(t) -> (2, nx, nz) traction perturbation

    def strength(self, v, theta):
        s = self.sigma0
        return self.a*s*np.arcsinh(v/(2*self.V0)*np.exp((self.f0 + self.b*np.log(self.V0*theta/self.Lc))/self.a))

    def dstrength(self, v, theta):
        s = self.sigma0
        e = np.exp((self.f0 + self.b*np.log(self.V0*theta/self.Lc))/self.a)/(2*self.V0)
        z = v*e
        return self.a*s*e/np.sqrt(1 + z*z)

    def theta_adv(self, theta, Va, Vb, dt):
        V = 0.5*(Va + Vb); V = np.maximum(V, 1e-30)
        return self.Lc/V + (theta - self.Lc/V)*np.exp(-V*dt/self.Lc)

    def local(self, Hp, theta, k11, k22):
        """v >= 0 with sum_c (Hp_c/(F(v) + |k_cc| v))^2 = 1; returns V vector"""
        hm = np.sqrt(Hp[0]**2 + Hp[1]**2)
        A1, A2 = -k11, -k22
        lo = np.full(hm.shape, -80.0); hi = np.log(np.maximum(hm/min(A1, A2), 1e-300)) + 1.0
        u = np.minimum(np.log(np.maximum(np.abs(self.V[0]) + np.abs(self.V[1]), 1e-30)), hi - 0.5)
        u = np.maximum(u, lo + 1)
        for _ in range(60):
            v = np.exp(u)
            F = self.strength(v, theta); dF = self.dstrength(v, theta)
            d1 = F + A1*v; d2 = F + A2*v
            f = (Hp[0]/d1)**2 + (Hp[1]/d2)**2 - 1.0
            df = (-2*Hp[0]**2/d1**3*(dF + A1) - 2*Hp[1]**2/d2**3*(dF + A2))*v      # d f / d u
            lo = np.where(f > 0, u, lo); hi = np.where(f > 0, hi, u)
            un = u - f/np.where(df != 0, df, -1.0)
            bad = ~np.isfinite(un) | (un <= lo) | (un >= hi)
            un = np.where(bad, 0.5*(lo + hi), un)
            if np.max(np.abs(un - u)) < 1e-13: u = un; break
            u = un
        v = np.exp(u)
        F = self.strength(v, theta)
        return np.stack([v*Hp[0]/(F + A1*v), v*Hp[1]/(F + A2*v)])

    def solve(self, H, theta, Vg, maxit=400, tol=1e-13):
        k11 = self.k1[0, 0, 1, 1]; k22 = self.k1[1, 1, 1, 1]
        V = Vg.copy()
        for it in range(maxit):
            Hp = H + self.apply_k1(V)
            if getattr(self, "fixed_rake", False): Hp[1] = 0.0
            Vn = self.local(Hp, theta, k11, k22)
            ch = np.abs(Vn - V).max()
            V = Vn
            if ch <= tol*max(np.abs(V).max(), 1.0): break
        return V, it + 1

    def step(self):
        tab = self.t; dt = tab.dt; t1 = self.time + dt
        H = self.known()
        if self.load is not None: H = H + self.load(t1)
        vmag = np.sqrt(self.V[0]**2 + self.V[1]**2)
        th_p = self.theta_adv(self.theta, vmag, vmag, dt)
        Vs, i1 = self.solve(H, th_p, self.V)
        vs = np.sqrt(Vs[0]**2 + Vs[1]**2)
        th_c = self.theta_adv(self.theta, vmag, vs, dt)
        V, i2 = self.solve(H, th_c, Vs)
        self.theta = th_c; self.V = V
        k11 = self.k1[0, 0, 1, 1]; k22 = self.k1[1, 1, 1, 1]
        tau = H + self.apply_k1(V); tau[0] += k11*V[0]; tau[1] += k22*V[1]
        self.tau = tau
        self.hist.append(self._pad_hat(V))
        if len(self.hist) > tab.L: self.hist.pop(0)
        self.delta += V*dt; self.n += 1; self.time = t1
        return i1 + i2
