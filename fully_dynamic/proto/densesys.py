"""Dense 2-D boundary-element system for fault and virtual free-surface
elements (Hok and Fukuyama 2011 idea), for prototyping and tests; no
Toeplitz structure, so element counts of a few hundred.

Elements: centre X, tangent s, half-length a; kind 0 = fault (one unknown:
slip along s), kind 1 = surface (antiplane: one unknown; in-plane: slip
along s and opening, with zero shear and normal traction).
Tables (end-of-step collocation, as in the rest of the prototype):
  K[l] = S(l dt) - S((l-1) dt), W[l] = K[l] - Ks dt, l = 1..L,
rows = traction components of receivers (shear for every element; normal
for surface elements in-plane, the unknown's partner; fault normal rows are
kept separately for normal-stress changes), columns = unknowns.
Marching: first-order (end of step) or the second-order NL scheme of Noda
et al. (2020); the instantaneous coefficients are mu/(2 beta) for shear
and rho alpha/2 for opening."""
import numpy as np
import inplane as ip
import antiplane as ap
import inplane_open as io

def _open_static(Xs, ss, a, Xr, sr, mu, nu):
    import fs_static2d as fs
    g = fs.stress_global(Xs, ss, a, Xr, 1)
    return fs.traction_rows(g, sr)

class DenseSys:
    def __init__(self, X, s, a, kind, dt, alpha, beta, rho, L, mode="ip", verbose=False):
        self.X = np.asarray(X, float); self.s = np.asarray(s, float); self.a = np.asarray(a, float)
        self.kind = np.asarray(kind, int); self.dt = dt; self.mode = mode; self.L = L
        self.alpha, self.beta, self.rho = alpha, beta, rho
        mu = rho*beta**2; nu = (alpha**2 - 2*beta**2)/(2*(alpha**2 - beta**2)); self.mu = mu
        N = len(self.X)
        # unknown components: (element, 0 = slip, 1 = opening); traction rows likewise (0 shear, 1 normal)
        comps = [(i, 0) for i in range(N)]
        if mode == "ip": comps += [(i, 1) for i in range(N) if self.kind[i] == 1]
        self.comps = comps; Nc = len(comps)
        idx = {c: k for k, c in enumerate(comps)}
        self.fault = np.array([k for k, (i, c) in enumerate(comps) if self.kind[i] == 0])
        self.surf = np.array([k for k, (i, c) in enumerate(comps) if self.kind[i] == 1], int)
        self.W = np.zeros((L, Nc, Nc)); self.Ks = np.zeros((Nc, Nc))
        self.Wn = np.zeros((L, N, Nc)); self.Kn = np.zeros((N, Nc))     # normal traction on all elements (fault normal stress)
        for (j, cj) in comps:
            col = idx[(j, cj)]
            for sgrp in np.unique(self.s, axis=0):                    # receivers grouped by orientation
                rec = np.where(np.all(np.isclose(self.s, sgrp), axis=1))[0]
                Xr = self.X[rec]
                if mode == "ap":
                    Ws, Wn, Ks_s, Ks_n = ap.W_table_ap(self.X[j], self.s[j], self.a[j], Xr, sgrp, L, dt, alpha, beta, rho)
                elif cj == 0:
                    Ws, Wn, Ks_s, Ks_n = ip.W_table(self.X[j], self.s[j], self.a[j], Xr, sgrp, L, dt, alpha, beta, rho)
                else:
                    Ws, Wn, Ks_s, Ks_n = self._W_open(j, Xr, sgrp, L, mu, nu)
                for r, i in enumerate(rec):
                    self.W[:, idx[(i, 0)], col] = Ws[r]; self.Ks[idx[(i, 0)], col] = Ks_s[r]
                    if (i, 1) in idx:
                        self.W[:, idx[(i, 1)], col] = Wn[r]; self.Ks[idx[(i, 1)], col] = Ks_n[r]
                    self.Wn[:, i, col] = Wn[r]; self.Kn[i, col] = Ks_n[r]
        self.K1 = self.W[0] + self.Ks*dt
        self.Kn1 = self.Wn[0] + self.Kn*dt
        self.E = np.array([mu/(2*beta) if c == 0 else rho*alpha/2 for (i, c) in comps])   # instantaneous coefficients
        self.idx = idx; self.Nc = Nc; self.N = N

    def _W_open(self, j, Xr, sr, L, mu, nu):
        xr, yr = ip.receiver_in_source_frame(self.X[j], self.s[j], Xr)
        tl = np.concatenate([[0.0], np.arange(1, L + 1)*self.dt])
        s = [np.zeros((len(xr), L + 1)) for _ in range(3)]
        if not hasattr(self, "_ocache"): self._ocache = {}
        for r in range(len(xr)):                       # row by row (quadrature memory), cached by relative position
            key = (round(float(xr[r]), 6), round(float(yr[r]), 6), round(float(self.a[j]), 6), L)
            if key not in self._ocache:
                o = io.segment_ramp_stress_open(np.full(L + 1, xr[r]), np.full(L + 1, yr[r]), tl, self.a[j], self.alpha, self.beta, self.rho, tshift=1e-4*self.dt)
                for k in range(3): o[k][0] = 0.0                     # t = 0
                self._ocache[key] = o
            o = self._ocache[key]
            for k in range(3): s[k][r] = o[k]
        d = [np.diff(c, axis=1) for c in s]
        ts, tn = ip.project(d[0], d[1], d[2], self.s[j], sr)
        import fs_static2d as fs
        g = fs.stress_global(self.X[j], self.s[j], self.a[j], Xr, 1)
        Ks_s, Ks_n = fs.traction(g, sr)
        return ts - Ks_s[:, None]*self.dt, tn - Ks_n[:, None]*self.dt, Ks_s, Ks_n

class Marcher:
    """V: unknown rates (Nc). Fault unknowns either prescribed (kinematic, vf_fun(t)) or from
    slip weakening (mu_s, mu_d, d0, sigma0, tau0); surface unknowns from zero traction."""
    def __init__(self, sysd, nl=True, tau0=None, sigma0=None, sw=None, vf_fun=None, normal_coupling=True):
        self.d = sysd; self.nl = nl; Nc = sysd.Nc
        self.V = np.zeros(Nc); self.delta = np.zeros(Nc); self.t = 0.0; self.hist = []
        self.tau0 = np.zeros(Nc) if tau0 is None else tau0
        self.sigma0 = sigma0; self.sw = sw; self.vf_fun = vf_fun; self.nc = normal_coupling
        E = np.diag(sysd.E)
        if nl: self.B = 0.5*(sysd.K1 + E) - E; self.C = 0.5*(sysd.K1 + E)
        else: self.B = sysd.K1; self.C = np.zeros_like(sysd.K1)
        self.Bn = 0.5*sysd.Kn1 if nl else sysd.Kn1; self.Cn = 0.5*sysd.Kn1 if nl else 0.0*sysd.Kn1
        self.tau = np.zeros(Nc); self.sigma = None

    def _known(self):
        d = self.d
        H = self.tau0 + d.Ks @ self.delta; Hn = d.Kn @ self.delta
        for l in range(2, min(d.L, len(self.hist) + 1) + 1):
            Vh = self.hist[-(l - 1)]
            H = H + d.W[l - 1] @ Vh; Hn = Hn + d.Wn[l - 1] @ Vh
        return H + self.C @ self.V, Hn + self.Cn @ self.V

    def step(self):
        d = self.d; dt = d.dt; t1 = self.t + dt
        H, Hn = self._known()
        f, s = d.fault, d.surf
        Vold = self.V.copy(); V = self.V.copy()
        if self.vf_fun is not None:
            V[f] = self.vf_fun(t1, f)
            V[s] = np.linalg.solve(self.B[np.ix_(s, s)], -(H[s] + self.B[np.ix_(s, f)] @ V[f]))
        else:
            # slip weakening on the fault, zero traction on the surface: eliminate the surface
            Bss_inv = np.linalg.inv(self.B[np.ix_(s, s)])
            Hf = H[f] - self.B[np.ix_(f, s)] @ Bss_inv @ H[s]
            Af = self.B[np.ix_(f, f)] - self.B[np.ix_(f, s)] @ Bss_inv @ self.B[np.ix_(s, f)]
            mu_s, mu_d, d0 = self.sw
            hw = 0.5 if self.nl else 1.0
            dslip = self.delta[f] + (0.5*dt*Vold[f] if self.nl else 0.0)
            el = np.array([d.comps[k][0] for k in f])
            vf = Vold[f].copy()
            for it in range(200):
                # normal stress at the end of the step (compression positive)
                Vt = V.copy(); Vt[f] = vf
                Vt[s] = -Bss_inv @ (H[s] + self.B[np.ix_(s, f)] @ vf)
                sig = self.sigma0 - ((Hn + self.Bn @ Vt)[el] if self.nc else 0.0)
                sig = np.maximum(sig, 0.0)
                r = np.minimum(np.abs(dslip)/d0, 1.0)
                str0 = sig*(mu_s - (mu_s - mu_d)*r)
                slope = np.where(np.abs(dslip) < d0, sig*(mu_s - mu_d)*hw*dt/d0, 0.0)
                # Jacobi on the fault system with the eliminated operator
                diag = np.diag(Af)
                Hp = Hf + Af @ vf - diag*vf
                vn = np.where(np.abs(Hp) > str0, np.sign(Hp)*(np.abs(Hp) - str0)/(-diag - slope), 0.0)
                if np.max(np.abs(vn - vf)) < 1e-12*max(1.0, np.abs(vn).max()):
                    vf = vn; break
                vf = vf + 0.7*(vn - vf)
            V[f] = vf
            V[s] = -Bss_inv @ (H[s] + self.B[np.ix_(s, f)] @ vf)
        self.tau = H + self.B @ V
        self.sigma_dyn = Hn + self.Bn @ V
        Vm = 0.5*(Vold + V) if self.nl else V
        self.hist.append(Vm)
        if len(self.hist) > d.L: self.hist.pop(0)
        self.delta = self.delta + Vm*dt
        self.V = V; self.t = t1
        return V
