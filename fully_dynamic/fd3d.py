"""
3-D fully dynamic boundary-integral prototype: one planar fault (x3 = 0),
uniform square elements dx, piecewise-constant slip rate (2 components) in
space and in time steps dt, end-of-step collocation at element centres;
kernels from the quadrantal-fault Green's functions of Tada (2006, App. A),
see tada_quad.py. Same formulation as the 2-D code:

  tau(t_n) = tau0 + Ks delta(t_n) + sum_{l>=1} W[l] V^(n+1-l),  W[l] = K[l] - Ks dt,
  K[l] = S(l dt) - S((l-1) dt), S = response to unit slip rate from t = 0.

In 3-D the medium response to an element becomes static once the S wave from
its farthest corner has passed (Huygens), so W[l] = 0 exactly beyond that lag.
Tables are Toeplitz in the 2-D element offset and stored as 2-D FFTs of the
full linear-convolution size (2 nx - 1) x (2 nz - 1).
Slip components: 0 along x1 (strike), 1 along x2 (dip); tractions likewise.
"""
import numpy as np
import tada_quad as tq

class Tables3D:
    def __init__(self, nx, nz, dx, dt, cL, cT, rho, L=None, verbose=True, cdtype=complex):
        self.nx, self.nz, self.dx, self.dt = nx, nz, dx, dt
        self.cL, self.cT, self.rho = cL, cT, rho
        a = dx/2
        ox = np.arange(-(nx - 1), nx)*dx; oz = np.arange(-(nz - 1), nz)*dx
        X1, X2 = np.meshgrid(ox, oz, indexing="ij")              # receiver - source centre
        rmax = np.sqrt((np.abs(X1) + a)**2 + (np.abs(X2) + a)**2)
        self.Loff = np.ceil(rmax/cT/dt).astype(int) + 1          # lags until W vanishes
        if L is None: L = int(self.Loff.max()) + 1
        self.L = L
        Mx, Mz = 2*nx - 1, 2*nz - 1
        self.Mx, self.Mz = Mx, Mz
        # static per unit slip (late-time slope, exact once all waves have passed)
        T1 = (rmax.max()/cT)*1.5 + 10*dt; T2 = 2*T1
        Ks = np.zeros((2, 2, Mx, Mz))
        for c in (0, 1):
            s1 = tq.rect_ramp_stress(X1, X2, 0.0, np.full_like(X1, T1), a, a, cL, cT, rho, comp=c + 1)
            s2 = tq.rect_ramp_stress(X1, X2, 0.0, np.full_like(X1, T2), a, a, cL, cT, rho, comp=c + 1)
            Ks[c, 0] = (s2["s31"] - s1["s31"])/(T2 - T1); Ks[c, 1] = (s2["s32"] - s1["s32"])/(T2 - T1)
        self.Ks_off = Ks
        self.Ks_hat = np.stack([[np.fft.fft2(self._wrap(Ks[c, b])) for b in (0, 1)] for c in (0, 1)])
        # W tables: S at levels 0, dt, ..., (L) dt
        Sprev = np.zeros((2, 2, Mx, Mz))
        self.W_hat = np.zeros((L, 2, 2, Mx, Mz), dtype=cdtype); self.cdtype = cdtype
        self.K1 = None
        for l in range(1, L + 1):
            tl = np.full_like(X1, l*dt)
            S = np.zeros((2, 2, Mx, Mz))
            for c in (0, 1):
                s = tq.rect_ramp_stress(X1, X2, 0.0, tl, a, a, cL, cT, rho, comp=c + 1)
                S[c, 0] = s["s31"]; S[c, 1] = s["s32"]
            K = S - Sprev
            if l == 1: self.K1 = K.copy()
            W = K - Ks*dt
            for c in (0, 1):
                for b in (0, 1):
                    self.W_hat[l - 1, c, b] = np.fft.fft2(self._wrap(W[c, b]))
            Sprev = S
        if verbose:
            print(f"Tables3D: {nx} x {nz} elements of {dx:g} m, dt {dt:.4g} s, window {L} lags, "
                  f"{self.W_hat.nbytes/1e6:.0f} MB")

    def _wrap(self, A):
        """offset array (index 0 = offset -(n-1)) to FFT order (offset 0 at index 0)"""
        return np.roll(np.roll(A, -(self.nx - 1), axis=0), -(self.nz - 1), axis=1)

    def conv(self, Ahat, Vfield_hat):
        """apply an offset kernel (FFT) to a source field (FFT, zero padded): receiver field"""
        out = np.fft.ifft2(Ahat*Vfield_hat).real
        return out[:self.nx, :self.nz]


class Fault3D:
    def __init__(self, tab, tau0, sigma0, mu_s, mu_d, d0, slip_mask=None):
        self.t = tab; nx, nz = tab.nx, tab.nz
        self.tau0 = tau0.copy()          # (2, nx, nz)
        self.sigma0 = sigma0             # (nx, nz)
        self.mu_s, self.mu_d, self.d0 = mu_s, mu_d, d0
        self.delta = np.zeros((2, nx, nz)); self.V = np.zeros((2, nx, nz))
        self.hist = []                    # FFTs of past V (most recent last), each (2, Mx, Mz)
        self.n = 0; self.time = 0.0
        self.mask = np.ones((nx, nz), bool) if slip_mask is None else slip_mask
        # first-lag operator: 3 x 3 stencil (P reach within dt), offsets -1..1
        c0x, c0z = nx - 1, nz - 1
        self.k1 = tab.K1[:, :, c0x - 1:c0x + 2, c0z - 1:c0z + 2]     # (2, 2, 3, 3)
        far = np.abs(tab.K1).copy(); far[:, :, c0x - 1:c0x + 2, c0z - 1:c0z + 2] = 0
        assert far.max() < 1e-12*np.abs(tab.K1).max(), "first-lag kernel reaches beyond the 3 x 3 stencil"

    def _pad_hat(self, F):
        P = np.zeros((2, self.t.Mx, self.t.Mz)); P[:, :F.shape[1], :F.shape[2]] = F
        return np.fft.fft2(P).astype(self.t.cdtype)

    def known(self):
        """tau0 + Ks delta + history (lags >= 2) for the step about to be taken"""
        tab = self.t
        dh = self._pad_hat(self.delta)
        H = self.tau0.copy()
        acc = np.zeros((2, tab.Mx, tab.Mz), dtype=complex)
        for b in (0, 1):
            for c in (0, 1):
                acc[b] += tab.Ks_hat[c, b]*dh[c]
        nh = len(self.hist)
        for l in range(2, min(tab.L, nh + 1) + 1):        # lag l pairs with V^(n+1-l) = hist[-(l-1)]
            Vh = self.hist[-(l - 1)]; Wl = tab.W_hat[l - 1]
            for b in (0, 1):
                acc[b] += Wl[0, b]*Vh[0] + Wl[1, b]*Vh[1]
        out = np.fft.ifft2(acc).real[:, :tab.nx, :tab.nz]
        return H + out

    def apply_k1(self, V):
        """first-lag traction from V on the 3 x 3 stencil, excluding the self term"""
        out = np.zeros_like(V); nx, nz = V.shape[1:]
        for c in (0, 1):
            for b in (0, 1):
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        if di == 0 and dj == 0: continue
                        k = self.k1[c, b, di + 1, dj + 1]
                        if k == 0: continue
                        src = np.zeros((nx, nz))
                        src[max(di, 0):nx + min(di, 0), max(dj, 0):nz + min(dj, 0)] = \
                            V[c, max(-di, 0):nx - max(di, 0), max(-dj, 0):nz - max(dj, 0)]
                        out[b] += k*src
        return out

    def step(self, maxit=400, tol=1e-13):
        tab = self.t; dt = tab.dt
        H = self.known()
        k11 = self.k1[0, 0, 1, 1]; k22 = self.k1[1, 1, 1, 1]
        dmag = np.sqrt(self.delta[0]**2 + self.delta[1]**2)
        s = self.sigma0
        mu_now = np.where(dmag < self.d0, self.mu_s - (self.mu_s - self.mu_d)*dmag/self.d0, self.mu_d)
        slope = np.where(dmag < self.d0, s*(self.mu_s - self.mu_d)*dt/self.d0, 0.0)
        V = self.V.copy()
        for it in range(maxit):
            Hp = H + self.apply_k1(V)
            Vn = np.zeros_like(V)
            if getattr(self, "fixed_rake", False): Hp[1] = 0.0      # scalar slip along x1, friction on tau_1 only
            hmag = np.sqrt(Hp[0]**2 + Hp[1]**2)
            str0 = s*mu_now
            act = (hmag > str0) & self.mask
            # local vector solve: (Hp_c/(str0 + (|k_cc| - slope) v))^2 summed = 1 for v >= 0
            A1 = -k11 - slope; A2 = -k22 - slope
            vlo = np.zeros_like(hmag); vhi = np.where(act, (hmag - str0)/np.minimum(A1, A2) + 1.0, 0.0)
            for _ in range(80):
                vm = 0.5*(vlo + vhi)
                f = (Hp[0]/(str0 + A1*vm))**2 + (Hp[1]/(str0 + A2*vm))**2 - 1.0
                vlo = np.where(f > 0, vm, vlo); vhi = np.where(f > 0, vhi, vm)
            v = np.where(act, 0.5*(vlo + vhi), 0.0)
            # cap at d0 crossing: beyond d0 the strength is mu_d (handled next step); keep simple
            Vn[0] = np.where(act, v*Hp[0]/(str0 + A1*v), 0.0)
            Vn[1] = np.where(act, v*Hp[1]/(str0 + A2*v), 0.0)
            change = np.abs(Vn - V).max()
            V = Vn
            if change < tol*max(1.0, np.abs(V).max()): break
        self.V = V
        tau = H + self.apply_k1(V); tau[0] += k11*V[0]; tau[1] += k22*V[1]
        self.tau = tau
        self.hist.append(self._pad_hat(V))
        if len(self.hist) > tab.L: self.hist.pop(0)
        self.delta += V*dt
        self.n += 1; self.time += dt
        return it + 1
