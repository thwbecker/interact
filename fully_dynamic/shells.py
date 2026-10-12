"""
WP2: discretised fault system, static operators, windowed dynamic-correction
tables and the history convolution.

Formulation (velocity form, end-of-step collocation, e_t = 1):

    tau_i(t_{n+1}) = tau0_i + load_i(t_{n+1})
                     + sum_j Ks[i,j] delta_j(t_{n+1})
                     + sum_{l=1}^{L} sum_j W[l][i,j] V_j^{n+1-l}

with W[l] = K[l] - Ks dt, K[l] the lag kernel of inplane.lag_kernel (stress
at l*dt from unit slip rate during [0,dt)), delta_j(t_{n+1}) = delta_j(t_n)
+ V_j^n dt and V_j^m the slip rate during [t_m, t_{m+1}). The l = 1 term
contains the current step's slip rate; its diagonal is close to -mu/(2 beta)
and is what the stepper treats implicitly.

Window: per fault pair, L is the last lag at which |W[l]| exceeds tol times
the static interaction envelope mu dx^2 dt / (2 pi (1-nu) r^2) for any
element pair (an envelope rather than the static value itself, because the
radiation pattern has nodal directions). Beyond L the correction is
dropped (2-D kernels decay algebraically; see fd_plan.md section 3).

Storage: for parallel pairs with identical element size and tangent, W[l]
is Toeplitz in i - j and is stored as its FFT along the padded offset axis,
so one history-convolution step costs O(L M) complex multiply-adds plus one
inverse FFT per receiver fault. Non-parallel pairs use dense tables (small
problems only).
"""
import numpy as np
import inplane as ip


class Fault:
    """Straight segment from X0 along unit tangent s, length Lf, n elements."""
    def __init__(self, X0, s, Lf, n, name=""):
        self.X0 = np.asarray(X0, float)
        self.s = np.asarray(s, float) / np.linalg.norm(s)
        self.n = int(n)
        self.Lf = float(Lf)
        self.dx = self.Lf / self.n
        self.a = self.dx / 2
        self.normal = ip.rot90(self.s)
        self.centers = self.X0[None, :] + (np.arange(self.n) + 0.5)[:, None] * self.dx * self.s[None, :]
        self.name = name


def _parallel(fs, fr, tol=1e-12):
    return abs(fs.dx - fr.dx) < tol * fs.dx and np.dot(fs.s, fr.s) > 1 - tol


class PairTable:
    """Dynamic-correction table W[l] for one (source, receiver) fault pair."""

    def __init__(self, fs, fr, dt, alpha, beta, rho, tol=1e-2, lmax=None, M=None, mode="inplane", cap_ok=False, et=1.0):
        self.fs, self.fr, self.dt = fs, fr, dt
        mu = rho * beta**2
        nu = (alpha**2 - 2 * beta**2) / (2 * (alpha**2 - beta**2))
        self.toeplitz = _parallel(fs, fr)
        if lmax is None:
            # generous upper bound: S arrival at the farthest pair times (1 + 2/sqrt(tol))
            d = np.linalg.norm(fr.centers[:, None, :] - fs.centers[None, :, :], axis=2)
            lmax = int(np.ceil((d.max() + fs.dx) / beta / dt * (1 + 2 / np.sqrt(tol)))) + 2
        if self.toeplitz:
            # receivers at offsets m = i - j relative to source element 0
            m = np.arange(-(fs.n - 1), fr.n)
            d0 = fr.centers[0] - fs.centers[0]
            base = np.dot(d0, fs.s)                       # along-tangent offset of fr element 0 from fs element 0
            perp = d0 - base * fs.s
            Xr = fs.centers[0][None, :] + (base + m * fs.dx)[:, None] * fs.s[None, :] + perp[None, :]
            if mode == "inplane":
                Ws, Wn, Ss, Sn = ip.W_table(fs.centers[0], fs.s, fs.a, Xr, fr.s, lmax, dt, alpha, beta, rho, et=et)
            else:
                from antiplane import W_table_ap
                Ws, Wn, Ss, Sn = W_table_ap(fs.centers[0], fs.s, fs.a, Xr, fr.s, lmax, dt, alpha, beta, rho, et=et)
            self.m = m
        else:
            raise NotImplementedError("dense non-parallel pairs: use DensePairTable")
        # window: |W| compared with the static interaction envelope at the
        # pair distance, mu dx^2 / (2 pi (1-nu) r^2) * dt (the radiation
        # pattern has nodal directions where the static value itself is small
        # but the dynamic correction is not)
        rr = np.maximum(np.linalg.norm(Xr - fs.centers[0][None, :], axis=1), fs.a)
        env = mu * fs.dx**2 / (2 * np.pi * ((1 - nu) if mode == "inplane" else 1.0) * rr**2) * dt
        big = (np.abs(Ws) > tol * env[:, None]) | (np.abs(Wn) > tol * env[:, None])
        last = np.where(big.any(axis=0))[0]
        self.L = int(last.max()) + 1 if last.size else 1
        if self.L >= lmax:
            if not cap_ok:
                raise RuntimeError("window reached lmax; increase lmax")
            self.L = lmax - 1          # explicit cap (e.g. the run is shorter than the natural window)
        Ws = Ws[:, :self.L]; Wn = Wn[:, :self.L]
        self.Ks_row, self.Kn_row = Ss, Sn          # static rows over offsets m
        # FFT storage along offsets
        self.M = M if M is not None else int(2 ** np.ceil(np.log2(fr.n + fs.n - 1)))
        self.has_normal = mode == "inplane"
        self.Ws_hat = np.empty((self.L, self.M), complex)
        self.Wn_hat = np.empty((self.L, self.M), complex) if self.has_normal else None
        for arr, hat in ((Ws, self.Ws_hat), (Wn, self.Wn_hat)):
            if hat is None:
                continue
            c = np.zeros((self.L, self.M))
            # offsets m >= 0 -> index m ; m < 0 -> index M + m
            pos = m >= 0
            c[:, m[pos]] = arr[pos, :].T
            c[:, self.M + m[~pos]] = arr[~pos, :].T
            hat[:] = np.fft.fft(c, axis=1)
        del Ws, Wn
        self.nbytes = self.Ws_hat.nbytes + (self.Wn_hat.nbytes if self.has_normal else 0)

    def static_blocks(self):
        """Dense static blocks Ks[i,j], Kn[i,j] (receiver x source) for this pair."""
        idx = self.fr_index_minus_src()
        return self.Ks_row[idx], self.Kn_row[idx]

    def fr_index_minus_src(self):
        i = np.arange(self.fr.n)[:, None]; j = np.arange(self.fs.n)[None, :]
        return (i - j) + (self.fs.n - 1)

    def W_dense(self, l):
        """Dense W[l] blocks (for tests)."""
        idx = self.fr_index_minus_src()
        c = np.fft.ifft(self.Ws_hat[l - 1]).real
        cn = np.fft.ifft(self.Wn_hat[l - 1]).real if self.has_normal else np.zeros(self.M)
        m = self.m
        Ws = np.where(m >= 0, c[np.clip(m, 0, None)], c[(self.M + m) % self.M])
        Wn = np.where(m >= 0, cn[np.clip(m, 0, None)], cn[(self.M + m) % self.M])
        return Ws[idx], Wn[idx]


class History:
    """Ring buffer of FFT'd slip-rate slices for one source fault."""
    def __init__(self, n, M, L):
        self.n, self.M, self.L = n, M, L
        self.buf = np.zeros((L, M), complex)
        self.k = 0            # number of slices pushed
    def push(self, V):
        vp = np.zeros(self.M); vp[:self.n] = V
        self.buf[self.k % self.L] = np.fft.fft(vp)
        self.k += 1
    def lagged(self, l, nlag):
        """indices of the slices at lags 1..nlag relative to the next slice to be pushed:
        slice at lag l is the one pushed (k - l)."""
        idx = (self.k - np.arange(1, nlag + 1)) % self.L
        valid = (self.k - np.arange(1, nlag + 1)) >= 0
        return idx, valid


class System:
    """Fault system with static operators, pair tables and histories."""
    def __init__(self, faults, dt, alpha, beta, rho, tol=1e-2, lmax=None, mode="inplane", cap_ok=False, et=1.0):
        self.faults = faults
        self.mode = mode
        self.et = et
        self.dt = dt
        self.alpha, self.beta, self.rho = alpha, beta, rho
        self.mu = rho * beta**2
        self.nu = (alpha**2 - 2 * beta**2) / (2 * (alpha**2 - beta**2))
        self.eta = self.mu / (2 * beta)
        self.N = sum(f.n for f in faults)
        self.off = np.cumsum([0] + [f.n for f in faults])
        nmax = max(f.n for f in faults)
        self.M = int(2 ** np.ceil(np.log2(2 * nmax - 1)))
        self.pairs = {}
        for a_, fs in enumerate(faults):
            for b_, fr in enumerate(faults):
                self.pairs[(a_, b_)] = PairTable(fs, fr, dt, alpha, beta, rho, tol=tol, lmax=lmax, M=self.M, mode=mode, cap_ok=cap_ok, et=et)
        self.L = max(p.L for p in self.pairs.values())
        self.hist = [History(f.n, self.M, self.L + 1) for f in faults]
        # dense static operators (receiver x source), physical convention
        self.Ks = np.zeros((self.N, self.N)); self.Kn = np.zeros((self.N, self.N))
        for (a_, b_), p in self.pairs.items():
            ks, kn = p.static_blocks()
            self.Ks[self.off[b_]:self.off[b_ + 1], self.off[a_]:self.off[a_ + 1]] = ks
            self.Kn[self.off[b_]:self.off[b_ + 1], self.off[a_]:self.off[a_ + 1]] = kn
        # first-lag kernel K[1] = W[1] + Ks dt, dense (banded in practice): used by the stepper
        self.K1s = np.zeros((self.N, self.N)); self.K1n = np.zeros((self.N, self.N))
        for (a_, b_), p in self.pairs.items():
            w1s, w1n = p.W_dense(1)
            ks, kn = p.static_blocks()
            self.K1s[self.off[b_]:self.off[b_ + 1], self.off[a_]:self.off[a_ + 1]] = w1s + ks * dt
            self.K1n[self.off[b_]:self.off[b_ + 1], self.off[a_]:self.off[a_ + 1]] = w1n + kn * dt
        self.nbytes = sum(p.nbytes for p in self.pairs.values()) + sum(h.buf.nbytes for h in self.hist)
        self._Kk = {1: (self.K1s, self.K1n)}

    def cumulative_kernel(self, k):
        """K_k = sum_{l=1}^{k} K[l] = sum_{l<=k} W[l] + k Ks dt: stress at k dt from
        unit slip rate during [0, k dt). Implicit operator for a step of k dt_min.
        Cached; dense (banded in practice)."""
        if k in self._Kk:
            return self._Kk[k]
        Kk_s = np.zeros((self.N, self.N)); Kk_n = np.zeros((self.N, self.N))
        for (a_, b_), p in self.pairs.items():
            ks, kn = p.static_blocks()
            acc_s = ks * self.dt * k; acc_n = kn * self.dt * k
            for l in range(1, min(k, p.L) + 1):
                ws, wn = p.W_dense(l)
                acc_s = acc_s + ws; acc_n = acc_n + wn
            Kk_s[self.off[b_]:self.off[b_ + 1], self.off[a_]:self.off[a_ + 1]] = acc_s
            Kk_n[self.off[b_]:self.off[b_ + 1], self.off[a_]:self.off[a_ + 1]] = acc_n
        self._Kk[k] = (Kk_s, Kk_n)
        return self._Kk[k]

    def split(self, v):
        return [v[self.off[a_]:self.off[a_ + 1]] for a_ in range(len(self.faults))]

    def push(self, V):
        """Append the slip-rate slice of the step just completed."""
        for a_, Va in enumerate(self.split(V)):
            self.hist[a_].push(Va)

    def history_stress(self, lag0=1):
        """sum_{l=lag0}^{L} W[l] V^{(next - l)} for all receivers; returns (shear, normal).
        With lag0 = 2 the current-step term is excluded (the stepper adds it
        through K[1] after solving for V^n); with lag0 = 1 after pushing V^n
        the full history term is returned."""
        ts = np.zeros(self.N); tn = np.zeros(self.N)
        for (a_, b_), p in self.pairs.items():
            h = self.hist[a_]
            # lags l = lag0..min(L, k): slice pushed at step k - l sits at position (k - l) mod Lbuf;
            # walk the ring buffer in at most two contiguous (reversed) segments, no copies
            Lb = h.L; k = h.k
            lmax_ = min(p.L, k)
            if lmax_ < lag0:
                continue
            acc_s = np.zeros(p.M, complex); acc_n = np.zeros(p.M, complex) if p.has_normal else None
            l = lag0
            while l <= lmax_:
                pos = (k - l) % Lb                  # position of lag l
                nseg = min(lmax_ - l + 1, pos + 1)  # lags l..l+nseg-1 occupy positions pos..pos-nseg+1
                Vh = h.buf[pos - nseg + 1:pos + 1][::-1]
                Wsl = p.Ws_hat[l - 1:l - 1 + nseg]
                acc_s += np.einsum("lk,lk->k", Wsl, Vh)
                if p.has_normal:
                    acc_n += np.einsum("lk,lk->k", p.Wn_hat[l - 1:l - 1 + nseg], Vh)
                l += nseg
            ts[self.off[b_]:self.off[b_ + 1]] += np.fft.ifft(acc_s).real[:p.fr.n]
            if p.has_normal:
                tn[self.off[b_]:self.off[b_ + 1]] += np.fft.ifft(acc_n).real[:p.fr.n]
        return ts, tn

    def reset_history(self):
        for h in self.hist:
            h.buf[:] = 0; h.k = 0
