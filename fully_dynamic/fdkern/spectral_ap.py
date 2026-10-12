"""
Independent check: spectral boundary integral method for a planar antiplane
fault (Perrin, Rice and Zheng 1995; Lapusta et al. 2000 displacement form):

    tau(x,t) = tau0 - eta V(x,t) + f(x,t),
    f_k(t)   = -(mu |k| / 2) int_0^t [J1(|k| beta s) / s] D_k(t - s) ds,

with D_k the Fourier coefficient of slip on a periodic domain of length Lam
(replication; Lam is chosen large enough that image effects are small over
the simulated time). Time convolution by the trapezoid rule on the slip
history at step ends; the current-step term is kept implicit.

This shares nothing with the space-time code except the physics, so it is
used to check the dynamic factor of the self-similar (Kostrov) crack.
"""
import numpy as np
from scipy.special import j1


def self_similar_crack_spectral(vfrac, beta=3464.0, rho=2700.0, dx=100.0, a_final_elems=60,
                                Lam_over_a=32, cfl=0.25, dsigma=1e6):
    mu = rho * beta**2
    eta = mu / (2 * beta)
    v = vfrac * beta
    a_f = a_final_elems * dx
    N = int(2 ** np.ceil(np.log2(Lam_over_a * a_f / dx)))
    Lam = N * dx
    x = (np.arange(N) - N // 2 + 0.5) * dx
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    ak = np.abs(k)
    dt = cfl * dx / beta
    nsteps = int(a_f / v / dt)
    # kernel samples K_k(s_j) = J1(|k| beta s_j) / s_j, with K(0) = |k| beta / 2
    s = np.arange(nsteps + 2) * dt
    with np.errstate(invalid="ignore", divide="ignore"):
        Kk = j1(ak[:, None] * beta * s[None, :]) / s[None, :]
    Kk[:, 0] = ak * beta / 2
    pref = -(mu * ak / 2)
    c0 = pref * (dt / 2) * Kk[:, 0]                       # implicit coefficient on D^{n+1}
    # real matrix of the implicit operator  V -> IFFT(c0 dt * FFT(V))
    I = np.eye(N)
    Mx = np.fft.ifft(c0[:, None] * dt * np.fft.fft(I, axis=0), axis=0).real
    A = Mx - eta * I
    D = np.zeros(N); Dk_hist = []
    t = 0.0
    for n in range(nsteps):
        t1 = (n + 1) * dt
        # known part of f(t_{n+1}): trapezoid over j = 0..n (weights 1/2 at j = 0), using D^j
        fk = np.zeros(N, complex)
        if Dk_hist:
            H = np.array(Dk_hist)                       # (n+1, N) with D^0 .. D^n
            w = np.ones(len(Dk_hist)); w[0] = 0.5
            lags = np.arange(len(Dk_hist), 0, -1)       # t_{n+1} - t_j in steps: n+1 .. 1
            fk = pref * dt * np.einsum("jk,kj->k", (w[:, None] * H), Kk[:, lags])
        fk = fk + c0 * np.fft.fft(D)                      # c0 * D^n part of the implicit term
        f_known = np.fft.ifft(fk).real
        free = np.abs(x) < v * t1
        V = np.zeros(N)
        if free.any():
            idx = np.where(free)[0]
            V[idx] = np.linalg.solve(A[np.ix_(idx, idx)], -dsigma - f_known[idx])
        D = D + V * dt
        Dk_hist.append(np.fft.fft(D))
        t = t1
    mid = np.abs(x) < 0.25 * v * t
    C = np.mean(D[mid] / np.sqrt(v * v * t * t - x[mid] ** 2))
    return C / (2 / mu * dsigma), x, D, V, t, v


if __name__ == "__main__":
    for s_ in (0.5, 0.8):
        r = self_similar_crack_spectral(s_)
        print(s_, r[0])
