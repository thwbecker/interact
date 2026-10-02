"""
Explicit (accelerated pseudo-transient) solver for a SEAS BP1-type problem.

2D antiplane strain, vertical strike-slip fault at x = 0, free surface at z = 0.
Symmetry about the fault: only x >= 0 is modeled and fault slip = 2 * u(0, z).
Quasi-dynamic: radiation damping eta * V on the fault.
Rate-and-state friction (regularized), aging law.

The quasi-static elastic problem  d/dx(mu du/dx) + d/dz(mu du/dz) = 0
is solved at each physical time step by an accelerated pseudo-transient
(damped wave) iteration. The fault enters as a nonlinear Robin condition:
at every PT iteration the slip rate V(z) is found from the current traction
by a bracketed root solve, and the fault displacement is set to
u0_old + V dt / 2.

Units: SI. Grid: nodal displacement u[i, j], i in x (nx+1 nodes), j in z (nz+1 nodes).
"""

import argparse
import time

import numpy as np


# ---------------------------------------------------------------------------
# BP1 parameters (Erickson et al., 2020)
# ---------------------------------------------------------------------------
class BP1:
    rho = 2670.0            # kg/m^3
    cs = 3464.0             # m/s
    mu = rho * cs**2        # 32.04 GPa
    sigma_n = 50.0e6        # Pa
    a0 = 0.010
    amax = 0.025
    b = 0.015
    Dc = 0.008              # m
    V0 = 1.0e-6             # m/s
    f0 = 0.6
    Vp = 1.0e-9             # m/s (plate rate)
    Vinit = 1.0e-9          # m/s
    H = 15.0e3              # m, depth of VW region
    h = 3.0e3               # m, transition width
    Wf = 40.0e3             # m, depth of RSF fault; creeping at Vp below
    eta = mu / (2.0 * cs)   # radiation damping, Pa s / m

    @classmethod
    def a_profile(cls, z):
        a = np.full_like(z, cls.amax)
        a[z < cls.H] = cls.a0
        mask = (z >= cls.H) & (z < cls.H + cls.h)
        a[mask] = cls.a0 + (cls.amax - cls.a0) * (z[mask] - cls.H) / cls.h
        return a


def friction(V, theta, a, p):
    """Regularized rate-and-state friction stress (Pa)."""
    arg = (p.f0 + p.b * np.log(p.V0 * theta / p.Dc)) / a
    return p.sigma_n * a * np.arcsinh(V / (2.0 * p.V0) * np.exp(arg))


def theta_aging_exact(theta_old, V, dt, p):
    """Aging law integrated exactly for constant V over dt."""
    x = V * dt / p.Dc
    # for small x use series to avoid cancellation
    ex = np.where(x < 1e-8, 1.0 - x, np.exp(-x))
    return p.Dc / V + (theta_old - p.Dc / V) * ex


# ---------------------------------------------------------------------------
# Elastic residual with PT iteration
# ---------------------------------------------------------------------------
def residual(u, mu, dx, dz, R):
    """
    R[i, j] = mu * (u_xx + u_zz) at interior x nodes 1..nx-1 and z nodes 0..nz-1.
    Free surface at j = 0 by mirror (du/dz = 0). Dirichlet at i = 0, i = nx, j = nz.
    Fills R in place; R has shape (nx+1, nz+1), only interior entries are meaningful.
    """
    inv_dx2 = 1.0 / dx**2
    inv_dz2 = 1.0 / dz**2
    # interior in z: j = 1..nz-1
    R[1:-1, 1:-1] = mu * (
        (u[2:, 1:-1] - 2.0 * u[1:-1, 1:-1] + u[:-2, 1:-1]) * inv_dx2
        + (u[1:-1, 2:] - 2.0 * u[1:-1, 1:-1] + u[1:-1, :-2]) * inv_dz2
    )
    # free surface row j = 0: mirror u[:, -1] = u[:, 1]
    R[1:-1, 0] = mu * (
        (u[2:, 0] - 2.0 * u[1:-1, 0] + u[:-2, 0]) * inv_dx2
        + 2.0 * (u[1:-1, 1] - u[1:-1, 0]) * inv_dz2
    )
    return R


def fault_traction_elastic(u, mu, dx):
    """Second-order one-sided mu * du/dx at x = 0 for every z node."""
    return mu * (-3.0 * u[0, :] + 4.0 * u[1, :] - u[2, :]) / (2.0 * dx)


def friction_dV(V, theta, a, p):
    """dF/dV for the regularized law at fixed state."""
    arg = (p.f0 + p.b * np.log(p.V0 * theta / p.Dc)) / a
    q = np.exp(arg) / (2.0 * p.V0)
    return p.sigma_n * a * q / np.sqrt(1.0 + (V * q)**2)


def solve_fault_V(u1, u2, u0_old, theta_old, a_f, tau0, dt, mu, dx, p, V_guess,
                  Vlo=1e-25, Vhi=1e3, newton_max=100, rtol=1e-12, g_tol=1e-4):
    """
    Solve, per RSF fault node,
      g(V) = tau0 + mu*(-3*(u0_old + V dt/2) + 4 u1 - u2)/(2 dx)
             - friction(V, theta_old) - eta V = 0
    State is held at theta_old (explicit in state, implicit in V), which keeps
    g strictly decreasing in V. Safeguarded Newton in s = log V starting from V_guess,
    with a bisection fallback on [Vlo, Vhi] if a Newton step leaves the bracket.
    """
    c0 = tau0 + mu * (4.0 * u1 - u2) / (2.0 * dx)
    cu = -3.0 * mu / (2.0 * dx)

    def g_and_dg(V):
        gv = c0 + cu * (u0_old + 0.5 * V * dt) - friction(V, theta_old, a_f, p) - p.eta * V
        dg = 0.5 * cu * dt - friction_dV(V, theta_old, a_f, p) - p.eta
        return gv, dg

    lo = np.full_like(u1, np.log(Vlo))
    hi = np.full_like(u1, np.log(Vhi))
    s = np.log(np.clip(V_guess, Vlo, Vhi))
    for _ in range(newton_max):
        V = np.exp(s)
        gv, dg = g_and_dg(V)
        pos = gv > 0.0
        lo = np.where(pos, s, lo)
        hi = np.where(pos, hi, s)
        ds = -gv / (dg * V)              # Newton step in log V
        # cap the step: for V far below the root, friction is linear in V and the
        # log-space Newton step overshoots by many decades
        ds = np.clip(ds, -5.0, 5.0)
        s_new = s + ds
        # keep the iterate inside the bracket: clamp to the bracket end (a point already
        # evaluated, from which Newton converges), and bisect only if clamping makes no progress
        s_clamped = np.clip(s_new, lo, hi)
        stalled = (s_clamped != s_new) & (s_clamped == s)
        s_new = np.where(stalled, 0.5 * (lo + hi), s_clamped)
        # converged if the log-V step is tiny or the traction imbalance is below g_tol (Pa);
        # the second test matters for locked nodes, where g is insensitive to V and
        # floating point noise in g prevents the first test from being met
        if np.all((np.abs(s_new - s) < rtol) | (np.abs(gv) < g_tol) | (hi - lo < rtol)):
            s = s_new
            break
        s = s_new
    solve_fault_V.last_it = _ + 1
    return np.exp(s)



# ---------------------------------------------------------------------------
# Optional numba-fused PT step (same algorithm as the NumPy path)
# ---------------------------------------------------------------------------
try:
    import numba as _nb

    @_nb.njit(cache=True)
    def _friction_nb(V, theta, a, sigma_n, f0, b, V0, Dc):
        arg = (f0 + b * np.log(V0 * theta / Dc)) / a
        return sigma_n * a * np.arcsinh(V / (2.0 * V0) * np.exp(arg))

    @_nb.njit(cache=True)
    def _friction_dV_nb(V, theta, a, sigma_n, f0, b, V0, Dc):
        arg = (f0 + b * np.log(V0 * theta / Dc)) / a
        q = np.exp(arg) / (2.0 * V0)
        return sigma_n * a * q / np.sqrt(1.0 + (V * q) ** 2)

    @_nb.njit(cache=True,parallel=True)
    #@_nb.njit(cache=True)
    def pt_step_nb(u, dudtau, R, V, u0_old, theta_old, a_f, tau0, rsf_idx, dt,
                   mu, dx, dz, sigma_n, f0, b, V0, Dc, eta, dtau, damp, tol_pa, max_it):
        nx1, nz1 = u.shape
        nx, nz = nx1 - 1, nz1 - 1
        inv_dx2, inv_dz2, inv_mu = 1.0 / dx**2, 1.0 / dz**2, 1.0 / mu
        cu = -3.0 * mu / (2.0 * dx)
        lo0, hi0 = np.log(1e-25), np.log(1e3)
        it = 0
        err = 1e300
        while True:
            # fault Robin condition
            for k in range(rsf_idx.shape[0]):
                j = rsf_idx[k]
                c0 = tau0[k] + mu * (4.0 * u[1, j] - u[2, j]) / (2.0 * dx)
                lo, hi = lo0, hi0
                s = np.log(min(max(V[k], 1e-25), 1e3))
                for _ in range(100):
                    Vs = np.exp(s)
                    gv = c0 + cu * (u0_old[k] + 0.5 * Vs * dt) \
                        - _friction_nb(Vs, theta_old[k], a_f[k], sigma_n, f0, b, V0, Dc) - eta * Vs
                    dg = 0.5 * cu * dt - _friction_dV_nb(Vs, theta_old[k], a_f[k], sigma_n, f0, b, V0, Dc) - eta
                    if gv > 0.0:
                        lo = s
                    else:
                        hi = s
                    ds = -gv / (dg * Vs)
                    ds = min(max(ds, -5.0), 5.0)
                    s_raw = s + ds
                    s_cl = min(max(s_raw, lo), hi)
                    if s_cl != s_raw and s_cl == s:
                        s_new = 0.5 * (lo + hi)
                    else:
                        s_new = s_cl
                    done = abs(s_new - s) < 1e-12 or abs(gv) < 1e-4 or (hi - lo) < 1e-12
                    s = s_new
                    if done:
                        break
                V[k] = np.exp(s)
                u[0, j] = u0_old[k] + 0.5 * V[k] * dt
            # residual on interior nodes
            err = 0.0
            #for i in range(1, nx):
            for i in _nb.prange(1, nx):
                for j in range(0, nz):
                    uxx = (u[i + 1, j] - 2.0 * u[i, j] + u[i - 1, j]) * inv_dx2
                    if j == 0:
                        uzz = 2.0 * (u[i, 1] - u[i, 0]) * inv_dz2
                    else:
                        uzz = (u[i, j + 1] - 2.0 * u[i, j] + u[i, j - 1]) * inv_dz2
                    r = mu * (uxx + uzz)
                    R[i, j] = r
                    err = max(err, abs(r))
                    #ar = abs(r)
                    #if ar > err:
                    #    err = ar
            # PT update
            #for i in range(1, nx):
            for i in _nb.prange(1, nx):
                for j in range(0, nz):
                    dudtau[i, j] = damp * dudtau[i, j] + dtau * R[i, j] * inv_mu
                    u[i, j] += dtau * dudtau[i, j]
            err *= dx
            it += 1
            if err < tol_pa or it >= max_it:
                break
        return it, err


    @_nb.njit(cache=True, parallel=True)
    def pt_step_layer_nb(u, dudtau, R, V, tau_f, dp_old, theta_old, a_f, tau0, rsf_idx, dt,
                         mu, dx, dz, sigma_n, f0, b, V0, Dc, eta, dtau, damp, tol_pa, max_it):
        """
        Quasi-viscous fault layer (Herrendoerfer et al. 2018 style, D = dx).
        Nodes sit at x = (i + 1/2) dx; the layer occupies [-dx/2, dx/2] across the
        symmetry plane, so the displacement jump across it is 2 u[0, j]. The layer
        carries plastic slip dp; its elastic shear stress perturbation is
        mu (2 u0 - dp) / dx and its total traction tau_f = tau0 + that. Node 0 is an
        unknown of the PT iteration for RSF nodes (Dirichlet for creeping nodes).
        Slip rate V is solved per node from tau_f(V) = friction(V, theta_old) + eta V
        with dp = dp_old + V dt (implicit in V, explicit in state).
        """
        nx1, nz1 = u.shape
        nx, nz = nx1 - 1, nz1 - 1
        inv_dx2, inv_dz2, inv_mu = 1.0 / dx**2, 1.0 / dz**2, 1.0 / mu
        lo0, hi0 = np.log(1e-25), np.log(1e3)
        it = 0
        err = 1e300
        while True:
            # fault layer: solve V per RSF node from current u0
            for k in range(rsf_idx.shape[0]):
                j = rsf_idx[k]
                c0 = tau0[k] + mu * (2.0 * u[0, j] - dp_old[k]) / dx
                cu = -mu * dt / dx
                lo, hi = lo0, hi0
                s = np.log(min(max(V[k], 1e-25), 1e3))
                for _ in range(100):
                    Vs = np.exp(s)
                    gv = c0 + cu * Vs \
                        - _friction_nb(Vs, theta_old[k], a_f[k], sigma_n, f0, b, V0, Dc) - eta * Vs
                    dg = cu - _friction_dV_nb(Vs, theta_old[k], a_f[k], sigma_n, f0, b, V0, Dc) - eta
                    if gv > 0.0:
                        lo = s
                    else:
                        hi = s
                    ds = -gv / (dg * Vs)
                    ds = min(max(ds, -5.0), 5.0)
                    s_raw = s + ds
                    s_cl = min(max(s_raw, lo), hi)
                    if s_cl != s_raw and s_cl == s:
                        s_new = 0.5 * (lo + hi)
                    else:
                        s_new = s_cl
                    done = abs(s_new - s) < 1e-12 or abs(gv) < 1e-4 or (hi - lo) < 1e-12
                    s = s_new
                    if done:
                        break
                V[k] = np.exp(s)
                tau_f[k] = c0 + cu * V[k]
            # residual, interior nodes i = 1..nx-1
            err = 0.0
            for i in _nb.prange(1, nx):
                for j in range(0, nz):
                    uxx = (u[i + 1, j] - 2.0 * u[i, j] + u[i - 1, j]) * inv_dx2
                    if j == 0:
                        uzz = 2.0 * (u[i, 1] - u[i, 0]) * inv_dz2
                    else:
                        uzz = (u[i, j + 1] - 2.0 * u[i, j] + u[i, j - 1]) * inv_dz2
                    r = mu * (uxx + uzz)
                    R[i, j] = r
                    err = max(err, abs(r))
            # residual, fault-side node i = 0 for RSF nodes
            for k in range(rsf_idx.shape[0]):
                j = rsf_idx[k]
                if j < nz:
                    dtx = (mu * (u[1, j] - u[0, j]) / dx - (tau_f[k] - tau0[k])) / dx
                    if j == 0:
                        uzz = 2.0 * (u[0, 1] - u[0, 0]) * inv_dz2
                    else:
                        uzz = (u[0, j + 1] - 2.0 * u[0, j] + u[0, j - 1]) * inv_dz2
                    r = dtx + mu * uzz
                    R[0, j] = r
                    err = max(err, abs(r))
            # PT update, interior
            for i in _nb.prange(1, nx):
                for j in range(0, nz):
                    dudtau[i, j] = damp * dudtau[i, j] + dtau * R[i, j] * inv_mu
                    u[i, j] += dtau * dudtau[i, j]
            # PT update, fault-side RSF nodes
            for k in range(rsf_idx.shape[0]):
                j = rsf_idx[k]
                if j < nz:
                    dudtau[0, j] = damp * dudtau[0, j] + dtau * R[0, j] * inv_mu
                    u[0, j] += dtau * dudtau[0, j]
            err *= dx
            it += 1
            if err < tol_pa or it >= max_it:
                break
        return it, err

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    HAVE_NUMBA = False

class PTElasticity:
    def __init__(self, nx, nz, Lx, Lz, mu, cfl=0.9, nu=4.0):
        self.nx, self.nz = nx, nz
        self.dx, self.dz = Lx / nx, Lz / nz
        self.mu = mu
        self.x = np.linspace(0.0, Lx, nx + 1)
        self.z = np.linspace(0.0, Lz, nz + 1)
        # PT parameters: wave speed 1 in pseudo-time, rho_pt = mu
        self.dtau = cfl * min(self.dx, self.dz) / np.sqrt(2.0)
        self.damp = 1.0 - nu / max(nx, nz)
        self.R = np.zeros((nx + 1, nz + 1))
        self.dudtau = np.zeros((nx + 1, nz + 1))

    def reset_momentum(self):
        self.dudtau[:] = 0.0

    def step(self, u, source=None):
        """One accelerated PT iteration on interior nodes. Returns max |R| * dx (Pa)."""
        R = residual(u, self.mu, self.dx, self.dz, self.R)
        if source is not None:
            R[1:-1, :-1] += source[1:-1, :-1]
        inner = (slice(1, -1), slice(0, -1))
        self.dudtau[inner] = self.damp * self.dudtau[inner] + self.dtau * R[inner] / self.mu
        u[inner] += self.dtau * self.dudtau[inner]
        return np.abs(R[inner]).max() * self.dx




def solve_layer_V(u0, dp_old, theta_old, a_f, tau0, dt, mu, dx, p, V_guess,
                  Vlo=1e-25, Vhi=1e3, newton_max=100, rtol=1e-12, g_tol=1e-4):
    """Layer-fault counterpart of solve_fault_V: g(V) = tau0 + mu(2u0 - dp_old - V dt)/dx - F - eta V."""
    c0 = tau0 + mu * (2.0 * u0 - dp_old) / dx
    cu = -mu * dt / dx
    lo = np.full_like(u0, np.log(Vlo))
    hi = np.full_like(u0, np.log(Vhi))
    s = np.log(np.clip(V_guess, Vlo, Vhi))
    for _ in range(newton_max):
        V = np.exp(s)
        gv = c0 + cu * V - friction(V, theta_old, a_f, p) - p.eta * V
        dg = cu - friction_dV(V, theta_old, a_f, p) - p.eta
        pos = gv > 0.0
        lo = np.where(pos, s, lo)
        hi = np.where(pos, hi, s)
        ds = np.clip(-gv / (dg * V), -5.0, 5.0)
        s_raw = s + ds
        s_cl = np.clip(s_raw, lo, hi)
        stalled = (s_cl != s_raw) & (s_cl == s)
        s_new = np.where(stalled, 0.5 * (lo + hi), s_cl)
        if np.all((np.abs(s_new - s) < rtol) | (np.abs(gv) < g_tol) | (hi - lo < rtol)):
            s = s_new
            break
        s = s_new
    V = np.exp(s)
    return V, c0 + cu * V


def pt_step_layer_np(u, pt, V, dp_old, theta_old, a_f, tau0, rsf, dt, p, tol_pa, max_it):
    """NumPy implementation of one time step in layer mode; mirrors pt_step_layer_nb."""
    dx, dz = pt.dx, pt.dz
    inv_dx2, inv_dz2 = 1.0 / dx**2, 1.0 / dz**2
    nx, nz = pt.nx, pt.nz
    R, dudtau = pt.R, pt.dudtau
    rsf_i = rsf.copy()
    rsf_i[nz] = False                      # bottom row is Dirichlet
    it = 0
    while True:
        V, tau_f = solve_layer_V(u[0, rsf], dp_old, theta_old, a_f, tau0, dt, p.mu, dx, p, V)
        residual(u, p.mu, dx, dz, R)
        # fault-side node column i = 0 for RSF nodes
        tau_f_full = np.zeros(nz + 1)
        tau_f_full[rsf] = tau_f - tau0
        dtx = (p.mu * (u[1, :] - u[0, :]) / dx - tau_f_full) / dx
        uzz = np.empty(nz + 1)
        uzz[0] = 2.0 * (u[0, 1] - u[0, 0]) * inv_dz2
        uzz[1:nz] = (u[0, 2:] - 2.0 * u[0, 1:-1] + u[0, :-2]) * inv_dz2
        uzz[nz] = 0.0
        R[0, :] = 0.0
        R[0, rsf_i] = (dtx + p.mu * uzz)[rsf_i]
        inner = (slice(1, -1), slice(0, -1))
        dudtau[inner] = pt.damp * dudtau[inner] + pt.dtau * R[inner] / p.mu
        u[inner] += pt.dtau * dudtau[inner]
        dudtau[0, rsf_i] = pt.damp * dudtau[0, rsf_i] + pt.dtau * R[0, rsf_i] / p.mu
        u[0, rsf_i] += pt.dtau * dudtau[0, rsf_i]
        err = max(np.abs(R[inner]).max(), np.abs(R[0, rsf_i]).max()) * dx
        it += 1
        if err < tol_pa or it >= max_it:
            break
    return it, err, V, tau_f

# ---------------------------------------------------------------------------
# Manufactured-solution test of the elastic solver
# ---------------------------------------------------------------------------
def test_manufactured(nx, nz, L=1.0, mu=1.0, tol=1e-10, nu=2.0, verbose=False):
    """
    u = sin(k x) cos(k z) satisfies du/dz = 0 at z = 0.
    -mu (u_xx + u_zz) = f  with f = 2 mu k^2 u. Dirichlet from exact u elsewhere.
    Returns (L2 error, iterations).
    """
    k = np.pi
    pt = PTElasticity(nx, nz, L, L, mu, nu=nu)
    X, Z = np.meshgrid(pt.x, pt.z, indexing="ij")
    uex = np.sin(k * X) * np.cos(k * Z)
    f = 2.0 * mu * k**2 * uex
    u = np.zeros_like(uex)
    u[0, :] = uex[0, :]
    u[-1, :] = uex[-1, :]
    u[:, -1] = uex[:, -1]
    it = 0
    while True:
        err = pt.step(u, source=f)
        it += 1
        if err < tol or it > 100 * max(nx, nz):
            break
    l2 = np.sqrt(np.mean((u - uex)**2))
    if verbose:
        print(f"manufactured nx={nx} it={it} resid={err:.2e} L2={l2:.3e}")
    return l2, it


# ---------------------------------------------------------------------------
# BP1 driver
# ---------------------------------------------------------------------------
def run_bp1(nx, nz, Lx, Lz, t_end_yr, tol_pa=1.0, xi=0.1, dt_max=None,
            nu=4.0, max_it=None, out_every=1, log_path=None, wall_limit=None,
            verbose=True, checkpoint=None, stop_Vmax=None, use_numba=False,
            fault="robin", log_dense_Vmax=1e-4):
    """
    fault = "robin": fault as a displacement discontinuity at x = 0 (Robin condition,
                     second-order one-sided traction); the original scheme.
    fault = "layer": quasi-viscous fault layer of width dx across x = 0 carrying plastic
                     slip (Herrendoerfer et al. 2018 style, D = dx); nodes at (i + 1/2) dx.
    log_dense_Vmax: log every step while the maximum slip rate exceeds this value, and
                    every out_every steps otherwise; None disables dense logging.
    """
    if fault not in ("robin", "layer"):
        raise ValueError("fault must be 'robin' or 'layer'")
    layer = fault == "layer"
    p = BP1
    # the layer node column has a larger diagonal stencil weight (5 vs 4 for square
    # cells), so its PT time step is reduced accordingly
    pt = PTElasticity(nx, nz, Lx, Lz, p.mu, nu=nu, cfl=0.8 if layer else 0.9)
    dx, dz = pt.dx, pt.dz
    z = pt.z
    yr = 365.25 * 86400.0
    if dt_max is None:
        dt_max = 0.1 * yr
    if max_it is None:
        max_it = 50 * max(nx, nz)

    rsf = z <= p.Wf                 # RSF nodes on fault
    rsf_idx = np.where(rsf)[0].astype(np.int64)
    if use_numba and not HAVE_NUMBA:
        raise RuntimeError("numba requested but not installed")
    a_f = p.a_profile(z[rsf])
    # steady state at Vinit; with this theta0 the tau0 below equals the BP1 prestress formula
    theta = np.full(rsf.sum(), p.Dc / p.Vinit)
    tau0 = friction(p.Vinit, theta, a_f, p) + p.eta * p.Vinit

    u = np.zeros((nx + 1, nz + 1))
    u_prev = u.copy()
    V = np.full(rsf.sum(), p.Vinit)
    dp = np.zeros(rsf.sum())          # plastic slip in the fault layer (layer mode only)
    tau_f = tau0.copy()
    t = 0.0
    dt_prev = None
    step = 0
    log = []
    import os
    append_log = False
    if checkpoint is not None and os.path.exists(checkpoint):
        ck = np.load(checkpoint)
        u, u_prev, V, theta = ck["u"], ck["u_prev"], ck["V"], ck["theta"]
        if "dp" in ck.files:
            dp = ck["dp"]
        elif layer:
            raise RuntimeError("checkpoint lacks plastic slip; it was written in robin mode")
        t, step = float(ck["t"]), int(ck["step"])
        dt_prev = float(ck["dt_prev"]) if ck["dt_prev"] >= 0 else None
        append_log = True
        if verbose:
            print(f"restarted from step {step}, t = {t/yr:.4f} yr")
    t_wall0 = time.time()
    iz = {d: int(round(d * 1e3 / dz)) for d in (0.0, 7.5, 15.0, 25.0)}
    keys = ["step", "t_yr", "dt_s", "it", "resid_pa", "Vmax", "V_0km", "V_7p5km",
            "V_15km", "tau_7p5km_MPa", "slip_7p5km"]
    flog = None
    if log_path:
        flog = open(log_path, "a" if append_log else "w")
        if not append_log:
            flog.write(",".join(keys) + "\n")

    while t < t_end_yr * yr:
        Vmax = V.max()
        dt = min(dt_max, xi * p.Dc / Vmax)
        t_new = t + dt
        u_far = 0.5 * p.Vp * t_new

        # warm start by linear extrapolation
        if dt_prev is not None:
            u_guess = u + (u - u_prev) * (dt / dt_prev)
        else:
            u_guess = u.copy()
        u_prev = u.copy()
        u0_old = u[0, rsf].copy()
        dp_old = dp.copy()
        theta_old = theta.copy()

        unew = u_guess
        unew[-1, :] = u_far
        unew[:, -1] = u_far
        unew[0, ~rsf] = u_far
        pt.reset_momentum()

        it = 0
        if layer and use_numba:
            V = V.copy()
            tau_f = np.empty_like(tau0)
            it, err = pt_step_layer_nb(unew, pt.dudtau, pt.R, V, tau_f, dp_old, theta_old, a_f,
                                       tau0, rsf_idx, dt, p.mu, dx, dz, p.sigma_n, p.f0, p.b,
                                       p.V0, p.Dc, p.eta, pt.dtau, pt.damp, tol_pa, max_it)
        elif layer:
            it, err, V, tau_f = pt_step_layer_np(unew, pt, V, dp_old, theta_old, a_f, tau0, rsf,
                                                 dt, p, tol_pa, max_it)
        elif use_numba:
            V = V.copy()
            it, err = pt_step_nb(unew, pt.dudtau, pt.R, V, u0_old, theta_old, a_f, tau0,
                                 rsf_idx, dt, p.mu, dx, dz, p.sigma_n, p.f0, p.b, p.V0, p.Dc,
                                 p.eta, pt.dtau, pt.damp, tol_pa, max_it)
        else:
            while True:
                # fault Robin condition from current interior field
                V = solve_fault_V(unew[1, rsf], unew[2, rsf], u0_old, theta_old,
                                  a_f, tau0, dt, p.mu, dx, p, V)
                unew[0, rsf] = u0_old + 0.5 * V * dt
                err = pt.step(unew)
                it += 1
                if err < tol_pa or it >= max_it:
                    break

        theta = theta_aging_exact(theta_old, V, dt, p)
        if layer:
            dp = dp_old + V * dt
        u = unew
        t = t_new
        dt_prev = dt
        step += 1

        if layer:
            slip = dp                       # plastic slip; total jump 2 u0 differs by the
                                            # elastic layer strain, which is negligible
        else:
            tau_f = tau0 + fault_traction_elastic(u, p.mu, dx)[rsf]
            slip = 2.0 * u[0, rsf]
        dense = log_dense_Vmax is not None and V.max() > log_dense_Vmax
        if step % out_every == 0 or dense:
            rec = dict(step=step, t_yr=t / yr, dt_s=dt, it=it, resid_pa=err,
                       Vmax=V.max(), V_0km=V[iz[0.0]], V_7p5km=V[iz[7.5]],
                       V_15km=V[iz[15.0]], tau_7p5km_MPa=tau_f[iz[7.5]] / 1e6,
                       slip_7p5km=slip[iz[7.5]])
            log.append(rec)
            if flog:
                flog.write(",".join(f"{rec[k]:.10e}" if isinstance(rec[k], float) else str(rec[k])
                                    for k in keys) + "\n")
                flog.flush()
            if verbose and step % out_every == 0:
                print(f"step {step:6d} t={t/yr:9.4f} yr dt={dt:9.3e} s it={it:5d} "
                      f"res={err:8.2e} Pa Vmax={V.max():9.3e} tau7.5={tau_f[iz[7.5]]/1e6:8.4f} MPa",
                      flush=True)
        if stop_Vmax is not None and V.max() > stop_Vmax:
            if verbose:
                print(f"stop_Vmax reached at step {step}, t = {t/yr:.6f} yr")
            break
        if wall_limit is not None and time.time() - t_wall0 > wall_limit:
            if verbose:
                print("wall limit reached")
            break

    if flog:
        flog.close()
    if checkpoint is not None:
        np.savez(checkpoint, u=u, u_prev=u_prev, V=V, theta=theta, dp=dp, t=t, step=step,
                 dt_prev=dt_prev if dt_prev is not None else -1.0, fault=fault)
    return log


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true", help="run manufactured-solution test")
    ap.add_argument("--nx", type=int, default=160)
    ap.add_argument("--nz", type=int, default=160)
    ap.add_argument("--Lx", type=float, default=80e3)
    ap.add_argument("--Lz", type=float, default=80e3)
    ap.add_argument("--t_end_yr", type=float, default=200.0)
    ap.add_argument("--tol_pa", type=float, default=1.0)
    ap.add_argument("--xi", type=float, default=0.1)
    ap.add_argument("--nu", type=float, default=4.0)
    ap.add_argument("--out_every", type=int, default=1)
    ap.add_argument("--wall_limit", type=float, default=None)
    ap.add_argument("--log", type=str, default="bp1_log.csv")
    ap.add_argument("--numba", action="store_true", help="use the fused numba PT step")
    ap.add_argument("--threads", type=int, default=None, help="numba thread count")
    ap.add_argument("--log_dense_Vmax", type=float, default=1e-4,
                    help="log every step while Vmax exceeds this (m/s); negative disables")
    ap.add_argument("--fault", type=str, default="robin", choices=["robin", "layer"],
                    help="fault treatment: robin (displacement discontinuity) or layer (quasi-viscous, D = dx)")
    ap.add_argument("--stop_Vmax", type=float, default=None, help="stop once max slip rate exceeds this")
    ap.add_argument("--checkpoint", type=str, default=None,
                    help="npz file; if it exists the run restarts from it and overwrites it at exit")
    args = ap.parse_args()
    if args.threads is not None:
        import numba
        numba.set_num_threads(args.threads)
    if args.test:
        for n in (32, 64, 128):
            test_manufactured(n, n, verbose=True)
    else:
        run_bp1(args.nx, args.nz, args.Lx, args.Lz, args.t_end_yr, tol_pa=args.tol_pa,
                xi=args.xi, nu=args.nu, out_every=args.out_every,
                log_path=args.log, wall_limit=args.wall_limit, checkpoint=args.checkpoint,
                stop_Vmax=args.stop_Vmax, use_numba=args.numba,
                fault=args.fault,
                log_dense_Vmax=args.log_dense_Vmax if args.log_dense_Vmax > 0 else None)
