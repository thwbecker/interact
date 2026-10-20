"""3-D full-space point moment-tensor source: displacement (Aki and Richards
2002, eq. 4.29) for a general moment function M_pq(t) = Mhat_pq m(t), and its
verification against the Navier equation (away from the source, after both
arrivals) and against the static double-couple limit. Building block of the
quadrature route for the 3-D rectangle kernels (second, independent route
to the closed-form kernels)."""
import sympy as sp

x1, x2, x3, t = sp.symbols("x1 x2 x3 t", real=True)
al, be, rho = sp.symbols("alpha beta rho", positive=True)
X = [x1, x2, x3]
r = sp.sqrt(x1**2 + x2**2 + x3**2)
g = [xi / r for xi in X]
d = lambda i, j: 1 if i == j else 0


def displacement(Mhat, m, mint):
    """u_n for moment tensor Mhat_pq m(t); m: function of time (sympy Lambda),
    mint(t): antiderivative pair giving int_{r/al}^{r/be} tau m(t - tau) dtau"""
    u = []
    for n in range(3):
        s = 0
        for p in range(3):
            for q in range(3):
                if Mhat[p][q] == 0:
                    continue
                gn, gp, gq = g[n], g[p], g[q]
                NF = (15*gn*gp*gq - 3*gn*d(p, q) - 3*gp*d(n, q) - 3*gq*d(n, p)) / (4*sp.pi*rho*r**4) * mint
                IPa = (6*gn*gp*gq - gn*d(p, q) - gp*d(n, q) - gq*d(n, p)) / (4*sp.pi*rho*al**2*r**2) * m(t - r/al)
                ISb = -(6*gn*gp*gq - gn*d(p, q) - gp*d(n, q) - 2*gq*d(n, p)) / (4*sp.pi*rho*be**2*r**2) * m(t - r/be)
                FPa = gn*gp*gq / (4*sp.pi*rho*al**3*r) * sp.diff(m(t), t).subs(t, t - r/al)
                FSb = -(gn*gp - d(n, p))*gq / (4*sp.pi*rho*be**3*r) * sp.diff(m(t), t).subs(t, t - r/be)
                s += Mhat[p][q] * (NF + IPa + ISb + FPa + FSb)
        u.append(s)
    return u


if __name__ == "__main__":
    import random
    tau = sp.symbols("tau", positive=True)
    m = sp.Lambda(t, t**6)                       # smooth moment function (valid after both arrivals)
    mint = sp.integrate(tau * (t - tau)**6, (tau, r/al, r/be))
    Mhat = [[0, 1, 0], [1, 0, 0], [0, 0, 0]]     # double couple, slip x1 on plane normal x2
    u = displacement(Mhat, m, mint)
    lam = rho*(al**2 - 2*be**2); mu = rho*be**2
    div = sum(sp.diff(u[i], X[i]) for i in range(3))
    vals = {al: 6.0, be: 3.4641, rho: 2.67, x1: 0.7, x2: -0.4, x3: 0.9, t: 1.3}
    worst = 0
    for n in range(3):
        lap = sum(sp.diff(u[n], X[k], 2) for k in range(3))
        res = rho*sp.diff(u[n], t, 2) - ((lam + mu)*sp.diff(div, X[n]) + mu*lap)
        sc = abs(float((rho*sp.diff(u[n], t, 2)).subs(vals)))
        worst = max(worst, abs(float(res.subs(vals))) / sc)
    print(f"Navier residual / |rho u_tt| = {worst:.2e}")
    # static limit: m(t) = 1 (step) at late time -> static double couple, compare with the
    # known static point double couple u_n = (1/(8 pi mu r^2)) [ ... ] via Okada-type formula:
    # static Green's function displacement gradient (Kelvin) contracted with Mhat
    ms = sp.Lambda(t, sp.Integer(1)); mints = sp.integrate(tau, (tau, r/al, r/be))
    us = displacement(Mhat, ms, mints)
    nu = lam/(2*(lam + mu))
    # static: u_n = M_pq * dG_np/d(xi_q) = -M_pq dG_np/dx_q (derivative with respect to the
    # source position, Aki and Richards convention), Kelvin G_np = [(3-4nu) d_np + g_n g_p] / (16 pi mu (1-nu) r)
    G = [[((3 - 4*nu)*d(n, p) + g[n]*g[p]) / (16*sp.pi*mu*(1 - nu)*r) for p in range(3)] for n in range(3)]
    uk = [-sum(Mhat[p][q]*sp.diff(G[n][p], X[q]) for p in range(3) for q in range(3)) for n in range(3)]
    err = max(abs(float((us[n] - uk[n]).subs(vals))) / max(abs(float(uk[k].subs(vals))) for k in range(3)) for n in range(3))
    print(f"static limit against Kelvin-derived double couple: rel err {err:.2e}")
