"""Dynamic check of the quadrant formulas (tada_quad.py); CHECK_OPENING=1 checks the
opening source (L_ij/3) with the moment tensor lambda delta + 2 mu n n. Response to slip D(t) = t^3/6 on a
square, computed (a) by double time integration of the Tada ramp response
and (b) by quadrature over the square of the Aki-Richards point source
(eq. 4.29) with moment function mu A t^3/6, in the regime where all source
points have both arrivals (post-S) and where all are between P and S."""
import numpy as np, sympy as sp, sys
sys.path.insert(0, ".")
import tada_quad as tq
cL, cT, rho = 6000., 3464., 2670.; mu = rho*cT**2
a = 50.0
x1, x2, x3, t, tau = sp.symbols("x1 x2 x3 t tau", real=True)
al, be, rh = sp.Float(cL), sp.Float(cT), sp.Float(rho)
X = [x1, x2, x3]; r = sp.sqrt(x1**2 + x2**2 + x3**2); g = [xi/r for xi in X]
d = lambda i, j: 1 if i == j else 0
mfun = lambda s: s**3/6
import os
OPEN = bool(os.environ.get("CHECK_OPENING"))
LAMR = (cL**2 - 2*cT**2)/cT**2                       # lambda/mu
def u_point(regime):
    if OPEN:   # opening across x3 = 0: M = (lambda delta + 2 mu n n) A, in units of mu A
        Mhat = [[sp.Float(LAMR), 0, 0], [0, sp.Float(LAMR), 0], [0, 0, sp.Float(LAMR) + 2]]
    else:
        Mhat = [[0, 0, 1], [0, 0, 0], [1, 0, 0]]          # slip in x1 on plane x3 = 0, unit mu*A
    upper = r/be if regime == "postS" else t
    mint = sp.integrate(tau*mfun(t - tau), (tau, r/al, upper))
    u = []
    for n in range(3):
        s = 0
        for p in range(3):
            for q in range(3):
                if Mhat[p][q] == 0: continue
                gn, gp, gq = g[n], g[p], g[q]
                w = Mhat[p][q]                                   # moment tensor component (units mu A)
                s0 = s; s = 0
                s += (15*gn*gp*gq - 3*gn*d(p, q) - 3*gp*d(n, q) - 3*gq*d(n, p))/(4*sp.pi*rh*r**4)*mint
                s += (6*gn*gp*gq - gn*d(p, q) - gp*d(n, q) - gq*d(n, p))/(4*sp.pi*rh*al**2*r**2)*mfun(t - r/al)
                s += gn*gp*gq/(4*sp.pi*rh*al**3*r)*sp.diff(mfun(t), t).subs(t, t - r/al)
                if regime == "postS":
                    s += -(6*gn*gp*gq - gn*d(p, q) - gp*d(n, q) - 2*gq*d(n, p))/(4*sp.pi*rh*be**2*r**2)*mfun(t - r/be)
                    s += -(gn*gp - d(n, p))*gq/(4*sp.pi*rh*be**3*r)*sp.diff(mfun(t), t).subs(t, t - r/be)
                s = s0 + w*s
        u.append(s)
    return u
def stress_fun(regime):
    u = u_point(regime)
    lam = rh*(al**2 - 2*be**2); MU = rh*be**2
    div = sum(sp.diff(u[i], X[i]) for i in range(3))
    s31 = MU*(sp.diff(u[2], x1) + sp.diff(u[0], x3))
    s32 = MU*(sp.diff(u[2], x2) + sp.diff(u[1], x3))
    s11 = lam*div + 2*MU*sp.diff(u[0], x1)
    s33 = lam*div + 2*MU*sp.diff(u[2], x3)
    return [sp.lambdify((x1, x2, x3, t), MU*e, "numpy") for e in (s31, s32, s11, s33)]   # times mu: moment = mu A D
from numpy.polynomial.legendre import leggauss
xg, wg = leggauss(48)
def quad_point(fun, X1, X2, X3, T):
    XI, ETA = np.meshgrid(xg*a, xg*a, indexing="ij"); W = np.outer(wg, wg)*a*a
    return float(np.sum(W*fun(X1 - XI, X2 - ETA, X3, T)))
def tada_cubic(X1, X2, X3, T, key, n=200001):
    s = np.linspace(0, T, n)
    S = tq.rect_ramp_stress(np.full(n, X1), np.full(n, X2), X3, s, a, a, cL, cT, rho, comp=3 if OPEN else 1)[key]
    return float(np.trapezoid((T - s)*S, s))
for regime, recv, T in (("postS", (900., 400., 300.), 0.40), ("postS", (300., -200., 0.0), 0.30),
                        ("between", (1800., 600., 500.), None)):
    if T is None:
        R = np.sqrt(recv[0]**2 + recv[1]**2 + recv[2]**2)
        T = 0.5*((R + 80)/cL + (R - 80)/cT)
    funs = stress_fun(regime)
    for key, f in zip(("s31", "s32", "s11", "s33"), funs):
        p = quad_point(f, *recv, T); q = tada_cubic(*recv, T, key)
        print(f"{regime:8s} x={recv} t={T:.3f}: {key} point-source quadrature {p:+.6e}  Tada (twice integrated) {q:+.6e}  rel diff {abs(p-q)/max(abs(p),1e-30):.1e}")
