"""
Derivation of the 2-D (plane strain) elastodynamic kernels used by the
in-plane fully dynamic prototype.

Steps
1. G2_ij(x,y,t): displacement at (x,y), time t, from an impulsive line
   force in direction j at the origin (full space). Obtained by
   integrating the 3-D Green's function (Aki and Richards eq. 4.23)
   along the line direction z. The z-integrals were done by hand; the
   result is verified here (i) symbolically against the Navier
   equation away from the source and (ii) numerically against direct
   z-quadrature of the 3-D near-field term plus the analytic far-field
   contributions.
2. Ghat_ij = double time integral of G2_ij (ramp response). The time
   antiderivatives were done by hand and are verified symbolically by
   differentiating twice.
3. Point edge dislocation created at t=0 with slip rate 1 (slip = t):
   equivalent double couple m_xy = m_yx = mu * t. Displacement
   u_i = mu (d_y Ghat_ix + d_x Ghat_iy), stresses from Hooke's law.
   The overall sign is fixed by the static test in the test suite.
4. Code generation: writes inplane_kernel.py with numpy functions
   returning (sxx, sxy, syy) for the ramp source, region-aware (before
   P: zero; between P and S: alpha terms only; after S: all terms).

Conventions: x_i x_j terms use xx, xy, yy; R = sqrt(x^2 + y^2);
alpha, beta wave speeds; rho density; mu = rho beta^2;
lambda = rho (alpha^2 - 2 beta^2).
"""
import sympy as sp

x, y, t = sp.symbols("x y t", real=True)
al, be, rho = sp.symbols("alpha beta rho", positive=True)
mu = rho * be**2
lam = rho * (al**2 - 2 * be**2)
R = sp.sqrt(x**2 + y**2)
X = [x, y]


def Z(c):
    return sp.sqrt(c**2 * t**2 - R**2)


# ---- 1. G2 (both waves arrived; alpha-only region = drop beta terms) ----
def G2_terms(c):
    """Pieces of the 2-D impulse response associated with wave speed c."""
    Zc = Z(c)
    A = Zc * (2 * Zc**2 + 3 * R**2) / (c**3 * t**2)   # near-field, x_i x_j  (factor 2/R^4)
    B = Zc / c                                        # near-field, delta_ij (factor 2/R^2)
    C = 1 / (t**2 * Zc)                               # far-field,  x_i x_j  (factor 2/c^3)
    D = 1 / Zc                                        # far-field S, delta_ij (factor 2/beta)
    return A, B, C, D


def G2(i, j, with_beta=True):
    Aa, Ba, Ca, Da = G2_terms(al)
    Ab, Bb, Cb, Db = G2_terms(be)
    xixj = X[i] * X[j]
    dij = 1 if i == j else 0
    expr = xixj * (2 / R**4 * Aa + 2 / al**3 * Ca) - dij * (2 / R**2 * Ba)
    if with_beta:
        expr += xixj * (-2 / R**4 * Ab - 2 / be**3 * Cb) - dij * (-2 / R**2 * Bb - 2 / be * Db)
    return expr / (4 * sp.pi * rho)


# ---- 2. double time integrals (hand derived), u = c t ----
def Ghat_terms(c):
    u = c * t
    Zc = Z(c)
    acosh = sp.acosh(u / R)
    acos = sp.acos(R / u)
    Ah = (Zc**3 / 3 - R**2 * Zc + R**3 * acos) / c**3
    Bh = (Zc**3 / 3 - R**2 * u * acosh + R**2 * Zc) / (2 * c**3)
    Ch = (Zc - R * acos) / R**2
    Dh = (u * acosh - Zc) / c**2
    return Ah, Bh, Ch, Dh


def Ghat(i, j, with_beta=True):
    Aa, Ba, Ca, Da = Ghat_terms(al)
    Ab, Bb, Cb, Db = Ghat_terms(be)
    xixj = X[i] * X[j]
    dij = 1 if i == j else 0
    expr = xixj * (2 / R**4 * Aa + 2 / al**3 * Ca) - dij * (2 / R**2 * Ba)
    if with_beta:
        expr += xixj * (-2 / R**4 * Ab - 2 / be**3 * Cb) - dij * (-2 / R**2 * Bb - 2 / be * Db)
    return expr / (4 * sp.pi * rho)


def verify():
    import numpy as np
    from scipy.integrate import quad
    print("verifying antiderivatives: d2/dt2 Ghat_terms == G2_terms ...")
    for c in (al, be):
        for h, g in zip(Ghat_terms(c), G2_terms(c)):
            d2 = sp.diff(h, t, 2) - g
            for sub in ({x: 1.3, y: 0.7, t: 1.1, al: 1.8, be: 1.0, rho: 1.0},
                        {x: -0.2, y: 2.1, t: 3.0, al: 1.73, be: 1.0, rho: 2.0},
                        {x: 5.0, y: -0.1, t: 2.9, al: 1.73, be: 1.0, rho: 2.0}):
                val = complex(d2.subs(sub).evalf())
                ref = complex(g.subs(sub).evalf())
                assert abs(val) < 1e-8 * max(1.0, abs(ref)), (c, h, val, ref)
            # first integral must vanish at the arrival time (continuity)
            h1 = sp.diff(h, t)
            for sub in ({x: 1.3, y: 0.7, al: 1.8, be: 1.0, rho: 1.0},):
                Rv = float(sp.sqrt(x**2 + y**2).subs(sub))
                cv = float(c.subs(sub))
                v0 = complex(h.subs(sub).subs(t, Rv / cv * (1 + 1e-4)).evalf())
                v1 = complex(h1.subs(sub).subs(t, Rv / cv * (1 + 1e-4)).evalf())
                assert abs(v0) < 1e-5 and abs(v1) < 0.1, (c, v0, v1)  # sqrt-type onset of far-field terms
    print("  ok")

    print("verifying Navier equation for G2 (after both arrivals) ...")
    subs = {x: 1.3, y: 0.7, t: 2.0, al: 1.8, be: 1.0, rho: 1.0}
    for j in range(2):
        for i in range(2):
            lhs = rho * sp.diff(G2(i, j), t, 2)
            div = sum(sp.diff(G2(k, j), X[k]) for k in range(2))
            rhs = (lam + mu) * sp.diff(div, X[i]) + mu * (sp.diff(G2(i, j), x, 2) + sp.diff(G2(i, j), y, 2))
            val = float((lhs - rhs).subs(subs).evalf())
            ref = float(lhs.subs(subs).evalf())
            assert abs(val) < 1e-9 * max(1.0, abs(ref)), (i, j, val, ref)
    # also in the alpha-only region
    subs2 = dict(subs); subs2[t] = 1.0; subs2[al] = 2.0  # R=1.476, t_alpha=0.74, t_beta=1.48
    for j in range(2):
        for i in range(2):
            Gij = G2(i, j, with_beta=False)
            lhs = rho * sp.diff(Gij, t, 2)
            div = sum(sp.diff(G2(k, j, with_beta=False), X[k]) for k in range(2))
            rhs = (lam + mu) * sp.diff(div, X[i]) + mu * (sp.diff(Gij, x, 2) + sp.diff(Gij, y, 2))
            val = float((lhs - rhs).subs(subs2).evalf())
            assert abs(val) < 1e-9, (i, j, val)
    print("  ok")

    print("verifying G2 against direct z-integration of the 3-D Green's function ...")
    a_, b_, r_ = 1.8, 1.0, 1.0
    xv, yv, tv = 1.3, 0.7, 1.1
    Rv = np.hypot(xv, yv)
    for i in range(2):
        for j in range(2):
            xi = [xv, yv]
            dij = 1.0 if i == j else 0.0

            def near(z, c1, c2):
                r = np.sqrt(Rv**2 + z**2)
                g = (3 * xi[i] * xi[j] / r**2 - dij) / r**3
                return g * tv * (float(r < c1 * tv) - float(r < c2 * tv))
            za = np.sqrt(max(a_**2 * tv**2 - Rv**2, 0))
            zb = np.sqrt(max(b_**2 * tv**2 - Rv**2, 0))
            # near field: integrate over |z| < za, subtract |z| < zb
            n1 = 2 * quad(lambda z: (3 * xi[i] * xi[j] / (Rv**2 + z**2) - dij) / (Rv**2 + z**2) ** 1.5 * tv, 0, za)[0]
            n2 = 2 * quad(lambda z: (3 * xi[i] * xi[j] / (Rv**2 + z**2) - dij) / (Rv**2 + z**2) ** 1.5 * tv, 0, zb)[0] if zb > 0 else 0.0
            # far field (delta functions, analytic)
            fa = 2 * xi[i] * xi[j] / (a_**3 * tv**2 * za)
            fb = -2 * (xi[i] * xi[j] / (b_**2 * tv**2) - dij) / (b_ * zb) if zb > 0 else 0.0
            ref = (n1 - n2 + fa + fb) / (4 * np.pi * r_)
            val = float(G2(i, j, with_beta=zb > 0).subs({x: xv, y: yv, t: tv, al: a_, be: b_, rho: r_}).evalf())
            assert abs(val - ref) < 1e-8 * max(1, abs(ref)), (i, j, val, ref)
    print("  ok")


def generate(path="inplane_kernel.py"):
    """Point edge dislocation ramp source: slip rate 1 along x on the
    plane y=0 at the origin for t>0. Returns stresses sxx, sxy, syy."""
    import re
    out = {}
    for tag, wb in (("ps", True), ("p", False)):
        Gh = [[Ghat(i, j, wb) for j in range(2)] for i in range(2)]
        # double couple m_xy = m_yx = mu * t   ->  u_i = mu (d_y Ghat_ix + d_x Ghat_iy)
        u = [mu * (sp.diff(Gh[i][0], y) + sp.diff(Gh[i][1], x)) for i in range(2)]
        exx = sp.diff(u[0], x)
        eyy = sp.diff(u[1], y)
        exy = (sp.diff(u[0], y) + sp.diff(u[1], x)) / 2
        div = exx + eyy
        sxx = lam * div + 2 * mu * exx
        syy = lam * div + 2 * mu * eyy
        sxy = 2 * mu * exy
        out[tag] = (sxx, sxy, syy)
    # code generation with common subexpression elimination
    lines = ['"""Generated by gen_inplane.py; do not edit by hand.',
             'Stresses at (x,y), time t, from a point edge dislocation at the origin',
             '(slip along x across the plane y=0) created at t=0 with unit slip rate.',
             'Regions: before the P arrival all stresses are zero; between P and S',
             'arrivals use the _p functions; after the S arrival the _ps functions."""',
             'import numpy as np',
             'from numpy import sqrt, arccos, arccosh, pi', '']
    for tag in ("ps", "p"):
        reps, red = sp.cse(out[tag], optimizations="basic")
        lines.append(f"def _sig_{tag}(x, y, t, alpha, beta, rho):")
        for s, e in reps:
            lines.append(f"    {s} = {sp.pycode(e)}")
        lines.append("    return (" + ", ".join(sp.pycode(e) for e in red) + ")")
        lines.append("")
    body = "\n".join(lines)
    body = body.replace("math.sqrt", "sqrt").replace("math.acos", "arccos").replace("math.acosh", "arccosh").replace("math.pi", "pi")
    body += '''

def sigma_point_ramp(x, y, t, alpha, beta, rho):
    """Vectorized over x, y, t (broadcast). Returns (sxx, sxy, syy)."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float); t = np.asarray(t, dtype=float)
    x, y, t = np.broadcast_arrays(x, y, t)
    R = np.hypot(x, y)
    mp = alpha * t > R
    ms = beta * t > R
    sxx = np.zeros(x.shape); sxy = np.zeros(x.shape); syy = np.zeros(x.shape)
    with np.errstate(all="ignore"):
        if ms.any():
            a, b, c = _sig_ps(x[ms], y[ms], t[ms], alpha, beta, rho)
            sxx[ms] = a; sxy[ms] = b; syy[ms] = c
        mpo = mp & ~ms
        if mpo.any():
            a, b, c = _sig_p(x[mpo], y[mpo], t[mpo], alpha, beta, rho)
            sxx[mpo] = a; sxy[mpo] = b; syy[mpo] = c
    return sxx, sxy, syy
'''
    with open(path, "w") as f:
        f.write(body)
    print("wrote", path)


if __name__ == "__main__":
    verify()
    generate()
