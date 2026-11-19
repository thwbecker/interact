"""Reflected waves of a 2-D in-plane (plane strain) half plane by Cagniard-de Hoop (fd_plan.md M11):
the time-domain alternative to the DWN tables (halfplane_dwn.py, fd_hpdwn.c), whose cost grows with
the square of the window. Prototype, validated against halfplane_dwn.py.

Rays: the up-going wave of a line moment (unit slip on a straight segment) is P or S and reflects
at the traction-free surface as P or S: four rays (P,P) (P,S) (S,P) (S,S). In the Laplace domain
(k = i q s, omega = i s) each ray's amplitude F(q) (stress at the receiver per unit slip, without
the exponential exp(-s tau(q)), tau = q X + eta1 zs + eta2 zr, eta_c = sqrt(1/c^2 - q^2), Re >= 0)
is homogeneous of degree 1 in s, so the ramp response (unit slip rate) of a point source is
R(t) = (1/pi) Im[F(q(t)) dq/dt] on the Cagniard path tau(q) = t (upper half plane, from the
Snell-law saddle q_s at the geometric arrival t_s = tau(q_s); closed form for (P,P) and (S,S),
Newton continuation for the converted rays). A segment: point responses integrated along it, the
1/sqrt(t - t0(tau)) fronts removed by a cosine substitution between the points where the arrival
time equals t. Late times: the rays grow like t^5 and cancel in their sum, so direct evaluation
holds to about 40 t0 (t0 the earliest reflected arrival); beyond, R = K t + c - b / t with K the
reflected static rate and c, b fitted at 30 and 40 t0.
Frame: x, z (depth); stresses sxx, sxz, szz in that frame, as halfplane_dwn.py.
"""
import numpy as np


AL, BE, RHO = 6000.0, 3464.0, 2670.0


def eta(q, c):
    e = np.sqrt(1.0/c**2 - q**2 + 0j)
    return np.where(e.real < 0, -e, e)


def amplitudes(q, s_t, alpha=AL, beta=BE, rho=RHO, s=1.0):
    """dict (c1, c2) -> (sxx, sxz, szz) at the receiver, Laplace variable s, k = i q s."""
    lam = rho*(alpha**2 - 2*beta**2)
    mu = rho*beta**2
    k = 1j*q*s
    w = 1j*s
    na = s*eta(q, alpha)
    nb = s*eta(q, beta)
    kb2 = (w/beta)**2
    sx, sz = s_t
    nx, nz = sz, -sx
    m = np.array([[2*sx*nx, sx*nz + sz*nx], [sx*nz + sz*nx, 2*sz*nz]])*mu
    # incident upgoing traction at the surface per c1 (without exp(-nu zs))
    tr = {}
    for c1, nu in (("P", na), ("S", nb)):
        Dr = [1j*k, nu]
        Ds = [-1j*k, -nu]
        g = 1.0/(2*nu)
        u = []
        for i in range(2):
            acc = 0
            for p in range(2):
                if c1 == "S":
                    Gip = (kb2*(1.0 if i == p else 0.0) + Dr[i]*Dr[p])*g/(rho*w**2)
                else:
                    Gip = -(Dr[i]*Dr[p])*g/(rho*w**2)
                for qq in range(2):
                    acc = acc + m[p, qq]*Ds[qq]*Gip
            u.append(acc)
        exx = Dr[0]*u[0]
        ezz = Dr[1]*u[1]
        szz = lam*exx + (lam + 2*mu)*ezz
        sxz = mu*(Dr[1]*u[0] + Dr[0]*u[1])
        tr[c1] = (sxz, szz)
    PxzA = mu*(1j*k*(-na) + (-na)*1j*k)
    PzzA = lam*(1j*k*1j*k + na**2) + 2*mu*na**2
    PxzB = mu*(1j*k*1j*k - nb**2)
    PzzB = lam*(1j*k*nb - nb*1j*k) + 2*mu*(-nb*1j*k)
    det = PxzA*PzzB - PxzB*PzzA
    out = {}
    for c1 in ("P", "S"):
        txz, tzz = tr[c1]
        A = (-txz*PzzB + tzz*PxzB)/det
        B = (-tzz*PxzA + txz*PzzA)/det
        # reflected P (A) and S (B) stresses at the receiver per amplitude (without exp(-nu zr))
        uxA, uzA = 1j*k*A, -na*A
        uxB, uzB = nb*B, 1j*k*B
        out[(c1, "P")] = ((lam + 2*mu)*(1j*k*uxA) + lam*(-na*uzA), mu*(-na*uxA + 1j*k*uzA), lam*(1j*k*uxA) + (lam + 2*mu)*(-na*uzA))
        out[(c1, "S")] = ((lam + 2*mu)*(1j*k*uxB) + lam*(-nb*uzB), mu*(-nb*uxB + 1j*k*uzB), lam*(1j*k*uxB) + (lam + 2*mu)*(-nb*uzB))
    return out



def tau(q, X, h1, c1, h2, c2):
    return q*X + eta(q, c1)*h1 + eta(q, c2)*h2


def dtau(q, X, h1, c1, h2, c2):
    return X - q*h1/eta(q, c1) - q*h2/eta(q, c2)


def d2tau(q, X, h1, c1, h2, c2):
    e1 = eta(q, c1)
    e2 = eta(q, c2)
    return -h1/(c1**2*e1**3) - h2/(c2**2*e2**3)


def saddle(X, h1, c1, h2, c2):
    """the real stationary point q_s in (-1/cmax, 1/cmax) (Snell's law) and the arrival t_s"""
    b = 1.0/max(c1, c2)*(1.0 - 1e-12)
    lo, hi = -b, b
    for it in range(200):
        m = 0.5*(lo + hi)
        if dtau(m, X, h1, c1, h2, c2).real > 0:      # tau' decreases in q
            lo = m
        else:
            hi = m
    qs = 0.5*(lo + hi)
    return qs, tau(qs, X, h1, c1, h2, c2).real


def path(ts, X, h1, c1, h2, c2):
    """q(t), dq/dt for ascending times ts (zeros before the arrival)"""
    qs, t0 = saddle(X, h1, c1, h2, c2)
    q = np.zeros(len(ts), complex)
    dq = np.zeros(len(ts), complex)
    d2 = d2tau(qs, X, h1, c1, h2, c2).real
    prev = None
    for n, t in enumerate(ts):
        if t <= t0:
            continue
        if prev is None or t - t0 < 1e-6*t0:
            g = qs + 1j*np.sqrt(2.0*(t - t0)/(-d2))      # series near the saddle
        else:
            g = prev
        for it in range(60):
            f = tau(g, X, h1, c1, h2, c2) - t
            fp = dtau(g, X, h1, c1, h2, c2)
            dg = f/fp
            g = g - dg
            if g.imag < 0:
                g = g.real + 1e-30j
            if abs(dg) < 1e-14*abs(g):
                break
        q[n] = g
        dq[n] = 1.0/dtau(g, X, h1, c1, h2, c2)
        prev = g
    return q, dq, t0


KEYS = (("P", "P"), ("P", "S"), ("S", "P"), ("S", "S"))


def vel(c):
    return AL if c == "P" else BE


def point_q(t, X, h1, c1, h2, c2, nst=24):
    """q and dq/dt at one time for one source point (continuation from the saddle)"""
    qs, t0 = saddle(X, h1, c1, h2, c2)
    if t <= t0:
        return None, None, t0
    ts = t0 + (t - t0)*(np.arange(1, nst + 1)/nst)**2
    q, dq, _ = path(ts, X, h1, c1, h2, c2)
    return q[-1], dq[-1], t0


def arrival(tau, xr, zr, xs, zs, st, key):
    X = xr - (xs + tau*st[0])
    h1 = zs + tau*st[1]
    return saddle(X, h1, vel(key[0]), zr, vel(key[1]))[1]


def gl(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def _root(key, t, xr, zr, xs, zs, st, lo, hi):
    flo = arrival(lo, xr, zr, xs, zs, st, key) - t
    for it in range(60):
        m = 0.5*(lo + hi)
        fm = arrival(m, xr, zr, xs, zs, st, key) - t
        if (fm <= 0) == (flo <= 0):
            lo, flo = m, fm
        else:
            hi = m
    return 0.5*(lo + hi)


def _argmin(key, xr, zr, xs, zs, st, lo, hi):
    gr = 0.5*(np.sqrt(5.0) - 1.0)
    x1, x2 = hi - gr*(hi - lo), lo + gr*(hi - lo)
    f1, f2 = arrival(x1, xr, zr, xs, zs, st, key), arrival(x2, xr, zr, xs, zs, st, key)
    for it in range(80):
        if f1 < f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - gr*(hi - lo)
            f1 = arrival(x1, xr, zr, xs, zs, st, key)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + gr*(hi - lo)
            f2 = arrival(x2, xr, zr, xs, zs, st, key)
    return 0.5*(lo + hi)


def seg_ray(key, t, xr, zr, xs, zs, st, a, ngrid=64, nq=24):
    """sxx, sxz, szz of one ray at time t; breakpoints where the arrival time equals t (a grid
    value equal to t counts), also around an interior minimum of the arrival time between grid
    points (as fd_hpcdh.c)"""
    tg = np.linspace(-a, a, ngrid + 1)
    f = np.array([arrival(v, xr, zr, xs, zs, st, key) for v in tg]) - t
    pts = [-a]
    for i in range(ngrid):
        if (f[i] <= 0) != (f[i + 1] <= 0):
            pts.append(_root(key, t, xr, zr, xs, zs, st, tg[i], tg[i + 1]))
        elif i > 0 and f[i] < f[i - 1] and f[i] <= f[i + 1] and f[i] > 0:
            tm = _argmin(key, xr, zr, xs, zs, st, tg[i - 1], tg[i + 1])
            if arrival(tm, xr, zr, xs, zs, st, key) - t < 0:
                if tm < tg[i]:
                    pts += [_root(key, t, xr, zr, xs, zs, st, tg[i - 1], tm), _root(key, t, xr, zr, xs, zs, st, tm, tg[i])]
                else:
                    pts += [_root(key, t, xr, zr, xs, zs, st, tg[i], tm), _root(key, t, xr, zr, xs, zs, st, tm, tg[i + 1])]
    pts.append(a)
    x, w = gl(nq)
    out = np.zeros(3)
    for ta, tb in zip(pts[:-1], pts[1:]):
        if tb <= ta or arrival(0.5*(ta + tb), xr, zr, xs, zs, st, key) >= t:
            continue
        th = 0.5*np.pi*(x + 1.0)
        tau = ta + (tb - ta)*0.5*(1.0 - np.cos(th))
        jac = (tb - ta)*0.5*np.sin(th)*0.5*np.pi
        for v, wv, jv in zip(tau, w, jac):
            X = xr - (xs + v*st[0])
            h1 = zs + v*st[1]
            q, dq, t0 = point_q(t, X, h1, vel(key[0]), zr, vel(key[1]))
            if q is None:
                continue
            F = np.array(amplitudes(np.array([q]), st)[key])[:, 0]
            out += wv*jv*np.imag(F*dq)/np.pi
    return out


def seg_total(t, xr, zr, xs, zs, st, a, **kw):
    return sum(seg_ray(k, t, xr, zr, xs, zs, st, a, **kw) for k in KEYS)



def first_arrival(xr, zr, xs, zs, st, a):
    tg = np.linspace(-a, a, 33)
    return min(arrival(v, xr, zr, xs, zs, st, ("P", "P")) for v in tg)


def series(ts, xr, zr, xs, zs, st, a, K, mult=40.0, nq=24):
    t0 = first_arrival(xr, zr, xs, zs, st, a)
    tc = mult*t0
    out = np.zeros((len(ts), 3))
    t1, t2 = 0.75*tc, tc
    r1 = seg_total(t1, xr, zr, xs, zs, st, a, nq=nq) - K*t1
    r2 = seg_total(t2, xr, zr, xs, zs, st, a, nq=nq) - K*t2
    b = (r2 - r1)/(1.0/t1 - 1.0/t2)               # r = c - b / t
    c = r2 + b/t2
    for n, t in enumerate(ts):
        if t <= t0:
            continue
        if t <= tc:
            out[n] = seg_total(t, xr, zr, xs, zs, st, a, nq=nq)
        else:
            out[n] = K*t + c - b/t
    return out, t0, tc
