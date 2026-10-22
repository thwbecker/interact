"""Static test of the virtual free surface in 2-D plane strain (Hok and
Fukuyama 2011 idea): a dipping fault plus a traction-free surface line made
of displacement-discontinuity elements (shear and opening) in a full plane,
against the exact half-plane solution (Crouch and Starfield's half-plane
displacement discontinuity, interact's tdd_coeff, built as a shared
library from src/green/tdd_coeff.F).
Fault traction per unit slip, element i from element j, both shear and
normal, compared as a relative Frobenius norm.
Build the half-plane library first (from the repository root):
  gfortran -cpp -shared -fPIC -O2 -DUSE_DOUBLE_PRECISION -Isrc/includes \
      -o /tmp/tdd/libtdd.so src/green/tdd_coeff.F
(or set TDD_LIB to its path). usage: python3 fs_static2d.py"""
import ctypes, numpy as np, sympy as sp, os, sys
import inplane as ip

MU, NU = 3.2e10, 0.25

# full-plane static stress of unit slip (shear) and unit opening in the source frame
def _cs_funcs():
    x, y, a, G, nu = sp.symbols("x y a G nu", real=True)
    f = -1/(4*sp.pi*(1 - nu))*(y*(sp.atan(y/(x - a)) - sp.atan(y/(x + a))) - (x - a)*sp.log(sp.sqrt((x - a)**2 + y**2))
                              + (x + a)*sp.log(sp.sqrt((x + a)**2 + y**2)))
    d = lambda *v: sp.diff(f, *v)
    op = [2*G*(d(y, y) + y*d(y, y, y)), 2*G*(-y*d(x, y, y)), 2*G*(d(y, y) - y*d(y, y, y))]
    return [sp.lambdify((x, y, a, G, nu), -e, "numpy") for e in op]     # minus: our convention D = u(+) - u(-)
_OPEN = _cs_funcs()

def stress_global(Xs, ss, a, X, comp):
    """global (sxx, sxy, syy) at points X from unit slip (comp 0) or opening (comp 1) on segment (Xs, ss, a)"""
    d = X - Xs; n = np.array([-ss[1], ss[0]])
    xl = d @ ss; yl = d @ n
    if comp == 0: sl = ip.static_segment_stress(xl, yl, a, MU, NU)
    else:
        with np.errstate(all="ignore"):
            sl = [np.asarray(f(xl, yl, a, MU, NU), float)*np.ones_like(np.asarray(xl, float)) for f in _OPEN]
            bad = ~np.isfinite(sl[0]) | ~np.isfinite(sl[1]) | ~np.isfinite(sl[2])
            if np.any(bad):
                # the closed form has removable 0/0 where the receiver lies on the normal through
                # an element end (x = +-a, y != 0): average over a tiny offset
                e = 1e-7*a
                for k, f in enumerate(_OPEN):
                    v = np.asarray(sl[k], float)
                    vp = f(np.asarray(xl, float) + e, yl, a, MU, NU); vm = f(np.asarray(xl, float) - e, yl, a, MU, NU)
                    v[bad] = (0.5*(np.asarray(vp)*np.ones_like(v) + np.asarray(vm)*np.ones_like(v)))[bad]
                    sl[k] = v
    sxx, sxy, syy = [np.asarray(v, float) for v in sl]
    c, s = ss[0], ss[1]          # local -> global rotation
    gxx = c*c*sxx - 2*c*s*sxy + s*s*syy
    gyy = s*s*sxx + 2*c*s*sxy + c*c*syy
    gxy = c*s*(sxx - syy) + (c*c - s*s)*sxy
    return gxx, gxy, gyy

def traction(g, sr):
    """shear (along sr) and normal (along n_r = rot90(sr)) traction on a plane with tangent sr"""
    gxx, gxy, gyy = g; n = np.array([-sr[1], sr[0]])
    tx = gxx*n[0] + gxy*n[1]; ty = gxy*n[0] + gyy*n[1]
    return tx*sr[0] + ty*sr[1], tx*n[0] + ty*n[1]

traction_rows = traction

# exact half plane (y <= 0 solid) from tdd_coeff
_lib = ctypes.CDLL(os.environ.get("TDD_LIB", "/tmp/tdd/libtdd.so"))
def tdd(x, y, cx, cy, a, cosb, sinb, ihalf):
    D = ctypes.c_double; I = ctypes.c_int
    args = [D(x), D(y), D(cx), D(cy), D(a), D(cosb), D(sinb), I(ihalf), D(NU), D(1 - 2*NU), D(2*(1 - NU)),
            D(1/(4*np.pi*(1 - NU))), D(2*MU), D(3 - 4*NU)]
    outs = [D() for _ in range(10)]; iret = I()
    _lib.tdd_coeff_(*[ctypes.byref(v) for v in args], *[ctypes.byref(v) for v in outs], ctypes.byref(iret))
    v = [o.value for o in outs]          # sxxs sxxn syys syyn sxys sxyn uxs uxn uys uyn
    return (v[0], v[4], v[2]), (v[1], v[5], v[3])     # (sxx, sxy, syy) for strike slip, for normal opening

def segments(x0, y0, ang, L, n):
    s = np.array([np.cos(ang), np.sin(ang)]); h = L/n
    X = np.array([[x0, y0]]) + ((np.arange(n) + 0.5)*h)[:, None]*s[None, :]
    return X, s, h/2

def run(dip_deg=30.0, Lf=10e3, nf=20, dtop=1e3, Ls=60e3, hs=500.0, verbose=True):
    ang = -np.radians(dip_deg)                 # fault goes down to the right from depth dtop
    Xf, sf, af = segments(0.0, -dtop, ang, Lf, nf)
    ns = int(round(2*Ls/hs)); Xs, ss, as_ = segments(-Ls + 0.0, 0.0, 0.0, 2*Ls, ns)
    Xs[:, 0] += Lf*np.cos(ang)/2               # centre the surface patch over the fault
    # exact half plane (calibrate the convention on the full plane first)
    Tfull = np.zeros((2, nf, nf)); Thalf = np.zeros((2, nf, nf)); Tours = np.zeros((2, nf, nf))
    for j in range(nf):
        g = stress_global(Xf[j], sf, af, Xf, 0)
        ts, tn = traction(g, sf); Tours[0, :, j] = ts; Tours[1, :, j] = tn
        for i in range(nf):
            for ih, T in ((0, Tfull), (1, Thalf)):
                st, so = tdd(Xf[i, 0], Xf[i, 1], Xf[j, 0], Xf[j, 1], af, sf[0], sf[1], ih)
                a_, b_ = traction(st, sf); T[0, i, j] = a_; T[1, i, j] = b_
    sgn = np.sign(np.sum(Tours[0]*Tfull[0]))
    cal = np.linalg.norm(Tours - sgn*Tfull)/np.linalg.norm(Tours)
    # virtual surface: unknowns (slip, opening) per surface element
    A = np.zeros((2*ns, 2*ns)); B = np.zeros((2*ns, nf)); C = np.zeros((2*nf, 2*ns))
    for j in range(ns):
        for c in (0, 1):
            ts, tn = traction(stress_global(Xs[j], ss, as_, Xs, c), ss)
            A[0::2, 2*j + c] = ts; A[1::2, 2*j + c] = tn
            ts, tn = traction(stress_global(Xs[j], ss, as_, Xf, c), sf)
            C[0:nf, 2*j + c] = ts; C[nf:, 2*j + c] = tn
    for j in range(nf):
        ts, tn = traction(stress_global(Xf[j], sf, af, Xs, 0), ss)
        B[0::2, j] = ts; B[1::2, j] = tn
    D = np.linalg.solve(A, -B)                 # surface discontinuities cancelling the surface tractions
    Tvs = np.concatenate([Tours[0], Tours[1]]) + C @ D
    Tex = sgn*np.concatenate([Thalf[0], Thalf[1]])
    err = np.linalg.norm(Tvs - Tex)/np.linalg.norm(Tex)
    err_full = np.linalg.norm(np.concatenate([Tours[0], Tours[1]]) - Tex)/np.linalg.norm(Tex)
    if verbose:
        print(f"dip {dip_deg} fault {Lf/1e3:g} km x{nf} top {dtop/1e3:g} km | surface +-{Ls/1e3:g} km, h {hs:g} m ({ns} el.): "
              f"full-plane convention check {cal:.1e}, rel. error virtual surface {err:.3e} (full plane without surface {err_full:.3e})")
    return err

if __name__ == "__main__":
    # the scan recorded in fd_plan.md
    for dtop in (250.0, 1000.0):
        for Ls in (10e3, 20e3, 40e3):
            for hs in (1000.0, 500.0, 250.0):
                run(dip_deg=30.0, dtop=dtop, Ls=Ls, hs=hs)
    for dip in (30.0, 60.0):
        for hs in (500.0, 250.0):
            run(dip_deg=dip, dtop=0.0, Ls=20e3, hs=hs)
    run(dip_deg=30.0, dtop=0.0, Ls=20e3, hs=125.0)
