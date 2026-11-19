"""3-D prototype of the free surface for a vertical strike-slip fault: mirror
image plus a layer of opening-only surface elements (fd_plan.md, stage A).
Kernels: Tada quadrant formulas (tada_quad.py) for rectangles, rotated
between element frames; tables cached by relative position (the geometry
is translation invariant along strike, and within the surface plane).
Unknowns: fault strike-slip rates (the image carries the same slip) and
surface opening rates (zero normal traction at their centres). Marching:
densesys.Marcher (NL or first order), prescribed fault slip rates."""
import sys, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import tada_quad as tq

AL, BE, RHO = 6000.0, 3464.0, 2670.0
MU = RHO*BE**2

FAULT_FRAME = np.array([[1.0, 0, 0], [0, 0, 1.0], [0, -1.0, 0]])     # rows e1 (strike), e2 (up), e3 = e1 x e2
SURF_FRAME = np.eye(3)

def stress_global(src_c, frame, a1, a2, comp, X, times):
    """global stress tensors (len(times), 3, 3) at X from unit rate (comp 1 slip along e1, 3 opening)
    on the rectangle |xi1| < a1, |xi2| < a2 of the source frame, ramp since t = 0"""
    d = np.asarray(X, float) - src_c
    x1, x2, x3 = frame @ d
    n = len(times)
    s = tq.rect_ramp_stress(np.full(n, x1), np.full(n, x2), x3, np.asarray(times, float), a1, a2, AL, BE, RHO, comp=comp)
    S = np.zeros((n, 3, 3))
    S[:, 0, 0] = s["s11"]; S[:, 1, 1] = s["s22"]; S[:, 2, 2] = s["s33"]
    S[:, 0, 1] = S[:, 1, 0] = s["s12"]; S[:, 0, 2] = S[:, 2, 0] = s["s31"]; S[:, 1, 2] = S[:, 2, 1] = s["s32"]
    return np.einsum("ia,nab,bj->nij", frame.T, S, frame)          # R^T S R

class FS3D:
    def __init__(self, Lx=4e3, W=2e3, top=0.0, h=500.0, hs=250.0, ext=3e3, margin=2e3, dt=None, L=60, with_layer=True, verbose=True):
        self.dt = dt if dt is not None else 0.5*h/BE; self.L = L
        nx, nz = int(round(Lx/h)), int(round(W/h))
        self.fault = [np.array([-Lx/2 + (i + 0.5)*h, 0.0, -top - (j + 0.5)*h]) for i in range(nx) for j in range(nz)]
        self.surf = []
        if with_layer:
            xs = -Lx/2 - margin + (np.arange(int(round((Lx + 2*margin)/hs))) + 0.5)*hs
            ys = -ext + (np.arange(int(round(2*ext/hs))) + 0.5)*hs      # y = 0 on an element boundary
            self.surf = [np.array([x, y, 0.0]) for x in xs for y in ys]
        Nf, Ns = len(self.fault), len(self.surf); Nc = Nf + Ns
        self.Nf, self.Ns, self.Nc, self.h, self.hs = Nf, Ns, Nc, h, hs
        times = np.arange(L + 1)*self.dt
        Tl = 2.5*(np.hypot(Lx + 2*margin, 2*ext + 2*W) + 2*W)/BE
        tlate = np.array([Tl, 2*Tl])
        cache = {}
        def resp(kind_src, src_c, X, kind_rcv):
            """(K[1..L], Ks) of the receiver's traction for a unit-rate source"""
            key = (kind_src, kind_rcv) + tuple(np.round(np.asarray(X) - src_c, 6)) + ((round(src_c[2], 6),) if kind_src == "f" else ())
            if key in cache: return cache[key]
            if kind_src in ("f", "i"):
                S = stress_global(src_c, FAULT_FRAME, h/2, h/2, 1, X, np.concatenate([times, tlate]))
            else:
                S = stress_global(src_c, SURF_FRAME, hs/2, hs/2, 3, X, np.concatenate([times, tlate]))
            if kind_rcv == "f":
                tr = np.einsum("i,nij,j->n", FAULT_FRAME[0], S, FAULT_FRAME[2])    # strike traction on the fault plane
            else:
                tr = S[:, 2, 2]                                                    # normal traction on the surface
            ks = (tr[-1] - tr[-2])/(tlate[1] - tlate[0])
            K = np.diff(tr[:L + 1])
            cache[key] = (K, ks); return cache[key]
        self.W = np.zeros((L, Nc, Nc)); self.Ks = np.zeros((Nc, Nc))
        rcv = [("f", X) for X in self.fault] + [("s", X) for X in self.surf]
        for c in range(Nc):
            for r, (kr, X) in enumerate(rcv):
                if c < Nf:
                    src = self.fault[c]; img = src*np.array([1.0, 1.0, -1.0])
                    K1, k1 = resp("f", src, X, kr); K2, k2 = resp("i", img, X, kr)
                    K = K1 + K2; ks = k1 + k2
                else:
                    K, ks = resp("s", self.surf[c - Nf], X, kr)
                self.Ks[r, c] = ks; self.W[:, r, c] = K - ks*self.dt
        self.K1 = self.W[0] + self.Ks*self.dt
        self.E = np.r_[np.full(Nf, MU/(2*BE)), np.full(Ns, RHO*AL/2)]
        self.fault_idx = np.arange(Nf); self.surf_idx = np.arange(Nf, Nc)
        # attributes for densesys.Marcher
        self.Kn = np.zeros((Nc, Nc)); self.Wn = np.zeros((L, Nc, Nc)); self.Kn1 = np.zeros((Nc, Nc))
        self.comps = [(k, 0) for k in range(Nc)]
        if verbose: print(f"FS3D: fault {nx} x {nz} of {h:g} m, top {top:g} m, surface layer {Ns} elements of {hs:g} m, {L} lags, {len(cache)} distinct kernels")

    @property
    def fault_view(self): return self

import densesys as ds
class _Adapter:
    """make an FS3D look like a densesys.DenseSys for the marcher"""
    def __init__(self, f):
        for k in ("W", "Ks", "K1", "E", "Kn", "Wn", "Kn1", "L", "dt", "Nc", "comps"): setattr(self, k, getattr(f, k))
        self.fault = f.fault_idx; self.surf = f.surf_idx

def kinematic(f, T=None, trise=1.0, nl=True):
    """uniform prescribed slip rate 1 on the fault for 0 < t <= trise; returns the fault traction history"""
    a = _Adapter(f); n = f.L if T is None else int(round(T/f.dt))
    vf = lambda t, idx: np.full(len(idx), 1.0 if t <= trise + 1e-12 else 0.0)
    m = ds.Marcher(a, nl=nl, vf_fun=vf)
    if len(f.surf_idx) == 0:
        m.B = a.K1 if not nl else 0.5*(a.K1 + np.diag(a.E)) - np.diag(a.E)
    tau = []; vs = []
    for k in range(n):
        if len(f.surf_idx):
            m.step()
        else:
            # no surface unknowns: just evaluate the traction for the prescribed rates
            H, Hn = m._known(); V = m.V.copy(); V[:] = vf(m.t + f.dt, np.arange(f.Nc))
            m.tau = H + m.B @ V; Vm = 0.5*(m.V + V) if nl else V
            m.hist.append(Vm); m.delta = m.delta + Vm*f.dt; m.V = V; m.t += f.dt
        tau.append(m.tau[:f.Nf].copy()); vs.append(np.abs(m.V[f.Nf:]).max() if len(f.surf_idx) else 0.0)
    return np.array(tau), np.array(vs), m
