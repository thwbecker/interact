"""3-D prototype of a dipping fault with a full virtual free surface (fd_plan.md,
stage B): planar fault with strike along x, dip delta, slip along strike
(rake 'strike') or along dip ('dip'); surface elements on z = 0 with three
unknowns (slip along x, slip along y, opening) and zero traction (sxz, syz,
szz) at their centres; fault receivers carry the shear traction along the
slip direction (the unknown's partner) and the normal traction (for normal
stress changes). Kernels: Tada quadrant formulas (tada_quad.py), rotated
between element frames; tables cached by relative position.
Fault frame: e1 = x (strike), e2 = up-dip, e3 = e1 x e2 (fault normal);
the fault dips towards -y for 0 < delta < 90 (hanging wall on the -y side)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tada_quad as tq
AL, BE, RHO = 6000.0, 3464.0, 2670.0
MU = RHO*BE**2

def fault_frame(dip_deg):
    d = np.radians(dip_deg)
    e1 = np.array([1.0, 0.0, 0.0]); e2 = np.array([0.0, np.cos(d), np.sin(d)])   # up dip (towards +y)
    return np.array([e1, e2, np.cross(e1, e2)])

def stress_global(src_c, frame, a1, a2, comp, X, times):
    d = np.asarray(X, float) - src_c
    x1, x2, x3 = frame @ d
    s = tq.rect_ramp_stress(np.full(len(times), x1), np.full(len(times), x2), x3, np.asarray(times, float), a1, a2, AL, BE, RHO, comp=comp)
    S = np.zeros((len(times), 3, 3))
    S[:, 0, 0] = s["s11"]; S[:, 1, 1] = s["s22"]; S[:, 2, 2] = s["s33"]
    S[:, 0, 1] = S[:, 1, 0] = s["s12"]; S[:, 0, 2] = S[:, 2, 0] = s["s31"]; S[:, 1, 2] = S[:, 2, 1] = s["s32"]
    return np.einsum("ia,nab,bj->nij", frame.T, S, frame)

class FSDip:
    def __init__(self, Lx=2e3, W=1e3, dip=30.0, top=500.0, h=500.0, hs=250.0, ext=2e3, margin=1e3, rake="dip",
                 L=45, surface=True, verbose=True):
        self.dt = 0.5*h/BE; self.L = L; self.h, self.hs = h, hs
        F = fault_frame(dip); self.F = F
        nx, nz = int(round(Lx/h)), int(round(W/h))
        # fault element (i, j): j = 0 at the top, centres going down dip
        top_c = np.array([0.0, 0.0, -top])
        self.fault = [np.array([-Lx/2 + (i + 0.5)*h, 0.0, 0.0]) + top_c - (j + 0.5)*h*F[1] for i in range(nx) for j in range(nz)]
        self.slip_dir = F[0] if rake == "strike" else F[1]
        self.fcomp = 1 if rake == "strike" else 2
        self.surf = []
        if surface:
            # surface centred over the fault's surface projection (the fault dips towards -y here)
            ymid = -0.5*W*np.cos(np.radians(dip)) * 1.0
            xs = -Lx/2 - margin + (np.arange(int(round((Lx + 2*margin)/hs))) + 0.5)*hs
            ny = int(round((2*ext + W*np.cos(np.radians(dip)))/hs))
            ys = ymid - ny*hs/2 + (np.arange(ny) + 0.5)*hs
            self.surf = [np.array([x, y, 0.0]) for x in xs for y in ys]
        Nf, Ns = len(self.fault), len(self.surf)
        # unknowns: fault slip (Nf), then surface slip x, slip y, opening (3 Ns)
        comps = [("f", k, 0) for k in range(Nf)] + [("s", k, c) for k in range(Ns) for c in (1, 2, 3)]
        self.comps = comps; Nc = len(comps); self.Nf, self.Ns, self.Nc = Nf, Ns, Nc
        times = np.arange(L + 1)*self.dt
        rmax = np.hypot(Lx + 2*margin, 2*ext + W) + W
        tl = 3*rmax/BE + 10*self.dt; tlate = np.array([tl, 2*tl])
        tt = np.concatenate([times, tlate])
        cache = {}
        SF = np.eye(3)
        def traction(kind, S):
            """rows: fault -> (shear along slip, normal); surface -> (sxz, syz, szz)"""
            if kind == "f":
                n = F[2]
                return np.stack([np.einsum("i,nij,j->n", self.slip_dir, S, n), np.einsum("i,nij,j->n", n, S, n)])
            return np.stack([S[:, 0, 2], S[:, 1, 2], S[:, 2, 2]])
        def resp(src_kind, si, comp, rcv_kind, X):
            src = self.fault[si] if src_kind == "f" else self.surf[si]
            rel = np.round(np.asarray(X) - src, 6)
            key = (src_kind, comp, rcv_kind) + tuple(rel) + ((round(src[2], 6),) if src_kind == "f" else ())
            if key not in cache:
                if src_kind == "f": S = stress_global(src, F, h/2, h/2, self.fcomp, X, tt)
                else: S = stress_global(src, SF, hs/2, hs/2, comp, X, tt)
                tr = traction(rcv_kind, S)                         # (rows, times)
                ks = (tr[:, -1] - tr[:, -2])/(tlate[1] - tlate[0])
                K = np.diff(tr[:, :L + 1], axis=1)
                cache[key] = (K, ks)
            return cache[key]
        self.W = np.zeros((L, Nc, Nc)); self.Ks = np.zeros((Nc, Nc))
        self.Wn = np.zeros((L, Nf, Nc)); self.Kn = np.zeros((Nf, Nc))      # normal traction on the fault
        rows_f = list(range(Nf)); rows_s = {(k, c): Nf + 3*k + (c - 1) for k in range(Ns) for c in (1, 2, 3)}
        for col, (sk, si, sc) in enumerate(comps):
            for r in range(Nf):
                K, ks = resp(sk, si, sc, "f", self.fault[r])
                self.W[:, r, col] = K[0] - ks[0]*self.dt; self.Ks[r, col] = ks[0]
                self.Wn[:, r, col] = K[1] - ks[1]*self.dt; self.Kn[r, col] = ks[1]
            for k in range(Ns):
                K, ks = resp(sk, si, sc, "s", self.surf[k])
                for c in (1, 2, 3):
                    row = rows_s[(k, c)]
                    self.W[:, row, col] = K[c - 1] - ks[c - 1]*self.dt; self.Ks[row, col] = ks[c - 1]
        self.K1 = self.W[0] + self.Ks*self.dt; self.Kn1 = self.Wn[0] + self.Kn*self.dt
        self.E = np.array([MU/(2*BE) if (k == "f" or c in (1, 2)) else RHO*AL/2 for (k, i, c) in comps])
        self.fault_idx = np.arange(Nf); self.surf_idx = np.arange(Nf, Nc)
        if verbose: print(f"FSDip: dip {dip:g}, {nx} x {nz} fault elements of {h:g} m ({rake} slip), top {top:g} m, {Ns} surface elements of {hs:g} m (3 unknowns each), {L} lags, {len(cache)} distinct kernels")

import densesys as ds
class _Adapter:
    def __init__(self, f):
        for k in ("W", "Ks", "K1", "E", "L", "dt", "Nc", "comps", "Wn", "Kn", "Kn1"): setattr(self, k, getattr(f, k))
        self.fault = f.fault_idx; self.surf = f.surf_idx

def kinematic(f, T, trise=0.4, nl=True):
    a = _Adapter(f); n = int(round(T/f.dt))
    vf = lambda t, idx: np.full(len(idx), 1.0 if t <= trise + 1e-12 else 0.0)
    m = ds.Marcher(a, nl=nl, vf_fun=vf)
    tau = []; vs = []; sn = []
    for k in range(n):
        if len(f.surf_idx): m.step()
        else:
            H, Hn = m._known(); V = np.ones(f.Nc)*vf(m.t + f.dt, [0])[0]
            m.tau = H + m.B @ V; Vm = 0.5*(m.V + V) if nl else V
            m.hist.append(Vm); m.delta = m.delta + Vm*f.dt; m.V = V; m.t += f.dt
        tau.append(m.tau[:f.Nf].copy()); vs.append(np.abs(m.V[f.Nf:]).max() if len(f.surf_idx) else 0.0)
        sn.append(m.sigma_dyn.copy() if hasattr(m, "sigma_dyn") else np.zeros(f.Nf))
    m.sn_hist = np.array(sn)
    return np.array(tau), np.array(vs), m
