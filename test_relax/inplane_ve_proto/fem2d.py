#!/usr/bin/env python3
"""
fem2d: independent plane-strain finite-element reference for the
inplane_ve_proto gates.  Bilinear quads on a regular grid, split
nodes along a fault segment (duplicate nodes with a prescribed
displacement jump, master-slave elimination), free surface on top,
zero displacement on the far bottom/side boundaries (half-space
approximation; keep the domain large and compare in an interior
window).  Optional Maxwell viscoelasticity below a given depth
(deviatoric relaxation, backward-Euler per-element with the exact
single-step relaxation factor, global re-solve each step).

Only vertical faults (x = xf, z in [d1, d2]) are supported, slip in
the z direction (dip slip on a 90-degree dipping fault); that is all
the gates need.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def build_mesh(Lx, Lz, nx, nz, xf, d1, d2):
    """regular grid on [-Lx,Lx]x[0,Lz]; returns node coords, elements,
    and the split-node bookkeeping.  Fault nodes at x=xf, d1<=z<=d2
    are duplicated; elements with centers x > xf reference the right
    copies."""
    xs = np.linspace(-Lx, Lx, nx + 1)
    zs = np.linspace(0.0, Lz, nz + 1)
    if not np.any(np.isclose(xs, xf)):
        raise SystemExit("fem2d: fault x must lie on a grid line")
    nnx, nnz = nx + 1, nz + 1
    nid = lambda i, j: j*nnx + i          # i x-index, j z-index
    coords = np.array([[xs[i], zs[j]] for j in range(nnz)
                       for i in range(nnx)])
    ifx = int(np.argmin(np.abs(xs - xf)))
    fault_j = [j for j in range(nnz) if d1 - 1e-9 <= zs[j] <= d2 + 1e-9]
    # duplicate ids for fault nodes
    dup = {}
    extra = []
    for j in fault_j:
        dup[nid(ifx, j)] = coords.shape[0] + len(extra)
        extra.append(coords[nid(ifx, j)])
    coords = np.vstack([coords, np.array(extra)]) if extra else coords
    elems = []
    for j in range(nz):
        for i in range(nx):
            n = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
            xc = 0.5*(xs[i] + xs[i + 1])
            if xc > xf:               # right side uses duplicates
                n = [dup.get(q, q) for q in n]
            elems.append(n)
    return coords, np.array(elems), dup, xs, zs

def quad_ke(dx, dz, D):
    """8x8 stiffness of a rectangular bilinear quad, 2x2 Gauss, and
    the strain-displacement matrices at the Gauss points."""
    gp = np.array([-1.0, 1.0])/np.sqrt(3.0)
    Ke = np.zeros((8, 8))
    Bs = []
    for xi in gp:
        for et in gp:
            dN = 0.25*np.array([[-(1 - et), (1 - et), (1 + et), -(1 + et)],
                                [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]])
            dNdx = dN[0]*2.0/dx
            dNdz = dN[1]*2.0/dz
            B = np.zeros((3, 8))
            B[0, 0::2] = dNdx
            B[1, 1::2] = dNdz
            B[2, 0::2] = dNdz
            B[2, 1::2] = dNdx
            w = (dx/2.0)*(dz/2.0)
            Ke += w*(B.T @ D @ B)
            Bs.append((B, w))
    return Ke, Bs

def Dmat(lam, mu):
    return np.array([[lam + 2*mu, lam, 0.0],
                     [lam, lam + 2*mu, 0.0],
                     [0.0, 0.0, mu]])

class Fem:
    def __init__(self, Lx, Lz, nx, nz, xf, d1, d2, lam1, mu1,
                 lam2=None, mu2=None, H=None):
        self.coords, self.elems, self.dup, xs, zs = \
            build_mesh(Lx, Lz, nx, nz, xf, d1, d2)
        self.dx, self.dz = xs[1] - xs[0], zs[1] - zs[0]
        self.nn = self.coords.shape[0]
        self.Lx, self.Lz = Lx, Lz
        # per-element moduli (substrate below H if given)
        zc = np.array([self.coords[e, 1].mean() for e in self.elems])
        self.sub = (zc > H) if H is not None else np.zeros(zc.size, bool)
        self.lam = np.where(self.sub, lam2 if lam2 is not None else lam1, lam1)
        self.mu = np.where(self.sub, mu2 if mu2 is not None else mu1, mu1)
        self.Ke = {}
        self.Bs = None
        # fixed dofs: bottom, left, right boundaries
        c = self.coords
        fixed = (np.isclose(c[:, 1], Lz) | np.isclose(c[:, 0], -Lx)
                 | np.isclose(c[:, 0], Lx))
        self.fixed = np.repeat(fixed, 2)
    def assemble(self, fac=None):
        """fac: per-element stiffness multiplier matrix (3x3) applied
        to D, used by the viscoelastic stepper; None = elastic."""
        rows, cols, vals = [], [], []
        self.gpB = []
        for ie, e in enumerate(self.elems):
            D = Dmat(self.lam[ie], self.mu[ie])
            if fac is not None:
                D = fac[ie] @ D
            Ke, Bs = quad_ke(self.dx, self.dz, D)
            self.gpB.append(Bs)
            dof = np.empty(8, dtype=int)
            dof[0::2] = 2*e
            dof[1::2] = 2*e + 1
            for a in range(8):
                for b in range(8):
                    rows.append(dof[a]); cols.append(dof[b])
                    vals.append(Ke[a, b])
        K = sp.csr_matrix((vals, (rows, cols)),
                          shape=(2*self.nn, 2*self.nn))
        return K
    def solve_elastic(self, buz):
        """prescribed jump [uz] = buz (right minus left) across the
        fault; returns nodal displacements."""
        K = self.assemble()
        u = np.zeros(2*self.nn)
        # master-slave: u[dup] = u[orig] + (0, buz): build transform
        # T u_red = u_full + g
        n = 2*self.nn
        keep = np.ones(n, bool)
        g = np.zeros(n)
        T = sp.identity(n, format="lil")
        for orig, d_ in self.dup.items():
            for c_ in range(2):
                keep[2*d_ + c_] = False
                T[2*d_ + c_, 2*d_ + c_] = 0.0
                T[2*d_ + c_, 2*orig + c_] = 1.0
            g[2*d_ + 1] = buz
        T = T.tocsr()[:, keep]
        Kr = (T.T @ K @ T).tocsr()
        fr = -T.T @ (K @ g)
        # dirichlet
        fixed_r = self.fixed[keep]
        free = ~fixed_r
        ur = np.zeros(keep.sum())
        ur[free] = spla.spsolve(Kr[free][:, free], fr[free])
        u = T @ ur + g
        self.T, self.g, self.keep = T, g, keep
        return u
    def stress(self, u, ie):
        """mean stress tensor (sxx, szz, sxz) of element ie"""
        e = self.elems[ie]
        dof = np.empty(8, dtype=int)
        dof[0::2] = 2*e
        dof[1::2] = 2*e + 1
        D = Dmat(self.lam[ie], self.mu[ie])
        s = np.zeros(3)
        wtot = 0.0
        for B, w in self.gpB[ie]:
            s += w*(D @ (B @ u[dof]))
            wtot += w
        return s/wtot
    def stress_field(self, u):
        return np.array([self.stress(u, ie)
                         for ie in range(len(self.elems))])
    def centers(self):
        return np.array([self.coords[e].mean(axis=0) for e in self.elems])
