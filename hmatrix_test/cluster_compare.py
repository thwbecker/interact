#!/usr/bin/env python3
"""
cluster_compare.py: does splitting an interaction operator by fault before
geometric clustering compress better than a single geometric tree over all
patches?

It reads the dense interaction matrix and the patch coordinates dumped by
compress_interaction_matrix (-dump_matrix / -dump_coords), builds a simple
H-matrix two ways on the SAME operator at the SAME accuracy, and reports the
storage of each:

  joint  : one geometric cluster tree over all patches (what the HTOOL /
           HACApK / hmmvp backends do, since they cluster on centroids only
           and are blind to which fault a patch belongs to).
  split  : partition patches by their group id (fault) first, then build a
           geometric tree within each fault. This realizes the block
           H-matrix F1F1, F2F2, F1F2, F2F1, each with its own tree on its
           own geometry.

The admissible (far) blocks are stored low-rank at a fixed relative tolerance
(rank from the singular values of the true subblock), inadmissible (near)
blocks are stored dense. Accuracy is therefore held fixed by `tol`; the only
thing that changes between joint and split is the partition.

This is a deliberately simple, backend-independent model meant to expose the
structural effect, not to reproduce any one library's exact storage. Results
are specific to the geometry, tolerance and admissibility settings tried; a
different operator, accuracy band, or admissibility constant may shift them,
and a production H-matrix library may differ in absolute numbers. See
rsf_solve_compression.md for the surrounding discussion.

Usage:
  cluster_compare.py -m A.bin -c coords.txt [-tol 1e-5] [-eta 1.0]
                     [-leaf 16] [-v]
"""
import sys
import argparse
import numpy as np


def load_matrix(binfile):
    info = open(binfile + ".info").readline().split()
    m, n = int(info[0]), int(info[1])
    A = np.fromfile(binfile, dtype="<f8")
    if A.size != m * n:
        raise SystemExit("matrix size %d != %d*%d from %s.info" %
                         (A.size, m, n, binfile))
    return A.reshape(m, n)


class Node:
    __slots__ = ("idx", "lo", "hi", "kids")

    def __init__(self, idx, lo, hi, kids):
        self.idx = idx        # patch indices in this cluster
        self.lo = lo          # bounding-box min corner
        self.hi = hi          # bounding-box max corner
        self.kids = kids      # list of child Nodes, or [] for a leaf

    def is_leaf(self):
        return not self.kids


def build_geometric(idx, xyz, leafsize):
    """Binary tree: split the widest box axis at the median, recursively."""
    pts = xyz[idx]
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    if len(idx) <= leafsize:
        return Node(idx, lo, hi, [])
    ax = int(np.argmax(hi - lo))
    med = np.median(pts[:, ax])
    left = idx[pts[:, ax] <= med]
    right = idx[pts[:, ax] > med]
    if len(left) == 0 or len(right) == 0:   # degenerate split, force halve
        order = idx[np.argsort(pts[:, ax])]
        h = len(order) // 2
        left, right = order[:h], order[h:]
    kids = [build_geometric(left, xyz, leafsize),
            build_geometric(right, xyz, leafsize)]
    return Node(idx, lo, hi, kids)


def build_root(idx, xyz, groups, leafsize, split_by_fault):
    if not split_by_fault:
        return build_geometric(idx, xyz, leafsize)
    # split: one geometric subtree per group/fault, joined under one root
    kids = []
    for g in np.unique(groups[idx]):
        gi = idx[groups[idx] == g]
        kids.append(build_geometric(gi, xyz, leafsize))
    los = np.array([k.lo for k in kids])
    his = np.array([k.hi for k in kids])
    return Node(idx, los.min(axis=0), his.max(axis=0), kids)


def box_dist(a, b):
    """Gap between two axis-aligned boxes (0 if they overlap/touch)."""
    gap = np.maximum.reduce([a.lo - b.hi, b.lo - a.hi,
                             np.zeros_like(a.lo)])
    return float(np.sqrt((gap * gap).sum()))


def diam(a):
    return float(np.sqrt(((a.hi - a.lo) ** 2).sum()))


def admissible(a, b, eta):
    d = box_dist(a, b)
    if d <= 0.0:
        return False
    return min(diam(a), diam(b)) <= eta * d


def collect_blocks(t, s, eta, out):
    """Partition the operator into admissible (far) and dense (near) leaves."""
    if admissible(t, s, eta):
        out.append((t.idx, s.idx, True))
    elif t.is_leaf() and s.is_leaf():
        out.append((t.idx, s.idx, False))
    elif t.is_leaf():
        for sk in s.kids:
            collect_blocks(t, sk, eta, out)
    elif s.is_leaf():
        for tk in t.kids:
            collect_blocks(tk, s, eta, out)
    else:
        for tk in t.kids:
            for sk in s.kids:
                collect_blocks(tk, sk, eta, out)


def block_rank(sub, tol):
    """Numerical rank at relative tolerance tol (largest singular value)."""
    sv = np.linalg.svd(sub, compute_uv=False)
    if sv[0] == 0.0:
        return 0
    return int(np.count_nonzero(sv > tol * sv[0]))


def storage(A, groups, root, eta, tol):
    blocks = []
    collect_blocks(root, root, eta, blocks)
    tot = 0          # stored scalars
    near = 0         # dense (near-field) scalars
    far = 0          # low-rank (far-field) scalars
    cross = 0        # scalars in blocks coupling two different faults
    nfar_dense = 0   # admissible blocks that did not actually compress
    for R, C, adm in blocks:
        full = len(R) * len(C)
        is_cross = (len(np.unique(groups[R])) == 1 and
                    len(np.unique(groups[C])) == 1 and
                    groups[R[0]] != groups[C[0]])
        if adm:
            r = block_rank(A[np.ix_(R, C)], tol)
            cost = r * (len(R) + len(C))
            if cost < full:
                tot += cost
                far += cost
                if is_cross:
                    cross += cost
                continue
            nfar_dense += 1   # admissible but low rank did not pay off
        tot += full
        near += full
        if is_cross:
            cross += full
    return dict(total=tot, near=near, far=far, cross=cross,
                nblocks=len(blocks), nfar_dense=nfar_dense)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-m", "--matrix", required=True, help="dense matrix .bin (with .info)")
    ap.add_argument("-c", "--coords", required=True, help="coords .txt from -dump_coords")
    ap.add_argument("-tol", type=float, default=1e-5, help="relative low-rank tolerance")
    ap.add_argument("-eta", type=float, default=1.0, help="admissibility constant")
    ap.add_argument("-leaf", type=int, default=16, help="leaf cluster size")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args()

    A = load_matrix(a.matrix)
    m, n = A.shape
    co = np.loadtxt(a.coords)
    if co.ndim == 1:
        co = co[None, :]
    xyz = co[:, 1:4].astype(float)
    groups = co[:, 12].astype(int)
    idx = np.arange(m)
    dense = m * n

    res = {}
    for name, split in (("joint", False), ("split", True)):
        root = build_root(idx, xyz, groups, a.leaf, split)
        res[name] = storage(A, groups, root, a.eta, a.tol)

    print("operator: %d x %d, dense scalars = %d" % (m, n, dense))
    ngroups = len(np.unique(groups))
    print("faults/groups: %d  (sizes %s)" %
          (ngroups, ", ".join(str(int((groups == g).sum())) for g in np.unique(groups))))
    print("settings: tol=%g  eta=%g  leaf=%d" % (a.tol, a.eta, a.leaf))
    print()
    hdr = "%-8s %12s %10s %12s %12s %10s" % (
        "tree", "stored", "compr_x", "near(dense)", "far(lowrank)", "cross")
    print(hdr)
    print("-" * len(hdr))
    for name in ("joint", "split"):
        r = res[name]
        print("%-8s %12d %9.2fx %12d %12d %10d" %
              (name, r["total"], dense / r["total"], r["near"], r["far"], r["cross"]))
    print()
    sj, ss = res["joint"]["total"], res["split"]["total"]
    if ss < sj:
        print("=> split stores %.1f%% less than joint (%.2fx vs %.2fx compression);"
              % (100.0 * (sj - ss) / sj, dense / ss, dense / sj))
        print("   isolating the faults helps here.")
    elif ss > sj:
        print("=> split stores %.1f%% MORE than joint; the joint geometric tree"
              % (100.0 * (ss - sj) / sj))
        print("   already separates the faults, so forcing a split does not pay.")
    else:
        print("=> joint and split store the same; no measurable clustering effect.")

    if a.verbose:
        print("\nblocks: joint %d (%d admissible-but-dense), split %d (%d)" %
              (res["joint"]["nblocks"], res["joint"]["nfar_dense"],
               res["split"]["nblocks"], res["split"]["nfar_dense"]))


if __name__ == "__main__":
    main()
