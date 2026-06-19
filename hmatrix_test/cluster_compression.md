# Fault-aware vs joint clustering, and a compressibility check

This note documents two things added under `scaling_tests/`: a way to export
an interact interaction operator for external clustering studies, and a small
experiment on whether isolating faults from each other before clustering
improves H-matrix compression. The numbers below are specific to the
geometries, accuracy band, and admissibility settings tried here; a different
operator, problem size, tolerance, or library version may shift them, and the
simple in-house H-matrix model used for the clustering comparison is meant to
expose the structural effect rather than to reproduce any one backend's exact
storage.

## Exporting the operator: `-dump_matrix` and `-dump_coords`

`compress_interaction_matrix` now accepts two optional flags that write the
operator and its point cloud to disk, so alternative cluster trees and
admissibility choices can be explored outside interact:

    -dump_matrix <file>   dense interaction matrix, raw row-major float64,
                          A[i*n+j] = stress at receiver i from unit slip at
                          source j (the same operator the H-matrix backends
                          approximate). A companion <file>.info records
                          "m n" and the layout. This path needs a single MPI
                          rank so the whole matrix is local; the coordinate
                          dump works at any rank count.

    -dump_coords <file>   one ASCII row per patch:
                          i x y z strike dip l w area nx ny nz group
                          i.e. centroid, orientation in degrees, half-sizes,
                          area, the unit normal, and the group/fault id. The
                          normal is included on purpose so a clustering
                          routine can use orientation, not just centroid
                          position.

Loading the matrix elsewhere is then a one-liner, for example in numpy:

    import numpy as np
    m, n = map(int, open("A.bin.info").read().split()[:2])
    A = np.fromfile("A.bin", dtype="<f8").reshape(m, n)

`scaling_tests/cluster_compare.py` is a worked example of using these two
files together.

## The clustering question

The H-matrix backends here (HTOOL, HACApK, hmmvp) cluster on patch centroids
and are blind to which fault a patch belongs to, which is a sensible and
common default. The question is whether, for two faults that are close to
each other but geometrically different, one could do better by splitting the
patches by fault first and clustering geometrically within each fault. That
is the block H-matrix F1F1, F2F2, F1F2, F2F1, each block with its own tree on
its own geometry.

The mechanism that could make splitting help: a far-field source cluster that
straddles two close, differently oriented faults radiates a superposition of
two orientation families, so its admissible far block tends toward rank
r1+r2 rather than max(r1,r2). Splitting by fault keeps each source cluster to
a single orientation family. The effect is confined to admissible far blocks;
genuine near-field inter-fault coupling stays dense either way, that is
physics, not a clustering choice.

`cluster_compare.py` tests this directly on a dumped operator. It builds a
simple H-matrix two ways at matched accuracy and reports the stored scalars
of each: `joint`, a single geometric binary tree over all patches; and
`split`, a partition by group id first, then a geometric tree within each
group. Admissible blocks are stored low-rank at a fixed relative tolerance
(rank taken from the singular values of the true subblock), inadmissible
blocks dense, so accuracy is held fixed and only the partition changes.

## What the experiment shows

`scaling_tests/cluster_compression_test.sh` runs a panel of makefault cases
through this comparison and, separately, through the real backends. Two things
govern whether splitting by fault helps: the orientation contrast between the
two faults, which sets the sign, and the resolution (patch count), which sets
the size of the effect.

Sweeping resolution at eta = 2, tol = 1e-5, leaf = 8, the storage saving from
splitting (positive means split stores fewer scalars than joint) behaves as:

    N        cross (perpendicular)    divergent (dip 90 / dip 45)
    192            +4 %                     -1 %
    384            +7 %                     -1 %
    640           +10 %                     -1 %
    1040           +9 %                     -1 %
    1536          +11 %                     +2 %
    6144          +16 %                     +6 %

    parallel (same orientation): about -2 to -4 % at all N
    separated (far apart):       0 % at all N

The intersecting case, two perpendicular faults sharing one box so that every
spatial cluster contains both orientations, gains the most, and the gain grows
with resolution, from a few percent at a few hundred patches to roughly fifteen
percent at several thousand. Two different-orientation faults that are close but
spatially divergent (dip 90 next to dip 45) are a wash at low resolution, where
the geometric tree separates them on its own, but turn positive as resolution
rises (about six percent at N = 6144): once the region where they pass close is
resolved into many fine patches of two orientations, the joint tree begins
mixing them in shared clusters, and splitting recovers that. Two parallel
same-orientation faults never benefit (splitting is slightly worse at all N:
there is no orientation contrast to exploit and the split only fragments the
far blocks), and two faults far apart are unaffected (both trees separate them).

The practical reading is that the intuition is mechanically correct, that the
gain is a bounded constant factor on the admissible far blocks rather than an
order of magnitude, and that it is larger than the smallest tests suggested
because it grows with resolution. It requires an orientation contrast between
the faults; same-orientation faults gain nothing. The regimes where it could
pay are therefore differently oriented faults at production resolution: most
strongly genuine spatial interleaving (intersecting or conjugate faults sharing
a volume, rough or non-planar faults whose adjacent patches differ in
orientation, or dense networks of varied orientation), and to a smaller but no
longer negligible degree even differently oriented faults that merely come
close at high resolution.

These figures are from the simplified in-house model and should be read as a
band rather than a precise curve; the resolution trend is not perfectly
monotonic (the cross case dips slightly near N = 1040). They are also the
headroom a fault-aware tree could add on top of the joint result: the real
backends, which cluster on centroids and cannot split, compress the same joint
operators to ratios close to the model's joint column (see the survey below),
so a production library that already compresses harder may realize only part of
the split saving.

One implementation caveat: none of the three backends can be told to split by
fault, since they take centroids only. HTOOL's native library does support a
user-supplied cluster tree, but the PETSc MATHTOOL wrapper and the HACApK and
hmmvp shims used here expose coordinate-based clustering only. So the split
numbers above are what a fault-aware tree could add on top of the joint
result, not something the current backends produce.

## Compressibility survey

The second table in the driver runs HTOOL, HACApK, and hmmvp on each joint
operator and reports their compression ratio (dense divided by H-matrix
storage). At the small sizes used by default for the clustering panel every
backend compresses little, which is expected: the far field is a small fraction
of a few-hundred-patch operator, and the separated case, with a large
well-isolated cross block, is the only one that compresses appreciably. Raising
the panel scale pushes N up and the backends into a more representative regime:
at SCALE = 48 (N = 1536) the three compress to roughly 1.6 to 3.6x (HTOOL and
HACApK) and 2.3 to 5.3x (hmmvp), and at SCALE = 96 (N = 6144) to roughly 3.4 to
8.2x and 6.3 to 13x, with hmmvp highest and the separated case compressing most.
HACApK in particular moves from near-dense at the few-hundred-patch sizes (its
reported memory sits slightly above the dense matrix there) to a few-fold
compression at these larger N, which confirms the small-N readings reflected the
operator size and the conservative default settings rather than a limitation of
the library. The broader, larger-N backend scaling and accuracy comparison lives
in `rsf_solve_compression.md`, `compress_interaction_matrix.md`, and
`hmat_scaling_test.sh`, and the recommendations there (which backend for
assembly-bound vs matvec-bound work) are the ones to use for production sizing.

## Files

  - `cluster_compare.py`            joint vs split storage on a dumped operator
  - `cluster_compression_test.sh`   makefault panel driver for both tables
  - dumps are produced by `compress_interaction_matrix -dump_matrix/-dump_coords`
