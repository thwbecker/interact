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
through this comparison and, separately, through the real backends. The
clustering result is consistent across the sizes tried:

  - Splitting helps only when the two faults genuinely interleave in space.
    The clear case is two perpendicular faults intersecting through the same
    box (vertical and horizontal, perpendicular normals): every spatial
    cluster contains both faults, so the joint tree cannot avoid mixing
    orientations, and splitting recovers a bounded amount. In the cases tried
    the saving was in the single digits to roughly ten percent of stored
    scalars, largest at a moderate admissibility constant (around eta = 2,
    where there are many far blocks but they are still large enough to carry
    the mixed rank) and smaller at the extremes.

  - For two planar faults of different orientation that are close but
    spatially divergent (for example dip 90 next to dip 45), splitting gives
    no gain and is slightly worse, because different-orientation planar
    faults also separate in space, so the geometric tree already isolates
    them and a forced split only fragments the far blocks.

  - For two parallel same-orientation faults, or two faults far apart,
    joint and split are effectively equal, as expected: the geometric tree
    separates them on its own.

The practical reading is that the intuition is mechanically correct but has
limited bite for two planar faults, because the configurations where
orientation mixing would hurt (different orientation) are largely the same
configurations where the faults are spatially separable and the geometric
tree already does the right thing. The regime where fault-aware clustering
could pay is genuine spatial interleaving of differently oriented patches:
intersecting or conjugate faults sharing a volume, rough or non-planar
single faults whose adjacent patches differ in orientation, or dense networks
of many small faults of varied orientation. Even there the gain is a bounded
constant factor on the admissible far blocks, not an order of magnitude.

One implementation caveat: none of the three backends can be told to split by
fault, since they take centroids only. HTOOL's native library does support a
user-supplied cluster tree, but the PETSc MATHTOOL wrapper and the HACApK and
hmmvp shims used here expose coordinate-based clustering only. So the split
numbers above are what a fault-aware tree could add on top of the joint
result, not something the current backends produce.

## Compressibility survey

The second table in the driver runs HTOOL, HACApK, and hmmvp on each joint
operator and reports their compression ratio (dense divided by H-matrix
storage). At the small sizes used for the clustering panel every backend
compresses little, which is expected: the far field is a small fraction of a
few-hundred-patch operator, and the separated case, with a large well-isolated
cross block, is the only one that compresses appreciably. Raising the panel
scale pushes N up and the backends into a more representative regime; the
broader, larger-N backend scaling and accuracy comparison lives in
`rsf_solve_compression.md`, `compress_interaction_matrix.md`, and
`hmat_scaling_test.sh`, and the recommendations there (which backend for
assembly-bound vs matvec-bound work) are the ones to use for production sizing.
The HACApK percentages above sit above 100 percent at these small N, meaning
its H-matrix is slightly larger than dense; that reflects the small operator
and the conservative default settings, not a limitation of the library, and it
compresses well at the larger N reported in the other notes.

## Files

  - `cluster_compare.py`            joint vs split storage on a dumped operator
  - `cluster_compression_test.sh`   makefault panel driver for both tables
  - dumps are produced by `compress_interaction_matrix -dump_matrix/-dump_coords`
