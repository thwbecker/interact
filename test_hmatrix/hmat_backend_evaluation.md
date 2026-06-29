# Evaluating H-matrix backends at comparable error (`test_hmatrix`)

How to compare HTOOL, HACApK, and hmmvp for the `interact` fault-interaction
operator on a common footing, and what settings to use for each. This
supersedes the per-tolerance comparisons in
`scaling_tests/compress_interaction_matrix.md`: the central point here is that
the nominal tolerance is not a comparable quantity across libraries or norm
modes, so the comparison must be made against a measured error.

These results are from the BP5 strike-slip geometry at 1km, 0.5km, and 0.25km
(N = 4000, 16000, 64000), 24 MPI ranks, June 2026, produced with
`sweep_hmat_bAx.sh` and `sweep_coherent_error.sh`. They reflect this kernel,
geometry, tolerance range, and machine; another configuration, or a newer
version of any of these externally maintained libraries, may rank differently.
The recommendations are starting points.

## The comparison problem: eps is not a common axis

Each library's tolerance controls a different quantity, so two backends at the
same nominal eps are not at the same accuracy. Two effects compound this.

First, the norm mode. A library can stop each low-rank block either when its
error is small relative to that block's own norm (block-local) or relative to
the whole-matrix norm (matrix-relative, also called global). For the generic
operator error the two are honest in different ways: in matrix-relative mode the
measured error comes out approximately equal to eps (the matrix-relative
Frobenius bound doing exactly what it promises), while in block-local mode the
measured error sits far below eps, because holding every small far-field block
to its own relative accuracy over-resolves relative to the global bound. In our
data hmmvp's matrix-relative mode gives a generic error of about eps (gen/eps
roughly 1.4 to 5.4), and its block-local mode gives a generic error roughly
1000x below eps.

Second, the libraries differ from each other even at matched mode. At the same
block-local tolerance hmmvp's compressor reaches an operator error tens of times
smaller than HTOOL's or HACApK's, while storing fewer scalars.

The consequence is that nominal eps spans about two orders of magnitude in
actual error across packages and modes at a fixed value. Comparing "at the same
eps" compares nothing. The fix is to measure an error and use it as the axis.

## Two error metrics, and which to use

`compress_interaction_matrix` reports both, against the dense operator built
from the same kernel:

- generic error, `|A_h x - A x| / |A x|` for random x, averaged over several
  draws. This estimates the operator error in a typical direction and is the
  fundamental `b = A x` accuracy. It is the axis to compare on.
- coherent error, the same quantity for the uniform-slip vector x = 1. This is
  the response to steady, spatially coherent loading, which is the direction
  that sets earthquake-cycle recurrence timing. It is the application-specific
  metric for `rsf_solve`.

The coherent error is a roughly fixed multiple of the generic error, about 7x at
1km and about 10x at 0.5km, growing with resolution because the coherent
far-field sums over more elements. So the uniform-slip direction is not a
different kind of error, it is a fixed amplification of the generic error, and
the only thing that changes its size dramatically is the norm mode (below).

## Result 1: norm mode is decisive for the coherent direction

Block-local compression preserves the coherent (uniform-slip) response far
better than matrix-relative, because matrix-relative stops a small far-field
block once its error is below the tolerance times the whole-matrix norm, which
for a block whose own norm is far below the matrix norm means almost no rank, so
the block is effectively dropped. Each such block contributes a small same-sign
amount to the row sums that set the loading, and dropping enough of them
corrupts the coherent response even though every block is within tolerance.
Block-local holds each block to its own relative accuracy and keeps them.

For hmmvp the effect is dramatic: block-local (BREM, `-hmmvp_inorm 1`) gives a
coherent error 800x to 25000x smaller than matrix-relative (MREM,
`-hmmvp_inorm 3`) at matched eps, at every resolution. For HACApK the two modes
are within about 1x to 80x, and its block-local mode (`-hacapk_inorm 1`, an
absolute threshold) is actually worse than its matrix-relative mode at loose
tolerance on fine meshes, because the absolute threshold is scale-dependent. See
`hmat_coherent_modes.png`.

## Result 2: the fundamental `b = A x` frontier

Plotting the two costs against the generic error (`hmat_size_vs_error.png`,
`hmat_matvec_vs_error.png`):

- size: hmmvp stores the least at every accuracy and resolution. For loose to
  moderate generic accuracy (about 1e-3 to 1e-5) its matrix-relative mode is
  cheapest, since dropping the far-field is exactly what matrix-relative does;
  its block-local mode takes over only for tight accuracy (about 1e-6 and
  below), where matrix-relative would need an impractically small eps. HTOOL and
  HACApK sit above hmmvp throughout.
- matvec: the ranking changes with N because matvec time tracks stored scalars.
  At 1km HTOOL has the fastest raw matvec (about 0.2 ms) with hmmvp close; by
  0.25km hmmvp is fastest (about 5 ms vs HTOOL 9 ms), because it stores far less.
  This settles an earlier caveat: hmmvp's MPI matvec scales fine, the
  non-scaling seen before was the separate OpenMP in-memory path.

At the production resolution (0.25km), to reach a generic error near 2e-6 (a
coherent error near 3e-5, adequate for cycle timing): hmmvp block-local at eps
1e-3 costs 73M scalars, 557 MB, 5.1 ms; HTOOL needs eps 1e-4 at 154M, 1.17 GB,
9.3 ms; HACApK block-local needs eps 1e-4 at 113M, 860 MB, 11.1 ms. hmmvp
delivers it at about half the memory and half the matvec time.

## Result 3: do not over-tighten hmmvp block-local on fine meshes

Block-local hmmvp over-resolves as eps tightens on a fine mesh. At 0.25km its
storage runs 73M at 1e-3, 155M at 1e-4, 704M at 1e-5, and 955M (7.3 GB) at 1e-6,
with the matvec climbing 5, 14, 52, 119 ms in step; 955M is 23 percent of the
4.1-billion-entry dense operator. The waste is concentrated in the 1e-4 to 1e-5
step, where storage jumps 4.5x (155M to 704M) for essentially no accuracy gain
(generic 1.8e-7 to 1.2e-7), so the 1e-5 setting is strictly dominated. Tightening
to 1e-6 does finally resolve much better (generic 6e-10, coherent 8e-9) but at a
cost (7.3 GB, 119 ms) that is rarely worth paying. Its matrix-relative mode does
not blow up but pays in coherent accuracy. The practical operating point at fine
resolution is block-local at a loose eps (1e-3 gives coherent 2.7e-5 at 73M and
5.1 ms), reserving 1e-6 only when a near-dense operator is genuinely required.

## Per-package recommended settings

hmmvp (`-use_hmatrix 4`): the strongest overall, and the clear choice at
production scale, where it is smallest and has the fastest matvec.

- norm mode: use block-local (`-hmmvp_inorm 1`, BREM) for `rsf_solve`
  earthquake-cycle work, where the coherent loading direction sets the timing,
  and whenever a tight generic accuracy is needed. Use matrix-relative
  (`-hmmvp_inorm 3`, MREM, hmmvp's own default) only when loose-to-moderate
  generic `b = A x` accuracy at minimum storage is all that is wanted; there eps
  reads out directly as the operator error.
- eps: on fine meshes operate block-local at a loose tolerance (1e-3, perhaps
  5e-4). Do not tighten it further: storage explodes for negligible accuracy
  gain (Result 3). At 0.25km, `-hmmvp_inorm 1 -hmmvp_tol 1e-3` is the sweet spot
  (generic 2e-6, coherent 2.7e-5, 73M).
- eta: `-hmmvp_eta 3` (the value tested). BREM costs more than MREM, factors of
  roughly 1.5 to 4 in storage, which is the price of the coherent accuracy.

HTOOL (`-use_hmatrix 1`): block-local only (its epsilon is block-relative; there
is no matrix-relative mode). Fastest raw matvec at small and moderate N; loses
to hmmvp on both size and matvec as N grows because it stores more.

- eta: `-mat_htool_eta 3` was used here and is reliable. Larger eta compresses
  more but the partial-pivoting ACA can mis-converge on large
  marginally-separated blocks for this sign-structured kernel; eta around 100 is
  not recommended (see `scaling_tests/compress_interaction_matrix.md`).
- eps: set a few times below the target generic error; the achieved error tracks
  eps closely (roughly eps/40 here).

HACApK (`-use_hmatrix 3`): assembly-friendly and deterministic
(rank-invariant), but stores more and has a slower matvec than hmmvp at these
sizes.

- norm mode: matrix-relative (`-hacapk_inorm 3`, `param(61)=3`) is the robust
  default and is HACApK's own default. Its block-local mode (`-hacapk_inorm 1`)
  is an absolute threshold rather than a block-relative one, so it is
  scale-dependent and can be worse than matrix-relative at loose tolerance on
  fine meshes; prefer it only at tight tolerance, or treat the absolute mode as
  needing a block-relative reimplementation before relying on it.
- eta: `-hacapk_eta 2` (`param(51)`, HACApK's default). ztol is conservative for
  this smooth Okada operator, so a loose ztol already gives high accuracy.

## How to reproduce

`sweep_hmat_bAx.sh` sweeps both modes and writes
`hmat_bAx.<res>.dat` with columns
`backend imode eps generic_err coherent_err stored_scalars mbytes matvec_ms`.
Set `res` and `ncore` at the top, run from `test_hmatrix/` with the PETSc
environment loaded as for `run_new_tests`, then plot the two costs against
`generic_err`. The generic-error probe is built into
`compress_interaction_matrix` (it prints `random-x rel err: mean = ...`
alongside the x=1 coherent error). `sweep_coherent_error.sh` is the
coherent-only variant kept for the cycle-timing question.

eta was held fixed per package in these sweeps; a dedicated eta sweep would
refine the size/accuracy trade but the norm-mode and per-library differences
dominate what eta does over the 2 to 3 range.
