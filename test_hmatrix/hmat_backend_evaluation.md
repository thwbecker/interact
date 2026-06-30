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

## Result 4: a two-fault operator confirms the ranking and the scaling

The single-fault tests above use one BP5 plane. To check that the ranking is not
specific to a single compact cluster, `make_two_fault.py` builds a controlled
two-fault geometry (two parallel 2:1 L:W faults, separated 0.5 W in the
fault-normal direction and offset 0.05 W along strike, same strike), and
`sweep_hmat_bAx_2fault.sh` runs the same metrics at total sizes N = 3600, 14400,
57600, 90000 (chosen to bracket the single-fault 0.5km and 0.25km counts). This
operator has two separated clusters and a coherent cross-fault far-field, which
is the structure a matrix-relative tolerance tends to discard. Results are for
this geometry on a 24-rank node; another configuration may differ.

Scaling at eps=1e-3 (`hmat_2f_scaling.png`). Storage is near-linear for every
backend and mode (about N^1.17 to N^1.31 over this range), with hmmvp
matrix-relative smallest in absolute terms, the block-local modes of hmmvp and
HACApK together in the middle, and HTOOL largest. The matvec is where they
separate: hmmvp scales near-linearly (block-local about N^1.1, matrix-relative
about N^1.0), HACApK about N^1.2, and HTOOL steepest at roughly N^1.5. HTOOL has
the fastest matvec at the smallest size but the worst scaling, so it crosses
above hmmvp by N of a few times 10^4; at production size hmmvp has both the
smallest matrix and the fastest matvec. The fitted exponents are approximate
(four sizes, with parallel-overhead noise at N=3600), but the ordering is clear.

Cost versus error per size is in `hmat_2f_size_vs_error.png` and
`hmat_2f_matvec_vs_error.png`. The picture matches the single-fault frontier:
hmmvp owns the envelope, block-local is needed for accuracy and the coherent
direction, and block-local over-resolving at tight eps worsens with N (stored
scalars at N=90000 run 129M at 1e-3, 296M at 1e-4, 821M at 1e-5, 1050M at 1e-6).

### Single versus two-fault at matched total N

Comparing the matched-size pairs directly (single 0.5km against two-fault j60 at
N near 16000, and single 0.25km against two-fault j120 at N near 64000;
`hmat_1f_vs_2f_size_per_dof.png`) isolates what the second fault changes. Three
points stand out.

Compressibility. In the block-local modes the two-fault operator stores about 10
to 20 percent more scalars per DOF than the single fault at matched generic
error (for example hmmvp block-local at eps=1e-3 holds about 875 scalars/DOF on
the single 0.5km fault versus about 1050 on two-fault j60, and about 1140 versus
1315 at the larger pair). That extra content is the genuine cross-fault coupling,
the admissible-but-nonzero blocks between the two clusters. In matrix-relative
mode (hmmvp i3) the per-DOF storage is essentially identical between the two
geometries (about 333 versus 334, and 421 versus 427), because matrix-relative
drops exactly those small cross-fault blocks. The storage scaling exponent is
unchanged by the second fault; only the prefactor shifts up slightly.

Matvec. At matched N the per-matvec times are comparable and track total stored
scalars in both geometries; the residual differences are within the spread set
by N not being exactly equal between the pairs and by differing cluster
locality, and should not be over-read.

Coherent direction. The uniform-slip (coherent) direction is amplified relative
to a generic direction in both geometries, and the amplification grows with
resolution: the median coherent-to-generic ratio runs about 7, 10, 14 on the
single fault over 1km, 0.5km, 0.25km, and about 5, 8, 11, 13 on the two-fault
over N = 3600 to 90000. At matched N the two-fault ratio is in fact slightly
lower than the single fault, not higher: spreading the slip over two partially
separated planes makes the uniform-slip loading marginally less coherent than on
one compact plane. Either way the coherent direction is amplified by roughly an
order of magnitude over a generic one, so matrix-relative mode, which targets the
generic Frobenius error, under-resolves it badly in both cases (mode 3 coherent
error about 2e-2 at eps=1e-3).

The operating point is unchanged from the single-fault case: hmmvp block-local at
eps=1e-3 at every size. The only regime where another package wins is generic
accuracy near 1e-7 at the largest N (HTOOL fastest matvec, HACApK smallest),
where hmmvp block-local has over-resolved; the earthquake-cycle application does
not need that accuracy.

## Per-package recommended settings

hmmvp (`-use_hmatrix 4`): the strongest overall, and the clear choice at
production scale, where it is smallest and has the fastest matvec.

- norm mode: use block-local (`-hmmvp_inorm 1`, BREM) for `rsf_solve`
  earthquake-cycle work, where the coherent loading direction sets the timing,
  and whenever a tight generic accuracy is needed. Use matrix-relative
  (`-hmmvp_inorm 3`, MREM, the hmmvp library's own default but no longer
  interact's) only when loose-to-moderate
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

These recommendations are now the defaults in `set_hmat_defaults_and_options`
(`src/petsc_interact.c`), shared by `compress_interaction_matrix`, the PETSc
matvec path, and `rsf_solve` (via `rsf_init`): hmmvp `inorm=1` (block-local) at
`tol=1e-3`, `eta=3`; HTOOL `eta=3` (previously 100) at `epsilon=1e-4`; HACApK
`inorm=3` at `ztol=1e-4`, `eta=2`. Any of these can still be overridden on the
command line. `bp5/run_new_tests` applies the same per-package modes for cycle
comparisons (it now uses HACApK `inorm=3` rather than the absolute mode 1).

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

The two-fault case (Result 4) is generated by `make_two_fault.py` (the 2:1,
0.5 W separation, 0.05 W offset geometry; resolution set by `jmax`, total
N = 4*jmax^2) and swept by `sweep_hmat_bAx_2fault.sh`, which writes
`hmat_bAx_2f.j<jmax>.dat` in the same column layout. `jmax = 60` and `120`
bracket the single-fault 0.5km and 0.25km counts for the matched-N comparison.

eta was held fixed per package in these sweeps; a dedicated eta sweep would
refine the size/accuracy trade but the norm-mode and per-library differences
dominate what eta does over the 2 to 3 range.
