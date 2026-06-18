# BigWham full-space backend in interact: implementation and initial scaling notes

## What we did

We wired the BigWham hierarchical-matrix library (Lecampion et al. 2025,
version 0.2.0) into interact's PETSc path as a full-space backend, selectable
with `use_hmatrix == 5`. BigWham owns its own full-space (infinite-medium)
kernel, so the integration assembles a BigWham hierarchical operator directly
from interact's patch geometry and wraps it as a PETSc `MATSHELL`:

- Each interact rectangular patch is mapped to a BigWham `3DR0` element, with the
  four corner nodes ordered so that BigWham's local frame aligns with interact's
  strike and dip directions.
- Elastic constants are taken from the same compile-time macros the Okada kernel
  uses, so the two operators are built from identical moduli.
- The shell `MatMult` applies BigWham's matvec and projects between interact's
  slip/traction ordering and BigWham's element-local degrees of freedom.

Because BigWham is an infinite-medium kernel, the backend is used only with
`-full_space 1`, and the implementation warns if it is invoked otherwise.

## Correctness

Against interact's native full-space Okada operator (the dense reference built
by `compress_interaction_matrix` with `-full_space 1`), the BigWham operator
agrees to machine precision for isolated patches and to the requested ACA
tolerance for compressed problems. Agreement holds across a range of strike and
dip values, for mixed-orientation and randomized-orientation meshes, and the
result is invariant to the matvec thread count. So the backend is correct.

## Scaling observations (this configuration)

Using OpenMP on our machine, across the problem sizes (a few hundred up to
20,000 patches), ACA tolerances, and single- and two-fault geometries we tested:

- Assembly parallelizes well. BigWham's OpenMP assembly showed good speedup with
  thread count and was, at the larger sizes we tried, competitive with the
  MPI-parallel Okada-based backends (HTOOL, HACApK, hmmvp).
- The matvec did not scale with thread count in our tests. It was fastest at a
  single thread and grew slower as the thread count increased, and changing the
  ACA tolerance or using a better-compressing (two-fault) geometry did not change
  that behaviour. As wired today, the per-matvec cost is therefore well above the
  MPI Okada-based backends for this full-space comparison.

## Caveats and status

These observations are specific to BigWham 0.2.0, our build of it, our machine's
OpenMP, and the test configurations above; a different version, build, parallel
environment, or geometry regime may behave differently, and a newer BigWham
release may have changed the matvec path. The assembly results in particular are
encouraging. We are following up with the BigWham authors on the matvec
behaviour. For now we are pausing further work on the BigWham backend; the
existing Okada-based backends remain the path for matvec-bound runs.
