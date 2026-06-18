# rsf_solve compression vs speed (single BP5 fault)

Companion note to the MPI scaling test. This sweep fixes the rank count and
varies each H-matrix backend's accuracy tolerance, mapping compression against
speed on the single SEAS BP5 fault. All numbers below are for one configuration
(1 km resolution, np = 8, a 50 yr interseismic run, rtol 1e-4) and one build of
each library; other resolutions, rank counts, builds, or library versions may
behave differently.

## Setup

The single BP5 fault is deterministic, so compression, memory, assembly time,
and the steady per-matvec cost are clean functions of tolerance with no
chaotic-recurrence confound. Tolerance is each package's own knob and is not
comparable across packages at equal nominal value: the same number means a
different operator accuracy for each (the tolerance-to-error mapping is in
compress_interaction_matrix.md). compr_x is dense memory divided by H-matrix
memory; matvec_ms is the mean per-MatMult cost.

## Results (1 km, np = 8, 50 yr, rtol 1e-4)

| backend | tol  | compr_x | mem_MB | assembly_s | matvec_ms | total_s | nsteps |
|---------|------|--------:|-------:|-----------:|----------:|--------:|-------:|
| dense   | -    |   1.0   | 122.1  |  3.12      | 1.897     | 7.37    | 341    |
| htool   | 1e-2 |  13.7   |   8.9  |  3.49      | 0.183     | 4.02    | 342    |
| htool   | 1e-3 |   8.6   |  14.2  |  3.49      | 0.189     | 4.01    | 341    |
| htool   | 1e-4 |   6.3   |  19.5  |  3.56      | 0.188     | 4.08    | 342    |
| htool   | 1e-5 |   4.9   |  24.7  |  3.54      | 0.202     | 4.09    | 342    |
| htool   | 1e-6 |   4.0   |  30.4  |  3.62      | 0.238     | 4.24    | 341    |
| hacapk  | 1e-1 |  10.1   |  12.1  |  0.30      | 0.253     | 0.95    | 337    |
| hacapk  | 3e-2 |   8.8   |  13.8  |  0.36      | 0.286     | 1.10    | 339    |
| hacapk  | 1e-2 |   7.9   |  15.4  |  0.43      | 0.306     | 1.19    | 338    |
| hacapk  | 3e-3 |   6.7   |  18.2  |  0.46      | 0.344     | 1.31    | 340    |
| hacapk  | 1e-3 |   6.0   |  20.2  |  0.48      | 0.371     | 1.38    | 341    |
| hmmvp   | 1e-3 |  16.3   |   7.5  |  0.47      | 0.312     | 1.24    | 336    |
| hmmvp   | 1e-4 |  11.6   |  10.5  |  0.58      | 0.320     | 1.38    | 340    |
| hmmvp   | 1e-5 |   8.3   |  14.7  |  0.67      | 0.333     | 1.50    | 342    |
| hmmvp   | 1e-6 |   6.2   |  19.8  |  0.80      | 0.365     | 1.70    | 341    |
| hmmvp   | 1e-7 |   4.8   |  25.5  |  0.91      | 0.380     | 1.83    | 341    |

## Findings

1. Deterministic and dynamics-neutral across the sweep. nsteps is 337 to 342 in
   every row including dense, so at 1 km none of the swept tolerances, down to
   the loosest, perturbs the integration. The cost differences are in the
   operator, not the time-stepping. This is specific to 1 km and these
   tolerances; coarser meshes or looser settings could change it.

2. htool: assembly is the expensive part (about 3.5 s, comparable to dense), but
   the matvec is the fastest of the three (about 0.18 to 0.24 ms versus 1.9 ms
   for dense) and is nearly insensitive to tolerance. Tightening epsilon from
   1e-2 to 1e-6 grows stored memory about 3.4x (8.9 to 30 MB) while the matvec
   moves only from 0.18 to 0.24 ms. So with htool one can tighten for accuracy
   at little matvec cost, paying in memory and in the one-time assembly.

3. hacapk: the opposite balance, very cheap assembly (about 0.3 to 0.5 s) and a
   matvec that does scale with tolerance (about 0.25 to 0.37 ms).

4. hmmvp: intermediate, with assembly growing as tolerance tightens (about 0.5
   to 0.9 s) and matvec about 0.31 to 0.38 ms.

5. Cross-backend comparison should be made at matched accuracy, not matched
   nominal tolerance. At the matched roughly 1e-6 band derived in
   compress_interaction_matrix.md (approximately htool epsilon 3e-5, hacapk ztol
   1e-1, hmmvp tol 1e-7) for this case: hacapk compresses best (about 10x, about
   12 MB), because its ztol is conservative for the smooth Okada kernel so a
   loose 1e-1 already reaches the band; htool is around 5 to 6x (about 20 MB)
   but with the fastest matvec; hmmvp is about 4.8x (about 25 MB) with the
   slowest matvec, the weakest of the three at matched accuracy here. The raw
   table can mislead on this point: hmmvp's 16x at tol 1e-3 is its low-accuracy
   end (its Frobenius-norm error estimate floors near 1e-6 at this N), not a
   matched-accuracy compression number.

## Backend choice for this configuration

- Long, matvec-dominated cycles: htool, for the cheapest matvec, which can be
  tightened for accuracy at little matvec cost. The price is an expensive
  one-time assembly and more memory at matched accuracy.
- Memory- or assembly-constrained workflows, or many short runs: hacapk at ztol
  near 1e-1, for the best compression at target accuracy with near-free assembly
  and a competitive matvec.
- hmmvp showed no regime where it was strictly best at 1 km in this test. This
  is specific to the size and settings and could differ at other sizes or with
  other builds; hmmvp is designed and tuned for large problems.

## Caveats

- These numbers are specific to the tested configuration (1 km, np 8, 50 yr,
  rtol 1e-4, this build of PETSc with HTOOL, HACApK, and hmmvp, and these BP5
  parameters). Other resolutions, rank counts, kernels, or library versions may
  give different compression and timing, and a newer release of any of these
  libraries may behave differently.
- The test measures compression and speed, not accuracy. The flat nsteps
  indicates the dynamics are unperturbed at these tolerances, but the
  operator-error-versus-tolerance mapping that justifies the matched-band
  reading lives in compress_interaction_matrix.md.
- Worth confirming at 0.5 km, where compression ratios rise and the matvec
  advantage over dense widens, whether the matched-accuracy ranking holds.
