# rsf_solve compression vs speed (single BP5 fault)

Companion note to the MPI scaling test. This sweep fixes the rank count and
varies each H-matrix backend's accuracy tolerance, mapping compression against
speed on the single SEAS BP5 fault. Numbers below are for two configurations
(1 km at np = 8 and 0.5 km at np = 16, both 50 yr interseismic runs at rtol
1e-4) and one build of each library; other resolutions, rank counts, builds, or
library versions may behave differently.

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

## Results (0.5 km, np = 16, 50 yr, rtol 1e-4)

N is about 16000, so the dense operator is roughly 2 GB and the dense baseline
is skipped. compr_x is relative to that dense size.

| backend | tol  | compr_x | mem_MB | assembly_s | matvec_ms | total_s | nsteps |
|---------|------|--------:|-------:|-----------:|----------:|--------:|-------:|
| htool   | 1e-2 |  35.4   |  55.2  | 42.57      | 0.631     | 45.98   | 722    |
| htool   | 1e-3 |  21.7   |  89.9  | 41.69      | 0.838     | 46.02   | 718    |
| htool   | 1e-4 |  15.4   | 127.2  | 41.36      | 1.337     | 48.14   | 717    |
| htool   | 1e-5 |   9.4   | 206.9  | 42.33      | 3.205     | 57.78   | 718    |
| htool   | 1e-6 |   4.4   | 442.1  | 47.27      | 7.682     | 83.64   | 716    |
| hacapk  | 1e-1 |  26.8   |  73.0  |  0.81      | 0.833     | 5.06    | 713    |
| hacapk  | 3e-2 |  22.0   |  88.6  |  0.94      | 0.975     | 5.92    | 723    |
| hacapk  | 1e-2 |  19.1   | 102.1  |  1.09      | 1.133     | 6.82    | 718    |
| hacapk  | 3e-3 |  16.3   | 119.5  |  1.28      | 1.414     | 8.27    | 715    |
| hacapk  | 1e-3 |  14.7   | 132.7  |  1.32      | 1.565     | 9.04    | 716    |
| hmmvp   | 1e-3 |  47.9   |  40.8  |  1.24      | 0.887     | 5.84    | 723    |
| hmmvp   | 1e-4 |  34.6   |  56.5  |  1.74      | 0.919     | 6.46    | 718    |
| hmmvp   | 1e-5 |  25.0   |  78.2  |  1.87      | 0.968     | 6.82    | 718    |
| hmmvp   | 1e-6 |  17.8   | 110.0  |  2.40      | 1.605     | 10.31   | 717    |
| hmmvp   | 1e-7 |  13.6   | 144.0  |  3.46      | 2.225     | 14.24   | 717    |

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
   1e-1, hmmvp tol 1e-7): at 1 km, hacapk compresses best (about 10x, about
   12 MB), htool is around 5 to 6x (about 20 MB) but with the fastest matvec
   (about 0.19 ms), and hmmvp is about 4.8x (about 25 MB) with the slowest
   matvec. The raw table can mislead on this point: hmmvp's headline 16x at tol
   1e-3 is its low-accuracy end (its Frobenius-norm error estimate floors near
   1e-6 at this N), not a matched-accuracy compression number.

6. The picture is size-dependent, and two of the 1 km readings on htool do not
   survive to 0.5 km (N about 16000):
   - Compression improves with size for all three, roughly 2 to 3x better at
     0.5 km than at 1 km (loose-end compr_x about 35x htool, 27x hacapk, 48x
     hmmvp).
   - htool's matvec is no longer tolerance-insensitive at large N. At 1 km it
     was flat near 0.19 ms across the tolerance range; at 0.5 km it scales from
     0.63 ms at 1e-2 to 7.68 ms at 1e-6, because the tight settings store a
     large fraction of the operator (442 MB at 1e-6, about 22 percent of dense).
   - htool's assembly grows steeply at large N: about 42 s at 0.5 km np 16,
     versus about 1 to 3 s for hacapk and hmmvp, and it dominates htool's total
     on a short run. Part of this is rank count (htool assembly scales with
     ranks, so more ranks would reduce it), but it remains the most expensive
     assembler here.
   - As a result the matched-accuracy ranking flips at 0.5 km: hacapk at ztol
     1e-1 is best on every axis (about 27x, 73 MB, 0.8 s assembly, 0.83 ms
     matvec), while htool at the matched eps needs roughly 2 ms matvec plus the
     large assembly, and hmmvp at 1e-7 is about 13.6x at 144 MB with a 2.2 ms
     matvec. So hacapk leads even on the matvec that dominates a long cycle,
     because its conservative ztol reaches the accuracy band while staying well
     compressed.

## Backend choice

The balance depends on problem size, so the choice is size-aware here.

- At 1 km, htool has the cheapest matvec (about 0.19 ms, nearly independent of
  tolerance) and is attractive for long, matvec-dominated cycles, at the cost of
  an expensive one-time assembly and more memory at matched accuracy. hacapk at
  ztol near 1e-1 is the choice when memory or assembly time matters, or for many
  short runs.
- At 0.5 km, hacapk at ztol near 1e-1 is the all-around pick at matched
  accuracy: best compression, near-free assembly, and the fastest matvec, so it
  wins even for long cycles. htool's matvec advantage does not hold at this size
  and its assembly is a liability; it becomes attractive mainly if its assembly
  is amortized over very long runs and given enough ranks to bring the assembly
  down, and even then its matvec at matched accuracy trails hacapk here.
- hmmvp showed no regime where it was strictly best in this test. Its raw
  loose-tolerance compression is the highest of the three, but that is below the
  matched accuracy band; at matched accuracy it trails. This is specific to the
  sizes and settings tested; hmmvp is designed and tuned for large problems and
  may compare differently elsewhere.

## Caveats

- These numbers are specific to the tested configurations (1 km at np 8 and
  0.5 km at np 16, 50 yr, rtol 1e-4, this build of PETSc with HTOOL, HACApK, and
  hmmvp, and these BP5 parameters). Other resolutions, rank counts, kernels, or
  library versions may give different compression and timing, and a newer
  release of any of these libraries may behave differently. The size dependence
  seen between 1 km and 0.5 km is itself a reason not to extrapolate any single
  size to another without measuring.
- The test measures compression and speed, not accuracy. The flat nsteps
  indicates the dynamics are unperturbed at these tolerances, but the
  operator-error-versus-tolerance mapping that justifies the matched-band
  reading lives in compress_interaction_matrix.md, and that mapping was
  characterized near N = 14400, so the exact matched-band values may shift
  somewhat at other sizes.
