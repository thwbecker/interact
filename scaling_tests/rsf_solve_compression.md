# rsf_solve compression vs speed (single BP5 fault)

Companion note to the MPI scaling test. This sweep fixes the rank count and
varies each H-matrix backend's accuracy tolerance, mapping compression against
speed on the single SEAS BP5 fault. Numbers below are from a tolerance sweep at
three resolutions (2 km, 1 km, 0.5 km, i.e. N about 1000, 4000, 16000) at 50 yr
and rtol 1e-4, run on walter, and one build of each library; other rank counts,
builds, or library versions may behave differently. The compression columns
(compr_x, mem) are independent of rank count; the timing columns (assembly,
matvec, total) depend on it, and the per-rank count was not recorded in the
CSVs, so cross-resolution timing should be read with that caveat while
within-resolution comparisons are clean.

## Setup

The single BP5 fault is deterministic, so compression, memory, assembly time,
and the steady per-matvec cost are clean functions of tolerance with no
chaotic-recurrence confound. Tolerance is each package's own knob and is not
comparable across packages at equal nominal value: the same number means a
different operator accuracy for each (the tolerance-to-error mapping is in
compress_interaction_matrix.md). compr_x is dense memory divided by H-matrix
memory; matvec_ms is the mean per-MatMult cost.

## Results (tolerance sweep, 50 yr, rtol 1e-4, walter)

### 2 km (N about 1000)

| backend | tol  | compr_x | mem_MB | assembly_s | matvec_ms | total_s | nsteps |
|---------|------|--------:|-------:|-----------:|----------:|--------:|-------:|
| dense   | -    |   1.0   |   7.6  |  0.16      | 0.040     | 0.277   | 226    |
| htool   | 1e-2 |   3.8   |   2.0  |  0.17      | 0.063     | 0.313   | 224    |
| htool   | 1e-3 |   2.5   |   3.0  |  0.18      | 0.066     | 0.323   | 226    |
| htool   | 1e-4 |   2.0   |   3.8  |  0.19      | 0.065     | 0.329   | 226    |
| htool   | 1e-5 |   1.6   |   4.7  |  0.20      | 0.061     | 0.336   | 226    |
| htool   | 1e-6 |   1.4   |   5.4  |  0.20      | 0.062     | 0.339   | 226    |
| hacapk  | 1e-1 |   3.8   |   2.0  |  0.07      | 0.116     | 0.307   | 223    |
| hacapk  | 3e-2 |   3.6   |   2.1  |  0.07      | 0.121     | 0.323   | 226    |
| hacapk  | 1e-2 |   3.3   |   2.3  |  0.07      | 0.121     | 0.322   | 226    |
| hacapk  | 3e-3 |   3.0   |   2.5  |  0.08      | 0.124     | 0.327   | 226    |
| hacapk  | 1e-3 |   2.7   |   2.8  |  0.08      | 0.128     | 0.336   | 226    |
| hmmvp   | 1e-3 |   5.8   |   1.3  |  0.14      | 0.115     | 0.362   | 224    |
| hmmvp   | 1e-4 |   4.0   |   1.9  |  0.15      | 0.117     | 0.382   | 226    |
| hmmvp   | 1e-5 |   3.0   |   2.5  |  0.16      | 0.118     | 0.392   | 226    |
| hmmvp   | 1e-6 |   2.5   |   3.1  |  0.18      | 0.122     | 0.421   | 226    |
| hmmvp   | 1e-7 |   1.9   |   3.9  |  0.20      | 0.126     | 0.439   | 226    |

### 1 km (N about 4000)

| backend | tol  | compr_x | mem_MB | assembly_s | matvec_ms | total_s | nsteps |
|---------|------|--------:|-------:|-----------:|----------:|--------:|-------:|
| dense   | -    |   1.0   | 122.1  |  1.70      | 0.506     | 2.88    | 341    |
| htool   | 1e-2 |  11.7   |  10.4  |  1.52      | 0.123     | 1.85    | 338    |
| htool   | 1e-3 |   7.4   |  16.6  |  1.50      | 0.131     | 1.86    | 342    |
| htool   | 1e-4 |   5.3   |  22.9  |  1.54      | 0.133     | 1.90    | 342    |
| htool   | 1e-5 |   4.2   |  29.4  |  1.54      | 0.142     | 1.92    | 341    |
| htool   | 1e-6 |   3.3   |  36.8  |  1.57      | 0.148     | 1.95    | 341    |
| hacapk  | 1e-1 |  10.1   |  12.1  |  0.19      | 0.214     | 0.72    | 337    |
| hacapk  | 3e-2 |   8.8   |  13.8  |  0.22      | 0.238     | 0.81    | 339    |
| hacapk  | 1e-2 |   7.9   |  15.4  |  0.25      | 0.259     | 0.87    | 338    |
| hacapk  | 3e-3 |   6.7   |  18.2  |  0.28      | 0.266     | 0.93    | 340    |
| hacapk  | 1e-3 |   6.0   |  20.2  |  0.31      | 0.288     | 1.00    | 341    |
| hmmvp   | 1e-3 |  16.3   |   7.5  |  0.35      | 0.240     | 0.93    | 336    |
| hmmvp   | 1e-4 |  11.6   |  10.5  |  0.39      | 0.242     | 0.98    | 341    |
| hmmvp   | 1e-5 |   8.3   |  14.7  |  0.43      | 0.242     | 1.02    | 342    |
| hmmvp   | 1e-6 |   6.2   |  19.8  |  0.50      | 0.261     | 1.13    | 341    |
| hmmvp   | 1e-7 |   4.8   |  25.5  |  0.59      | 0.267     | 1.23    | 341    |

### 0.5 km (N about 16000)

The dense operator is roughly 2 GB, so the dense baseline is skipped; compr_x is
relative to that dense size.

| backend | tol  | compr_x | mem_MB | assembly_s | matvec_ms | total_s | nsteps |
|---------|------|--------:|-------:|-----------:|----------:|--------:|-------:|
| htool   | 1e-2 |  35.4   |  55.2  | 42.76      | 0.669     | 46.33   | 722    |
| htool   | 1e-3 |  21.7   |  89.9  | 42.11      | 0.855     | 46.50   | 718    |
| htool   | 1e-4 |  15.4   | 127.2  | 41.98      | 1.362     | 48.86   | 717    |
| htool   | 1e-5 |   9.4   | 206.9  | 42.56      | 3.169     | 58.09   | 718    |
| htool   | 1e-6 |   4.4   | 442.1  | 47.02      | 7.726     | 83.40   | 716    |
| hacapk  | 1e-1 |  26.8   |  73.0  |  0.86      | 0.841     | 5.15    | 713    |
| hacapk  | 3e-2 |  22.0   |  88.6  |  1.00      | 0.982     | 6.02    | 723    |
| hacapk  | 1e-2 |  19.1   | 102.1  |  1.14      | 1.143     | 6.94    | 718    |
| hacapk  | 3e-3 |  16.3   | 119.5  |  1.26      | 1.373     | 8.07    | 715    |
| hacapk  | 1e-3 |  14.7   | 132.7  |  1.40      | 1.554     | 9.08    | 716    |
| hmmvp   | 1e-3 |  47.9   |  40.8  |  1.28      | 0.892     | 5.82    | 720    |
| hmmvp   | 1e-4 |  34.6   |  56.5  |  1.56      | 0.903     | 6.27    | 722    |
| hmmvp   | 1e-5 |  25.0   |  78.2  |  1.89      | 0.981     | 6.93    | 721    |
| hmmvp   | 1e-6 |  17.8   | 110.0  |  2.45      | 1.589     | 10.29   | 717    |
| hmmvp   | 1e-7 |  13.6   | 144.0  |  2.87      | 2.194     | 13.51   | 717    |

### 0.25 km (N about 64000)

The dense operator is roughly 31 GB. A dense baseline run is included here for
timing context (compr_x is relative to it); the H-matrix rows are what matter.

| backend | tol  | compr_x | mem_MB | assembly_s | matvec_ms | total_s | nsteps |
|---------|------|--------:|-------:|-----------:|----------:|--------:|-------:|
| dense   | -    |   1.0   | 31250  |  319.90    | 302.509   | 3535    | 1581   |
| htool   | 1e-2 | 122.6   |  254.9 | 4026.65    | 29.452    | 4088    | 1619   |
| htool   | 1e-3 |  71.4   |  437.4 | 3996.27    | 28.657    | 4072    | 1586   |
| htool   | 1e-4 |  16.2   | 1934.3 | 4024.99    | 63.657    | 4559    | 1588   |
| htool   | 1e-5 |   5.2   | 5986.6 | 4020.10    | 131.035   | 5362    | 1587   |
| htool   | 1e-6 |   3.7   | 8462.1 | 4013.50    | 158.807   | 5708    | 1586   |
| hacapk  | 1e-1 |  68.5   |  456.1 |    5.16    | 6.194     | 73.8    | 1572   |
| hacapk  | 3e-2 |  57.7   |  541.4 |    5.93    | 7.024     | 86.1    | 1580   |
| hacapk  | 1e-2 |  50.3   |  621.4 |    6.66    | 7.887     | 95.3    | 1572   |
| hacapk  | 3e-3 |  43.2   |  723.8 |    7.57    | 8.946     | 106.5   | 1577   |
| hacapk  | 1e-3 |  38.6   |  810.4 |    8.23    | 9.770     | 115.6   | 1579   |
| hmmvp   | 1e-3 | 149.8   |  208.6 |    6.29    | 5.083     | 66.6    | 1636   |
| hmmvp   | 1e-4 | 111.1   |  281.2 |    7.91    | 5.553     | 71.2    | 1583   |
| hmmvp   | 1e-5 |  79.7   |  392.3 |    9.88    | 6.095     | 79.0    | 1583   |
| hmmvp   | 1e-6 |  56.0   |  558.1 |   12.77    | 9.883     | 125.3   | 1587   |
| hmmvp   | 1e-7 |  42.3   |  739.4 |   14.92    | 11.612    | 144.4   | 1587   |

At this resolution hmmvp and HACApK take the 31 GB operator down to 200 to 810 MB
and assemble it in 5 to 15 s, for a 25 to 50x end-to-end speedup over dense.
HTOOL compresses comparably well (up to 123x at 1e-2) and its matvec stays 2 to
10x faster than dense, but its assembly is about 4000 s, roughly 12x the dense
assembly and 500 to 800x the HACApK and hmmvp assembly at the same size and rank
count, so it is net slower than dense end to end here. The assembly is also
nearly flat across tolerance (4013 to 4027 s), so it is not the ACA work that
dominates. nsteps is consistent across all backends and tolerances (1572 to
1636), so the compression is not corrupting the trajectory. The HTOOL assembly
cost is examined in its own section below.

## Compression versus N (matched ~1e-6 accuracy)

Reading the compression ratio at the matched band across the three sizes (htool
epsilon 3e-5 interpolated between its 1e-4 and 1e-5 rows; hacapk ztol 1e-1 and
hmmvp tol 1e-7 are exact):

| N (resolution) | htool | hacapk | hmmvp |
|----------------|------:|-------:|------:|
| 1000 (2 km)    |  1.8  |   3.8  |  1.9  |
| 4000 (1 km)    |  4.7  |  10.1  |  4.8  |
| 16000 (0.5 km) | 11.9  |  26.8  | 13.6  |

See compression_vs_N.png (compression ratio versus N, log-log) and
compression_speed_tradeoff.png (per-matvec cost versus compression along each
tolerance sweep, one panel per size).

- Compression grows with N for all three, a bit slower than the O(N log N)
  memory scaling reference (compr_x roughly proportional to N / log N), going
  from single-digit at N = 1000 to 12 to 27x at N = 16000. So the larger the
  fault discretization, the more there is to gain from compression, which is the
  expected H-matrix behaviour and means these strategies matter most exactly
  where the dense operator becomes unaffordable.
- At matched accuracy the three are close at every size, with HACApK
  consistently highest and pulling further ahead as N grows (about 2x htool and
  hmmvp at 0.5 km). Its ztol is conservative for the smooth Okada kernel, so it
  reaches the band while admitting more low-rank blocks. htool and hmmvp track
  each other closely.
- The per-matvec cost frontier (second figure) is flat in compression at 2 km
  and 1 km but opens up at 0.5 km, where htool's matvec runs from 7.7 ms at its
  least-compressed (tight tol) setting down to 0.67 ms at its most-compressed,
  while hacapk and hmmvp stay in the 0.8 to 2.2 ms range. So at large N htool is
  only matvec-competitive at its high-compression end, and its cheap-and-flat
  matvec from the smaller sizes does not carry over.

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
   12 MB), htool is around 4.7x (about 25 MB) but with the fastest matvec
   (about 0.13 ms), and hmmvp is about 4.8x (about 25 MB) with the slowest
   matvec. The raw table can mislead on this point: hmmvp's headline 16x at tol
   1e-3 is its low-accuracy end (its Frobenius-norm error estimate floors near
   1e-6 at this N), not a matched-accuracy compression number.

6. The picture is size-dependent, and two of the 1 km readings on htool do not
   survive to 0.5 km (N about 16000):
   - Compression improves with size for all three, roughly 2 to 3x better at
     0.5 km than at 1 km (loose-end compr_x about 35x htool, 27x hacapk, 48x
     hmmvp).
   - htool's matvec is no longer tolerance-insensitive at large N. At 1 km it
     was flat near 0.13 ms across the tolerance range; at 0.5 km it scales from
     0.63 ms at 1e-2 to 7.68 ms at 1e-6, because the tight settings store a
     large fraction of the operator (442 MB at 1e-6, about 22 percent of dense).
   - htool's assembly grows steeply at large N: about 0.3 s at 1 km, 42 s at
     0.5 km, and about 4000 s at 0.25 km, versus about 1 to 15 s for hacapk and
     hmmvp at the same sizes, and it dominates htool's total. At 0.25 km it is
     net slower than dense end to end. This is examined in "HTOOL assembly at
     large N" below; it is not a serial-assembly bug (the work does distribute
     across ranks), but its cost at production sizes is the open question.
   - As a result the matched-accuracy ranking flips at 0.5 km: hacapk at ztol
     1e-1 is best on every axis (about 27x, 73 MB, 0.8 s assembly, 0.83 ms
     matvec), while htool at the matched eps needs roughly 2 ms matvec plus the
     large assembly, and hmmvp at 1e-7 is about 13.6x at 144 MB with a 2.2 ms
     matvec. So hacapk leads even on the matvec that dominates a long cycle,
     because its conservative ztol reaches the accuracy band while staying well
     compressed.

## Backend choice

The balance depends on problem size, so the choice is size-aware here.

- At 1 km, htool has the cheapest matvec (about 0.13 ms, nearly independent of
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

## Relation to the forward-operator test (hmat_scaling_test.sh)

scaling_tests/hmat_scaling_test.sh is the compress_interaction_matrix analog: a
two-fault 14400-patch geometry swept over MPI ranks at the same matched band,
reporting forward error, assembly, and matvec. It is consistent with this sweep
on the points that matter most, and it validates the matched band directly,
which this compression test only assumes.

- Matched band confirmed. hmat_scaling measures the forward error we do not
  measure here: eps 3e-5 gives about 6.6e-7, ztol 1e-1 about 2.2e-7, and tol
  1e-7 about 1.6e-6. So those three settings really are a matched roughly 1e-6
  band, which underlies every cross-backend statement above. It also confirms
  the hmmvp floor: its error stays near 1.5 to 1.6e-6 from tol 1e-7 out to 1e-3,
  so hmmvp's loose-tolerance high compression is below the accuracy band.
- Assembly agrees in shape and ranking: htool dominates and is strongly
  rank-dependent (about 568 s at np 1 down to about 7 s at np 48 for 14400),
  while hacapk and hmmvp assemble cheaply (about 25 s at np 1 to about 1 s at np
  48). Same ordering as here.
- The matvec cross-backend ranking differs between the two tests, and this is a
  measurement distinction rather than a contradiction. hmat_scaling does a block
  forward apply of 100 vectors, which amortizes per-call overhead and reflects
  raw operator throughput; there hmmvp's matvec looks fastest. rsf_solve applies
  the operator one vector at a time inside the ODE loop through a PETSc MATSHELL,
  paying the full per-call cost every step (including the gather and scatter,
  and for hmmvp a scatter to root on every apply), which is why hmmvp ranks
  among the slowest here while hacapk and htool lead. The two harnesses also use
  different geometry (two makefault faults versus the single planar BP5 fault),
  which contributes to the per-backend differences as well. For choosing a
  backend for actual cycle runs the per-call rsf_solve numbers are the relevant
  ones; for raw forward throughput hmat_scaling is the right reference. If
  hmmvp's throughput is attractive, the lever is its per-call MATSHELL overhead,
  not its compression.

## HTOOL assembly at large N

The HTOOL assembly cost climbs sharply with resolution: about 0.3 s at 1 km
(N about 4000), 42 s at 0.5 km (N about 16000), and about 4000 s at 0.25 km
(N about 64000), the last roughly 12x the dense assembly and 500 to 800x the
HACApK and hmmvp assembly at the same size and rank count. Because that looked
like it could be an effectively serial assembly, the wiring and scaling were
checked directly with this build, mostly on a single core. The short answer is
that it is not serial by construction, and the usual suspects are ruled out, but
the production cost is not reproduced by the small-N behavior and remains open.

What was checked:

  - Not serial. The PETSc layout is an even row split (MatSetSizes with
    PETSC_DECIDE), and an np=1 versus np=2 run on one physical core shows the
    assembly work distributing: wall time stays flat as ranks double, which is
    the signature of constant total work split across ranks, matching the dense
    control, and the per-rank row ranges split evenly. A redundant or serial
    assembly would instead have roughly doubled the one-core wall time at np=2.

  - Not the admissibility constant. Sweeping -mat_htool_eta over 2, 10, 30, 100
    at N=6400 left eta=100 (the default) the fastest assembler and the best
    compressed, because a tighter eta forces more dense near-field blocks through
    the per-entry Okada kernel. Whether this tradeoff inverts at much larger N,
    where loose eta makes the top-level admissible blocks very large and the ACA
    on them expensive, was not reachable on a single core.

  - Not OpenMP oversubscription in this build. htool assembly time at N=6400 was
    flat across OMP_NUM_THREADS of 1, 2, 4, 8 on one core, so the assembly is not
    meaningfully OpenMP threaded here (unlike hmmvp, which pins itself to one
    thread internally).

  - Benign in the reachable range. Single-core assembly scaled about as N^1.35
    from N=1600 to 12780 (about 1.0, 2.5, 6.5, 16.6 s) and stayed faster than the
    dense assembler throughout.

The gap is that extrapolating that benign single-core trend to N=64000 predicts
on the order of 100 to 200 core-seconds, a few seconds on 16 ranks, two to three
orders of magnitude below the observed 4000 s. The 0.5 km point (42 s on 16
ranks) is already about 30x above the single-core extrapolation. So the cost at
production sizes is explained neither by serial wiring nor by the small-N
scaling. The two candidates that remain, both untested here, are a sharp
super-linear breakdown above N about 16000 (the 42 s to 4000 s jump from 0.5 km
to 0.25 km is about N^3 over that interval) and poor parallel efficiency at high
rank count, where adding ranks stops reducing the assembly.

The decisive check needs the cluster: hold N fixed at 0.5 km (N about 16000) and
vary ranks (1, 2, 4, 8, 16), timing the htool assembly alone. If it is roughly
flat in rank count the assembly is effectively serial at scale despite the
correct row split, which would point at a serial section in htool's parallel
assembly path or a load imbalance in its block distribution, and would be worth
raising with the htool or PETSc MATHTOOL maintainers. If it falls with rank count
the assembly is simply expensive and super-linear at this size, and the practical
answer is to prefer the cheaper assemblers. A second, independent check is
-mat_htool_eta at 0.25 km (for example 100 versus 10 versus 2), since the eta
tradeoff measured at small N may invert once the admissible blocks become very
large. All of this is specific to this build of PETSc and htool and to the BP5
geometry; a different htool version or configuration may behave differently.

Until that is resolved, HACApK at ztol near 1e-1 and hmmvp at 1e-3 to 1e-5 are
the backends to use at 0.25 km: both assemble in seconds, compress 40 to 150x,
and run the cycle 25 to 50x faster than dense, while HTOOL's otherwise strong
compression and matvec are buried under its assembly cost at this resolution.

## Caveats

- These numbers are specific to the tested configurations (2 km, 1 km, and
  0.5 km on walter, 50 yr, rtol 1e-4, this build of PETSc with HTOOL, HACApK, and
  hmmvp, and these BP5 parameters). The per-rank count was not recorded in the
  CSVs, so the timing columns should be compared within a resolution rather than
  across resolutions; the compression columns do not have this issue. Other rank
  counts, kernels, or library versions may give different compression and
  timing, and a newer release of any of these libraries may behave differently.
  The size dependence seen across 2 km, 1 km, and 0.5 km is itself a reason not
  to extrapolate any single size to another without measuring.
- The test measures compression and speed, not accuracy. The flat nsteps
  indicates the dynamics are unperturbed at these tolerances, but the
  operator-error-versus-tolerance mapping that justifies the matched-band
  reading lives in compress_interaction_matrix.md, and that mapping was
  characterized near N = 14400, so the exact matched-band values may shift
  somewhat at other sizes.
