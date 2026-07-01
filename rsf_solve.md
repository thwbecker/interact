# `rsf_solve`: H-matrix backends and performance notes

`rsf_solve` integrates a rate-and-state quasi-dynamic earthquake cycle on an
interact patch geometry. The dense elastic interaction operator is applied once
per Runge-Kutta stage, so the dominant cost at scale is the **matrix-vector
product** with the stress-interaction matrix. `rsf_solve` can apply that operator
through several hierarchical-matrix (H-matrix) backends, selected with
`-use_hmatrix`:

| `-use_hmatrix` | backend | parallelism | tolerance knob (default) | notes |
|:--:|----|----|----|----|
| 0 | dense | MPI | n/a | exact reference; O(N^2) memory and matvec |
| 1 | **HTOOL** (PETSc `MATHTOOL`) | MPI | `-mat_htool_epsilon` (3e-5) | also `-mat_htool_eta` (100), `-mat_htool_compressor` |
| 2 | H2OPUS | serial/GPU | `-h2opus_eta` | |
| 3 | **HACApK** (lattice H-matrix) | MPI | `-hacapk_ztol` (1e-1) | shell matvec into HACApK |
| 4 | hmmvp | MPI | `-hmmvp_tol` (1e-5), `-hmmvp_eta` (3) | run under MPI; matched-accuracy comparisons use tol 1e-7 |

This note focuses on the two production MPI backends, **HTOOL** and **HACApK**,
benchmarked against the **dense** ground truth.

## Source layout

`rsf_solve` was split (2026-06) from a single `rsf_solve.c` into four translation
units plus a shared header, separating the driver, the physics, the option/IC
parsing, and the output. The split is mechanical: the RHS, friction inversion,
initial conditions, backslip loading, and defaults are unchanged from the
single-file version (verified by reproducing the BP5 cycle, see below). All
rate-and-state parameters that used to live in file-scope statics now hang off
`medium->rsf` (a `struct rsf_vars`, allocated in `init_medium_rsf`), so no mutable
rsf state is global.

| file | role | key routines |
|----|----|----|
| `src/rsf_solve.c` | driver | `main`; `rsf_solve_run` (geometry, interaction matrices, backslip loading, initial conditions, TS setup, solve, teardown); `rsf_event_function` (slip-rate threshold crossing for event onset/arrest) |
| `src/rsf_init.c` | setup | `init_medium_rsf` (allocate `medium->rsf`, set `dim`); `rsf_get_settings` (all option parsing, unit conversions, defaults, H-matrix backend selection) |
| `src/rsf_engine.c` | physics | `vel_from_rsf` (invert the regularized friction law for slip rate); `rsf_ODE_RHSFunction` (the quasi-dynamic aging-law RHS, where the `Is` matvec and the backslip `sinc` loading enter); `rsf_domain_check` (per-step validity gate) |
| `src/rsf_output.c` | output | `rsf_init_monitor_and_event` and `rsf_finalize_monitor_and_event` (set up and tear down the TS monitor, event detector, gather, and output files); `rsf_TS_Monitor` (state-change logging to `rsf_monitor.dat`, periodic stats, optional per-group GMT slip-rate fields); `rsf_post_event`; and the field-output group helpers (`rsf_build_groups`, `rsf_group_coords`, `rsf_write_group_geometry`, `rsf_free_groups`, `rsf_dcmp`) |
| `src/includes/rsf.h` | header | the output/grid/settings structs (`rsf_out_ctx`, `rsf_group_grid`, `rsf_solve_settings`), the rsf prototypes, and the single shared `rsf_par_static` (the geometry pointer the domain check reads) |

The build links all four objects via `RSF_SOLVE_OBJS` in the makefile:

    RSF_SOLVE_OBJS = $(ODIR)/rsf_solve.o $(ODIR)/rsf_engine.o $(ODIR)/rsf_init.o $(ODIR)/rsf_output.o

One build caveat worth flagging: this split also moved three `Vec` fields out of
`struct med` and added the `medium->rsf` pointer, which shifts the offsets of
later `struct med` members. After pulling that change, do a clean rebuild
(`rm objects/*.o objects/*.a`, or `make clean`) so every object sees the same
`structures.h` layout. A partial rebuild that reuses a stale `petsc_interact.o`
(which writes `medium->rs/re/rn`) against freshly built rsf objects will read the
ownership range at the wrong offset and segfault in the monitor on the first
output step. This is a stale-object hazard, not a code defect: a clean build
fixes it.

## BP5 accuracy and refactor validation

The split is validated by reproducing the SEAS BP5-QD benchmark and comparing
against the single-file predecessor and against HBI (`sozawa94/hbi`).

Setup: `bp5/` 1 km inputs (`geom_bp5_1km.in`, `rsf_bp5_1km.dat`, `ic_bp5_1km.in`,
`dc_bp5_1km.in`), dense backend (`-use_hmatrix 0`), `-rtol 1e-4`, run to 250 yr
(just past the first spontaneous system event). Accuracy is read from the first
system-event onset in `rsf_events.dat` and the slip-rate trace in
`rsf_monitor.dat`.

| code | first spontaneous event | note |
|----|----|----|
| interact `rsf_solve`, split (`rsf_*.c`) | 234.2937 yr | this refactor |
| interact `rsf_solve`, single-file predecessor | 234.2938 yr | matches to 1e-4 yr |
| HBI (lattice HACApK, near-dense) | about 234.2 to 234.4 yr | exact value depends on the onset threshold |

The split matches the predecessor to about 1e-4 yr over a 234 yr interval, which
is time-step accumulation at the ULP level rather than a physics difference, and
tracks HBI to within roughly 0.14 yr (about 0.06 percent), the same agreement
documented for the pre-split code (`bp5/README.md`, `bp5_physics_crosscheck.md`).
The slip-rate traces overlay across the seeded event near t = 0, the interseismic
plateau at the plate rate (log10 of max slip rate near -9), and the first
spontaneous event near 234 yr (`bp5_rsf_vs_hbi_split.png`). This is specific to
the 1 km dense BP5 configuration; the per-backend H-matrix accuracy at this and
other resolutions is in the benchmark table below.

## Recommended defaults for a single-fault cycle

Based on the scaling and compression studies (`scaling_tests/rsf_solve_scaling_test.sh`,
`scaling_tests/hmat_scaling_test.sh`, `scaling_tests/rsf_solve_compression.md`):

- Resolution: 1 km is the official BP5 resolution. The cohesive zone
  `L_b = G D_RS / (b sigma)` is about 6 km, so 2 km resolves it with only ~3
  cells: it still runs (recurrence ~236.8 yr versus 234.3 yr at 1 km, see
  `bp5/README.md`) but is under-resolved, so use 1 km or finer for production
  accuracy. (A strongly-coupled two-fault test did stall at 2 km when stepping
  through its first event, but that reflects that more aggressive coupled
  nucleation, not the single BP5 fault.)
- Integrator: `5dp` (the default) when matching HBI or minimizing phase error,
  `-ts_rk_type 3bs` for production (roughly 24 to 41 percent fewer matvecs, more
  at finer resolution, at a recurrence-interval error well under 0.02 yr). ODE
  `-rtol 1e-4`.
- Backend: HACApK (`-use_hmatrix 3`) is the robust default (trivial assembly,
  best compression at matched accuracy, competitive-to-best matvec, flat low
  memory), strongest at the np = 16 to 24 sweet spot. HTOOL and hmmvp have
  headroom at higher rank counts on cleanly-scaling hosts. Run hmmvp under MPI.
- H-matrix tolerance at the matched ~1e-6 accuracy band: HACApK `ztol` 1e-1,
  HTOOL `epsilon` 3e-5, hmmvp `tol` 1e-7. These sit about a decade below the
  `rtol 1e-4` integration floor that actually limits the cycle accuracy, so the
  dynamics are unperturbed; looser is possible once `nsteps` is confirmed flat.
  These are now the code defaults for HACApK and HTOOL (see below); hmmvp's
  default is left at 1e-5, which is still below the integration floor.

All of the above is specific to the tested builds and BP5 geometry; other
kernels, sizes, or library versions may shift the picture, so measure rather
than extrapolate.

---

## Default change (2026-06): H-matrix tolerance defaults to the matched band

The HACApK and HTOOL tolerance defaults in `petsc_interact.c` were loosened to the
matched ~1e-6 accuracy band:

```
HACApK  -hacapk_ztol        1e-4  ->  1e-1     # set in set_hacapk_defaults_and_options
HTOOL   -mat_htool_epsilon  1e-6  ->  3e-5     # set in set_htools_defaults_and_options
```

The earlier values were tighter than the cycle needs. For the smooth Okada
kernel HACApK's `ztol` is conservative, so 1e-1 already reaches a forward error
~2.2e-7 (well below the `rtol 1e-4` cycle floor) while compressing far more than
1e-4, which is effectively near-dense for this kernel. HTOOL `epsilon` 3e-5 gives
~6.6e-7 and, unlike 1e-6, does not land HTOOL at its least-compressed,
slowest-matvec point at large N. Both are overridable on the command line as
before. See `scaling_tests/rsf_solve_compression.md` for the compression-versus-N
evidence and `compress_interaction_matrix.md` for the tolerance-to-error mapping.

---

## Default change (2026-06): HTOOL compressor → `sympartialACA`

PETSc's `MATHTOOL` defaults its low-rank compressor to **SVD**, which is the most
accurate but by far the most expensive to assemble. For `rsf_solve` the accuracy
of the compressor is irrelevant in practice (see below, the error floor is set
by the ODE tolerance, not the H-matrix), so `rsf_solve` now sets

```
-mat_htool_compressor sympartialACA      # set in rsf_solve.c when -use_hmatrix 1
```

unless the user supplies `-mat_htool_compressor` on the command line. Override
with `SVD`, `fullACA`, or `partialACA` as desired. This change cut HTOOL assembly
from ~48 s to ~8 s on the BP5 4000-cell test, with no measurable accuracy loss.

---

## Benchmark: SEAS BP5-QD, 1 km / 4000 cells, single core

Setup: `bp5/` 1 km inputs, run to just past the first natural event (~234 yr),
`-rtol 1e-4`. "matvec" is the mean wallclock per `MatMult` from `-log_view`;
"setup" is total minus time-stepping (dominated by H-matrix assembly); accuracy
is measured against the dense solution.

| backend (tolerance) | recurrence | Δ vs dense | max Δtrace | memory | setup | step | matvec | **total** |
|----|----|----|----|----|----|----|----|----|
| dense | 234.294 yr | n/a | n/a | 128 MB | 19.5 s | 199.5 s | 20.14 ms | **219 s** |
| **HACApK** ztol 1e-4 | 234.293 | −0.0008 | 0.0036 | 25.2 MB (20%) | 6.5 s | 19.8 s | **1.75 ms** | **26 s** |
| HTOOL SVD eps 1e-4 | 234.294 | +0.0002 | 0.0041 | 18.9 MB (15%) | 48.3 s | 27.2 s | 2.46 ms | 75 s |
| HTOOL ACA eps 1e-4 | 234.285 | −0.009 | 0.0050 | 23.5 MB (18%) | 8.0 s | 32.0 s | 2.94 ms | 40 s |
| HTOOL SVD eps 1e-6 | 234.294 | −0.0002 | 0.0039 | 28.9 MB (23%) | 49.1 s | 33.8 s | 3.16 ms | 83 s |

(The HACApK row here is at `ztol` 1e-4, which is now known to be near-dense for
this kernel and is no longer the default. At the recommended `ztol` 1e-1 it
compresses far more: roughly 10x at 1 km and 27x at 0.5 km versus the ~5x of
the 1e-4 setting, with a faster matvec, at the same cycle accuracy. See
`scaling_tests/rsf_solve_compression.md`. The relative ranking of the backends
in this single-core table is unaffected; only HACApK's memory and matvec improve
at the recommended setting.)

### Findings

1. **Accuracy is a wash.** Every H-matrix variant reproduces the dense
   recurrence to <0.01 yr (~3e-6 relative) and the entire `log10 max|V|` trace to
   <=0.005 log units. Tightening HTOOL from `epsilon`=1e-4 to 1e-6 does *not*
   improve accuracy, the error floor is set by the ODE `rtol` (1e-4), not the
   H-matrix tolerance. So at `ztol`/`epsilon` <= 1e-4 both backends are
   effectively exact for this problem; there is no accuracy reason to prefer one.

2. **HACApK is fastest on a single core**, by a wide margin: fastest assembly
   (6.5 s) *and* fastest matvec (1.75 ms, ~1.4-1.8x quicker than any HTOOL
   variant and ~11x faster than the dense matvec). Since the matvec is paid on
   every one of ~9,700 RK-stage applications, the per-matvec edge dominates the
   total wallclock (26 s vs 40-83 s).

3. **HTOOL's SVD default is the trap.** SVD assembly (~48 s) is the single
   biggest cost in the HTOOL runs; `sympartialACA` removes it at no accuracy
   cost (hence the default change). Even so, HTOOL's matvec stays slower than
   HACApK's here, so its total remains above HACApK on one core.

4. **Memory**: at matched tolerance HTOOL compresses slightly more (15-18% of
   dense vs HACApK's 20%), but the gap is small and is outweighed by HACApK's
   faster matvec in wallclock terms.

---

## Why does this differ from the `compress_interaction_matrix` microbenchmark?

The standalone `compress_interaction_matrix` tests have indicated HTOOL
outperforming HACApK, which is the opposite of the BP5 result above. The
benchmarks are not measuring the same thing, and several factors plausibly flip
the ranking. Candidate explanations, roughly in order of suspected importance:

1. **Compressor default.** The microbenchmark may exercise HTOOL with ACA (or
   weight assembly differently), whereas the earlier BP5 runs used HTOOL's SVD
   default. SVD assembly alone is ~6x the ACA cost; if assembly dominates the
   microbenchmark, SVD-vs-ACA explains a large swing.
2. **What is timed.** The microbenchmark times assembly plus a handful of
   matvecs; `rsf_solve` times ~10^4 matvecs inside an adaptive ODE integrator.
   The two reward different things, assembly throughput vs steady-state matvec
   latency. HACApK's shell matvec is a direct call into its optimized routine;
   HTOOL's goes through the PETSc `MATHTOOL` layer.
3. **Admissibility (`eta`) and geometry.** HTOOL `eta`=100 is a loose
   admissibility condition tuned for general 3-D point clouds; the BP5 fault is a
   thin planar rectangle (100x40 km, one cell thick), whose low-rank block
   structure may favor HACApK's clustering. The microbenchmark geometry may be
   more 3-D/isotropic.
4. **Core count.** Both results above are single core. HTOOL is built into
   PETSc's distributed `Mat` layer and may scale differently from HACApK's
   lattice-H decomposition; the microbenchmark may have run multi-core. **This is
   the most important open question and is exactly what the scaling tests below
   target.**

These are hypotheses to test, not conclusions. The single-core verdict, 
"HACApK is the better default for serial / small-core BP5-type problems", should
not be read as universal.

---

## Parallel scaling (MPI)

Measured with `scripts/bench_hmatrix.sh` on a 24-core node (np=48 uses
hyperthreads), BP5 1 km / 4000 cells, short timing run (~357 steps, no event).
`step` is the matvec-bound time-stepping cost; `total` includes assembly.

> **Scope and caveats, read the numbers below as indicative, not definitive.**
> They come from a *single* machine, a *single* problem size (4000 cells), and a
> *single* timing run per point (no repetition, so expect ~10-20% run-to-run
> noise, especially at high rank counts where wallclocks are ~1-2 s). Absolute
> times and the exact crossover points depend on CPU, memory bandwidth, MPI
> implementation, and pinning, and will differ on other hardware. The np=48 point
> uses hyperthreads on a 24-core node, so 24 -> 48 is not a fair core-doubling.
> Treat the conclusions as "what happened in this configuration"; the robust,
> transferable observations are the *trends* (HACApK cheap assembly; HACApK
> matvec saturating as cells-per-rank drops; dense matvec scaling near-ideally),
> not the precise numbers.

**Total wallclock (s):**

| np | dense | HTOOL (ACA) | HACApK | fastest |
|---:|---:|---:|---:|:--|
| 1 | 30.92 | 42.21 | **10.91** | HACApK |
| 2 | 15.78 | 18.32 | **5.98** | HACApK |
| 4 | 10.02 | 10.60 | **3.76** | HACApK |
| 8 | 5.51 | 4.77 | **2.34** | HACApK |
| 16 | 2.63 | 3.10 | **1.59** | HACApK |
| 24 | 1.99 | 2.22 | **1.56** | HACApK |
| 48 | **1.23** | 1.49 | 1.64 | dense |

**Component speedup, np=1 → 48:** assembly, HTOOL 35x, dense 23x, HACApK 10x;
matvec/step, dense 28x, HTOOL 18x, **HACApK only 5.7x**; total, HTOOL 28x,
dense 25x, HACApK 6.6x.

**Memory (MB):** dense flat 122; HACApK flat 25.2; HTOOL grows with rank count
(17.9 at np=1 → 39.3 at np=48, per-rank clustering duplicates near-field blocks),
crossing above HACApK around np=24.

### Findings: the picture appears core-count dependent

These observations are specific to this machine and problem size (see caveats
above); the directional trends should transfer, the exact crossovers may not.

1. **HACApK was fastest at all physical-core counts here (np <= 24)**, fastest
   assembly at *every* rank count tested (3.8 s → 0.4 s) and the fastest serial
   matvec. On this 24-core box it looks like a sensible default for BP5-scale
   problems (1.3-4x ahead). On other hardware the margin will vary.

2. **HACApK's matvec appears not to strong-scale at this size.** Its step time
   flattens by ~np=16 (5.7x at 48 vs dense's near-ideal 28x). The likely cause is
   granularity: at 4000 cells / 48 ranks that is only ~83 cells/rank, where
   communication tends to dominate a hierarchical matvec while the dense matvec
   stays embarrassingly parallel (BLAS gemv + one reduction). In this run the
   ranking flips at np=48, dense (1.23 s) < HTOOL (1.49 s) < HACApK (1.64 s), 
   though those three are within ~30% of each other and of the same order as the
   run-to-run noise, so read it as "HACApK has lost its lead by 48," not a precise
   ordering.

   This should *not* be read as a limitation of HACApK. The lattice-H HACApK is
   designed explicitly for large-scale distributed-memory machines (Ozawa et al.,
   2021; Ida et al.), a regime this 4000-cell single-node test does not exercise, 
   we are simply starving it of per-rank work. Other or newer HACApK
   builds/configurations may well scale considerably better at high np, and the
   saturation we see is most plausibly a small-problem artifact rather than
   anything inherent. The same caution applies to HTOOL and to PETSc's
   integration layer: all three are capable libraries kindly shared by their
   authors, and these numbers reflect one build, one tuning, one problem.

3. **HTOOL showed the best scaling overall** (28x total) but from the most
   expensive assembly; it only caught HACApK at the highest rank counts, and its
   memory advantage eroded there too. ACA (now the default) is what makes this
   viable, SVD would inflate the already-dominant assembly.

4. **Any crossover is expected to move with problem size.** HACApK likely
   saturates here because the per-rank work is small; a larger model (0.5 km /
   16k cells, ~4x the cells-per-rank) was predicted to push the saturation point, 
   and the dense/HTOOL crossover, to higher rank counts. This has since been
   measured and the prediction held (see "Second resolution" below). The np=48
   point also crosses into hyperthreading on this 24-core node, which tends to
   penalize bandwidth-bound H-matrix matvecs (HACApK step time *rose* 24 → 48 at
   1 km).

### On the `compress_interaction_matrix` discrepancy: a plausible partial explanation

The scaling data is *consistent with* two of the earlier hypotheses, without
proving either. (i) HTOOL compresses better and scales better, so a microbenchmark
weighted toward **memory/compression or run at high core count** would plausibly
favor HTOOL, matching that test's ranking. (ii) The earthquake-cycle run instead
stresses **serial-to-moderate-core matvec latency**, where HACApK did better here.
The two benchmarks may simply probe different regimes rather than disagreeing.
This is a working explanation, not a settled one, confirming it would need the
microbenchmark re-run with matched compressor, tolerance, and core count.

### Second resolution: 0.5 km / 16000 cells (measured)

Same machine and protocol, finer grid (4x the cells, ~333 cells/rank at np=48
instead of ~83). Same caveats as above, one machine, one run per point, np=48 is
hyperthreaded. Total wallclock (s):

| np | dense | HTOOL (ACA) | HACApK | fastest |
|---:|---:|---:|---:|:--|
| 1 | 786 | 2428 | **120** | HACApK |
| 2 | 390 | 1269 | **62** | HACApK |
| 4 | 190 | 771 | **36** | HACApK |
| 8 | 113 | 160 | **21** | HACApK |
| 16 | 75.6 | 47.8 | **12.6** | HACApK |
| 24 | 57.1 | 26.3 | **9.8** | HACApK |
| 48 | 47.8 | 14.3 | **9.3** | HACApK |

Observations (specific to this configuration):

1. **The predicted shift happened.** At 1 km, HACApK lost the lead at np=48; at
   0.5 km it is fastest at *every* rank count tested. Its total speedup improved
   to ~13x (from ~6.6x at 1 km), consistent with the small-cells-per-rank
   saturation being the 1 km culprit. (It still flattens 24 → 48, but that step is
   hyperthreaded here.)

2. **HTOOL assembly becomes the dominant cost at this size.** Single-core ACA
   assembly was ~2299 s (~38 min) at 16000 cells, vs ~21 s for HACApK and ~263 s
   for dense, roughly a 73x jump for 4x cells (~quadratic), far worse than
   HACApK's ~5.7x. It parallelizes very strongly (down to ~10 s at np=48), but for
   serial or few-core use at this size HTOOL assembly is impractical. This may be
   tunable (admissibility `eta`, clustering) and could differ with another PETSc
   or Htool build; it should not be taken as a fixed property of the library.

3. **But HTOOL has the cheapest matvec at high rank count** (~5.0 ms/step at np=48
   vs HACApK's ~10.1 ms/step). The benchmark's short ~755-step run is
   assembly-dominated, which favors HACApK; for a *long* production run the matvec
   dominates and HTOOL's per-step edge amortizes its assembly. A rough linear
   extrapolation puts the run-length crossover near ~1700 steps at np=48 (~3700 at
   np=24), i.e. a full multi-event BP5 cycle (~10^4 steps) could plausibly favor
   HTOOL at high core count, if one accepts the large one-time assembly. This is an
   extrapolation from a short run, not a measured production comparison, and
   assumes a roughly constant per-step matvec cost.

4. **Memory.** HACApK flat at 167 MB (8.5% of the 1953 MB dense); HTOOL best at
   low np (102 MB, ~5%) but growing to ~204 MB by np=48 (per-rank duplication);
   dense's ~2 GB footprint is itself a constraint and its matvec appears
   bandwidth-saturated by np≈24.

### Multi-host update (six hosts): the saturation is real and host-dependent

The single-node 0.5 km numbers above have since been reproduced on six shared-memory
hosts (committed scaling CSVs and summary plots). They confirm the picture and sharpen
the caveats, so the single-node claim "HACApK fastest at *every* rank count" should be
read as a property of *that* node rather than a general result:

- **np ≈ 16-24 is the host-robust sweet spot.** At np=24, `hacapk` was the fastest
  matvec on *every* host tested, the one ordering that holds everywhere.
- **np=48 is host-dependent.** On cleanly-scaling nodes `htool` and `hmmvp` keep scaling
  and **overtake `hacapk` by np=48** (`hacapk`, the leanest and most bandwidth-bound
  matvec, plateaus near np≈16-24); bandwidth-limited nodes instead *degrade* past np≈24;
  and one node's np=48 collapsed across all backends, reproducibly, an oversubscription
  artifact (fewer usable cores than ranks), not a backend property. The node in the table
  above, where `hacapk` stayed fastest to np=48, sits at the favorable end of that spread.
- The transferable statements are therefore: `hacapk` is the best low-to-mid-rank choice
  (and the cheapest serial matvec), while `htool`/`hmmvp` are the backends with headroom
  left at high rank. The cross-code version of the same conclusion (vs HBI's lattice) is
  in `hbi_tests/rsf_solve_vs_hbi_scaling.md`.

### Practical default (provisional)

On the evidence so far, keeping HACApK (`-use_hmatrix 3`) is a reasonable default:
it had the cheapest assembly at every size and rank count tested and won on total
wallclock for short-to-moderate runs, comfortably so for parameter sweeps,
ensembles, and single-node jobs. Three caveats temper that, all pointing the same
way, measure for your actual workload:

- **Run length matters as much as core count.** The benchmarks are short
  (~755 steps) and therefore assembly-weighted. For long production cycles
  (~10^4 steps) the matvec dominates, and at high rank counts HTOOL's cheaper
  matvec may amortize its (large) assembly and overtake HACApK, extrapolated,
  not yet measured.
- **Problem size moved the crossover.** At 0.5 km HACApK led to np=48, whereas at
  1 km it did not, so conclusions at one resolution don't transfer.
- **HTOOL assembly scales poorly with N on few cores** here (~38 min at 16k cells,
  np=1), which may be tunable and build-dependent.

So treat HACApK as the starting point, but for a single large model on a many-core
or multi-node allocation, especially a long one, benchmark dense and HTOOL(ACA)
at the target size, rank count, *and* representative step count with
`scripts/bench_hmatrix.sh` rather than assuming. It is also worth trying the
MPI+OpenMP hybrid: HACApK's matvec is OpenMP-threaded, and the rank counts where it
flattened are exactly where threads-per-rank may help.

---

## Integrator order: an orthogonal speedup lever

Everything above concerns the *matvec* (which operator backend, how many ranks). A
second, independent lever is the **Runge-Kutta order**, set with `-ts_rk_type`. The
default is Dormand-Prince `5dp` (5(4), 6 stage matvecs per accepted step), chosen partly
so its order matches HBI's Cash-Karp for cross-code comparison. A lower-order embedded
pair does fewer stage matvecs per step: `3bs` (Bogacki-Shampine RK3(2)) costs 3, and
although it takes more (and more frequently rejected) steps, it still nets fewer matvecs.

Judged by the physically meaningful metric, the event **recurrence interval**, which is
far better converged than the absolute event phase (the phase drifts by ~`rtol`
run-to-run; the interval does not), `-ts_rk_type 3bs` reproduced the `5dp` recurrence to
within ~0.005-0.02 yr (<~0.01% of the ~230 yr interval) at every resolution tried (BP5
2 km dense; 1 km dense and HACApK; 0.5 km HACApK), while reducing matvecs by very roughly
**24% (2 km), 35% (1 km), 41% (0.5 km)**, the saving grows with resolution because finer
meshes take more steps, so the cheaper-per-step method compounds. These are single-host
serial runs over one event sequence, so confirm for your own case.

Practical guidance (also captured in comments in `src/rsf_solve.c`):

- Keep `5dp` (the default) for anything compared against HBI, and where you want the most
  converged recurrence.
- Consider `-ts_rk_type 3bs` for production where a recurrence error of order ~0.02 yr is
  acceptable; it stacks with the step controller `-ts_adapt_type dsp` (which trims
  rejected steps at equal accuracy).
- Do **not** loosen `-rtol` below ~1e-4: at 1e-3 the rejection rate climbs enough that the
  run is both less accurate and not faster, and over long runs can make very slow progress.
- `2a` (RK2(1)) needed far more matvecs here; `5bs` behaved poorly on the stiff coseismic
  phase; higher order (`8vr`) only repaid its per-step cost for very tight *absolute* event
  times, not long multi-cycle runs.

---

## Quick reference: running the benchmark

```bash
# HACApK (recommended serial default; ztol defaults to 1e-1)
mpirun -np 1 bin/rsf_solve -use_hmatrix 3   <bp5 flags>

# HTOOL (now defaults to sympartialACA; relax epsilon for speed)
mpirun -np 1 bin/rsf_solve -use_hmatrix 1 -mat_htool_epsilon 1e-4   <bp5 flags>

# dense ground truth
mpirun -np 1 bin/rsf_solve -use_hmatrix 0   <bp5 flags>
```

Add `-log_view` to get per-event timings (`MatMult`, `TSStep`).

Velocity-field snapshots (the `tmp_rsf/vel-*-gmt` files used for slip-evolution plots)
now default **off**; enable them with `-slip_line_dt_yr <interval>` (e.g. `1` for the fine
SEAS-style interseismic cadence, or `10`-`20` for a cycle-scale view). They are written
once per interval and are *not* purged between runs, so for long, fine-grid runs (e.g.
0.5 km / N=16000, ~1 MB per frame) they accumulate quickly and can fill the working disk;
clean `tmp_rsf` between runs, and note that a full disk makes the snapshot and
`rsf_monitor.dat` writes fail, which can look like a solver stall.

---

## Task 1 (2026-07): SEAS-style event catalog, rupture-time field, slip budget

Three optional, default-off outputs that complement the existing `rsf_events.dat`
event tracker rather than replace it. They ride on the same slip-rate threshold
crossings that `-track_events` already detects (so `-track_events`, on by
default, is required), and add only the per-event quantities that
`rsf_events.dat` does not carry. All three are inert unless explicitly enabled,
so a run without the new flags is unchanged (verified: no new files, identical
behavior).

| flag | output | contents |
|----|----|----|
| `-rsf_catalog` | `rsf_catalog.dat` | one row per completed event: onset/arrest time, duration, ruptured-cell count and area, mean and max coseismic slip, mean and max static stress drop, peak slip rate, seismic moment `M0` and `Mw` |
| `-rsf_rupture_time` | `rsf_rupture_time.dat` | per-cell time of first crossing of the rupture-front threshold during event 1, referenced to the earliest crossing (the initiation time), with a ruptured flag; the Jiang et al. (2022) Fig 4 / Fig 7 diagnostic |
| `-slip_budget` | `rsf_slip_budget.dat` | area-integrated on-fault slip versus the plate-rate reference `vpl*t*area` at each stats interval, a long-term loading-consistency check |
| `-rupture_vth <v>` | | rupture-front threshold [m/s] for the catalog and rupture-time field; defaults to the event onset threshold `vel_event`. The benchmark uses a distinct, larger value (0.1 or 0.03 m/s), passed here |

Definitions, so the values are unambiguous: coseismic slip and stress drop are
per-cell differences between the onset snapshot and the arrest state, with the
drop signed so a stress decrease is positive; the mean is over ruptured cells
only (cells that exceeded `rupture_vth` at any accepted step during the event),
the max is over the same set; `M0 = G * sum(area_i * coseismic_slip_i)` over
ruptured cells, `Mw = (2/3)(log10 M0 - 9.05)`; peak slip rate is the maximum over
cells and over accepted steps within the event. These are threshold-dependent
and, for stress drop, domain-dependent, so they are reported relative to the
chosen thresholds rather than as resolution-independent quantities.

Implementation is confined to `rsf_output.c` (setup, per-step tracking in the
monitor, and the onset-snapshot / arrest-reduction in the event handler),
`rsf.h` (context and settings fields, prototypes), `rsf_init.c` (option
parsing), and `rsf_solve.c` (two calls). The per-cell tracking runs every
accepted step while an event is in progress, outside the change-triggered
monitor gate, so it follows the solver's own step density, which collapses
through the coseismic phase and resolves the rupture front. The rupture-time
field is gathered to rank 0 in geometry order through a row-layout `Vec` scatter,
so it is correct in parallel.

One associated fix went in with this work: the event tracker previously
initialized its slipping state from the scalar `vel_init`, which for a per-cell
initial-condition file (the BP5 nucleation patch is seeded at 3e-2 m/s while
`vel_init` is the background `vpl`) mislabeled the decaying seed as a spurious
onset then arrest. The slipping state is now initialized from the true initial
maximum slip rate over the fault, so the seeded first event is treated as an
event already in progress at t=0, which is what both `rsf_events.dat` and the
new catalog want. This changes only the labeling of the initial seeded
transient; the spontaneous recurrence times used for the HBI cross-check are
unaffected.

### Verification

On the BP5 2 km dense test to a few years (long enough for the seeded event 1 to
arrest), the standalone reducer unit test (`test_hmatrix/catalog_unit_test.c`,
compiled and run outside the interact build) reproduces the catalog arithmetic
(mean/max coseismic slip, signed stress drop, ruptured area, `M0`/`Mw`, peak slip
rate, and the initiation-referenced rupture time) against hand-computed values,
including the degenerate no-rupture case. The end-to-end run produces a coherent
outward-propagating rupture-time field (initiation at 0 s at the nucleation
patch, growing with distance), and the HACApK backend reproduces the catalog to
within the expected sub-percent compression-level differences (for example 307
versus 309 ruptured cells, `Mw` 7.486 versus 7.491 at `ztol` 1e-1). All of this
is specific to the tested 2 km configuration; the absolute event size at 2 km is
coarse and under-resolved, so the numbers are for validating the machinery, not
for benchmark comparison, which needs 1 km or finer.

---

## Task 2 (2026-07): dip-slip mode and the normal-stress path

`rsf_solve` was strike-slip only: the source slip and the shear interaction
matrix `Is` were fixed to the strike direction, so the normal-stress evolution
`-calc_sigma_dot` (the `In` matrix) was only ever "normal stress due to strike
slip", which is small on a planar fault and never exercised by the BP5 setup.
A thrust or normal fault is dip slip, and that is the case where slip genuinely
changes the normal stress, so a slip-direction option was added.

| flag | meaning |
|----|----|
| `-rsf_slip_mode <0\|1>` | slip direction the rate-and-state solve operates on: 0 strike (default, unchanged), 1 dip. Sets the source slip direction and makes the `Is` shear matrix resolve onto `t_dip`; `In` is unaffected |
| `-calc_sigma_dot` | evolve the normal stress (`dsigma/dt` from `In`). Was previously settable only in code; now a command-line boolean |
| `-rsf_sigma_file FILE` | optional per-cell initial normal stress `sigma0` [Pa], geometry order, one value per patch. Overrides the uniform `-sigma_init`. Companion to a dipping fault where `sigma0` varies with depth (item 6a) |

Mechanically the change is confined to the slip-direction argument now threaded
through `calc_petsc_Isn_matrices`: for the shear matrix (`mode` 0) the source
slip and the resolved stress component both follow `-rsf_slip_mode`, while the
normal matrix (`mode` 1) always resolves onto the fault normal. The backslip
loading (`sinc`) and the RHS need no separate change, since they are just `Is`
and `In` applied to the plate rate, so they inherit the dip direction
automatically. At the default `slip_mode = strike` every call reduces to the
original strike-slip path, so the strike-slip results are unchanged (verified:
BP5 2 km is byte-identical in `rsf_events.dat` and the event catalog).

### Thrust test and validation

`thrust/make_thrust.py` generates a surface-breaking, listric thrust in the
spirit of Ozawa et al. (2023): about 50 km along strike and 20 km down the
(curved) dip path, with the dip tapering from 30 degrees at the surface to 10
degrees at depth, following the `calc_quad_base_vecs` convention (for strike 0
the down-dip direction is `(cos t, 0, -sin t)` and the normal is
`(sin t, 0, cos t)`). It writes the geometry, `a`/`b`, initial conditions, and an
optional depth-dependent `sigma0` file, with a seeded interior nucleation patch
so that slip, and therefore a normal-stress change, actually occurs.

The geometry was checked numerically: the shallow row breaks the surface (top
edge at `x = 0`, `z = 0`), the rows tile the curve continuously (row-to-row edge
gaps at the millimetre level from rounding), the dip-path length is 20 km, the
dip runs 30 to 10 degrees, and the along-strike extent is 50 km centred on the
origin.

The normal-stress path was validated as an on/off contrast on the 2 km thrust
with `-rsf_slip_mode 1` (`thrust/test_thrust.sh`): with `-calc_sigma_dot` the
`In` matrix is built and the normal stress moves away from the uniform 58 MPa as
the seed slips (a spread of about 1.2 MPa on this coarse mesh, with both clamping
and unclamping, consistent with dip slip on a bending, surface-breaking fault);
without the flag the `In` matrix is not built and the normal stress stays pinned
at 58 MPa for the whole run. The `-rsf_sigma_file` path was checked separately by
confirming that a depth-dependent `sigma0` file is honoured at `t = 0`.

These results are specific to the tested configuration. At 1 km cells with the
Ozawa-like `D_c` the cohesive zone is only marginally resolved, so this is a
validation of the dip-slip and normal-stress machinery, not a converged
benchmark; quantitative comparison to a published thrust would need a finer mesh
and a matched parameter set, and any difference from another code may reflect
resolution or kernel choices rather than a defect in either code.
