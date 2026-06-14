# `rsf_solve` — H-matrix backends and performance notes

`rsf_solve` integrates a rate-and-state quasi-dynamic earthquake cycle on an
interact patch geometry. The dense elastic interaction operator is applied once
per Runge-Kutta stage, so the dominant cost at scale is the **matrix-vector
product** with the stress-interaction matrix. `rsf_solve` can apply that operator
through several hierarchical-matrix (H-matrix) backends, selected with
`-use_hmatrix`:

| `-use_hmatrix` | backend | parallelism | tolerance knob (default) | notes |
|:--:|----|----|----|----|
| 0 | dense | MPI | — | exact reference; O(N^2) memory and matvec |
| 1 | **HTOOL** (PETSc `MATHTOOL`) | MPI | `-mat_htool_epsilon` (1e-6) | also `-mat_htool_eta` (100), `-mat_htool_compressor` |
| 2 | H2OPUS | serial/GPU | `-h2opus_eta` | |
| 3 | **HACApK** (lattice H-matrix) | MPI | `-hacapk_ztol` (1e-4) | shell matvec into HACApK |
| 4 | hmmvp | OpenMP | `-hmmvp_tol` (1e-5), `-hmmvp_eta` (3) | |

This note focuses on the two production MPI backends, **HTOOL** and **HACApK**,
benchmarked against the **dense** ground truth.

---

## Default change (2026-06): HTOOL compressor → `sympartialACA`

PETSc's `MATHTOOL` defaults its low-rank compressor to **SVD**, which is the most
accurate but by far the most expensive to assemble. For `rsf_solve` the accuracy
of the compressor is irrelevant in practice (see below — the error floor is set
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
| dense | 234.294 yr | — | — | 128 MB | 19.5 s | 199.5 s | 20.14 ms | **219 s** |
| **HACApK** ztol 1e-4 | 234.293 | −0.0008 | 0.0036 | 25.2 MB (20%) | 6.5 s | 19.8 s | **1.75 ms** | **26 s** |
| HTOOL SVD eps 1e-4 | 234.294 | +0.0002 | 0.0041 | 18.9 MB (15%) | 48.3 s | 27.2 s | 2.46 ms | 75 s |
| HTOOL ACA eps 1e-4 | 234.285 | −0.009 | 0.0050 | 23.5 MB (18%) | 8.0 s | 32.0 s | 2.94 ms | 40 s |
| HTOOL SVD eps 1e-6 | 234.294 | −0.0002 | 0.0039 | 28.9 MB (23%) | 49.1 s | 33.8 s | 3.16 ms | 83 s |

### Findings

1. **Accuracy is a wash.** Every H-matrix variant reproduces the dense
   recurrence to <0.01 yr (~3e-6 relative) and the entire `log10 max|V|` trace to
   <=0.005 log units. Tightening HTOOL from `epsilon`=1e-4 to 1e-6 does *not*
   improve accuracy — the error floor is set by the ODE `rtol` (1e-4), not the
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
   The two reward different things — assembly throughput vs steady-state matvec
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

These are hypotheses to test, not conclusions. The single-core verdict —
"HACApK is the better default for serial / small-core BP5-type problems" — should
not be read as universal.

---

## Parallel scaling (MPI)

Measured with `scripts/bench_hmatrix.sh` on a 24-core node (np=48 uses
hyperthreads), BP5 1 km / 4000 cells, short timing run (~357 steps, no event).
`step` is the matvec-bound time-stepping cost; `total` includes assembly.

> **Scope and caveats — read the numbers below as indicative, not definitive.**
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

**Component speedup, np=1 → 48:** assembly — HTOOL 35x, dense 23x, HACApK 10x;
matvec/step — dense 28x, HTOOL 18x, **HACApK only 5.7x**; total — HTOOL 28x,
dense 25x, HACApK 6.6x.

**Memory (MB):** dense flat 122; HACApK flat 25.2; HTOOL grows with rank count
(17.9 at np=1 → 39.3 at np=48 — per-rank clustering duplicates near-field blocks),
crossing above HACApK around np=24.

### Findings — the picture appears core-count dependent

These observations are specific to this machine and problem size (see caveats
above); the directional trends should transfer, the exact crossovers may not.

1. **HACApK was fastest at all physical-core counts here (np <= 24)** — fastest
   assembly at *every* rank count tested (3.8 s → 0.4 s) and the fastest serial
   matvec. On this 24-core box it looks like a sensible default for BP5-scale
   problems (1.3-4x ahead). On other hardware the margin will vary.

2. **HACApK's matvec appears not to strong-scale at this size.** Its step time
   flattens by ~np=16 (5.7x at 48 vs dense's near-ideal 28x). The likely cause is
   granularity: at 4000 cells / 48 ranks that is only ~83 cells/rank, where
   communication tends to dominate a hierarchical matvec while the dense matvec
   stays embarrassingly parallel (BLAS gemv + one reduction). In this run the
   ranking flips at np=48 — dense (1.23 s) < HTOOL (1.49 s) < HACApK (1.64 s) —
   though those three are within ~30% of each other and of the same order as the
   run-to-run noise, so read it as "HACApK has lost its lead by 48," not a precise
   ordering.

3. **HTOOL showed the best scaling overall** (28x total) but from the most
   expensive assembly; it only caught HACApK at the highest rank counts, and its
   memory advantage eroded there too. ACA (now the default) is what makes this
   viable — SVD would inflate the already-dominant assembly.

4. **Any crossover is expected to move with problem size.** HACApK likely
   saturates here because the per-rank work is small; a larger model (e.g. 0.5 km
   / 16k cells, ~4x the cells-per-rank) should push the saturation point — and the
   dense/HTOOL crossover — to higher rank counts, but this has not yet been
   measured. The np=48 point also crosses into hyperthreading, which tends to
   penalize bandwidth-bound H-matrix matvecs (HACApK step time *rose* 24 → 48).

### On the `compress_interaction_matrix` discrepancy — a plausible partial explanation

The scaling data is *consistent with* two of the earlier hypotheses, without
proving either. (i) HTOOL compresses better and scales better, so a microbenchmark
weighted toward **memory/compression or run at high core count** would plausibly
favor HTOOL — matching that test's ranking. (ii) The earthquake-cycle run instead
stresses **serial-to-moderate-core matvec latency**, where HACApK did better here.
The two benchmarks may simply probe different regimes rather than disagreeing.
This is a working explanation, not a settled one — confirming it would need the
microbenchmark re-run with matched compressor, tolerance, and core count.

### Practical default (provisional)

On the evidence so far, keeping HACApK (`-use_hmatrix 3`) is reasonable for
typical use — parameter sweeps, ensembles, and single-node jobs at BP5-like
sizes. For a single large model on a many-core or multi-node allocation, don't
assume any backend wins: benchmark dense and HTOOL(ACA) at the target rank count
and resolution first with `scripts/bench_hmatrix.sh`, since the crossover is
expected to move with cells-per-rank and is hardware-dependent. The default is a
starting point, not a recommendation to skip measuring on your own system.

---

## Quick reference — running the benchmark

```bash
# HACApK (recommended serial default)
mpirun -np 1 bin/rsf_solve -use_hmatrix 3 -hacapk_ztol 1e-4   <bp5 flags>

# HTOOL (now defaults to sympartialACA; relax epsilon for speed)
mpirun -np 1 bin/rsf_solve -use_hmatrix 1 -mat_htool_epsilon 1e-4   <bp5 flags>

# dense ground truth
mpirun -np 1 bin/rsf_solve -use_hmatrix 0   <bp5 flags>
```

Add `-log_view` to get per-event timings (`MatMult`, `TSStep`).
