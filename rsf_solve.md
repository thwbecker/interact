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
| `src/rsf_engine.c` | physics | `vel_from_rsf` (invert the regularized friction law for slip rate); `rsf_state_rate` (the state evolution rate for all five laws, with optional analytic partial derivatives; the single place a new law is added, shared by both integrator paths); `rsf_compute_vel_and_stressing` (fill the slip-rate vector from the state and apply the `Is` and, if enabled, `In` matvecs; shared by both paths); `rsf_ODE_RHSFunction` (the monolithic quasi-dynamic RHS used by the default explicit integrator, where the backslip `sinc` loading enters); `rsf_domain_check` (per-step validity gate) |
| `src/rsf_imex.c` | IMEX variant (`-imex`) | `rsf_IMEX_RHSFunction` (explicit part: matvec stressing rates and the non-stiff local terms); `rsf_IMEX_IFunction` (implicit part: state evolution and the radiation-damping correction, purely local, no matvec); `rsf_IMEX_IJacobian` (exact analytic per-cell block Jacobian); `rsf_IMEX_setup` (switch a TS to ARKIMEX with these functions, create the block-diagonal Jacobian matrix, set the `imex_` options prefix and stage-solver defaults); `rsf_IMEX_set_stage_solver` (assert the exact block stage solve after the options pass) |
| `src/rsf_output.c` | output | `rsf_init_monitor_and_event` and `rsf_finalize_monitor_and_event` (set up and tear down the TS monitor, event detector, gather, and output files); `rsf_TS_Monitor` (state-change logging to `rsf_monitor.dat`, per-frame field output and its `rsf_vel.times` index, optional per-group GMT slip-rate fields); `rsf_post_event`; and the field-output group helpers (`rsf_build_groups`, `rsf_group_coords`, `rsf_write_group_geometry`, `rsf_free_groups`, `rsf_dcmp`) |
| `src/includes/rsf.h` | header | the output/grid/settings structs (`rsf_out_ctx`, `rsf_group_grid`, `rsf_solve_settings`), the rsf prototypes, and the single shared `rsf_par_static` (the geometry pointer the domain check reads) |

The build links all four objects via `RSF_SOLVE_OBJS` in the makefile:

    RSF_SOLVE_OBJS = $(ODIR)/rsf_solve.o $(ODIR)/rsf_engine.o $(ODIR)/rsf_imex.o $(ODIR)/rsf_init.o $(ODIR)/rsf_output.o

(`rsf_imex.o` was added with the IMEX option, 2026-07, see the IMEX section
below.  The state-law formulas that used to be inlined in the RHS moved to the
shared `rsf_state_rate` at the same time; the default explicit path calls the
same helper.  The refactor was verified to be bit-identical for the aging law
over a 0.2 yr BP5 test; the PRZ trajectory was identical for 297 monitor rows
before diverging by one unit in the ninth digit of dt, consistent with a
floating-point contraction difference from the function extraction rather than
a formula change.)

One build caveat worth flagging: this split also moved three `Vec` fields out of
`struct med` and added the `medium->rsf` pointer, which shifts the offsets of
later `struct med` members. After pulling that change, do a clean rebuild
(`rm objects/*.o objects/*.a`, or `make clean`) so every object sees the same
`structures.h` layout. A partial rebuild that reuses a stale `petsc_interact.o`
(which writes `medium->rs/re/rn`) against freshly built rsf objects will read the
ownership range at the wrong offset and segfault in the monitor on the first
output step. This is a stale-object hazard, not a code defect: a clean build
fixes it.

## Command-line options and the `-h` help

`rsf_solve -h` prints a grouped summary of every option the code parses, with
defaults, plus a short list of the PETSc options that matter here and two worked
examples. It is the quick reference; run it whenever the flag set feels
unwieldy. PETSc's own `-help` still works and dumps the full registered-option
list, which is much longer. Options can also be placed in `petsc_settings.yaml`,
which is read automatically at startup, so anything on the command line can live
there instead.

The groups, in the order `-h` prints them, are:

- geometry and input files: `-geom_file`, `-rsf_file`, the optional per-cell
  files `-rsf_ic_file`, `-rsf_dc_file`, `-rsf_sigma_file`, plus `-full_space`
  and `-tv`.
- elastic and rate-and-state parameters: `-shear_modulus`, `-s_wave_speed`
  (these two set the radiation damping coefficient `eta = G/(2 c_s)`; note
  that a SMALLER wave speed means MORE damping, so `c_s = 0` is not a way to
  switch damping off), `-rd_fac` (scales the damping coefficient: 1 is the
  standard quasi-dynamic value and the default, 0 switches damping off, the
  `c_s -> infinity` limit, values above 1 give enhanced damping; caution:
  with damping off the quasi-dynamic coseismic phase is unregularized and
  the slip rate is unbounded at instabilities, so integrators fail there,
  and the intended use is stable sliding, nucleation, or slow-slip studies,
  verified on the slider: creep runs cleanly with `-rd_fac 0`, stick-slip
  runs stall at the first event),
  `-f0`, `-dc`, `-v0`, `-vpl`, and `-rsf_slip_mode` (0 strike, 1 dip).
- initial conditions: `-sigma_init`, `-tau_init`, `-vel_init`, `-rand_amp`.
- normal-stress evolution and limiter: `-calc_sigma_dot`, `-limit_sigma`,
  `-min_sigma`, `-max_sigma`.
- state evolution: `-state_law` (1 aging, the default; 2 slip; 3 PRZ; 4 Sato;
  5 Kato and Tullis; see the law notes at the end of this file) and
  `-vmin_state`.  The gated-law parameters (sato_beta, kt_vc) currently have
  compile-time defaults only (1e-2 and 1e-2 v0) and are not yet runtime
  options.
- time stepping: `-rtol`, `-atol_slip`, `-dt_init`, `-dt_max`, `-stop_time_yr`,
  `-imex` to switch from the default explicit Runge-Kutta to the IMEX
  (ARKIMEX) formulation (see the IMEX section below; its stage solvers are
  tunable under the `imex_` options prefix, e.g. `-imex_ksp_type`), and
  `-domain_check_max_reject` (abort after that many CONSECUTIVE
  out-of-domain trial steps, default 1000, `<= 0` disables; such rejection
  storms bypass the PETSc rejection limits and previously could grind a run
  indefinitely).  One solver caution from the slider benchmark, specific to
  the configurations tested there: the 5dp method showed a failure band at
  loose tolerances (rtol 1e-4 to 1e-5), dying at nucleation in exactly such
  storms, and generally runs at a much higher rejection rate than 3bs, which
  completed every benchmark run; where robustness matters, 3bs has been the
  safer choice, consistent with experience on BP5.
- monitor and event detection: `-print_interval_yr`, `-dt_monitor_yr`,
  `-rdx_monitor`, `-adx_monitor`, `-monitor_tmin_yr`, `-track_events`,
  `-vel_event`, `-vel_event_hyst`, `-event_tmin_yr`.
- optional outputs (all default off): `-rsf_catalog`, `-rsf_rupture_time`,
  `-slip_budget`, `-rupture_vth`, `-field_step_interval`, `-field_tmin_yr`,
  `-slip_line_dt_yr`.
- interaction-matrix backend: `-use_hmatrix` (0 dense, 1 HTOOL, 2 H2OPUS,
  3 HACApK, 4 hmmvp, 5 BigWham) and the per-backend tolerance and admissibility
  knobs (`-hacapk_ztol`, `-hmmvp_tol`, and so on).

Relevant PETSc options include `-ts_rk_type`, `-ts_adapt_type`, `-ts_max_steps`,
`-ts_monitor`, `-log_view`, `-mat_htool_compressor`, and `-options_file`. The
help text is generated by `rsf_print_help` in `src/rsf_init.c`, kept next to the
option parser so the two stay in step; when an option is added to the parser, its
line should be added to the help.

### State evolution law (`-state_law`)

The state is carried in the regularized variable `psi`, with
`v = 2 v0 exp(-psi/a) sinh(tau/(sigma a))`. Five evolution laws are available,
all implemented in one place, `rsf_state_rate` in `src/rsf_engine.c` (together
with their analytic partial derivatives, which the IMEX path uses), so both
integrator paths see identical formulas. Throughout, `Omega = |v| theta/dc`
with `theta = (dc/v0) exp((psi-f0)/b)`, so `Omega = 1` is the classical steady
state.

- `-state_law 1`, the aging (Dieterich) law, and the default. In `theta` form
  `d theta/dt = 1 - Omega`; in `psi`:
  `d psi/dt = b/dc ( v0 exp((f0-psi)/b) - |v| )`.
  Heals at stationary contact; the state rate is bounded by `b|v|/dc` above
  steady state.
- `-state_law 2`, the slip (Ruina) law. In `theta` form
  `d theta/dt = -Omega ln Omega`; in `psi`:
  `d psi/dt = -(|v|/dc) (psi - psi_ss)`.
  Heals only while slipping; the rate grows linearly in the `psi` lag.
- `-state_law 3`, the PRZ law (Perrin, Rice and Zheng, 1995), normalized so
  that `Omega_ss = 1`:
  `d theta/dt = (1 - Omega^2)/2`; in `psi`:
  `d psi/dt = b/(2 dc) ( v0 exp((f0-psi)/b) - v^2/v0 exp((psi-f0)/b) )`.
  Heals at rest, but at HALF the aging rate under this normalization (see the
  normalization note at the end of this file for the trade-off).  Above
  steady state the rate grows like `Omega^2`, which is what makes PRZ stiff
  for explicit integrators at rupture fronts (see the IMEX section) and slow
  to run at fine resolution.
- `-state_law 4`, the Sato-type gated composite law:
  `d theta/dt = exp(-Omega/beta) - Omega ln Omega`, the slip law plus aging
  healing shut off by a gate keyed on `Omega` (`beta = 1e-2`, compile-time).
- `-state_law 5`, the Kato and Tullis (2001) composite law:
  `d theta/dt = exp(-|v|/Vc) - Omega ln Omega`, the same construction with the
  gate keyed on `|v|` (`Vc = v0/100`, compile-time); its low-speed steady
  state sits at `Omega_ss = 1.763` rather than 1, see the law notes at the
  end of this file.

Laws 1 to 4 share the same steady state, `psi_ss = f0 - b ln(|v|/v0)`,
corresponding to the usual `f_ss = f0 + (a-b) ln(|v|/v0)`; the laws differ in
how `psi` relaxes toward it. `|v|` (not `v`) is used throughout, as in HBI,
which matters for the sinh regularization that permits `v < 0`. In the slip,
Sato, and Kato and Tullis laws `|v|` is floored at `-vmin_state` (default
1e-16 m/s) so that `ln(|v|/v0)` stays finite as `v` approaches zero, where
`d psi/dt` tends to zero in any case; `b = 0` cells return a zero state rate.

In a test on the SEAS BP4-QD setup at 2 km, run in this sandbox with an otherwise
identical configuration, the slip law gave shorter recurrence and more frequent
events than the aging law. This is the expected direction, since the slip law heals
less between events, but the comparison was a single coarse configuration and
should not be read as a general or converged result.

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

### Event detection and spatial resolution

Events are bracketed by crossings of the global maximum slip rate through
`-vel_event` (arrest at `vel_event * -vel_event_hyst`), and a bracket is only
written to `rsf_catalog.dat` if at least one cell reached `-rupture_vth`. That
last filter matters more than it first appears: the maximum slip rate routinely
dips back below the arrest threshold and re-crosses it while a rupture is still
nucleating, but as long as those excursions stay below `rupture_vth` they contain
no ruptured cell and are correctly discarded.

Whether that holds depends on spatial resolution, and the slip law is the more
demanding case. In tests on the SEAS BP4-QD setup in this sandbox (single
configuration, not a general result), with the process zone
`Lb = G dc / (b sigma)` at about 2 km:

- at `dx` = 2 km, so roughly one cell per `Lb`, individual cells reached seismic
  slip rates and ruptured independently. These are genuine ruptures in the
  discretized model, so they are logged as events, and the catalog filled with
  one-cell rows. No event-detection rule can filter these without also discarding
  legitimate small events; the problem is the mesh, not the logging.
- at `dx` = 1 km, about two cells per `Lb`, a single rupture appeared as several
  catalog rows: the nucleation front advanced cell by cell, with individual cells
  exceeding `rupture_vth` in turn and the maximum slip rate dipping in between
  (the observed arrest-to-onset gaps were roughly 20 to 500 s).
- at `dx` = 0.5 km, about four cells per `Lb`, the nucleation crossings persisted
  but peaked around 1e-3 m/s, far below a `rupture_vth` of 0.1 m/s, so no cell
  ruptured in them and the existing filter discarded them. The catalog was clean.

The practical guidance is therefore to check cells per `Lb` before trusting a
catalog, and to be aware that a mesh that looks adequate under the aging law may
not be under the slip law, which localizes more strongly. Fragmented or one-cell
events in the catalog should be read as a resolution diagnostic rather than as
something to be filtered away in post-processing.

---

## Task 3 (2026-07): IMEX time integration (`-imex`)

`-imex` switches the time integrator from the default explicit embedded
Runge-Kutta (TSRK) to a PETSc additive Runge-Kutta IMEX method (TSARKIMEX,
tableau `l2`: L-stable, second order, stage order 2; `-ts_arkimex_type`
overrides).  The tableau default matters for CORRECTNESS here, not just
efficiency: the nominally more accurate ARKIMEX3 suffers stiff stage-order
reduction on these stage problems, its embedded error estimator
underestimates the state-variable error through fast transitions, and the
controller then accepts steps whose true error is far above the requested
tolerance.  On the single-patch slider benchmark at rtol 1e-6 this cost a
factor of about 50 in event-time accuracy relative to `l2` (about 8e-2
versus 2e-3 yr RMS, uniformly across the evolution laws), ARKIMEX4 was worse
still, and ars443 and bpr3 numerically quenched the stick-slip cycle
outright (seven events became zero); see slider/README.md.  A cheaper
stage-order-2 sub-family (`2c`/`2d`/`2e`, about 3.4x fewer steps than `l2`
at equal tolerance) passed every slider test, including event-time accuracy,
but was REJECTED as a default after a fault-scale check: on the BP5 0.5 km
PRZ problem at rtol 1e-4, `2d` produced numerous spurious small events
absent from the explicit and `l2` solutions, i.e. it altered the event
statistics at the operating tolerance by exciting marginally stable
partial-rupture modes that the single-patch benchmark cannot represent.
`l2` therefore remains the default (slower, statistics-faithful); the cheap
tableaus stay available via `-ts_arkimex_type` but should only be used with
a fault-scale event-statistics verification (for example against an
explicit or `l2` run, or at demonstrably converged tolerance).  The general
lesson is recorded in slider/README.md: the slider validates phase accuracy
and robustness, not the stability of marginal spatial modes, so tableau or
solver changes need one fault-scale statistics check before adoption.  The price of `l2` is second
order: at rtol 1e-4 on the BP5 2 km aging problem it needs about ten times
the steps of ARKIMEX3 for the same (correct) first event, which further
reinforces that the default explicit path remains the production tool and
`-imex` is infrastructure for stiff local physics.  The ODE system is
unchanged: the implicit and
explicit parts sum to the same right hand side as `rsf_ODE_RHSFunction`, so
the option changes the integrator, not the physics.  The motivation is the
state-evolution stiffness of some laws, most severely PRZ (`-state_law 3`),
whose theta rate grows like Omega^2 above steady state: at rupture fronts the
state lags steady state by many e-folds regardless of law, and the stable
EXPLICIT step then collapses like 2 dc/(|v| Omega) (about 3e-8 s in the BP5
tests at 1 to 2 km), while the nonlocal stress transfer is not stiff on those
scales.

### The split

Per cell, with X = (psi, tau, sigma, slip), v the regularized flow-law slip
rate, eta = G/(2 cs), and denom = 1 + eta dv/dtau:

- explicit `G` (in `rsf_IMEX_RHSFunction`): the `Is`/`In` matvec stressing
  rates plus backslip, the dv/dsigma damping term, and the slip rate;
- implicit `F = Xdot - Fimpl` (in `rsf_IMEX_IFunction`): the state rate
  S(psi, |v|) and the radiation-damping correction -eta (dv/dpsi) S / denom.

The implicit part is purely local: the IFunction performs no matvec, and the
IJacobian (`rsf_IMEX_IJacobian`) is analytic, exact (verified against finite
differences to about six digits, including at states where the stage solves
are difficult; `-imex_check_jacobian` re-runs that verification), and block
diagonal with one small block per cell.  Block Jacobi with ILU(0) is an exact
factorization for this sparsity pattern, so the default stage linear solve
(`preonly` + `bjacobi`) is direct, costs per-cell work only, and adds no
matvecs.

### Options and the `imex_` prefix

- `-imex` enables the formulation (default off).
- The stage solvers (SNES, KSP, PC) live under the options prefix `imex_`:
  deliberate overrides are `-imex_snes_type`, `-imex_ksp_type`,
  `-imex_pc_type`, `-imex_snes_max_it`, and so on.  UNPREFIXED solver options
  do not reach them.  This is intentional: `petsc_settings.yaml` is read
  automatically and its entries are indistinguishable from command-line
  options, and a yaml tuned for the dense solves of the interact main program
  (`ksp_type fgmres`, `pc_type jacobi`, `ksp_rtol 1e-8`) was observed to
  silently replace the exact block stage solve, producing about 780 linear
  iterations per Newton solve, 594 nonlinear failures, and a factor of order
  30 to 100 slowdown of an aging BP5 2 km run on 32 cores (8.1e6 linear
  iterations, 146 s, versus seconds when configured as intended).  With the
  prefix, the same yaml is harmless.  A correctly configured run reports
  `total number of linear solver iterations` EQUAL to the nonlinear count in
  the TS summary; a larger ratio means the stage solve has been overridden.
- ARKIMEX flavor and adaptivity remain tunable through the usual TS options
  (`-ts_arkimex_type`, `-ts_adapt_type`, ...).

### Validation (single core, dense operator, PETSc 3.19.6, rtol 1e-4, BP5)

The following was observed in one configuration and other systems may differ
in detail; the qualitative picture is expected to carry over.

- Aging, 2 km: 1 yr in 358 IMEX steps versus 427 explicit (3bs), with the
  first event's onset, arrest, and peak slip rate agreeing to all printed
  digits.  100 yr completes in 468 steps in under a second of wall time.
- Slip law, 2 km: 2693 IMEX steps versus 3420 explicit for 1 yr, event
  identical to printed precision.
- PRZ (ssvinit IC), 2 km: 1 yr in about 18000 steps for BOTH integrators, with
  the same double event (peak rates 3.48/5.99 versus 3.48/6.01 m/s explicit)
  and the event phase shifted by about 0.0007 yr, comparable to the
  documented 3bs-versus-5dp phase sensitivity.  Pre-event trajectories agree
  to 4 to 6 digits at matched times.

### When to use it, honestly

At 2 km the in-event step is ACCURACY-limited, not stability-limited, so IMEX
cannot reduce the step count there, and each IMEX step costs several explicit
steps (four explicit stage evaluations with matvecs plus four local implicit
solves); for production 2 km runs of any law the default explicit path is the
faster tool.  The intended regime is PRZ at about 1 km and finer, where the
explicit integrator is stability-limited: there `-imex` improved the in-event
step from about 3e-8 s to about 2.6e-5 s in the test configuration, roughly a
factor 1000, but remains about a factor 50 short of the estimated accuracy
limit because the Newton stage solves fail during coseismic phases and cap dt.
The cause is characterized: at large stage dt the damping-only implicit
per-cell problem undergoes a saddle-node (a nucleating cell has no local
elastic feedback to arrest it within a stage), Newton converges to the fold
with a nonzero residual, and TS retries with a smaller step.  Fine-resolution
PRZ production runs are therefore NOT yet practical with either integrator; a
per-cell stage solver that can cross the instability jump (with the cell
self-stiffness moved into the implicit part) exists as an experimental patch
and is parked pending further work, as is the alternative of an operator-split
analytic state update outside the TS framework.

### File and interface notes

The IMEX code lives in `src/rsf_imex.c` (see the source layout table above);
the state laws themselves are in the shared `rsf_state_rate` in
`src/rsf_engine.c`, together with their analytic partial derivatives, so a new
evolution law is added in exactly one place and both integrator paths pick it
up.  The `use_imex` flag lives in `rsf_solve_settings`, the per-run fields
(`imex_cap_efolds`, reserved for the parked stage-solver work) in
`struct rsf_vars`.  The driver branches once, at TS setup in `rsf_solve.c`;
monitor, event detection, domain check, tolerances, and adaptivity are shared
with the explicit path unchanged.

---

## Output files

| file | cadence | contents |
|---|---|---|
| `rsf_monitor.dat` | adaptive (`-dt_monitor`, `-adx_monitor`, `-rdx_monitor`) | `step time[s] time[yr] dt[s] log10(max|v|) mean_slip mean_mu max_sigma min_sigma`. The main time series: dense through an event, sparse in the interseismic. Computed from the distributed solution with two reductions, so it is cheap and needs no gather. Flushed on every write, so a running job's progress is visible |
| `rsf_vel.times` | one row per field frame (`-field_step_interval`) | `frame step time[yr] time[s] log10(max|v|) mean|v| std|v| min|v| mean_slip`. The index for the `tmp_rsf/rsf_vel.gGGG.NNNNNN.bin` frames: every row corresponds to a frame on disk |
| `rsf_events.dat` | one row per slip-rate threshold crossing | `time[s] time[yr] onset(1)/arrest(-1) log10(max|v|) mean_slip mean_mu`. A raw crossing log; see the note on resolution above before reading events off it directly |
| `rsf_catalog.dat` | one row per completed event (`-rsf_catalog`) | onset and arrest time, duration, ruptured-cell count and area, mean and max coseismic slip, mean and max static stress drop, peak slip rate, `M0`, `Mw` |
| `rsf_geom.gGGG.dat` | once | per-group patch geometry for the field frames |

The velocity statistics in `rsf_vel.times` (`mean|v|`, `std|v|`, `min|v|`) were
previously written to a separate `rsf_stats.dat`. That file has been removed. It
duplicated the time and maximum slip rate already present elsewhere, and its
statistics could only be formed on the gathered solution, which is exactly what
the field-frame path already does, so the columns now ride along with the frame
they describe and every row indexes a frame that exists on disk.

Note that these are slip *speeds*. `vel_from_rsf` returns a signed velocity, since
the sinh regularization admits `v < 0`, and a signed mean or maximum would not
summarize how fast the fault is moving; the statistics are therefore formed on
`|v|`. The rate-and-state solve carries a single slip component whose direction is
set by `-rsf_slip_mode`, so `|v|` is the speed in that direction. Should the solve
ever be generalized to carry strike and dip slip simultaneously, the natural
generalization is `sqrt(v_strike^2 + v_dip^2)`, and the statistics block is where
that would go. For the same reason the gathered field now stores the solved
velocity in `fault[].u[slip_mode]` rather than always in `fault[].u[STRIKE]`, so a
dip-slip run is not mislabelled for anything downstream that reads `fault[].u`.

#### The gated (composite) laws: Sato-type and Kato and Tullis

Two further laws are available, both of which are the slip law plus the aging
law's healing term multiplied by a gate that shuts healing off away from
stationary contact. In `theta` form,

    d theta/dt = gate - Omega ln Omega,     Omega = |v| theta/dc

- `-state_law 4`, Sato-type: `gate = exp(-Omega/beta)`, keyed on `Omega`
- `-state_law 5`, Kato and Tullis (2001) composite law: `gate = exp(-|v|/Vc)`,
  keyed on `|v|`

In `psi` both become

    d psi/dt = b/dc ( v0 exp((f0-psi)/b) gate - |v| ln Omega )

whose second term is exactly the slip law and whose first is the aging law's
healing term times the gate. The two constants are held fixed for now, at the
values of the MATLAB reference implementation: `beta = 1e-2` and `Vc = v0/100`.
They are set in `rsf_init.c` and are not exposed as command-line options.

The steady states differ, and this matters when comparing runs. The Sato gate is
negligible near `Omega = 1`, so `Omega_ss = 1` and its steady state coincides with
aging, slip and PRZ; reproducing the steady-state canon is the point of that form.
The Kato and Tullis gate, keyed on `|v|` rather than `Omega`, does not vanish for
`|v| << Vc`, so there `Omega ln Omega = 1` and `Omega_ss = e^W(1) = 1.763`, raising
`psi_ss` by `b W(1) = 0.567 b` (about 0.007 for `b = 0.013`). This is the composite
law's well known steady-state offset, not an artifact of this implementation, and
it means a Kato and Tullis run is not directly comparable at low slip rate to the
other laws at equal `dc`. For `|v| >> Vc` the gate vanishes and the offset goes
away.

The `psi` forms of all five laws were checked against an independent MATLAB
implementation in the Gu et al. (1984) nondimensionalization, and agree to machine
precision. Note the numbering differs by one: MATLAB `evol` 1 to 5 correspond to
`-state_law` 0 to 4.

#### Normalization of PRZ, and comparability of the laws

All laws share the steady state `psi_ss = f0 - b ln(|v|/v0)`, i.e. the usual
`f_ss = f0 + (a-b) ln(|v|/v0)`, and in the normalization used here,
`d theta/dt = (1 - Omega^2)/2`, PRZ also shares the linearization about it
(`d theta/dt ~ -(Omega - 1)`), so the laws agree on linear stability and `h*`
at equal `dc`.

That choice has a consequence away from steady state that matters for cycle
comparisons: the rest healing rate.  At `v -> 0`, aging heals with
`d theta/dt -> 1` while `(1 - Omega^2)/2 -> 1/2`, so PRZ restrengthens a
locked fault at HALF the aging rate.  No constant rescaling of PRZ can match
aging in all three properties at once (steady state, linearized relaxation,
rest healing); the convention here matches the first two.  The alternative,
`d theta/dt = 1 - Omega^2`, matches steady state and rest healing but DOUBLES
the linearized relaxation rate, i.e. behaves like `dc/2` for stability and
nucleation and thus halves `h*`.

The healing difference appears to be a modest effect for cycle recurrence.
The onset strength is lowered by about `sigma b ln 2` (0.52 MPa at BP5
parameters), of order 10 percent of recurrence against BP5-like reloading
rates, and a direct test of the alternative normalization on the
single-patch slider benchmark (see slider/README.md) shifted the slider
recurrence by 2.3 yr, under 2 percent; a runtime option for the alternative
was tried and removed as not useful.  In a 2 km, 600 yr, dense, serial BP5
test (one configuration; other setups may differ), the measured full-rupture
(Mw >= 7) recurrence was 234.6 yr for aging versus 172.1 yr for PRZ (a
factor 0.73), while the mean coseismic stress drop was 7.60 MPa versus
3.38 MPa (a factor 0.44), with PRZ also producing about 190 small partial
events (median Mw about 5.3) between mainshocks.  The shorter PRZ recurrence
is therefore dominated by the law's coseismic behavior (strong weakening,
roughly half the stress drop, arrest at higher residual stress, plus stress
pre-release by the partial-event activity) rather than by the healing
normalization.  Caveat: the partial-event population is exactly the feature
most sensitive to under-resolution at 2 km (the process zone is about one
cell there), so the recurrence ratio should be checked at adequate
resolution before being treated as a property of the law.

PRZ is also written `d theta/dt = 1 - (|v| theta/(2 dc))^2`, the same law
under `theta -> theta/2`, which places `theta_ss` at `2 dc/|v|` and shifts
`psi_ss` by `b ln 2`; anyone porting results in or out should check which
convention the other code uses.
