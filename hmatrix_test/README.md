# Two-fault rsf_solve accuracy-vs-speed sweep (HTOOLS, HACAPK, HMMVP)

A test harness for choosing H-matrix backend settings in `rsf_solve` on an
interacting two-fault problem, scoring each setting on both accuracy and speed
against the dense solve. BigWham is deliberately excluded here; the focus is the
three Okada-fed backends.

## What it does

1. Builds an interacting two-fault problem with `make_two_fault.py`: two
   BP5-QD-style vertical strike-slip faults laid down as parallel planes a few
   km apart, with their nucleation patches on opposite along-strike ends so the
   two systems couple through elastic stress transfer rather than rupturing in
   lockstep. The per-fault friction/medium/initial-condition recipe is copied
   unchanged from `bp5/make_bp5.py`.
2. Runs a DENSE reference once (`-use_hmatrix 0`), then runs HTOOLS, HACAPK, and
   HMMVP each across a range of its accuracy tolerance.
3. Scores every run on accuracy versus the dense solve and on speed, and writes
   a table plus a CSV.

## Files

- `make_two_fault.py` two-fault geometry + RSF input generator
- `sweep_rsf_hmatrix.sh` the sweep driver (runs everything, builds the table)
- `rsf_accuracy.py` compares one run to the dense reference
- `README.md` this file

## Quick start

On a multi-core node, with `rsf_solve` already built and `PETSC_DIR`/`PETSC_ARCH`
set (and PETSc libs on `LD_LIBRARY_PATH`):

```
./sweep_rsf_hmatrix.sh 1.0 260 8
#                       ds  yr  procs
```

This generates a 1 km two-fault problem, runs to 260 yr (past the first event),
and uses 8-way parallelism per run. Use ds=1.0 or finer: a 2 km grid is too
coarse to integrate through the coseismic rupture (the adaptive step size
collapses at the first event and the run appears to stall). Edit the CONFIG block at the
top of `sweep_rsf_hmatrix.sh` for paths, fault separation, the per-backend
tolerance lists, and `MPIRUN` (set it to `ibrun` or `srun` on a cluster). On a
single-core box, pass `procs=1` (or set `MPIRUN` to oversubscribe) for a
functional check; the timings are only meaningful with real cores.

## Accuracy metric

This follows the methodology already in `rsf_solve.md`. Both the reference and
each H-run write `rsf_monitor.dat` (a `log10(max|V|)` time trace) and, with
`-track_events`, `rsf_events.dat` (event onset/arrest times). For each run we
report, relative to the dense solve:

- `maxdV`, `rmsdV`: the maximum and RMS deviation of the `log10(max|V|)` trace,
  interpolated onto the reference time samples over the overlapping span (log10
  units).
- `dEv1_yr`: the first-event onset-time difference (yr).
- `dRec_yr`: the mean-recurrence difference when at least two events are present
  (yr).

These are differences relative to the dense run, not absolute errors, and they
are specific to the configuration (resolution, `rtol`, monitor cadence, friction
setup). A small deviation means the backend reproduces the dense sequence for
this problem.

## Speed metric

Parsed from `-log_view`: `total_s` (wallclock), `asm_s` (total minus
time-stepping, dominated by H-matrix assembly), `mv_ms` (mean per-`MatMult`),
`nstep` (TS steps). `mem_MB` is a best-effort read of each package's reported
H-matrix memory (dense is the exact `N^2` array; HTOOLS/HACAPK from their own
logs; HMMVP left blank).

## Reading the table

At matched tolerance, the deviation columns say how closely each backend tracks
the dense sequence and the speed columns say what it costs. On the single BP5
fault, `rsf_solve.md` found accuracy is a wash at tolerances at or below 1e-4,
because the ODE `rtol` (1e-4) sets the error floor, not the H-matrix tolerance;
in that regime the loosest tolerance whose deviation is still in the `rtol`
noise is the efficient choice, since tightening further buys accuracy you cannot
see while still costing assembly time. The reason this two-fault test exists is
to check whether that conclusion survives once cross-fault stress transfer is in
play, where event timing could in principle be more sensitive to the elastic
operator. Treat the single-fault result as a hypothesis to confirm here, not a
given.

## Caveats

- Resolution must carry the rupture. At 2 km the BP5 process zone is only ~3
  cells and the solve cannot step through the first coseismic event: the step
  size collapses to zero and every backend (including the exact dense solve)
  stalls. Use 1 km, where the single-fault BP5 cycle was validated, or 0.5 km if
  the two-fault coupling still stalls. The progress poll prints the latest dt; a
  dt collapsing toward 0 with the time stuck is the signature of an
  under-resolved rupture, distinct from a backend simply being slow.
- The dense reference is the cost driver at these sizes (no compression, O(N^2)
  matvec). To keep it tractable you can limit `stop_yr` to just past the first
  event rather than several recurrence cycles, since the trace and first-event
  timing already give the accuracy signal; or, if the dense reference is
  impractical at your chosen resolution, use a tight-tolerance H-run (for
  example hacapk at `ztol` 1e-6) as the reference instead and compare the looser
  settings against it.
- HMMVP needs a tighter tolerance than the others for clean cycle stepping. Its
  matvec is in-memory and cheap (the MPI path loads the operator into RAM once;
  per-matvec cost is comparable to HACAPK), but at a loose tolerance its forward
  operator is the least accurate of the three, and the adaptive rate-and-state
  integrator responds to that operator noise by collapsing its step size. In
  testing at a coarse resolution, HMMVP at `-hmmvp_tol` 1e-2 took roughly 200x
  smaller steps than HACAPK and barely advanced in simulated time, while at
  `-hmmvp_tol` 1e-7 it stepped normally and finished as fast as HACAPK. The
  effect is the cost showing up as step count rather than per-step cost, so watch
  the `nstep` column: a backend/tolerance that inflates `nstep` is paying for
  operator error through the integrator. The `RUN_TIMEOUT` config guards against
  any backend/tolerance combination whose step count blows up from stalling the
  whole sweep; set it generously on a real run, or to 0 to disable. This is the
  accuracy-vs-speed coupling the sweep is meant to expose, so include tight
  HMMVP tolerances in `HMTOL_LIST` rather than only loose ones.
- Parallelism: all four backends (dense, HTOOLS, HACAPK, HMMVP) run on `procs`
  MPI ranks, so each gets `procs`-way parallelism in its native mode. HMMVP now
  runs under MPI; set `HMMVP_OMP=1` to instead run it on one rank with `procs`
  OpenMP threads. Cross-backend timing comparisons are fairest at a fixed
  `procs` understood this way.
- Default tolerance lists are a starting point; widen or narrow them in CONFIG.
  Secondary knobs (`-mat_htool_eta`, `-mat_htool_compressor`, `-hmmvp_eta`) are
  left at their `rsf_solve` defaults and can be added to the sweep if of
  interest.
