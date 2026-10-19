# Fully dynamic extension of `rsf_solve` in C: design

Written 2026-10-17 after reading `src/rsf_solve.c`, `rsf_engine.c`,
`rsf_imex.c`, `rsf_ve.c` and the test scripts. Companion to fd_plan.md (what
to compute) and the Python prototype in this directory (how it was
verified). Nothing here is implemented yet.

## 1. Constraints, in order

1. Existing runs stay bit-identical. Elastic and visco-elastic `rsf_solve`
   runs without `-fd_mode` must produce byte-identical output files before
   and after every change, and run at the same speed (within timing noise,
   a 2 percent budget). This is the same contract the VE extension kept
   ("option-absent runs are bit-identical to the elastic code path") and it
   is checked mechanically, not by inspection (section 4).
2. Clean separation. All FD code lives in new files; existing functions are
   not edited except for option registration and a struct pointer. The FD
   machinery is attached to the TS only when `-fd_mode` is on, so the
   quasi-dynamic code path, including `rsf_compute_vel_and_stressing` and
   the IMEX functions, contains no FD branch at all.
3. Behaviour matches the prototype. Every milestone is accepted against a
   Python result produced earlier in this directory (tables to 1e-12,
   TPV205-2D and BP1-FD sequences to the printed precision), before any
   new physics is added in C.

## 2. How the FD phase attaches to the existing solver

`rsf_solve` integrates x = (psi, tau, sigma, slip) per patch with a PETSc
TS, tau being a state variable with the quasi-dynamic strength = stress
condition built into dtau/dt (radiation damping through
`shear_mod_over_2cs_si_for_damping`). The prototype instead switches
integrators at a slip-rate threshold and marches the dynamic phase with its
own fixed-step scheme. The C version does the same, and the switch is the
only contact with the TS:

- A TS monitor, `rsf_fd_monitor`, registered only when `-fd_mode > 0`,
  checks the maximum slip rate after each accepted step. Above `-fd_vsw`
  (default 1e-3 m/s) it calls the FD marcher with the current X and time.
- The marcher (`rsf_fd_run`) seeds the fine-grid slip-rate history from the
  last accepted steps (the prototype uses the dense output of the QD
  integrator; the C version keeps a short ring of accepted (t, vel) pairs
  for the same purpose), marches with variable steps k dt_min until the
  maximum slip rate has stayed below `-fd_vback` for one window, and
  writes psi, tau, sigma, slip back into X at the new time. tau is set
  from the strength = stress condition at the final slip rate, which is the
  quantity the QD ODE integrates, so the hand-back is consistent.
- The monitor then calls `TSSetTime`, `TSSetTimeStep` (a small step) and
  `TSRestartStep`, which is what the existing `rsf_post_event` already
  does after a discontinuous change; the next accepted step continues the
  QD integration from the new state. Whether modifying the solution inside
  a monitor is accepted by the PETSc version in use is the first thing M1
  tests; the fallback is to split `TSSolve` at FD events, which the driver
  comment warns against for output continuity but is acceptable for a
  discontinuous hand-over (the state is exact, not interpolated).

Nothing in the RHS, IMEX or VE code is touched. Restrictions in the first
version, refused at option parsing with a clear message: one MPI rank,
`-ve_mode 0`, 2-D segment geometries (types 30 to 33) made of straight runs
of equal, collinear patches (the Toeplitz structure the tables rely on),
no H-matrix requirement for the dynamic part (the static part uses the
assembled `Is`/`In` through `MatMult`, so any backend works for it).

## 3. Files and what moves from the prototype

| C file (new, `src/`) | prototype origin | content |
|---|---|---|
| `fd_kernels.c` | inplane.py, inplane_kernel.py, inplane_series.py, antiplane.py, ref_ando2007.py | closed-form Tada-Madariaga primitives, antiplane Ando primitives, large-time series switch, segment ramp stress S(x, y, t) for both modes |
| `fd_tables.c` | shells.py (PairTable, W_table) | per segment pair: lag table W[l] = S((l+e_t)dt) - S((l-1+e_t)dt) - Ks dt on the Toeplitz offset axis, FFT along offsets, window from the 1e-2 envelope criterion, cumulative kernels K_k |
| `fd_fft.c` | numpy.fft | complex power-of-two FFT (no new dependency); FFTW optional later |
| `fd_history.c` | shells.py (History, history_stress, push) | fine-grid ring buffer of FFT'd slip-rate slices, copy-free two-segment convolution |
| `fd_friction.c` | stepper.py (RateState, SlipWeakening, FlashRateState) | banded Newton in log V on the coupled first-lag system with the bracketed elementwise fallback, exact aging-law and slip-law state updates, slip weakening for the benchmarks |
| `fd_band.c` | scipy.linalg.solve_banded | banded LU for the first-lag system (bandwidth grows with k) |
| `rsf_fd.c` | cycle.py | options, setup (segment detection, table build, memory report), monitor, marcher, hand-over, restart records |
| `includes/fd.h` | | `struct fd_vars` hung off `struct rsf_vars` as a pointer (NULL when off) |

Geometry conventions to carry over exactly: `Is` is assembled in the
stress-drop convention, the prototype's kernels are stress-from-slip
(global sign flip, gated at setup against `Is` as the VE code gates its
K0 block); slip sign and normal direction as in the interact 2-D elements
(to be read from `read_geometry.c` and checked against the static Ks of
the prototype on a two-segment geometry, M2).

## 4. The regression harness (built first, run at every step)

`fully_dynamic/c_port/regress_rsf.sh`:

1. Builds a reference binary from a named commit into `bin/rsf_solve.ref`
   (once) and the working-tree binary into `bin/rsf_solve`.
2. Runs a fixed case set with both, identical options, in scratch
   directories:
   - 2-D antiplane elastic, BP5-style depth-dependent a-b, 300 yr
     (`test_twod/anti_plane/gen_bp2d.py`, as in `run_bp2d_cycle_test`
     stage 1);
   - the same with `-ve_mode 2` (plate over Maxwell) and `-ve_mode 1`;
   - 2-D in-plane with `-ve_mode 3` if its kernel generator runs in
     reasonable time, otherwise the elastic in-plane case;
   - one small 3-D case from `interact_examples`;
   - each with `-track_events -rsf_stations` so that `rsf_events.dat`,
     `rsf_monitor.dat` and the station files exist.
3. Compares `md5sum` of every output file between reference and current;
   any difference fails the run.
4. Times each case three times with each binary, reports the ratio of the
   minima; a ratio outside 0.98 to 1.02 is flagged.

The harness is the acceptance test for every patch to `src/`, including
the ones that only add files, and its log is attached to each patch.

## 5. Milestones

M1. Harness and skeleton. `regress_rsf.sh`; `rsf_fd.c` with option parsing
(`-fd_mode`, `-fd_vsw`, `-fd_vback`, `-fd_tol`, `-fd_et`, `-fd_kmax`,
`-fd_xi`), struct allocation, a monitor that only logs the maximum slip
rate, and a `TSRestartStep` exercise behind `-fd_mode -1` (debug) to settle
the hand-over mechanism. Acceptance: harness identical and same speed.

M2. Kernels and tables. `fd_kernels.c`, `fd_fft.c`, `fd_tables.c` with a
standalone test program `test_fd_tables` that writes S, W, K_1 and K_k for
a given segment pair and compares against reference tables dumped by the
prototype (`shells.py`, text) to 1e-12; segment detection on an interact
2-D geometry; static Ks of the detected segments gated against `Is`.
Acceptance: table comparison passes for in-plane and antiplane, parallel
and offset pairs; harness identical.

M3. Marcher with slip weakening. `fd_history.c`, `fd_band.c`,
`fd_friction.c` (slip weakening only), `rsf_fd_run` driven by a standalone
driver `fd_dynamic` that reads an interact 2-D geometry and a prestress and
runs one dynamic rupture. Acceptance: TPV205-2D at 100 and 50 m matches
the prototype's station series to the printed precision (arrivals, final
slip, peak slip rates), and the Kostrov factors to 1e-3.

M4. Rate-and-state, hand-over, sequences. Banded Newton with fallback,
exact state updates, history seeding from accepted QD steps, monitor and
hand-back. Acceptance: BP1-FD at 50 m reproduces the prototype's event
times (196.37, 301.80, 420.54 yr), slip per event and station peaks to the
printed precision; BP1-QD unchanged; harness identical.

M5. Two segments. Segment pair tables, normal-stress coupling, the
stepover sequences of WP6 reproduced for two configurations (delayed and
immediate triggering). Acceptance: delay distributions match the prototype
within the threshold-event ambiguity documented in fd_plan.md.

M6. Performance. Partitioned FFT convolution over the lag axis (the
O(log L) amortized scheme noted in fd_plan.md), then OpenMP over Fourier
modes if the machine warrants it. Acceptance: identical results to M5 to
1e-10, speed-up reported.

M7. Beyond the prototype. Non-parallel pairs (dense lag tables with the
same marcher), VE memory carried through the FD phase, MPI distribution of
the tables, 3-D elements. Each with its own validation case; not planned in
detail here.

## 6. Decisions taken for M1 (2026-10-17)

- Options to `rsf_solve`, not a separate binary. The FD phase reuses the
  whole of `rsf_solve` (geometry, operator assembly and H-matrix
  backends, friction parameters, initial conditions, loading, events,
  catalog, stations, checkpoints); a standalone `rsf_fd_solve` would
  duplicate all of it for the sake of a marcher that needs none of it
  changed. Standalone binaries are used only for unit tests and the
  slip-weakening benchmarks (`src/testing/rsf_solve/`), which the
  rate-and-state framework of `rsf_solve` cannot express.
- Monitor, not TS event. Agreed: the switch threshold is immaterial
  within an order of magnitude (prototype, BP1), so root-finding on the
  crossing buys nothing, and the monitor keeps the existing event
  handler (`-track_events`) untouched. PETSc 3.25.5; the `TSRestartStep`
  inside a monitor is exercised by `-fd_mode -1` before the marcher
  depends on it.
- Source layout. The `rsf_*` sources moved to `src/rsf_solve/` (VPATH),
  `src/rsf_solve.c` stays as the driver; new FD sources go to
  `src/rsf_solve/`, their tests to `src/testing/rsf_solve/`. The move is
  pure (no code changes) and the harness checks it.
- FFT. Built-in radix-2 (`fd_fft.c`) behind a plan-based interface
  (`fd_fft.h`) that is a subset of FFTW's; `-DUSE_FFTW` plus `FFTW_LIBS`
  in the machine config switches the backend without touching callers.
  `bin/test_fd_fft` (target `rsf_testing`) checks both.
- Reference commit for the harness: the master at the time of each
  patch; the harness takes it as an argument and builds it in a git
  worktree.

## 7. The VE framework under the same contract

`rsf_ve_apply_sink` is called inside `rsf_compute_vel_and_stressing` on
every RHS evaluation and returns immediately when `ve_np == 0`. That is
one function call and one integer test per RHS evaluation against a
dense or H-matrix matvec of O(N^2) or O(N log N): not measurable, and
with no floating-point operation on the path it cannot change results.
So the elastic QD code already satisfies the bit-identity and speed
contract with respect to VE, and the harness confirms it on every run
(the elastic case is in the set). The only change worth making, if one
wanted the RHS to contain no VE reference at all, is a function pointer
`rsf->stressing_hook` set at setup and called unconditionally; that
replaces a branch by an indirect call and gains nothing. Left as is.
The FD extension needs no RHS hook at all, which is the stronger form
of separation. VE and FD together (memory states carried through a
dynamic phase) is M7 work.

## 8. Questions that remain

- Which machine the 2 percent timing budget refers to (the harness prints
  ratios; the budget is applied by whoever reads them).
- A quick 3-D case for the harness (`EXTRA_CASES`), so that the 3-D
  elastic path is covered as well.
- Whether `-fd_mode -1` behaves on PETSc 3.25.5 as the design assumes
  (the first thing to run after building M1).

## 11. M5 (2026-10-18)

Two parallel in-plane segments with `-calc_sigma_dot 1`: the marcher's
normal coupling (nsign = +1, section 9) is active through the phase,
`Hn_dot = sinc[1]`, sigma handed back. The prototype's stepover case (two
24 km segments, 1 km step, 2 km overlap, 0.5 km velocity-strengthening
margins, 100 m) gives segment 0 at 1.0549 yr with 1.173 m mean slip over
all 480 elements (prototype 1.175), segment 1 after 53 min (prototype 53
min) with 1.550 m (prototype 1.546). Gates 4e-16 on both operators. The
post-event callback now runs the tracker's logging and catalog capture on
the pre-phase state and the phase afterwards. Sandbox setup for a hand-off
is in c_port/SANDBOX.md; usage and theory for users in rsf_solve_fd.md.

## 12. State after the first production use (2026-10-19)

- The workstation sweep with run_step_fd.sh (16 configurations, 16 cycles
  each, QD and FD) ran without intervention; results in fd_plan.md (WP6,
  C sweep). Per FD cycle on two 24 km segments at 100 m: a few minutes.
- BP1-FD at 25 m: reproduces the prototype (fd_plan.md); about 1 h per
  phase, no faster than the Python because the history convolution (same
  algorithm, L M per step and source) dominates; M6 addresses it.
- Unequal segments are in gen_step.py/run_step_fd.sh (no solver change;
  the defaults reproduce the previous inputs byte for byte).
- Known limits unchanged (section 3 of rsf_solve_fd.md). Next code items,
  in order: station output of the QD stations during phases (the
  coseismic part is currently only in fd_phase_*_st*.dat); M6 partitioned
  convolution; M7 scope list.

## 13. History alignment fix (2026-10-20)

Found while writing down the exact convolution that M6 must reproduce.
The history term evaluated before a step of k fine substeps must be

    sum_{l >= k+1} W[l] V^(next + k - l),

lag l counted from the end of the step (the step's own lags 1..k enter
through K_k). The prototype (shells.history_stress) and the C port
(fd_sys_history) used V^(next - l), i.e. the kernel was misaligned by k
slices: by one step for k = 1 and by up to 64 steps (kmax) during slow
phases. Checked against a direct space-time sum of the closed-form kernels
(n = 12, k = 1, 2, 4): corrected version to 3e-14, old version off by 40
to 70 percent. Both codes are fixed identically; the C TPV205 test
reference was regenerated from the corrected prototype and the C code
again matches it to 1e-10.

Effect of the fix, all with the corrected code:
- self-similar crack factors: within 2e-3 of the exact values in modes II
  and III at 100 and 50 m, as before (the old alignment was also within
  3e-3; this test does not discriminate);
- TPV205-2D against Dunham's MDSBI 6.25 m: main-front arrivals at -4.5 and
  -7.5 km are now about 1.5 steps early instead of 1.5 steps late (100 m:
  -18 and -31 ms against +25 and +13 ms; 50 m: -11 and -16 ms against +11
  and +5 ms); final slip within 0.015 m (100 m) and 0.009 m (50 m). The
  scheme is first order in time either way; the corrected version is the
  one consistent with the stated discretization;
- single-fault sequence: recurrence 82.23 yr (82.06 before), peak slip
  rate 2.475 m/s (2.430), mean slip unchanged;
- BP1-FD at 50 m: events at 196.360, 301.695, 420.482 yr (301.778,
  420.520 before); recurrence 105.33 and 118.79 yr (sem2dpack 106.5 and
  118.7); station peak slip rates 4.08, 3.29, 1.43 m/s in event 1 (4.13,
  3.35, 1.27 before).
Conclusions drawn from earlier runs are unaffected within these margins.
The time accuracy itself (first order, arrivals within about 1.5 steps)
is a property of end-of-step collocation with piecewise-constant slip
rate; a second-order variant (midpoint collocation or a corrected
first-lag operator) is a candidate for later work, not part of M6.

## 14. M6: speed (2026-10-20)

Measured first: on a BP1 phase at 50 m (one core) the history convolution
was 90 percent of the time, and two thirds of the steps were long
substeps (k = 64) in the afterslip tail. That step mix limits what a
block-based convolution can do, because a step that jumps past a whole
block of output times needs a new block each time. What was done, each
change checked against the direct code (TPV205 and both Harris and Day
cases to 1e-8, identical step counts, peak slip rates and slips on the
sequence and BP1 phases):

1. Friction solves use the band of the first-lag operator (it is pruned to
   its band, so the band product sums the same nonzero terms in the same
   order as the dense one): solve time on a BP1 phase 35 s to 10 s.
2. A step of k substeps pushes its slice with one spatial FFT, copied into
   k ring slots: push time 8.6 s to 1.3 s.
3. Partitioned convolution (fd_sys_set_conv): lags 1..D summed directly,
   lags D+1..L in partitions of B lags by overlap-save FFTs of length 2B
   along time, per spatial mode, with a frequency-domain delay line of
   input windows; D = q0 B - 1 with (q0 - 1) B >= kmax, so a far block never
   needs inputs that a variable step has not pushed yet. Time FFTs are
   batched over the spatial modes (contiguous butterflies). Mode 1
   (partitioned only) releases the time-domain tables beyond D and the
   ring beyond D + 2B + kmax (about twice the direct memory) but is slower
   than direct for long substeps; mode 2 (hybrid, default) keeps the direct
   representation (about three times the memory), uses the far blocks for
   steps with 2k <= B and the direct sum otherwise, and forms input windows
   only when a far block needs them. B is chosen automatically (64 for
   BP1 at 50 m); 256 was slower.
4. Quasi-steady hand-back (-fd_back_eps, off by default): the hereditary
   term of a slip rate that is constant over one window equals the QD
   radiation damping, so a phase may end once every slip rate has changed
   by less than back_eps over a window (and the maximum is below half the
   switch threshold). With 0.05, BP1 event times move by about 1e-5
   relative and peak slip rates by under 1 percent; the single-fault
   sequence by 2e-5 in event times.

BP1 50 m, first phase, one core:

| configuration | wall |
|---|---|
| before M6 (direct, -fd_back_eps 0) | 312 s |
| direct, items 1 and 2 | 262 s |
| hybrid (default), -fd_back_eps 0 | 183 s |
| hybrid, -fd_back_eps 0.05 | 110 s (history 95, solve 9, static 4, push 1) |

So 1.7x at unchanged results and 2.8x with the quasi-steady hand-back. The
remaining time is the direct sum for long substeps. Options not taken: a
larger maximum substep (-fd_kmax 256/1024: fewer steps but slower per step
and 0.004 yr different event times; the operator cache was raised to 16
sizes so they run), larger partitions (slower), a non-uniform partition
(small blocks for early lags, larger for late ones; the structural fix for
long substeps, more complex, deferred until the 25 m runs show the
convolution still dominating). FD_PROFILE=1 in the environment prints the
per-phase split (history, static, solve, push).

## 15. Coarse history for long substeps (2026-10-21)

The profile after M6 (BP1 at 25 m: history 97 percent of the phase time)
showed the long substeps of the afterslip tail as the remaining cost. Such
a step of K = kmax substeps pushes K identical slices. If the steps start
at multiples of K, the input is constant on aligned blocks of K slices, and
once a pair's whole window consists of such blocks the history term is
exactly

    sum_{l >= K+1} W[l] V^(T - l) = sum_{j >= 1} A[j] V_(block n - j),
    A[j] = sum_{l = jK+1}^{(j+1)K} W[l],

with K times fewer terms. Implemented as fd_sys_set_coarse (aggregated
tables built once from the full time-domain tables, so not with -fd_conv 1)
and a step selector that only grows k when the current fine index is a
multiple of the new k (-fd_coarse 1, default).

Exactness: with the alignment rule in place, switching the coarse sum on or
off changes the single-fault sequence by roundoff (8e-8 yr in the second
event time, one step). The alignment rule itself changes the step sequence
during nucleation and arrest: single-fault event times by 1e-4 relative,
BP1 50 m event times by 0.003 yr and peak slip rates by about 1 percent
(the quantity most sensitive to the time discretization).

Speed, BP1 50 m, one core: 50, 43, 47 s per phase against 183 s (M6
hybrid) and 312 s (before M6); history 20 s, solve 17 s, static 7 s,
push 3 s per phase. Single-fault sequence: 8 s per phase against 13 to
16 s.

BP1 25 m with -fd_back_eps 0.05 (Thor, before the coarse history): event
times 196.3638, 302.0828, 420.8594 yr against 196.3638, 302.0832, 420.8598
with the default criterion (1e-6 relative), peak slip rates within 0.04
percent, but only 12 percent faster, because the slip-rate floor
(-fd_back_vmax) binds for most of the tail at this resolution. With the
coarse history the tail is cheap anyway; -fd_back_eps stays opt-in.
