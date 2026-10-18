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

## 9. Conventions established on the real operators (M2, 2026-10-18)

Measured with the static gate on rsf_solve's assembled operators (2-D
full plane, antiplane and in-plane, one and two segments):
- `Is` is assembled in SI (the scale is applied inside
  `calc_petsc_Isn_matrices`) with the stress-from-slip sign: `Is` equals
  `+Ks` of the tables to 3e-16. No scaling, no sign flip.
- `In` (with `-calc_sigma_dot`) equals `+Kn` of the tables to 2e-16,
  including the cross-segment blocks. interact's positive slip is
  left-lateral and `In` is compression positive; the prototype's positive
  slip is right-lateral and its `Kn` is tension positive; the two flips
  cancel. The marcher therefore uses the tables as they are and treats
  `Kn . delta` as a compression-positive change, consistent with the QD
  path `sigma_dot = In . v`.
- The element size enters as 2 l (interact's l is the half-length) and
  `t_strike` gives the tangent; dt_min = CFL dx / beta with CFL 0.5.

## 10. The hand-over as built (M4, 2026-10-18)

- PETSc 3.25 locks the solution vector during monitors, and replacing the
  solution in a TSPostStep hook left the event machinery with a stale
  previous value that it then tried to bracket by re-stepping from a state
  that no longer existed (a divergent second phase 250 s after the first
  hand-back, absent without -track_events). The switch is therefore a
  second TS event in rsf_solve's existing handler (slot
  uc->fd_event_index, upward crossing of -fd_vsw), and the phase runs in
  the post-event callback, which PETSc reinitialises the event state
  after. With -track_events off the handler is registered with the FD
  event alone. The monitor only records accepted (t, slip rate) pairs for
  the history seeding.
- The marcher works in time relative to the switch: tau0 is the elastic
  stress tau + eta V at the switch, the static term uses slip since the
  switch, and the backslip loading enters as sinc (t - t_switch). The
  prototype could use absolute time because its stress origin was t = 0;
  with absolute time here the deepest BP1 element received 1.3 GPa in the
  first step.
- Hand-back: psi, tau = strength at the final slip rate (= elastic stress
  at the end of the last step, the quantity the QD ODE integrates),
  sigma, slip += slip of the phase; TSSetTime to the end of the phase
  (the phase lasts up to half an hour of afterslip) and TSSetTimeStep to
  -fd_dt_restart.
- Known cosmetic issue: when the tracker's 1e-3 m/s crossing and the FD
  switch fire together, rsf_events.dat logs the onset with the post-phase
  state (peak slip rate below 1e-3). To be fixed by logging before the
  phase.
- Verified: BP1-FD three events and the single-fault sequence against the
  prototype (fd_plan.md); src/testing/rsf_solve/test_fd_sequence.sh is the
  one-minute regression of the latter.
