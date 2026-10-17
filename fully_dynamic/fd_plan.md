# Fully dynamic extension of `rsf_solve`: plan and status

Updated 2026-10-12. Replaces the earlier dated log; history is kept only
where it explains a current choice. Science target: rupture jump
probabilities across stepovers, comparing quasi-dynamic (QD) and fully
dynamic (FD) sequences on identical fault systems.

## 1. Approach (settled)

Space-time boundary integral formulation with the static operators of
`rsf_solve` kept and a hereditary dynamic correction added:

    tau_i(t) = tau0 + loading + sum_j Ks[i,j] delta_j(t)
               + sum_{l=1}^{L} sum_j W[l][i,j] V_j(t - l dt) - eta V_i(t)

with W[l] = K[l] - Ks dt, K[l] the pre-integrated kernel for unit slip rate
during one step observed l steps later, and eta = mu/(2 beta) the exact
instantaneous term. Radiation damping is therefore not an approximation;
QD drops the W term. The same holds for normal tractions (Kn, Wn).

Consequences that shape everything else:
- Interseismic phases need no history (W times a slow V is negligible), so
  the existing adaptive QD integrator stays; FD runs only during events.
- In 2-D the correction decays algebraically, so a window L is needed; in
  3-D it is exactly zero once the S wave has crossed the pair.
- The first-lag operator K[1] couples neighbours reached by the P wave
  within one step and must be solved as a coupled (banded) system; a
  lagged or predictor-corrector treatment produces a first-order delay of
  the rupture front.

First target is 2-D plane strain (map view, mode II), because stepover
jumping is controlled by dynamic normal-stress change and the stopping
phase, neither of which exists in antiplane. Known bias: infinitely deep
faults; 2-D results are read as QD-versus-FD differences, not absolute
probabilities.

## 2. What exists and what it has been checked against

Python prototype `fdkern` (its README has module-level detail). All items
below are covered by tests unless marked otherwise.

| component | status | checked against |
|---|---|---|
| in-plane segment kernels (closed form, Tada and Madariaga 2001 / Ando et al. 2007) | done | independent derivation from the 3-D Green's function, 4e-14 median; static dislocation limits; plane-wave identity |
| antiplane kernels | done | same route, 6e-15 |
| large-time series of the kernels (needed for the window) | done | closed form at moderate beta t / r |
| windowed W tables, Toeplitz/FFT for parallel pairs, history buffers | done | dense convolution; kinematic self-similar crack |
| fixed-step FD stepper, end-of-step collocation, exact banded K[1] solve | done | exact Kostrov/Broberg self-similar factors (mode II 0.9631/0.8651/0.7662 at v/beta 0.2/0.5/0.8, mode III 0.8257/0.7052; reproduced to 1e-3); independent spectral mode III solver |
| linear slip weakening, forced nucleation | done | SCEC TPV205-2D against Dunham MDSBI 6.25 m: arrivals within 1 to 3 steps at 100/50/25 m, final slip within 0.4 percent, rms slip 0.011 m at 25 m |
| rate-and-state (aging, regularized), banded Newton in log V, exact aging-law state update | done | SEAS BP1-QD and BP1-FD at 50 and 25 m (event times, slip per event, station histories; deep-station peak slip rate still resolution-limited) |
| normal-stress coupling across segments | done | static sign check; Harris and Day cases D and A: trigger delays and locations within 0.2 s and 0.4 km |
| adaptive QD integrator (DOP853 on slip and state) | done | consistency with FD at slow slip rates; periodic sequences |
| QD/FD hand-off both ways, history rebuilt from QD dense output | done | switch threshold 1e-3 versus 1e-2 m/s: 1 percent in slip |
| variable FD steps (k dt_min, k <= 64, from xi dc / Vmax, cumulative kernels) | done | k = 1 versus k = 2 on the self-similar crack; event cost about 2000 steps |
| FD versus QD sequence on one fault | done | recurrence +24 percent, slip +30 percent, peak slip rate x2.5, extent +9 percent, healing by stopping phases visible (compare_qd_fd.py); unchanged at 50 m |
| strong rate weakening (flash heating, slip law) | built | steady-state checks only; low-prestress setup not yet run |
| two-segment sequence driver with jump classification (run_step.py) | running | 1 km step, 2 km overlap, 100 m: QD and FD both give a delayed jump (segment 2 ruptures 28 d after segment 1 in QD, 16 d dilational and 6 d compressional in FD), no dynamic triggering |
| non-parallel segment pairs (dense tables) | not built | needed for branches, TPV14/15-2D |
| free surface, antiplane (image method on a doubled fault) | done | surface-breaking crack relaxes to the static half-space solution; BP1 |
| free surface, in-plane or 3-D | not built | needed for in-plane half-space problems and the 3-D stage |
| C implementation in rsf_solve | not started | WP7 |

Open discrepancy, documented rather than tuned: Harris and Day case A
(S = 0.49) runs near the P-wave speed in our model (2.7 s to 14 km against
their 3.35 s), insensitive to nucleation radius, forced speed and dx, and a
compressional 4 km step is triggered in ours but not in theirs. Their 250 m
grid is coarser than L_c (about 200 m) and they used artificial damping.
Andrews (1976b) or Dunham (2007) transition lengths would settle the speed
question; Hu et al. (2016) or Lozos et al. (2012) the stepover one.

## 3. Numerical choices in force

- beta dt_min = dx / 2, end-of-step collocation (e_t = 1). This lies on
  the stability island of Tada and Madariaga (2001) for mode II and has
  shown no instability through supershear transitions. A systematic
  stability scan over (h_T, e_t) with rate-and-state friction has not been
  done.
- Window tolerance 1e-2 relative to the static interaction envelope
  mu dx^2 dt / (2 pi (1-nu) r^2); gives L of about 1.3 to 5 S-travel times
  across a pair depending on geometry. Tested on BP1-FD at 50 m: 1e-3
  (window three times longer) changes nothing to the printed precision
  over two events. A spot check on one stepover case with a dynamic jump
  remains worthwhile but is no longer a prerequisite.
- Switch to FD at 1e-3 m/s; return to QD after max V has stayed below
  1e-4 m/s for one full window.
- Resolution: cohesive zone resolved by at least 5 elements. For
  TPV5-like slip weakening that is 50 m; for the rate-and-state test fault
  (dc = 16 mm, b = 0.015, sigma = 50 MPa) 100 m.
- Cost (one core, Python): FD step dominated by the history convolution,
  O(L M) per pair. 30 km fault at 50 m, 12 s event: about 2 minutes; two
  28 km segments at 100 m: 1 minute per 12 s; rate-and-state event at
  N = 300: 30 to 40 s. Memory: tables 8 L M doubles per segment pair.

## 4. Plan

### WP6, stepover sequences: interim results (2026-10-16)
Two parallel 24 km segments, aging law, sigma 50 MPa, b - a = 0.005 on
the weakening parts, 2 km overlap, 100 m elements, QD and FD on identical
setups (run_step.py, sweeps step/step2/step3). Every sequence has the same
structure: segment 0 ruptures, segment 1 follows after a delay, a full
cycle passes (87 to 96 yr QD, 104 to 117 yr FD). With identical segments
the receiver is always within hours to weeks of failure when the source
ruptures, so the jump probability within a month is 1 in both modes and
the QD/FD difference is in the delay:
- 4 km margins (weakening patch 4 km from the step): delayed jumps of
  days to months at all widths 0.5 to 3 km; FD 1.5 to 2.5 times faster.
- 0.5 km margins, 0.5 km width: joint ruptures in both modes.
- 1 km width, dilational, margins 0.25 to 0.5 km: FD immediate (0.4 s
  after the source front reaches the step) or 10 to 30 min; QD 3 to 8 h.
  Immediate triggering occurs only in FD and only here; the case sits on
  the dynamic threshold (a threshold event flips between 0.4 s and 12 min
  at 50 m and with tol 1e-3, everything else reproduces).
- 1 km, margins 0.75 to 1 km: FD 0.5 to 5 h, QD 5 to 24 h. Compressional
  1 km: FD 1 to 3 h, QD 2 to 36 h.
- 1.5 km: FD 3 to 10 h / QD 30 to 40 h (dilational), FD 17 h to 6 d / QD
  2 to 43 d (compressional); 2 km: FD 2 to 3 d / QD 6 to 8 d and 12 to 14
  d / 27 to 30 d (margins 0.5 and 1 km).
- FD events carry about 25 percent more slip and 3 times the peak slip
  rate; dilational triggers faster than compressional throughout.
Design consequence: a probability needs variation of the receiver's
cycle phase at the time of the source event, which identical segments
cannot give. Next sweep: unequal segments (24 and 16 km, or different
a - b on the receiver) so that successive source events find the receiver
at different loading states; the jump probability within a window is then
the fraction of source events followed by a receiver rupture, per mode.

### WP6, stepover sequences (plan)
1. Driver for two parallel segments under backslip loading with
   rate-and-state friction, built from cycle.py (run_step.py exists; its
   first QD sequence predates the clock fix and the state-update fix and
   is rerun first); runs at 100 m first.
2. Event bookkeeping: per event, moment per segment, rupture extent,
   peak slip rate, receiver prestress at source nucleation, time of
   receiver nucleation relative to source arrest, receiver stress series at
   the stopping-phase arrival.
3. Jump classification: dynamic (within wave travel plus a few seconds of
   source arrest), postseismic (hours to days), independent. The QD/FD
   difference is expected mostly in the first class.
4. Parameters: step width 0.5 to 5 km, overlap including gaps,
   extensional versus compressional sense, a/b contrast, optional
   velocity-strengthening arrest zones at segment ends. Paired QD/FD runs
   share loading and spin-up; compare events from the spun-up cycle.
5. Ensemble size: tens of events per configuration. At 100 m this is
   hours per configuration in Python; the C port (WP7) is needed for 50 m
   ensembles.
Deliverable: jump probability per class and sense versus step width, QD
and FD, with window and resolution sensitivities on a subset.

### WP5, external checks (in progress)
- Supershear transition length, against Andrews (1976b) and Dunham (2007).
  Andrews: linear slip weakening, forced nucleation at 0.5 beta (the same
  procedure as run_transition.py), grid 0.1 L_c(A) where L_c(A) =
  (8/pi) mu (lambda+mu)/(lambda+2mu) G/(tau_0 - tau_f)^2, which is twice the
  Griffith half-length used in our scripts; his Fig. 9 band at S = 0.75 is
  L_c(A)/L = 0.11 to 0.17, i.e. 6 to 9 L_c(A), with picks "uncertain by a
  factor of 2". Ours at S = 0.75, same nucleation, 38 elements per static
  cohesive zone: 6.0 to 6.4 L_c(A) by the breakdown front (slip reaching
  d0) and 9.2 L_c(A) by first slip (daughter-crack formation), inside his
  band. Dunham: Dugdale friction, static BCS crack plus 20 percent
  perturbation, 40 elements per cohesive zone; at S = 0.5 we get 25.8
  L_fric by the breakdown front against his points read as 18 to 22, and
  his linear-to-Dugdale ratio of about one third is consistent with our
  linear results. The "Andrews fit" quoted earlier from memory is not in
  either paper and is dropped.
- Resolution is the controlling issue for transition problems: with 8 to
  20 elements per cohesive zone the transition occurs 4 to 6 times too
  early (S = 0.75 linear: 1.0 to 1.5 L_c(A) at 8 elements against 6.4 at
  38; Dugdale S = 0.5: 6.4 L_fric at 20 against 25.8 at 40), and the
  insensitivity to dx seen between 100 and 50 m was only because both were
  far from converged. Rule adopted: at least 40 elements per (static)
  cohesive zone wherever supershear transition or near-threshold dynamic
  triggering matters; for slip, stress drop and arrival times (TPV205-2D)
  10 to 20 suffice. The rate-and-state sequence runs at 8 elements per
  cohesive zone should be checked at 50 m for the FD/QD differences.
- Consequence for Harris and Day case A: L_c about 200 m, so our 50 to 100
  m runs had 1 to 2 elements per cohesive zone and are not meaningful for
  the rupture speed or the 4 km compressional jump; their 250 m grid is in
  the same regime. The case must be rerun rescaled (larger d0 at fixed S
  and stress ratios) with 40 elements per cohesive zone. The earlier
  attribution to their grid and damping is withdrawn as unsupported either
  way.
- Transition detection: a slip-rate threshold marks the daughter-crack
  formation; the breakdown front (slip reaching d0) marks the instability.
  Both are now reported.
- Time stepping: CFL 0.25 at e_t = 1 is violently unstable in mode II
  (Tada and Madariaga's map); CFL 0.5 is stable with grid-scale slip-rate
  noise of 10 to 25 percent of the peak behind the front. The
  Peirce-Siebrits epsilon scheme (collocation at e_t dt, implemented as
  the `et` argument of System) halves that noise at e_t = 2 at the cost of
  a 0.5 to 0.7 percent bias in the self-similar crack factors (a time
  delay of order (e_t - 1) dt). Default stays e_t = 1; e_t = 2 is the
  option when noise matters. A systematic (h_T, e_t) scan with
  rate-and-state friction is still open.
- SEAS BP1-QD/FD: done at 50 m (sandbox) and 25 m (workstation, about 2.5
  h for the FD sequence). QD matches group 4 to 0.01 yr and 0.01 m at both
  resolutions; FD matches sem2dpack in event times (302.11 against 302.21
  yr for the second event at 25 m), recurrence and slip per event, and
  Jiang's 25 m SBIM in station peak slip rates at 2.5 and 7.5 km (within 3
  percent at 25 m). Converged with respect to element size: timing, slip,
  recurrence, station histories. Not converged: peak slip rate at the deep
  station (12.5 km, in the nucleation zone), 1.50 m/s at 25 m against 2.27
  (Jiang) and 2.39 (sem2dpack), moving by 18 percent per halving. A 12.5 m
  run would settle it but is not needed for the stepover work. A comparison
  with the rsf_solve BP1-QD output is still worth doing.
- Integration rule adopted from BP1: the aging-law state must be advanced
  with its exact solution over the step; explicit Heun steps produced early
  healing and 40 percent too little shallow slip. The same holds for the
  slip law (exponential relaxation, already in place).
- Enhanced weakening: to be added as strong rate weakening in the
  Noda-Dunham-Rice form (flash heating, slip law), which needs no extra
  state and is what Lambert and Lapusta use; thermal pressurization later.

### Performance options for the prototype (noted 2026-10-16, not started)
About 80 percent of an FD step is the history convolution, which streams
the whole W table (L x M complex, 280 MB at N = 800) once per step and is
memory-bound. Consequences:
- CPU threading buys at most a factor of 2 and is not worth doing;
  parallelism across configurations is already process-level.
- The worthwhile change is algorithmic: a uniformly partitioned FFT
  convolution over the lag axis (per Fourier mode) turns the O(L) work per
  step into O(log L) amortized, a factor of 20 to 50 at L of a few
  thousand, in Python or C alike. Compatible with the variable-step scheme
  because the history stays on the fine grid. Should be prototyped in
  Python before the C port adopts it.
- A GPU version is cheap: the convolution is an elementwise product and a
  reduction, which CuPy runs 20 to 50 times faster than one core, with
  the tables of a 25 m BP1 run (2.3 GB) fitting an 8 GB card. A backend
  switch (numpy or cupy) in shells.py and the stepper's matvecs would do
  it, QD phases staying on the CPU. Cannot be tested in the sandbox;
  would be delivered as a patch with a --gpu flag and a CPU/GPU identity
  test on the single-fault sequence.
Order if pursued: partitioned convolution first (reduces the work), then
the GPU switch for the remainder.

### WP7, C implementation in rsf_solve
Design in fd_c_design.md (2026-10-17): monitor-attached FD phase, no
changes to the QD code path, bit-identity and timing harness run at every
step, milestones M1 to M7 with their acceptance tests. M1 (source
reorganisation, skeleton, FFT, harness) delivered 2026-10-17. M2 (kernels,
tables, segment detection, static gate against Is; C tables agree with
the prototype dumps to 1e-12) delivered 2026-10-17. M3 (marcher core:
fine-grid history with copy-free convolution, cumulative kernels, LAPACK
banded solve, slip-weakening active-set solve, standalone fd_dynamic
driver) delivered 2026-10-18: TPV205-2D at 100 and 50 m reproduces the
prototype station series to 1e-10 relative, 25 times faster than the
Python; harness identity PASS.
- Kernels: Tada-Madariaga primitives for 2-D segments (four primitive
  evaluations per pair and lag) plus the large-time series; Tada (2006)
  triangles or space-time quadrature for 3-D.
- Tables: lag shells per pair; Toeplitz only when geometry allows; memory
  and cost as in section 3.
- Integration with the IMEX split: local damping and state implicit,
  nonlocal history explicit; the banded K[1] solve replaces the local
  friction inversion during FD phases.
- History restart and output of dynamic fields.

### Later
- Non-parallel pairs (dense tables) for branches and TPV14/15-2D.
- 2.5-D or 3-D stage before absolute jump probabilities; half-space
  static kernels with full-space dynamic correction as the first 3-D
  version, auxiliary traction-free surface as the consistent alternative.
- H-matrix compression of the dynamic shells: per-lag compression of the
  front is poor (factor 2); the arrival-shifted frame compresses to rank 4
  to 5 but needs a directional (PWTD or FDP=H) matvec. Stage-3 item for N
  above about 1e4.

## 5. Decisions open

1. Resolution for the first stepover ensembles: 100 m in Python now, or
   wait for the C port and start at 50 m.
2. Whether to spend a day on the Andrews/Dunham transition-length check
   before the ensembles, or carry the Harris and Day case A difference as a
   caveat.
3. Enhanced weakening (thermal pressurization, flash heating): the
   literature puts the largest QD/FD differences there; not in scope so
   far.
