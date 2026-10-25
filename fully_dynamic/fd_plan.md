# Fully dynamic extension of `rsf_solve`: plan and status

Updated 2026-10-25. Section 0 is the current summary; sections 2 and 3
are the state of 2026-10-12 and are kept for the record; section 4 holds
the dated work notes (newest topics first within the 3-D block). Science
target: rupture jump probabilities across stepovers, comparing
quasi-dynamic (QD) and fully dynamic (FD) sequences on identical fault
systems. For how to resume the work see HANDOFF.md.

## 0. Status summary (2026-10-25)

### 0.1 What the code does

rsf_solve runs seismic-cycle simulations as QD between events and fully
dynamic (space-time boundary integral, inertia and waves) above a slip-rate
threshold (-fd_mode 1), handing stress, state, slip (and normal stress
where coupled) over in both directions. The QD path is untouched when
-fd_mode is not set (bit identity against the master is checked by the
harness, fully_dynamic/c_port/regress_rsf.sh, after every change).

| | 2-D (in-plane and antiplane) | 3-D (planar fault, squares) |
|---|---|---|
| geometry | one fault or parallel segments (stepovers), one plane strain or antiplane problem | one planar fault, strike or dip slip |
| friction | rate and state (aging law, regularized), slip weakening (tests) | the same |
| time scheme | second-order (Noda, Sato and Kurihara 2020) default, first order optional | the same |
| steps | variable substeps (powers of two to 64) with an exact coarse history | fixed fine steps; quasi-steady hand-back (-fd_back_eps) for speed; long steps experimental |
| normal-stress coupling | parallel segments | from the free surface (In plus the surface's dynamic part), with -calc_sigma_dot |
| free surface | antiplane: mirror (exact); in-plane: prototype only | vertical strike slip: mirror (exact for shear) plus optional opening layer (exact); any dip: full virtual surface |
| statics | interact's exact operators (Is, In) | the same, gated against the tables |

### 0.2 What it has been checked against

- 2-D: Kostrov and self-similar cracks (exact), TPV205-2D against Dunham,
  Harris-Day stepover cases, BP1-QD and BP1-FD against the SEAS
  submissions, bit identity of the C code with the Python prototype.
- Time scheme: on BP1-FD the second-order scheme converges 3 to 5 times
  faster than the first-order one; on BP5-FD at 1 km it matches Lambert's
  BICyclE arrivals to 0.1 s and peak slip rates to 2 percent at the same
  step; random-perturbation stability maps for both schemes in 2-D.
- 3-D kernels: Tada's quadrant formulas (slip along x1 and x2, opening)
  against an independent point-source quadrature (1e-8 to 1e-11), Kelvin
  statics, the plane-wave limits; C against Python at 600 points.
- 3-D free surface: mirror plus opening layer converges to Okada's half
  space (4.5e-4 at dx/4); the full virtual surface for buried dipping
  faults against Okada statically (1e-3 to 1e-4) and against the
  independent 2-D plane-strain prototype dynamically (0.3 to 0.6 percent
  in slip, identical rupture times, with and without normal-stress
  coupling, thrust and normal sense).
- SCEC TPV101 (whole space, aging law) at 100 m with the second-order
  scheme: arrivals within 0.08 to 0.13 s, final slip within 1 percent and
  peak slip rates within 2 percent of Dunham's MDSBI and 8 percent of
  Kaneko's SPECFEM3D; second-order convergence from 250 to 100 m. SCEC
  TPV10 (60 degree normal fault reaching the surface, full virtual surface)
  converging towards SPECFEM3D at 1 km and 500 m (final surface slip within
  2 percent at 500 m; arrivals and peaks need finer elements).
- BP4-QD and BP5 setups reproduce the published first-event timing once
  the participants' initial conditions are used (BP4) and the second-order
  scheme is on (BP5).

### 0.3 Findings worth remembering

1. The first-order scheme with end-of-step collocation is stable only
   near cS dt/dx = 0.5 in mode II; the second-order scheme is stable from
   0.4 to at least 2 (2-D). It is the default; the only measured drawback is
   a too-early triggering right at a stepover's jump threshold at coarse
   resolution (Harris-Day, 100 m), where the first-order scheme happened
   to be closer.
2. BP4-QD was 80 yr late against the published runs because the benchmark
   description's uniform pre-stress is not what the participants used
   (steady state at V_init, seed 0.001 m/s); with their setup it matches.
3. For BP5 (vertical strike slip, constant normal stress) the mirror's
   normal-traction error changes the first event by less than 0.5 percent:
   the mirror is adequate and 2 to 12 times cheaper than a surface.
4. The 3-D virtual surface has an unstable checkerboard mode for
   cS dt/hs between about 0.85 and 1.05 (a mixed-mode grid-scale mode that
   does not exist in 2-D; per-mode analysis reproduces the band; Peirce and
   Siebrits' "intermittent instabilities"); the code keeps the ratio at 2.
5. For a fault whose top edge lies in the surface, the full surface
   converges slowly at the trace: the top element's edge carries a
   singular field that only fault refinement removes (surface refinement,
   graded rows and a tied extension row above the surface do not); the
   mirror has no such edge. Buried dipping faults are fine.
6. Low-rank compression of the fault-surface tables after the along-strike
   FFT gains little (thin blocks); the surface is kept small instead
   (narrow strips, graded rows), 0.1 to 1.5 GB at 1 to 2 km.
7. The 3-D long-step (variable substep) solve settles on spurious fast
   solutions near arrested rupture edges with a damped Jacobi iteration;
   it needs a proper nonlinear solver before 3-D sequences become cheap.

### 0.4 Open items, in the order I would take them

1. Surface-breaking dipping faults: a fault-resolution study at the trace
   (quasi-2-D thrust with the top at the surface at 1 km, 500 m, 250 m
   against the 2-D prototype), then TPV10 when the SCEC data return.
2. A robust long-step solver in 3-D (Newton-Krylov with the FFT operator),
   then the coarse history, so sequences on dipping faults and BP5-FD
   sequences become affordable.
3. Stepover ensembles in 2-D and 3-D with the second-order scheme (the
   original goal; the WP6 work below predates the scheme change).
4. 2-D: non-parallel segments in the C code (bends, branches), slip law
   and flash heating; TPV101 station comparison when the data return.
5. Speed: the surface's outer iteration (2 to 12 times a mirror run),
   OpenMP over wavenumbers and lags.

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

## 2. What exists and what it has been checked against (state of 2026-10-12, kept for the record)

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

## 3. Numerical choices in force (2026-10-12; superseded where section 0 says so)

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

### WP6, C sweep with rsf_solve (2026-10-19, Thor's workstation)
run_step_fd.sh / summarize_rsf_step.py, two 24 km segments, 100 m, 2 km
overlap, 16 cycles per configuration, QD and FD on identical setups
(+W is the dilational side in interact's left-lateral convention):

| width, margin | QD delay | FD delay |
|---|---|---|
| 1.0 km, 0.5 km, dilational | 3 h (15/15) | 6 joint ruptures and 7 within 24 to 53 min |
| 1.0 km, 1.0 km, dilational | 21 h | 2 to 5 h |
| 1.0 km, 0.5 km, compressional | 2 to 36 h (median 21 h) | 2 to 3 h |
| 1.0 km, 1.0 km, compressional | 8 to 14 h | 43 min to 2 h |
| 1.5 km, 0.5 km, dilational | 31 h | 2 to 10 h (median 5 h) |
| 1.5 km, 1.0 km, dilational | 7 d | 2 d |
| 1.5 km, 0.5 km, compressional | 11 to 42 d (median 13 d) | 8 h to 6 d (median 3 d) |
| 1.5 km, 1.0 km, compressional | 7 to 15 d | 23 h to 15 d (median 42 h) |

Every source event is followed by a jump within a year in all 16
configurations and both modes (identical segments: the receiver is always
near failure). Inertia shortens the delay by 3 to 10 and lengthens the
recurrence by about 20 percent (102 to 115 yr against 87 to 97). New with
16 cycles: at 1 km width, 0.5 km margins, dilational, FD settles into a
two-cycle rhythm, segment 0 then segment 1 after 25 to 30 min in one
cycle, segment 1 nucleating and both rupturing in one dynamic phase in the
next; QD gives 3 h in every cycle. FD delays are also much less variable
than QD delays at the same configuration. The Python results (step2/step3)
are confirmed throughout.

Next sweep (to turn delays into a probability): unequal segments, so that
the receiver's cycle phase at the source event varies. Plan: segment 2 of
16 km (or b - a reduced on segment 2), both senses, widths 1 and 1.5 km,
margins 0.5 and 1 km, 30 cycles; report the fraction of source events
followed by a receiver rupture within one day and within one month, QD
against FD. Tooling in place: run_step_fd.sh W ov vs tend mode L2 dba2
(gen_step.py writes segments.txt; stations follow the segment centres),
summarize_rsf_step.py prints P(jump within 1 d) and within 30 d.

### WP6, unequal segments (2026-10-20, Thor's workstation) -- pinned
Second segment 16 km with b - a = 0.004 (first: 24 km, 0.005), 2 km
overlap, 2500 yr, run_step_fd.sh ... 2500 both 16 0.004. Fraction of source
events followed by a rupture of the other segment within 1 d / 30 d:

| width, margin, sense | QD | FD |
|---|---|---|
| 1.0, 0.5, dilational | 0.10 / 0.20 | 0.35 / 0.68 |
| 1.0, 1.0, dilational | 0.03 / 0.31 | 0.20 / 0.51 |
| 1.0, 0.5, compressional | 0.35 / 0.67 | 0.34 / 0.68 |
| 1.0, 1.0, compressional | 0.32 / 0.62 | 0.33 / 0.67 |
| 1.5, 0.5, dilational | 0.00 / 0.09 | 0.08 / 0.30 |
| 1.5, 1.0, dilational | 0.02 / 0.08 | 0.00 / 0.18 |
| 1.5, 0.5, compressional | 0.00 / 0.10 | 0.02 / 0.22 |
| 1.5, 1.0, compressional | 0.00 / 0.17 | 0.08 / 0.26 |

Caveats that limit the meaning of these fractions: most runs are limit
cycles (periods of about 176 to 300 yr with two or three source events per
period; 0.33 and 0.67 are 1/3 and 2/3 of one period), so the effective
sample is a few events, not 40 to 60; jumps in both directions are pooled;
partial re-ruptures within months are counted as source events. A
probability needs the periodicity broken (ensembles over initial
conditions, or along-strike heterogeneity) and per-direction counts. The
exploratory sweeps are paused here; these runs predate the history
alignment fix (fd_c_design.md section 13), whose effect on sequences is
0.1 to 0.2 percent in recurrence.

### WP6, stepover sequences (original plan)
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

### WP8, 3-D (started 2026-10-21)
Decisions: rectangular (square) elements on planar faults, 2-D Toeplitz
tables and 2-D FFTs over the fault plane, otherwise the same scheme as in
2-D; full space first (TPV3, TPV101, BP4-FD), the free surface afterwards.

Kernels (prototype, fully_dynamic/):
- gen3d_point.py: point moment-tensor displacement (Aki and Richards eq.
  4.29) for a general moment function; Navier residual 3e-14, static
  Kelvin double couple 1e-15.
- tada_quad.py: stress response to a constant slip rate on a quadrantal
  fault, transcribed from Tada (2006, GJI 164, Appendix A, eqs A12-A23; slip
  along x2 by the x1 <-> x2 symmetry), and the square element as four
  quadrants. Checks (check_tada_quad.py and inline): early-time self term
  equals -mu/(2 cT) exactly; late-time slope equals the static stress from
  a Gauss quadrature of the Kelvin double couple to 9 digits (in-plane,
  off-element); response to a cubic slip history, computed by twice
  integrating the Tada ramp response, equals the quadrature of the
  point-source solution over the element to 1e-9 (off-plane and in-plane
  receivers; s31, s32, s11; after both waves and between P and S).
- The response becomes exactly static once the S wave from the farthest
  corner of an element has passed (Huygens), so W[l] vanishes beyond a
  finite, offset-dependent lag.

Marcher (fd3d.py): one planar fault, two slip components, vector slip
weakening with the exact local solve along the traction direction, the
first-lag operator on the 3 x 3 stencil (the P wave reaches the eight
neighbours within dt = 0.5 dx/cS, checked) solved by Jacobi iteration,
history by 2-D FFTs, end-of-step collocation.

SCEC TPV3 (run_tpv3.py; Day et al. 2005, whole space, 30 x 15 km, slip
weakening), station values (PI at 7.5 km along strike, PA at 6 km along
dip):

| dx | PI: V > 0.1 m/s | PI final slip | PA: V > 0.1 m/s | PA final slip |
|---|---|---|---|---|
| 500 m | 2.742 s | 3.177 m | 2.237 s | 2.123 m |
| 250 m | 2.851 s | 3.123 m | 2.670 s | 2.070 m |
| 150 m | 2.901 s | 3.125 m | 2.750 s | 2.061 m |
| Day et al. 2005, Fig. 7 (read off) | about 2.97 s | about 3.12 m | about 2.86 s | about 2.06 m |

Final slip agrees within 0.5 percent at 250 and 150 m. The rupture arrival
converges at first order towards the reference: extrapolating 250 and
150 m gives 2.976 s at PI and 2.87 s at PA. The stress peak at PI is
81.24 MPa (the static strength) at 2.97 s with a single-element station
at 200 m, as in the reference. (The 200 m run is not in the table: the
nucleation patch of 3 km is not a whole number of 200 m elements and came
out 2.8 km wide.) Run times on one core: 10 s at 500 m, 1 min at 250 m,
2 min at 200 m, 7 min at 150 m (single-precision tables). 100 m needs
about 4 GB of tables in this form; the C port with real-to-complex FFTs
and the symmetries of the square element will cut that by about 4.

SCEC TPV101 (fd3d_rs.py, run_tpv101.py; rate and state, aging law, whole
space; computational fault x in [-21, 21] km, y in [-6, 21] km, locked
beyond): vector rate-and-state friction with an exact local solve along the
traction direction, exact aging-law update, predictor-corrector in the
state, nucleation load as specified. Rupture times (bilinear in the
rupture-time field) along strike at y = 7.5 km and along dip at x = 0:

| dx | x = 3, 6, 9, 12, 15 km | y = 4.5, 3, 1.5, 0 km |
|---|---|---|
| 300 m | 2.30, 4.22, 6.43, 8.16, 9.55 s | 2.53, 4.07, 5.78, 7.04 s |
| 250 m | 2.26, 4.02, 5.88, 7.45, 8.79 s | 2.47, 3.81, 5.20, 6.37 s |

The rupture arrests in the transition layer along dip and strike, as it
should, but is far from converged at these resolutions (9 percent faster at
9 km from 300 to 250 m): the process zone is about 400 m, so this problem
needs 100 m or less, which is beyond the Python prototype in memory and
time. TPV101 is therefore a test for the C version. (A first run had the
along-dip half-width of the velocity-weakening region wrong, W instead of
W/2, and ruptured into the velocity-strengthening region; fixed.)

Fixed rake: TPV3 at 250 m with slip and friction restricted to the strike
component (FIXED_RAKE=1) gives the same station values to the printed
digits (PI 2.851 s, 3.122 m; PA 2.670 s, 2.070 m) as the vector version.
rsf_solve's quasi-dynamic state has one slip component per patch (fixed
rake), so the C port keeps that: one table (strike slip to strike
traction) instead of four, and vector slip stays a later extension of both
the QD and the FD parts (BP4 asks for both slip components).

C port, standalone part (2026-10-22): src/rsf_solve/fd3d.c, fd3d.h,
src/testing/rsf_solve/fd3d_test.c. Tada quadrant kernels, square-element
tables as 2-D Toeplitz kernels with 2-D FFTs built from the 1-D FFT (only
the Hermitian half of each spectrum stored; -DFD3D_SINGLE stores it in
single precision, about 1e-5 relative differences in slip-weakening runs),
the 3 x 3 first-lag stencil with a Jacobi solve converged to 1e-13, slip
weakening and rate and state (aging law), slip along a fixed direction.
Against the prototype (fixed rake, 500 m; both codes' Jacobi tolerance
tightened to 1e-13): TPV3 over 12 s to 4e-10 in slip and 2e-9 in slip rate,
TPV101 over 6 s to 6e-9 and 4e-8; both cases are in run_tests.sh.

TPV101 with the C code (fault x in [-19.5, 19.5] km, y in [-4.5, 19.5] km,
FD3D_SMALL=1; at 300 m the station values equal those on the larger fault):

| dx | arrival at (9, 7.5) km | arrival at (0, 3) km | slip at (0, 7.5) km at 12 s |
|---|---|---|---|
| 300 m | 6.32 s | 3.90 s | 6.98 m |
| 250 m | 5.81 s | 3.68 s | 7.24 m |
| 200 m | 5.51 s | 3.61 s | 7.25 m |
| 150 m (single precision) | 5.22 s | 3.40 s | 7.29 m |
| 100 m (Thor) | 5.02 s | 3.31 s | 7.34 m |

With 75 m (Thor: 4.94 s at (9, 7.5) km, 3.26 s at (0, 3) km, hypocentral
slip 7.36 m, peak slip rate 4.77 m/s) the arrivals converge at about first
order (exponents of about 1.3 and 1.0 from 150, 100, 75 m; an earlier fit
without 75 m suggested 2.2 and was wrong), consistent with the first-order
time scheme at fixed CFL. Extrapolated arrivals: about 4.74 s at (9, 7.5) km
and 3.12 s at (0, 3) km, so errors of about 0.19 and 0.14 s (4 percent) at
75 m; hypocentral slip is converged to about 0.3 percent, the peak slip rate
at (9, 7.5) km is still rising (3.96, 4.49, 4.77 m/s from 150 to 75 m). The comparison with other codes is open: the SCEC pages for TPV101
show cross-code summaries but no station time series we could download
(2026-10-22); to return to, e.g. by asking the benchmark organisers or
from published comparisons. Cost: 200 m 7 min, 150 m 13 min on one
core in this sandbox.

rsf_solve 3-D setup and static gate (2026-10-22): with all elements
3-D rectangles, -fd_mode builds the fd3d tables from the geometry (squares
of one size on one plane forming a complete grid, -full_space 1, slip along
strike or dip) and compares their static limit with the assembled Is: a
16 x 8 vertical strike-slip fault of 500 m squares agrees to 3e-10 outside
the 3 x 3 self neighbourhood, with the same sign convention as in 2-D
(Is = +Ks). src/testing/rsf_solve/test_fd3d_gate.sh is the check.
-fd_mode 1 refuses 3-D until the phases are connected.

BP4: the SEAS BP4 description specifies a quasi-dynamic problem
(radiation damping); there is no fully dynamic BP4 reference (Jiang et al.
2022 compare QD only). fully_dynamic/c_port/run_bp4_qd.sh now wraps the
repository's seas_bp_tests/bp4/make_bp4.py (velocity seed Vi = 0.03 m/s in
the nucleation square, state unchanged, as the published runs; the
description's 1.1 tau0 prestress bump nucleates slowly, near 200 yr) and
uses the 3bs integrator as seas_bp_tests/bp4/run_bp4.sh does (without it a
2 km run stopped in a rejection storm at 136 yr). With the prestress bump
(the first version of this script) the events came at 213, 337, 465, ... yr
at 1 km and 192, 353, 508, ... yr at 500 m: recurrence 123 to 128 and 154
to 160 yr, as in Jiang et al. Fig. 5, but shifted (their first large events
near 110 to 125 yr and 155 yr). Dense and H-matrix runs agree to 0.02 yr.
With the velocity seed alone the seeded event died within a minute
(fault-wide mean slip 0.013 m at 1 km, 0.017 m at 500 m) and the next
event came at 193 or 222 yr. The cause is the initial state: the
participants (Jiang et al. 2022, eqs 15-16 and Table 3) started from steady
state at Vinit everywhere (theta = L/Vinit, tau0 the local steady-state
stress with the local a) and seeded BP4 with Vi = 0.001 m/s (0.03 m/s is
BP5's value). With the description's uniform tau0 from amax the
velocity-weakening region starts about 6 MPa below its steady-state stress
and a full rupture would have a dynamic stress drop of only about 0.4 MPa,
so the seed cannot spread. make_bp4.py now defaults to the Jiang et al.
version (init_state = "steady", Vi = 0.001). Results (Thor, 2026-10-23):

| dx | events (yr) | recurrence |
|---|---|---|
| 1 km | 0.00, 117.77, 246.64, 372.32, 499.71, 627.49, 755.21, 882.95 | 125.7 to 128.9 yr |
| 500 m | 0.00, 151.40, 306.03, 460.58, 615.11, 769.65, 924.18 | 154.5 to 154.6 yr |

Jiang et al. (2022), Fig. 5: first spontaneous event at about 110 to
125 yr (1 km) and about 155 yr (500 m), recurrence about 130 and 155 yr.
BP4-QD with rsf_solve agrees with the published runs at both resolutions.

Remaining for M9: the C port (M9): square elements on one planar fault (grid detection
like fd_segments), the quadrant kernels (tada_quad.py) in C, 2-D Toeplitz
tables with 2-D FFTs built from the 1-D FFT, the 3 x 3 first-lag stencil,
the coarse history; the static gate against interact's 3-D full-space
operator; a standalone driver for TPV3 and TPV101 against this prototype
at roundoff, then TPV101 at 100 m against the SCEC reference series, then
BP4 with rsf_solve.

### Time accuracy: second-order marching (2026-10-24)

Why: the BP5-FD first event (half the time step moves arrivals by 0.5 s at
1 km) and TPV101 (arrivals converge at first order) showed that the
first-order time scheme limits accuracy.

Tried and rejected: collocation inside the step (Stepper with sys.et < 1,
e.g. 0.5, kept as an option). TPV205-2D against Dunham: arrivals converge,
but final slip is off by 0.16 to 0.20 m at every resolution and peak slip
rates grow with refinement (+82 percent at 25 m): the oscillatory
instability of interior collocation with piecewise-constant rates that
Tada and Madariaga (2001) mapped.

Adopted: the "NL" scheme of Noda, Sato and Kurihara (2020, Earth Planets
Space 72:76), their adaptation of Noda and Lapusta (2010) to space-domain
boundary integrals. The tables stay as they are; the history stores
midpoint slip rates V^(n-1/2) = (V^(n-1) + V^n)/2, the first lag uses the
end-of-step kernel without its instantaneous part (K1 + eta I), and the
radiation damping acts on the slip rate at the collocation time:
  tau^n = H + (K1 + eta I) (V^(n-1) + V^n)/2 - eta V^n,
slip advances with the trapezoid rule, the slip-weakening strength uses
the slip at t_n, the aging-law state the mean rate over the step. We solve
this implicitly (their corrector iterated to convergence). They report
second-order convergence (1.5 in some mode II cases) and stability for all
modes for 0.40 <= cS dt/dx <= 1.14, which includes our 0.5.
Prototype: Stepper.enable_nl() (fixed steps, k = 1, no normal coupling);
FD_NL=1 for run_tpv205.py and run_rs_rupture.py.

TPV205-2D (slip weakening) against Dunham, arrival (V > 0.1 m/s) / peak
slip rate / final slip:

| scheme | dx | st-045 | st-075 |
|---|---|---|---|
| current | 100 m | -18.8 ms / -4.3 % / -0.014 m | -35.7 ms / -6.6 % / -0.012 m |
| current | 50 m | -11.6 ms / -2.8 % / -0.009 m | -14.1 ms / -4.8 % / -0.006 m |
| NL | 100 m | +39.0 ms / +1.3 % / +0.003 m | +22.0 ms / -4.2 % / -0.003 m |
| NL | 50 m | +24.5 ms / +2.2 % / +0.003 m | +14.8 ms / -2.4 % / -0.001 m |
| NL | 25 m | +10.1 ms / +0.6 % / +0.002 m | +4.0 ms / -0.8 % / -0.0001 m |

With slip weakening the arrival error is not reduced at 100 and 50 m
(cohesive-zone resolution in space dominates), but slip and peak slip
rate are much better.

Rate-and-state rupture (run_rs_rupture.py: in-plane, 20 km, aging law,
a = 0.008, b = 0.012, L = 0.01 m, sigma 50 MPa, steady state at 1e-9 m/s
with an 8 MPa overstressed 2 km patch), against NL at 12.5 m (arrival
interpolated in log V):

| scheme | dx | x = 4 km | x = 7 km |
|---|---|---|---|
| current | 100 m | +18.7 ms / -7.7 % / +0.006 m | +33.7 ms / -16.2 % / +0.017 m |
| current | 50 m | +6.3 ms / -5.9 % / +0.002 m | +11.7 ms / -11.2 % / +0.007 m |
| current | 25 m | +1.9 ms / -3.6 % / -0.000 m | +3.8 ms / -6.6 % / +0.002 m |
| NL | 100 m | +13.1 ms / +3.5 % / +0.019 m | +26.3 ms / +9.2 % / +0.023 m |
| NL | 50 m | +2.7 ms / +1.0 % / +0.008 m | +7.4 ms / +1.7 % / +0.009 m |
| NL | 25 m | +0.6 ms / +0.3 % / +0.002 m | +1.6 ms / +0.7 % / +0.003 m |

NL arrivals and peak slip rates converge at about second order (peak
slip rate error at 25 m 0.3 to 0.7 percent against 3.6 to 6.6 percent);
the reference is itself NL, so its own error is not in these numbers.
C, 3-D (2026-10-24): fd3d.c with nl = 1 (rsf_solve -fd_nl 1; fd3d_test
with FD3D_NL=1). The 3 x 3 first-lag stencil carries half the neighbour
terms, the V_old half goes into the known part, the self coefficient
becomes kdiag/2 - eta/2 (kdiag includes the mirror self-image), the
history receives the step's midpoint rate and slip advances with it. The
mirror check against the explicit doubled fault holds in this scheme too
(4.6e-11, in run_tests.sh). First-order remains the default.

TPV3 (slip weakening, stations PI and PA; reference about 2.97 / 2.86 s,
3.12 / 2.06 m): at 250 m first order 2.851 / 2.670 s, second order
3.067 / 2.959 s; final slip 3.111 / 2.062 m (first order 3.122 / 2.070).
As in 2-D, slip weakening arrivals are limited by the cohesive zone in
space; at 500 m the second-order run is far off (3.54 / 3.61 s).

BP5-FD first event at 1 km, CFL 0.5, against Lambert (BICyclE, 1 km):

| station | Lambert | first order | first order, dt / 2 | second order |
|---|---|---|---|---|
| strk+00dp+00, V > 0.1 m/s | 16.74 s | 17.47 s | 16.96 s | 16.74 s |
| peak slip rate | 2.492 m/s | 2.262 | 2.407 | 2.492 |
| slip at 39 s | 6.690 m | 6.550 | 6.612 | 6.646 |
| strk+16dp+00, V > 0.1 m/s | 20.02 s | 20.64 s | 20.14 s | 19.92 s |
| peak slip rate | 4.072 m/s | 3.421 | 3.755 | 3.988 |
| slip at 39 s | 6.380 m | 6.217 | 6.331 | 6.419 |
| strk-24dp+10, peak slip rate | 1.147 m/s | 1.071 | 1.113 | 1.137 |
| slip at 39 s | 6.569 m | 6.614 | 6.636 | 6.657 |

The second-order scheme at the original step matches Lambert's arrivals
to 0.1 s or better and peak slip rates to 2 percent, better than the
first-order scheme at half the step, at the same cost per step.

Variable substeps (prototype, 2026-10-24): with a step of k substeps the
slip rate varies linearly from V_old to V across the step, substep j
carries (1 - w_j) V_old + w_j V with w_j = (j - 1/2)/k, and
tau = H + C_k V_old + B_k V with B_k = sum_l w(l) K[l] + (w(1) - 1) eta I,
C_k = sum_l (1 - w(l)) K[l] + (1 - w(1)) eta I, w(l) = (k - l + 1/2)/k.
Checks: for k = 1 these are the NL operators exactly; B_k + C_k equals the
first-order cumulative operator K_k to 6e-17 (constant rate). Rate-and-
state rupture at 50 m with k fixed during the whole rupture (a stress
test; the step selector uses k = 1 at high slip rate): the first-order
scheme blows up for k = 2 and 4, NL stays stable (slip error -0.036 m at
k = 2, -0.24 m at k = 4; peak slip rate within 4 percent).

Stability (stab_random.py, after Noda et al. 2020): traction-free in-plane
fault of 64 elements with random initial traction, growth rate of ||V||
for c_p t/dx from 160 to 320 (negative is stable):

| cS dt/dx | 0.3 | 0.4 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|---|
| first order | +0.034 | +0.014 | -0.010 | +0.027 | +0.022 | +0.015 | +0.006 |
| NL | +0.032 | -0.002 | -0.010 | -0.001 | -0.009 | -0.003 | -0.010 |

The first-order scheme with end-of-step collocation is stable only near
our 0.5, an island in mode II (cf. Tada and Madariaga 2001); NL is stable
from 0.4 to at least 2.0, consistent with Noda et al. (0.40 to 1.14 in
their periodic setting). Our first-order runs were stable because they
sit on that island and because rate-and-state friction adds positive
instantaneous rate dependence.

C, 2-D (2026-10-24): fd_march.c with fd_state.nl = 1 (rsf_solve -fd_nl 1,
fd_dynamic with FD_NL=1). fd_sys_ops_nl builds B_k and C_k (band storage
for B_k, dense C_k), fd_sys_push_lin pushes the k substep midpoint rates
from the spatial FFTs of V_old and V (two FFTs per source and step),
fd_step_slipweak and fd_step_rsf branch to the NL variants. TPV205 at
100 m equals the prototype's NL run to 1e-10 in slip and 2e-10 in slip
rate (ref_tpv205_nl_dx100.txt, in run_tests.sh).
Not yet: the coarse history (an NL step suspends it, so slow phases cost
more: the sequence test's phases take 1.8 times the wall time), and
normal-stress coupling (refused with a message; needed for the stepover
work).
Sequence test (test_fd_sequence.sh setup, 100 m, variable substeps):
first order / NL phase 1 2018 / 2019 steps, peak slip rate 2.473 / 2.828
m/s, mean slip 0.9836 / 0.9899 m, recurrence 82.23 / 82.58 yr. The higher
peak slip rate is in line with the TPV205, rate-and-state rupture and
BP5-FD comparisons, where the first-order scheme underestimates it.

Coarse history for NL (2026-10-25): with a rate linear within each
aligned block of Kc substeps, the block's contribution is
A[j] s_1 + (B[j] - w_1 A[j]) (s_Kc - s_1) Kc/(Kc - 1), with the first and
last slices s_1, s_Kc of the block, w_1 = 1/(2 Kc), and the first moments
B[j] = sum_l W[l] ((j+1) Kc - l + 1/2)/Kc built next to A[j]
(fd_sys_set_coarse). On the sequence test the NL run with the coarse sum
equals the run with FD_COARSE_OFF=1 to all printed digits (phase times,
steps, peak slip rates, slip) and takes the first-order wall time again
(5.8 s against 11.5 s per phase).

Normal-stress coupling in NL (2026-10-25): sigma = sigma0 + nsign (Hn +
Cn V_old + Bn V) in both step functions; the refusal in rsf_solve is gone.
Harris and Day case D, W = +500 m (fd_dynamic hd, FD_NL=1), maximum final
slip on segment 2:

| dx | first order | NL |
|---|---|---|
| 100 m | 0.249 m | 1.244 m |
| 50 m | 0.234 m | 0.567 m |
| 25 m | 0.231 m | 0.293 m |

Both converge to an aborted jump (about 0.23 m), but NL at 100 m and 50 m
produces a sustained rupture on segment 2 that is not there at 25 m; the
first-order scheme is already within 8 percent at 100 m for this
quantity. Without normal coupling (FD_NOCOUP=1) both give no jump. The
traces agree until segment 2 starts slipping; the normal stress shows no
oscillation. Whether this is a remaining error in the NL normal coupling
or a genuine resolution sensitivity of near-threshold triggering is open;
it matters directly for the jump statistics. Until it is understood NL is
not the default, and stepover runs should use the first-order scheme or
check both.

NL normal coupling, checks (2026-10-25):
- Operators: Bn + Cn equals the first-order normal operator to 1e-16
  (k = 1 and 4), as B + C does for shear.
- Prototype: Stepper.enable_nl() now includes normal coupling
  (sigma = sigma0 - (Hn + Cn V_old + Bn V)); dump_hd.py takes FD_NL=1. The
  C run equals the prototype on Harris-Day W = 500 m at 100 m to 5e-11 in
  slip and 2e-10 in slip rate (ref_hd_nl_W500_dx100.txt, in run_tests.sh),
  so the over-triggering at 100 and 50 m is a property of the scheme as
  formulated, not an implementation error in C.
- W = -500 m (the other sense, a sustained jump in both schemes): maximum
  slip on segment 2 1.365 / 1.327 m (first order / NL) at 100 m,
  1.381 / 1.364 m at 50 m; the schemes agree.
- W = +500 m at 12.5 m does not fit in the sandbox (tables of several GB);
  to run on a workstation:
  fd_dynamic hd 500 12.5 and FD_NL=1 fd_dynamic hd 500 12.5.
Reading: in this near-threshold case (segment 2 nucleates or not) NL
needs finer elements than the first-order scheme to give the converged
answer; where a jump clearly happens or clearly does not, the schemes
agree. For jump statistics, resolution checks are needed with either
scheme; NL stays opt-in.

Harris-Day W = +500 m, maximum final slip on segment 2 (12.5 m: Thor):

| dx | first order | NL |
|---|---|---|
| 100 m | 0.249 m | 1.244 m |
| 50 m | 0.234 m | 0.567 m |
| 25 m | 0.231 m | 0.293 m |
| 12.5 m | 0.226 m | 0.225 m |

The converged answer is an aborted jump (about 0.225 m). The first-order
scheme is within 11 percent at 100 m; NL gives the wrong qualitative
answer at 100 and 50 m and converges only at 12.5 m. For near-threshold
triggering at the resolutions used for sequences, the first-order scheme
is the safer choice.

BP1-FD at 50 m (2026-10-25), switch times of the dynamic phases, fault-wide
peak slip rate and mean slip per event; first order at 25 m from Thor's
earlier run; group 5 is the densest SEAS FD submission in ref/
(onset at 7.5 km):

| | event 2 | event 3 | peak slip rate, events 1-3 | mean slip, events 1-3 |
|---|---|---|---|---|
| first order, 50 m | 301.69 yr | 420.48 yr | 12.54, 10.81, 11.69 m/s | 1.561, 1.303, 1.461 m |
| first order, 25 m | 302.08 yr | 420.86 yr | 14.17, 12.04, 13.01 m/s | |
| NL, 50 m | 302.55 yr | 421.42 yr | 15.68, 13.13, 14.28 m/s | 1.565, 1.316, 1.464 m |
| group 5 | 302.21 yr | 420.93 yr | | |

Extrapolating the first-order runs at first order (50 and 25 m) gives
about 302.47 and 421.24 yr and peak slip rates of about 15.8, 13.3 and
14.3 m/s: NL at 50 m sits at that limit. Peak slip rates at the stations
(2.5, 7.5, 12.5 km), events 1-3:

| | event 1 | event 2 | event 3 |
|---|---|---|---|
| first order, 50 m | 4.08, 3.30, 1.43 | 2.89, 2.58, 1.18 | 3.51, 3.38, 1.32 |
| NL, 50 m | 4.91, 3.88, 2.07 | 3.62, 2.84, 1.75 | 4.28, 3.99, 1.90 |
| group 5 | 6.08, 4.50, 2.39 | 4.25, 3.41, 2.09 | 4.83, 4.41, 2.24 |

NL is closer to group 5 in all nine values. Wall time per phase in the
sandbox: 38 to 44 s (NL) against 31 to 35 s (first order).

Summary of the two schemes: NL gives smooth quantities (event times,
peak slip rates, arrivals) at about the accuracy of the first-order
scheme at half the element size or better (BP1-FD, BP5-FD, TPV205,
rate-and-state rupture), but at coarse resolution it can mispredict
near-threshold triggering (Harris-Day W = +500 m), where the first-order
scheme is robust. Default stays first order; -fd_nl 1 with resolution
checks for rupture-front details.

BP1-FD with NL at 25 m (Thor): switch times 196.364, 302.628, 421.568 yr,
peak slip rates 16.24, 13.58, 14.73 m/s, mean slip 1.566, 1.316, 1.465 m;
from 50 to 25 m NL changes event times by +0.08 and +0.14 yr and peak slip
rates by 3.2 to 3.6 percent, against +0.39 and +0.38 yr and 11 to 13
percent for the first-order scheme. Wall time 300 to 340 s per phase,
about 25 percent more than first order.

Harris-Day sweep (2026-10-25, tend 10 s instead of 8 s; maximum slip on
segment 2 and triggering time, first V > 0.1 m/s on segment 2):

| W | dx | first order | NL |
|---|---|---|---|
| +200 m | 100 / 50 m | 2.219 m, 5.48 s / 2.231 m, 5.45 s | 2.176 m, 5.60 s / 2.205 m, 5.51 s |
| +300 m | 100 / 50 m | 2.101 m, 5.51 s / 2.109 m, 5.46 s | 2.059 m, 5.61 s / 2.087 m, 5.52 s |
| +400 m | 100 / 50 m | 2.006 m, 5.66 s / 2.017 m, 5.70 s | 1.960 m, 5.63 s / 1.987 m, 5.54 s |
| +500 m | 100 / 50 m | 1.423 m, 7.72 s / 1.407 m, 7.73 s | 1.877 m, 5.99 s / 1.600 m, 7.26 s |
| +700, +1000 m | 100 / 50 m | no jump | no jump |
| -500 m | 100 / 50 m | 1.974 m, 5.85 s / 1.980 m, 5.79 s | 1.958 m, 5.93 s / 1.966 m, 5.83 s |
| -1000 to -2500 m | 100 / 50 m | no jump | no jump |

Correction of the earlier reading: at W = +500 m both schemes jump; the
8 s runs ended just as the delayed first-order jump started (7.7 s), so
"0.23 m, aborted" was an artefact of the window. The schemes agree on jump
or no jump at every width and both resolutions; they differ in the delay
near the threshold, where NL at 100 m triggers 1.7 s early (6.0 against
about 7.7 s) and at 50 m 0.5 s early, while the first-order delay is
already stable at 100 m. Away from the threshold the triggering times
change by 0.03 to 0.10 s from 100 to 50 m in both schemes.

Default (2026-10-25): -fd_nl 1. NL is more accurate for event times, peak
slip rates and arrivals in every benchmark; the only measured
disadvantage is the near-threshold triggering delay at coarse
resolution, which needs a resolution check with either scheme.
test_fd_sequence.sh passes -fd_nl 0 to keep its comparison with the
first-order prototype run; the reference numbers in run_bp1_fd.sh are
first-order values.

### 3-D phases: speed (2026-10-26)

Variable substeps in 3-D (fd3d_step_k, -fd_kmax3 1, experimental): steps of
k = 2..64 fine substeps for both time schemes; known part for a step of k
substeps, operators B_k, C_k (or K_k) in the half spectrum with their self
coefficients, the history receives k slices per step (linear in time for
NL). The long-step operator reaches far beyond the 3 x 3 stencil, so the
implicit solve iterates (damped Jacobi, each sweep applies the operator by
FFT and solves locally with the self coefficient); a step that does not
converge or whose peak slip rate more than doubles is rejected and retried
with half the substeps. The step selector adds two 3-D limits to the 2-D
slip-distance criterion: at most 10 percent relative change of the peak
slip rate per step (from the previous step), and at most a doubling of k
per step. Check: forcing k = 1 through the iterative path reproduces the
stencil path on TPV101 (6 s) to the output precision (1e-5 m).
Result on the BP5-FD first event at 2 km: without the rejection test the
run took 475 instead of 12,607 steps (30 s instead of 360 s) with
identical arrivals, peak slip rates within 1.5 percent and slip at 60 s to
1e-4 m, but one long step in the afterslip produced a spurious solution
(0.46 m of slip at the hypocentral station, negative shear stress); with
the rejection test the iteration fails for k >= 2 near the arrested
rupture edge (cells near x2 = 27-29 km, 5 km depth: the damped iteration
settles on a spurious fast solution even at peak slip rates of 5e-4 m/s),
so the run falls back to single steps. Not usable yet; open: a robust
solver for the long-step system (Newton-Krylov with the FFT operator, or
an implicit near field with an explicit far field), and the coarse
history in 3-D.

Quasi-steady hand-back in 3-D (-fd_back_eps, as in 2-D): BP5-FD first
event at 2 km, -fd_back_eps 0.05 against the full phase (both NL): phase
806 s instead of 3639 s, wall 91 s instead of 399 s; slip at the stations
480 and 640 afterwards (from the QD continuation) 6.8625 / 6.8626,
6.8896 / 6.8898, 6.9173 / 6.9172 m and 6.5957 / 6.5958, 6.6135 / 6.6138,
6.6331 / 6.6333 m at 5000, 10000 and 20000 s. This is the speed-up to use
in 3-D for now.
The 3-D phase files fd_phase_<n>_st<i>.dat now have a fifth column, the
number of substeps of the step, and the phase summary lists how many steps
used each k and how many were retried.

### Free surface for dipping faults (started 2026-10-26)

Literature: Ando (2016, GJI 207, FDP=H) reviews the options (mirror image
for vertical strike-slip faults; numerically integrated half-space Green's
functions, Zhang and Chen 2006; a virtual traction-free fault representing
the ground surface, Ando and Okuyama 2010). Hok and Fukuyama (2011, GJI
184) give the method in detail: surface elements carry three displacement
components, zero traction in all components, the same kernels as the
fault; validated against Okada (static) and a discrete-wavenumber code
(waveforms); surface elements finer than the fault elements near the
fault (0.5 against 1 km in the static test); a surface patch about 2.5
times the fault length, graded; a growing checkerboard oscillation on the
surface elements, which they suppressed with an ad hoc spatial filter.
Plan: static part from interact's exact half-space operator, the virtual
surface for the dynamic correction only; second-order marching (stable in
our traction-free test, cS dt/dx 0.4 to 2) in the hope of avoiding the
filter; first in 2-D in-plane (dipping fault plus surface line), then 3-D.

Step 1 (2026-10-26): opening kernels for 2-D in-plane. gen_inplane_open.py
generates the point opening source from the same Green's function as the
shear source (moment density m_xx = lambda t, m_yy = (lambda + 2 mu) t);
inplane_open.py integrates it over a segment with the complement identity
on the segment (whole-line integral = plane P wave). test_inplane_open.py:
early-time self term equals the plane P wave (sxx = -lambda/(2 alpha),
syy = -rho alpha/2) to 1e-9; static limit against Crouch and Starfield's
constant displacement discontinuity to 2e-5 at four receivers (their
D = u(-) - u(+) is the opposite convention; the existing shear kernel
shows the same overall sign relation to 1e-9).
Step 2 (2026-10-26), static test of the virtual surface in 2-D
(fs_static2d.py): a fault of 10 km (20 elements of 500 m) dipping 30 or 60
degrees, plus a traction-free surface line of displacement discontinuities
(shear and opening) in a full plane, against the exact half-plane solution
(interact's tdd_coeff, Crouch and Starfield, built as a shared library;
the full-plane conventions agree to 1e-16). Relative error of the fault's
static traction operator (shear and normal):

| fault top | surface half-length | h_s = 1000 m | 500 m | 250 m | 125 m |
|---|---|---|---|---|---|
| 1 km | 10 km | 1.6e-3 | 8.6e-4 | 8.7e-4 | |
| 1 km | 20 km | 1.4e-3 | 2.6e-4 | 2.6e-4 | |
| 1 km | 40 km | 1.4e-3 | 6.8e-5 | 6.8e-5 | |
| 0.25 km | 10-40 km | 4.3e-2 | 2.6e-2 | 0.8e-3 to 1.3e-3 | |
| 0 (30 deg) | 20 km | | 1.5e-1 | 4.2e-2 | 6.4e-3 |
| 0 (60 deg) | 20 km | | 3.8e-2 | 5.3e-3 | |

(no surface at all: 1.7e-2, 5.7e-2, 1.3e-1, 9.7e-2 for the four fault
geometries). Rules: surface elements no larger than the depth of the fault
top (else no better than no surface); for faults that reach the surface,
surface elements much finer than the fault elements near the trace (a
quarter of the fault element for below 1 percent); a surface half-length
of two fault lengths gives about 1e-4. One degenerate case (60 degrees,
125 m) evaluates a receiver exactly on an element end and returns NaN: the
junction of fault trace and surface needs staggered element positions.
This matches Hok and Fukuyama's need for fine surface elements near the
fault and suggests a refined strip along the trace with a coarser outer
surface, as separate segments with their own tables.
Next: dense tables for non-parallel segment pairs with shear and opening
sources, a "surface" segment type (zero traction in both components), and
the dynamic test against the exact half-plane static limit and the mirror
image (vertical fault) before a dipping fault.

### Free surface, step 3: dynamics in 2-D (2026-10-27)

densesys.py: a dense fault-plus-surface system for the prototype (no
Toeplitz structure; a few hundred unknowns): fault elements with one slip
component, surface elements with slip and (in-plane) opening and zero
shear and normal traction; tables W[l] = K[l] - Ks dt for every pair
(opening sources by quadrature, cached by relative position), both time
schemes (instantaneous coefficients mu/(2 cS) for shear, rho cP/2 for
opening), the surface eliminated at each step (Schur complement), slip
weakening with normal-stress coupling on the fault, or prescribed fault
slip rates. test_fs_dyn2d.py runs the checks below (quick version about 2
min: python3 test_fs_dyn2d.py).
Fixes found by the tests: the opening kernel cancelled catastrophically
when a front crossed a segment end exactly at a sample time (spikes of
1e13 in the tables); tables now use the mean of the values at t -+ 1e-4 dt
(segment_ramp_stress_open(tshift=...)). The closed-form static opening
field has a removable 0/0 for receivers on the normal through an element
end; fs_static2d.stress_global averages over a tiny offset there.

1. Antiplane, vertical fault reaching the surface, prescribed rupture
   running up from the bottom at 0.8 cS, against the exact image solution:
   maximum fault traction error over 4 s 3.1e-3, 2.1e-3, 8.9e-4 for
   surface elements of 250, 125, 62.5 m (surface +-8 km); 9.3e-3 with
   +-4 km (edge reflections); stable over 8 s.
2. In-plane, buried 30 degree fault (top 500 m), prescribed slip for 1 s:
   before the surface reflection arrives the fault traction equals the
   full-space one to roundoff; at 8 s it is within 3.2e-3 (surface +-6 km)
   or 1.9e-3 (+-12 km) of the exact half-plane static traction of the
   final slip (full space: 7.8e-2). The first-order scheme blows up on the
   surface elements (surface slip rates of 1e4), the instability Hok and
   Fukuyama suppressed with a filter; NL stays stable without one.
3. Spontaneous rupture (slip weakening, normal coupling), 6 km fault
   dipping 30 degrees, top at 400 m:
   - thrust without / with surface: slip at the top 0.42 / 0.92 m, mid
     fault 1.32 / 2.18 m, arrival at the top 3.12 / 2.71 s, peak slip rate
     2.9 / 6.9 m/s; normal sense with surface: 0.42 / 1.57 m, 3.58 s,
     2.2 m/s (the full space is the same for both senses): the expected
     thrust-normal asymmetry of free-surface effects (Oglesby et al. 1998);
   - normal sense, 200 against 100 m elements: final slip within 0.086 m
     of 1.58 m (fine result averaged onto the coarse elements), arrival at
     the top 3.58 against 3.75 s.
4. Thrust reaching the surface (3 km, 200 m fault elements): blew up at
   first, for two reasons now fixed in the test setup: the time step was
   set by the finer surface elements, which put the fault elements at
   cS dt/dx = 0.125 to 0.25, below NL's stability range (0.4 to 2), and a
   surface collocation point sat on the fault tip. With dt from the fault
   elements (surface elements then at 0.5 hf/hs, 1 for hs = hf/2 and 2
   for hf/4, inside the range) and the trace on a surface element
   boundary: stable; slip at the trace 2.38 / 2.25 m for hs = 100 / 50 m;
   with 100 m fault elements (hs 50 m) final slip within 0.089 m of
   2.25 m, arrival at the top 1.27 against 1.33 s, peak slip rate 14 /
   19 m/s.
Rules so far: NL (no filter needed); time step from the fault elements;
surface elements between a quarter and a half of the fault element size
near the trace; the fault trace on a surface element boundary.
Next: the design for the C code and 3-D (a refined strip along the trace
plus a coarser outer surface, cross-plane tables, statics from interact's
exact half-space operator).

### Free surface in 3-D: design (2026-10-27)

Structure. Fault and surface are planes that share the strike direction
(x1), so every coupling is translation invariant along strike: tables are
Toeplitz along strike (1-D FFT, half spectrum) for all pairs, and Toeplitz
in the second direction only within a uniform grid (fault-fault, 2-D FFT as
now). Across strike the surface needs non-uniform rows (a fine strip along
the trace, coarser outside: the 2-D rules, h_s a quarter to a half of h_f
near the trace, a patch of about two fault widths either side), so
surface-surface and fault-surface couplings are dense in (fault row,
surface row) and (surface row, surface row) per along-strike wavenumber.
Surface elements are rectangles (h_f along strike, h_s across), which the
Tada quadrant kernels allow; the along-strike lattice is that of the fault.
Kernels: the full stress tensor at off-plane receivers for slip along x1
and x2 and for opening (Tada 2006, A18-A27; opening, L_ij/3, not yet
ported). Statics: interact's exact Okada half-space operator for the fault,
as with the mirror; the surface supplies only the dynamic correction.
Time step from the fault elements (2-D rule); NL marching (no filter).

Vertical strike-slip faults (BP5): the mirror image already makes both
shear tractions on the surface vanish; what remains is the normal traction.
Opening elements on the surface plane, in a full space, change no shear
traction on that plane and are symmetric under the mirror, so a surface
layer with opening only, added to the mirror, should give the exact free
surface at a third of the surface unknowns. To verify statically against
Okada before any dynamics.

Memory (fs3d_cost.py; single-precision half spectra, the window of the
farthest pair, which ends exactly in 3-D once the S wave has passed):

| case | fault | surface rows | tables |
|---|---|---|---|
| BP5-FD, mirror + opening-only layer | 100 x 40 km, 1 km | 56 (16 at 250 m) | about 8 GB |
| BP5-FD, full virtual surface | same | same | about 29 GB |
| dipping thrust 60 x 30 km, 30 deg | 1 km | 63 (24 at 250 m) | about 11 GB |
| same | 500 m | 127 (48 at 125 m) | about 180 GB |
| dipping thrust 20 x 10 km, 30 deg | 250 m | 85 (32 at 62.5 m) | about 27 GB |

The dense across-strike surface blocks dominate. 1 km problems fit a
workstation; finer dipping-fault problems need compression (candidates:
low-rank or H-matrix blocks per wavenumber and lag, coarser outer rows with
a taper, or shorter windows for far pairs).

Plan:
A. BP5 with mirror plus opening-only layer: port Tada's opening kernels
   (Python first, checked against a late-time slope and the 2-D analogue),
   static check of fault traction against interact's half-space Is (the
   mirror alone deviates 1.8 to 2.5 percent), then the dynamic tables in
   C and BP5-FD against Lambert and against the mirror-only run.
B. Dipping faults with the full three-component surface, at 1 km first,
   then compression if the science needs finer resolution.

### Free surface in 3-D, step A1: opening kernels and the static check (2026-10-27)

tada_quad.L_3: Tada's opening source (A24-A27; L_22/3 and L_32/3 by the
x1 <-> x2 symmetry); rect_ramp_stress(comp=3). Checks: early time at the
centre of a large rectangle gives the plane P wave (s33 = -rho cP/2,
s11 = -lambda/(2 cP)) exactly; against the independent point-source
quadrature (check_tada_quad.py with CHECK_OPENING=1, moment tensor
lambda delta + 2 mu n n) to 1e-8 to 1e-11 for s31, s32, s11, s33 after both
waves and between P and S, off-plane and in-plane. (The check script had
omitted the moment-tensor weights, harmless for the shear source whose
entries are 1; fixed, and the shear check still agrees to 5e-8.)

Static check (fs3d_static.py; all blocks from interact's Okada routine,
dc3dmod.F built as a shared library, half space as the reference): a
vertical strike-slip fault reaching the surface, 4 x 2 km of 1 km
elements, strike traction operator against the exact half space:

| construction | relative error |
|---|---|
| full space | 8.5e-2 |
| mirror image | 2.2e-2 |
| mirror + opening layer, surface 1000 x 250 m | 2.0e-2 |
| mirror + opening layer, surface 500 x 250 m | 4.8e-3 |
| mirror + opening layer, surface 500 x 125 m | 5.2e-3 |
| mirror + opening layer, surface 250 x 250 m | 4.5e-4 |

(surface +-4 km across strike, 2 km beyond the fault ends; doubling the
extent changes nothing at 1000 m along strike.) The mirror plus a layer of
opening-only surface elements converges to the exact free surface, as the
symmetry argument says, but the surface elements must be fine in both
directions near the trace: about a quarter of the fault element, along
strike as well. For the Toeplitz design this means an along-strike lattice
of h_f/4 for the surface (fault columns every fourth lattice point), which
multiplies the surface table memory by about four compared with the
estimate above, or a separate fine strip along the trace.
Next: refine the memory estimate for that lattice, then the dynamic
tables in the prototype (vertical fault, mirror plus opening layer) and a
comparison with the mirror-only BP5-FD first event at 2 km.

### Free surface: literature check (2026-10-27)

What exists, so that we do not reinvent it:
- Virtual traction-free surface elements in a full-space space-time BIEM:
  Ando and Okuyama (2010, GRL 37), Hok and Fukuyama (2011, GJI 184; three
  displacement components on the surface, finer surface elements near the
  fault, a large graded patch, an ad hoc spatial filter against a growing
  checkerboard instability), Ando et al. (2017), Ando (2016, GJI 207; FDP=H
  for the cost of the time convolution in such models).
- Exact half-space Green's functions in the BIEM: Zhang and Chen (2006a,b,
  GJI), Chen and Zhang (2006); for vertical strike-slip fault systems the
  kernel splits into the full-space part and a free-surface part (Chen
  and Zhang 2010, Earthquake Science 23); closed-form Lamb's-problem
  solutions (Feng and Zhang 2018 and later) make image-type methods
  practical for simple geometries.
- Mirror image for vertical strike-slip faults: Tse and Rice (1986), Quin
  (1990), Lapusta et al. (2000), Aochi and Fukuyama (2002), Lapusta and
  Liu (2009); exact in antiplane 2-D, approximate in 3-D (normal traction
  on the surface left over). BICyclE uses it for BP5 (Jiang et al. 2022).
- Hybrid treatments: spectral BIEM for planar self-interactions plus
  H-matrix-compressed space-time kernels for interactions between
  surfaces, with the free surface as an auxiliary planar boundary whose
  self-interaction is spectral and whose coupling to the faults is
  compressed (Ciardo and Romanet 2026, arXiv:2603.00383, 2-D antiplane so
  far; they also use exact causality-based window truncation, which in 3-D
  is the exact end of the kernel once the S wave has passed, as here);
  FDPM combined with H-matrices for log-linear cost (Sato and Ando 2021);
  coupled BIEM-FEM (Ma et al. 2008-2011 style SH/P-SV hybrids) and 3-D
  hybrid FE-spectral BI (Albertini et al. 2021) put the free surface in a
  volume region next to the fault.
What does not seem to be in the literature: completing the mirror for a
vertical strike-slip fault with an opening-only surface layer (step A1
above). It is a small specialisation of the virtual-surface idea, exact by
symmetry, and the static check confirms it; we keep it as stage A.
Consequences for our design:
1. Our Toeplitz/FFT tables for fault-fault and surface-surface play the
   role of the spectral self-interactions in Ciardo and Romanet; the
   fault-surface coupling is where compression belongs. Per lag and per
   along-strike wavenumber the (fault row) x (surface row) blocks are
   small dense matrices of smooth, far-field kernels: low-rank compression
   (truncated SVD or ACA) per block and lag is the natural first step,
   before FDPM.
2. The finding that the surface must be about four times finer than the
   fault near the trace (both directions) is consistent with Hok and
   Fukuyama; a fine strip with its own lattice plus a coarser outer
   surface, coupled through compressed blocks, follows Ciardo and
   Romanet's divide-and-conquer pattern.
3. NL time marching removes the instability that forced Hok and Fukuyama
   to filter; we have not found that combination reported.

### Free surface in 3-D, step A2: dynamic prototype of mirror + opening layer (2026-10-28)

fs3d_dyn.py: vertical strike-slip fault (image tied to the fault slip) plus
a layer of opening-only rectangles on the surface plane, Tada kernels
rotated between element frames, tables cached by relative position (3122
distinct kernels for a 4 x 2 fault of 500 m and a 256-element layer of
250 m instead of 70,000 pairs), marching with densesys.Marcher.
tada_quad.rect_ramp_stress now averages over a tiny symmetric offset when
a receiver lies on the extension of an element edge (removable 0/0 in the
closed forms).
test_fs3d_dyn.py (prescribed uniform slip rate on the fault for 0.4 s):
- static limit of the tables (late-time slopes through the layer) against
  Okada's half space: 1.2e-4 (mirror alone 2.9e-3; top at 500 m);
- causality: with the fault top at 500 m the run equals the mirror-only
  run exactly until 2 top / cP (0.167 s);
- static traction of the final slip, maximum relative error after the
  pulse (exact: Okada half space):

| fault top | layer | time | with layer | mirror only |
|---|---|---|---|---|
| 500 m | 250 m | 3.5 s | 6.1e-4 (NL), 4.2e-4 (first order) | 9.5e-3 |
| 500 m | 250 m | 8 s | 5.5e-4 (both) | 9.5e-3 |
| 0 (reaching the surface) | 250 m | 3.5 s | 1.05e-3 | 2.5e-2 |
| 0 | 125 m | 3.5 s | 7.3e-4 | 2.5e-2 |

- stable over 8 s in both schemes (surface opening rates decay to 2e-4 with
  NL, 1e-3 with first order; no filter).
So the opening layer reduces the free-surface error of the mirror by 15 to
35 times dynamically as well as statically, and is stable.
Cost outlook for BP5 at 1 km (Toeplitz along strike on a 250 m lattice,
fault columns every fourth point; layer rows of 250 m within +-2 km of the
trace and graded to +-20 km beyond): cross tables about 1025 wavenumbers
x 40 fault rows x 40 layer rows x 2 directions per lag, about 26 MB, over
about 400 lags, 10 GB, plus about 5 GB for the layer's self tables, so
about 15 GB before compression. Next: measure the compressibility of the
cross blocks, then the C implementation (the 3-D marcher with a second
plane of unknowns).

### Free surface in 3-D, step A3: compressibility, C kernels (2026-10-28)

Compressibility of the fault-to-layer tables (fs3d_compress.py): strike
slip on a vertical 20 x 10 km fault of 1 km (plus image) to the normal
traction on a layer with a 250 m along-strike lattice and 52 rows (250 m
within +-2 km of the trace, 1 km beyond, +-20 km), 94 lags, FFT along
strike: the (10 fault rows) x (52 layer rows) block of each wavenumber and
lag has numerical rank 5.8 / 7.8 / 8.8 / 9.0 on average (of at most 10) at
relative tolerances 1e-3 / 1e-4 / 1e-5 / 1e-6, i.e. storage ratios of
0.63 / 0.84 / 0.91 / 0.92. Per-block low-rank compression after the FFT
along strike gains little: the blocks are thin (the fault has few rows).
Real-space hierarchical compression per lag (Ciardo and Romanet 2026) or
causal partitioning (Ando 2016; Sato and Ando 2021) would be the
alternatives; for BP5 at 1 km the uncompressed estimate (about 15 GB)
fits a large workstation, so stage A proceeds uncompressed.

C kernels: fd3d_rect_ramp_tensor gives the full stress tensor (11, 22,
33, 12, 31, 32) for slip along x1 or opening on a rectangle a1 x a2 (port
of tada_quad L_1, L_3). fd3d_test kernels compares it with tada_quad at
400 random points (both sources, rectangles of 250 and 500 m, receivers
off and in the source plane, off-plane receivers on edge extensions):
agreement to 4e-11 relative to each point's largest component (in
run_tests.sh).
Conditioning, found while building that test: in the source plane the
quadrant fields are singular along each element edge line and cancel only
between quadrants. Receivers on such a line outside the element are
handled by a symmetric offset of 1e-3 a (1e-7 a off the plane), but
receivers within about 1e-2 a of a line, or within 1e-1 a of the source
plane near an edge extension, lose digits (Python and C, with or without
fused multiply-add, differ by up to 1e-4 there). Grid rule to avoid them:
coplanar elements on a common lattice, with layer row widths odd
multiples of the finest width, so receivers (odd multiples of half the
finest width) never fall on edge lines (multiples of the finest width).
Next: the layer in the C marcher (tables: layer-layer 2-D Toeplitz,
fault-layer Toeplitz along strike with dense rows; the layer's opening
rates solved each step together with the fault; fault statics from the
exact Is, the layer's static coupling to the fault left out to avoid
counting it twice).

### Free surface in 3-D, step A4: the layer in the C marcher, dense version (2026-10-28)

fd3d.c: the single fine step is split into step_begin / step_solve /
step_end (same arithmetic as before: the TPV3, TPV101 and mirror checks
and the BP5-FD 2 km first event are unchanged: 2791 steps over 805.7 s,
peak 3.814 m/s, mean slip 1.9378 m with -fd_back_eps 0.05). A layer
(struct fd3d_layer, fd3d_layer_create) attached to a mirrored 3-D fault
adds opening-only rectangles on the surface plane; each step iterates the
fault solve (with the layer's first-lag traction as an extra known part)
and a direct solve of the layer's own first-lag system (LU), both schemes,
until the rates of both change by less than 1e-12 relative. With the exact
static operator (static_fn) the layer's static coupling to the fault is
left out (it is part of Is). Tables are dense in this first version:
memory L x (layer elements)^2, so only small problems (tests).
fd3d_test layer: a 2 x 1 km vertical fault of 500 m reaching the surface
with a 16 x 16 layer of 250 m, prescribed slip rate 1 for 0.4 s, NL,
against fs3d_dyn.py: fault traction history equal to 8e-5 Pa of 1.6e7
(5e-12 relative), in run_tests.sh. QD bit-identity harness: PASS.
Next: the FFT-structured layer for production sizes (layer-layer 2-D
Toeplitz with an iterative first-lag solve; fault-layer Toeplitz along
strike with dense rows), checked against this dense version, then the
rsf_solve option and BP5-FD at 2 km with and without the layer.

### Free surface in 3-D, step A5: a narrow FFT-structured layer in rsf_solve (2026-10-28)

Shrinking the layer. Static check (fs3d_static.py, 8 x 4 km fault of 1 km
reaching the surface, layer elements of 500 m): a strip of only +-1 / +-2 /
+-4 km across the trace gives 1.7e-3 / 1.0e-3 / 8.1e-4 against Okada's half
space (mirror alone 1.5e-2). The correction lives near the trace, so the
layer is a narrow strip: a few fault elements either side, elements of
dx/2 along and across strike.
FFT structure (fd3d_layer_create_fft): the strip on a lattice of dx/m
along strike; layer-layer tables 2-D Toeplitz; fault-layer tables
Toeplitz along strike (the fault's columns every m-th lattice point, an
offset c0 between the two lattices) with dense (fault row) x (layer row)
blocks per wavenumber; histories as spectra; the layer's first-lag system
solved by Jacobi with its self coefficient inside the outer fault-layer
iteration. Equal to the dense version and to fs3d_dyn.py on the 2 x 1 km
test to 5e-12 (fd3d_test layer with FD3D_LAYER_FFT=1, in run_tests.sh), 11
times faster there.
rsf_solve: -fd_layer n (half-width in fault elements; 0, the default, is
off), -fd_layer_m (refinement, default 2), -fd_layer_margin (fault elements
beyond the ends, default n); only for a mirrored 3-D fault (half space,
vertical strike slip reaching the surface). The layer restarts from rest at
each phase.
Memory: BP5 at 2 km with -fd_layer 2: 108 x 8 layer elements of 1 km, 119
lags, about 90 MB; estimate for BP5 at 1 km: about 0.65 GB (the
fault-layer blocks dominate), instead of tens of GB for a full surface.
BP5-FD first event at 2 km, mirror only against mirror plus layer
(-fd_back_eps 0.05): phase 2791 steps over 805.7 s in both; peak slip rate
3.814 / 3.836 m/s, mean slip 1.9378 / 1.9386 m; stations 480 (x2 = 0,
near the surface), 244 (hypocentre), 640 (x2 = 16 km): arrivals
identical, peak slip rates 2.026 / 2.035, 1.236 / 1.236, 3.693 / 3.686
m/s, slip at 60 s 6.6883 / 6.6891, 6.5504 / 6.5522, 6.5160 / 6.5162 m.
So for BP5 (vertical strike slip, constant normal stress) the mirror's
normal-traction error changes the first event by less than 0.5 percent:
the mirror used by BICyclE (and by default here) is adequate for BP5; the
layer matters where the normal traction acts on friction (dipping faults,
normal-stress coupling) and is the basis for stage B. Wall time with the
layer 253 s against 76 s per phase (the outer iteration repeats FFTs;
not optimised yet). QD bit-identity harness: PASS.

### Free surface in 3-D, step A6: layer speed (2026-10-28)

Profile of the BP5-FD first event at 2 km with -fd_layer 2: 3.6 outer
(fault <-> layer) and 25 inner (layer Jacobi) iterations per step; time in
the known part 62 s, fault solves 160 s, layer solves 8 s, finish 1 s
(mirror only: 76 s in total). The repeated full fault solves dominate.
Changes: the outer iteration stops as soon as the new layer rates would
change the fault's known part by less than 1e-12 of the normal stress
(the fault solve just done is then the converged one), and later fault
solves start from the previous outer iterate. Result: 3.3 outer
iterations, fault solves 101 s, phase wall time 173 s instead of 236 s
(2.3 times the mirror-only run), station values unchanged (arrivals
identical, peak slip rates and slips to the printed digits). The phase
summary now lists the layer's iteration counts and times. (A restricted
re-solve of only the cells near the trace was tried and dropped: with rate
and state the state predictor couples all cells within a step, so it is
not exact.)

### Free surface for dipping faults in 3-D, step B1: prototype (2026-10-29)

fs3d_dip.py: a planar fault with strike along x and dip delta (frame e1 =
strike, e2 = up dip, e3 = normal; dipping towards -y), slip along strike or
dip, and a full virtual surface on z = 0 with three unknowns per element
(slip along x and y, opening) and zero traction (sxz, syz, szz) at the
element centres; fault receivers carry the shear traction along the slip
and the normal traction (normal-stress changes). Tada kernels for all
three source components (L_2 by the x1 <-> x2 symmetry), rotated between
the frames, tables cached by relative position (18,157 distinct kernels for
an 8-element fault and 840 surface elements), marching with
densesys.Marcher (NL), normal traction on the fault tracked.
test_fs3d_dip.py (2 x 1 km fault of 500 m, 30 degrees, top at 500 m, dip
slip, surface elements of 250 m over +-4 km around the fault's projection):
- consistency: the full-space operators agree with Okada's full space to
  1.8e-5;
- static operators through the surface against Okada's half space: shear
  1.5e-4 (full space 2.8e-2), normal 6.2e-3 (the full space has none:
  error 1); with +-2 km the errors are 6.3e-4 and 2.1e-2, with 125 m
  elements at +-2 km 6.3e-4 and 2.0e-2: the surface extent limits the
  normal operator, not the element size (for this buried fault);
- prescribed slip pulse, NL: equal to the full-space run before 2 top / cP;
  at 3.5 s the shear traction of the final slip is within 2.3e-4 and the
  normal traction within 4.8e-3 of Okada's half space (full space 0.12 for
  shear); surface rates decay (stable, no filter).
Next (C, stage B): the surface with three unknowns in the FFT-structured
layer code (coplanar shear-shear 2 x 2 plus opening-normal tables;
fault <-> surface tables for one fault component and two fault tractions),
L_2 in C, normal-stress changes in the 3-D marcher (exact static In plus the
surface's dynamic part), the rsf_solve half-space branch for dipping
faults (statics from interact's Is and In, gated), checks against this
prototype at roundoff, then spontaneous thrust and normal ruptures
(asymmetry, resolution) and a dipping-fault benchmark (SCEC TPV10 if the
station data are available).

### Free surface for dipping faults in 3-D, step B2: C and rsf_solve (2026-10-30)

C (fd3d.c): slip along the second grid axis in the kernels (fd3d_rect_ramp_tensor comp 2;
fd3d_test kernels: 600 points, all three sources, 1.5e-8 relative) and in the
fault-fault tables (fd3d_tab_create_c, dip slip with axis 1 along strike); the
full virtual surface (fd3d_surf_create, fd3d_step_surf): three unknowns per surface
element on a lattice of dx/m, surface-surface tables 2-D Toeplitz (shear-shear 2 x 2
plus opening-normal), fault <-> surface tables Toeplitz along strike with dense
(fault row) x (surface row) blocks for the fault's slip and its shear and normal
traction, the outer fault-surface iteration as for the layer, normal-stress changes
on the fault (sigma = sigma0 - tn; the exact static part from a callback), a sign
for the marcher's slip direction. fd3d_test surf: 2 x 1 km thrust of 500 m dipping
30 degrees, top at 500 m, 840 surface elements, prescribed slip pulse: equal to
fs3d_dip.py to 6e-12 (shear) and 5e-12 (normal traction) over 48 steps.
A latent bug found on the way and fixed (layer and surface): the fault->surface and
surface->fault Toeplitz tables were filled from one loop over the union of both
offset ranges, which can exceed the FFT length and alias; each direction is now
filled over its own range. The earlier layer results are unaffected (their union
fitted: 215 offsets in 256 for BP5 at 2 km; 31 in 32 for the layer test).
Stability (fd3d_test surfstab): the surface alone, fault at rest, random initial
rates, NL, 12 x 6 fault of 1 km dipping 30 degrees: with surface elements of dx/2
(cS dt/hs = 1) the two shear components grow (max rate 0.5 -> 1.2e-2 after 300 steps
here, divergent in a larger patch: a factor 10 per 100 steps); dx/1 (0.5) and dx/4
(2) decay (dx/4: 0.5 -> 1.1e-7), dx/3 (1.5) decays slowly; the first-order scheme
also grows at dx/2. This differs from 2-D, where NL was stable at 1 in the
traction-free test; not understood yet. The surface default is dx/4, which also
meets the static accuracy rule.
rsf_solve: a planar fault in a half space that is not vertical strike slip (or any
fault with -fd_surface 1) uses the full surface: grid axis 1 along strike, rows by
depth, slip along strike or dip, the surface frame from the fault geometry (strike,
up dip, sign of the slip direction relative to interact's t_strike and t_dip), the
surface over the fault's projection plus -fd_surface_ext fault elements (default the
fault width), the fault trace on a surface row boundary; statics from the exact Is,
and with -calc_sigma_dot the exact static normal traction -In delta plus the
surface's dynamic part, the normal stress handed back to QD at the phase end.
Test: 12 x 6 km thrust of 1 km, dipping 30 degrees, top at 500 m, rate and state
(velocity weakening interior, seeded patch), -calc_sigma_dot: with surface elements
of dx/2 the run diverged after about 100 s (slip of kilometres, the instability
above, amplified by the normal-stress coupling); with dx/4 (72 x 45 surface elements,
about 214 MB) the phase runs 731 steps over 105.5 s, peak slip rate 1.5 m/s, mean
slip 0.053 m, minimum normal stress 48.9 MPa, surface rates decaying to 4e-6 m/s.
Without the surface (full space) the same seed gives 1.23 m/s and 0.040 m.
QD bit-identity harness PASS; BP5-FD 2 km (mirror) unchanged.
Next: thrust against normal faulting and a resolution check in rsf_solve, a fault
reaching the surface, the stability question at cS dt/hs = 1, and a dipping-fault
benchmark (TPV10-type).

### Dipping faults in 3-D, step B3: surface stability band, first physics checks (2026-10-30)

Stability band of the surface (fd3d_test-style random test, surface alone, NL, 12 x 6
fault of 1 km, surface elements dx/2, varying the step): cS dt/hs = 0.8 decays (4e-5
after 300 steps), 0.94 grows strongly (rates of order 7 after 300 steps), 1.0 grows,
1.06 decays slowly (4e-4); 0.5 and 2.0 decay strongly, 1.5 slowly. So there is an
unstable band for the surface's shear components at roughly 0.85 to 1.05, a genuine
property of the scheme, not an artefact of wave fronts meeting sample times exactly.
rsf_solve now refines the surface automatically to keep cS dt/hs >= 2 (dx/4 at the
default CFL 0.5, dx/8 at 0.25) unless -fd_surface_m is given, and warns if a chosen
refinement falls in 0.8 to 1.2.
Physics checks with rsf_solve (12 x 6 km fault of 1 km dipping 30 degrees, rate and
state with a seeded patch, -calc_sigma_dot, first dynamic phase):

| case | peak slip rate | mean slip | phase length |
|---|---|---|---|
| full space (either sense) | 1.232 m/s | 0.040 m | 7.1 s |
| half space, thrust, top at 500 m | 1.499 m/s | 0.053 m | 105.5 s |
| half space, normal sense, top at 500 m | 1.254 m/s | 0.040 m | 7.9 s |
| half space, thrust reaching the surface | 1.573 m/s | 0.056 m | 31.5 s |

The free surface promotes thrust slip and hardly changes the normal sense, the
asymmetry of the 2-D prototype and of Oglesby et al. (1998); the thrust reaching the
surface runs stably (trace on a surface row boundary). test_fd3d_dip.sh (8 x 4 km,
about 15 s) checks the branch is chosen, the phase is stable and thrust slip exceeds
the full space (1.612 against 1.232 m/s, 0.1025 against 0.0901 m).
Still to do: matching published numbers (TPV10-type, when the SCEC site is back),
resolution studies, and understanding the band (3-D shear on a plane combines in-plane
and antiplane behaviour; a 2-D antiplane traction-free test would be the next check).

### Dipping faults in 3-D, step B4: why the surface is unstable near cS dt/hs = 1 (2026-10-30)

Literature. Space-time BIEMs are known for "intermittent instabilities" that appear
and disappear as the time step to element size ratio changes (Peirce and Siebrits
1997, IJNME 40, with a z-transform stability analysis of time-stepping schemes for
displacement-discontinuity boundary elements and the epsilon and half-step schemes
designed from it; Siebrits and Peirce 1997; Birgisson, Siebrits and Peirce 1999).
Tada and Madariaga (2001) map the 2-D stability of collocation schemes; Noda, Sato
and Kurihara (2020) give the 2-D ranges of the NL scheme (modes II and III
separately); Hok and Fukuyama (2011) report a growing checkerboard on their 3-D
free-surface elements with the conventional scheme and filter it.
Tests (all with the second-order NL scheme unless stated):
1. 2-D traction-free random tests (stab_random.py, FD_MODE=inplane or antiplane):
   stable at cS dt/dx = 0.5, 0.94, 1.0, 1.06 in both modes.
2. 3-D planar fault alone, one slip component, fd3d tables (3 x 3 first-lag stencil):
   stable at 0.5, 0.94, 1.0, 1.06, 2.0.
3. 3-D surface alone: with the shear cross coupling (x-slip <-> y-slip) switched off,
   a single shear component on the surface still grows at 0.94 (0.5 -> 7.6 in 400
   steps), so the coupling between the two shear components is not the cause.
4. The growing field is a pure checkerboard: sign changes between every neighbour
   along strike and across (dominant wavenumber (pi, pi); fd3d_test surfstab fields).
5. Per-mode analysis (fd3d_test surfsym writes the Fourier symbols Ks dt and W(l) of the
   x-slip -> xz-traction tables on a grid of wavenumbers; fs3d_surf_stability.py
   iterates the traction-free recursion of each mode and reports the growth per step):

| cS dt/hs | worst NL mode, growth per step | checkerboard NL | first-order scheme |
|---|---|---|---|
| 0.5 | (3/8, 0) 1.009 | 0.982 | 1.019 at the checkerboard |
| 0.8 | (3/8, 0) 1.014 | 0.991 | divergent |
| 0.94 | checkerboard 1.049 | 1.049 | divergent |
| 1.0 | checkerboard 1.040 | 1.040 | divergent |
| 1.06 | (1, 3/4) 1.015 | 1.012 | divergent |
| 1.5 | (7/8, 0) 1.031 | 0.926 | 1.050 |
| 2.0 | 0.997 (all decay) | 0.898 | 1.049 |

(periodic-lattice analysis with a 60-lag window: growth factors within about 1 percent
of 1 are not resolved, and a finite patch radiates energy, which is why the mild values
at 0.5 and 0.8 are not seen in the runs.)
Reading: in 2-D each grid-scale mode is purely in-plane (mode II, along the slip) or
purely antiplane (mode III, across it), and NL is stable for both near a ratio of 1.
In 3-D the surface carries the checkerboard mode (pi, pi), which has no 2-D
counterpart: its stress combines the mode II response of the neighbours along the
slip, the mode III response of those across, and the P-wave coupling of the diagonal
neighbours, which enters the first lag once cP dt exceeds about 1.4 hs (ratio about
0.8). In the band 0.85 to 1.05 that combination gives the second-order recursion a
root outside the unit circle; near 1.5 the along-slip sawtooth is marginal, and at 2
every mode decays. The first-order scheme is unstable over most of the range, which
is the instability Hok and Fukuyama filtered. Consequences: the default ratio of 2
(dx/4 at CFL 0.5), the warning near 1, and the two tools for checking other
settings. Possible later improvements, from Peirce and Siebrits: a half-step or
epsilon-type collocation for the surface elements, which they designed to remove such
intermittent instabilities.

### Dipping faults in 3-D, step B5: cross-check of the two free-surface treatments (2026-10-31)

BP5-FD first event, vertical strike-slip fault reaching the surface, the same fault
resolution, three treatments: mirror image, mirror plus opening layer (-fd_layer 2),
full virtual surface (-fd_surface 1, surface elements dx/4, -fd_surface_ext 2):

| | 4 km: hypocentre arrival / peak / slip at 100 s | shallow station slip at 100 s | 2 km: x2 = 0 arrival / peak / slip at 60 s | x2 = 16 km |
|---|---|---|---|---|
| mirror | 2.31 s / 1.511 m/s / 5.42 m | 0.234 m | 15.88 s / 2.03 m/s / 6.69 m | 19.34 s / 3.69 m/s / 6.52 m |
| mirror + layer | 2.31 / 1.510 / 5.43 | 0.225 | 15.88 / 2.04 / 6.69 | 19.34 / 3.69 / 6.52 |
| full surface | 1.73 / 1.503 / 5.51 | 0.407 | 15.30 / 2.23 / 6.52 | 18.48 / 4.13 / 6.39 |

Mirror and mirror plus layer agree (the layer completes the mirror exactly, so the
mirror's normal-traction error is small for BP5). The full surface ruptures faster
near the surface (0.6 to 0.9 s earlier at 2 km), with about 10 percent higher peak
slip rates and about twice the shallow slip at 4 km, and the gap is the same at 4 and
2 km. Tests of the causes (4 km, early part of the phase):
- surface margin 5 instead of 2 fault elements: no change (1.73 s, 1.504 m/s, 5.50 m);
- surface elements dx/8 instead of dx/4: hypocentre unchanged, shallow slip larger
  still (0.43 and 0.82 m against 0.41 and 0.65): not converging towards the mirror;
- the tables' own statics instead of Is (-fd_table_statics 1, a consistent but
  truncated system): the mirror's arrival and peak (2.31 s, 1.513 m/s) but slip badly
  off (5.94 m), so the exact statics are needed and are not the cause;
- cost: 915 s of wall time per 2 km phase against 76 s for the mirror.
Diagnosis: a fault element edge lying in the surface. In the full space that edge
carries the dislocation's singular field, which the piecewise-constant surface elements
must cancel; they converge slowly. Static check (fs3d_dip.py, vertical strike-slip
2 x 1 km fault of 500 m against Okada): reaching the surface, the shear operator error
is 6.4e-3 with 250 m surface elements and 4.0e-3 with 125 m; buried 500 m it is 1.7e-3
with 250 m. The mirror has no such edge (fault plus image is one dislocation through
the surface), and mirror plus layer converged to 4.5e-4 at dx/4. A tied mirror image
of a dipping fault does not help (fs_static2d.py, image=True: 30 degree fault reaching
the surface, 1.3e-2 at 125 m against 6.4e-3 without the image): the image makes a kink
at the surface, which is singular again.
Consequences: for vertical strike-slip faults use the mirror (plus layer if normal
traction matters). The full surface is accurate for buried faults (steps B1, B2: 1e-3
to 1e-4) and converges slowly for faults that reach the surface; its results there
need a resolution study of the surface rows next to the trace, which needs graded
rows (fine next to the trace only; the Toeplitz-along-strike structure allows dense
rows of unequal width). TPV10-type benchmarks (surface-breaking normal faults) will
test exactly this, so that is the next development before matching numbers.

### Dipping faults in 3-D, step B6: independent dynamic check against the 2-D prototype (2026-10-31)

Quasi-2-D comparison: a long thrust dipping 30 degrees (24 km along strike, 6 km down
dip, 1 km squares, top at 500 m), slip weakening (0.677 / 0.525 / 0.40 m), tau0 70 MPa
with 81.6 MPa on the bottom 2 km along the whole strike, sigma 120 MPa, full virtual
surface (elements of 250 m, margin 2 km), NL, 6 s; the central column of the 3-D run
against the 2-D plane-strain prototype (densesys, surface elements 250 m, same
friction and loading), element by element, top to bottom:

| | 3-D central column: slip (m) | 2-D: slip (m) | rupture times (s), identical in both |
|---|---|---|---|
| full space | 1.151 1.519 1.735 1.827 2.064 1.628 | 1.151 1.515 1.724 1.818 2.054 1.621 | 2.454 2.021 1.588 0.866 0.144 0.144 |
| surface, no normal coupling | 2.168 2.481 2.651 2.622 2.667 1.987 | 2.169 2.480 2.642 2.613 2.662 1.983 | 2.454 1.876 1.588 0.866 0.144 0.144 |
| surface, normal coupling, thrust | 2.907 3.104 3.112 2.956 2.891 2.125 | 2.907 3.104 3.108 2.955 2.884 2.122 | 2.165 1.876 1.443 0.866 0.144 0.144 |

Agreement to 0.3 to 0.6 percent in slip (the remainder is the 3-D along-strike edge
effect at 24 km) with identical rupture times, in all three cases: the free surface
doubles the shallow slip, and the normal-stress coupling adds another third for the
thrust. The 2-D prototype was itself checked against exact statics (tdd_coeff), the
image solution (antiplane) and causality, so this is an independent check of the 3-D
dipping-fault machinery (fault tables for dip slip, surface tables for three
components, normal-stress changes) with the free surface active. A first attempt with
the opposite slip sense in 3-D gave less slip with coupling (the normal-faulting
sense is clamped), which is the expected asymmetry, not an error.
fd3d_test q2d (16 km strike, 24 s) is in run_tests.sh: 1.3 percent in slip against the
2-D reference (ref_fd3d_q2d.txt, from test_fs_dyn2d.inplane_rupture), rupture times
equal.
So the 3-D dipping fault with the full surface works for buried faults; the open issue
is the convergence at the trace for faults that reach the surface (step B5), and the
TPV10-type comparison when the data are available.

### Dipping faults in 3-D, step B7: graded surface rows, fault-plane extension (2026-10-25)

Graded rows (fd3d_surf_create_rows): the surface-surface tables are Toeplitz along
strike with dense (row, row) blocks, so rows can have any widths; uniform rows are
detected and their blocks computed once per row offset. rsf_solve builds the rows
(-fd_surface_grade 1, default) with widths of at most half the distance from the surface
point to the nearest part of the fault's cross-section (vertical depth over the
projection, distance to the top or bottom edge line outside it), powers of two times
dx/m between dx/m and dx, at most doubling from row to row. A first version graded by
distance from the trace alone was wrong for dipping faults: on the quasi-2-D buried
thrust it put 1 km rows over a fault 0.5 to 3 km below the surface and the top element's
slip came out 5 percent low (2.763 against 2.906 m) with a shifted rupture time; with
the distance rule the regression thrust (test_fd3d_dip.sh) gives 1.609 m/s and
0.1024 m with 17 graded rows against 1.612 m/s and 0.1025 m with 30 uniform rows.
Extension rows (fd3d_create3, -fd_surface_extend, default 1 when the fault's top edge
lies in the surface): the fault plane continues above the surface by rows whose slip
rate copies the top row, with the mirror image's mechanism (table rows nz + next, the
extension in the first-lag stencil and self coefficient, and as a tied source in the
fault-to-surface tables), so the dislocation crosses the surface without an edge. 2-D
static test (fs_static2d.py, 30 degree fault reaching the surface): error 7.8e-2 /
2.4e-2 / 2.6e-3 for surface elements of 500 / 250 / 125 m against 1.5e-1 / 4.2e-2 /
6.4e-3 without it.
BP5 at 4 km (early phase, hypocentre slip at 100 s / shallow station slip at 100 s):
mirror 5.42 / 0.23 m; uniform surface rows dx/4 5.51 / 0.41; distance-graded rows with
the extension 5.51 / 0.43 (rows of 1 and 2 km, 139 MB instead of 187 MB). The trace
error of a surface-breaking fault is therefore not removed by the extension or by
finer surface rows at the fault's resolution of 4 km; it shrank between 4 and 2 km
(step B5), so it is governed by the fault element size at the trace.
State after this step: buried dipping faults validated; surface-breaking dipping
faults need a fault-resolution study at the trace (and the TPV10 data) before numbers
can be matched. QD harness, run_tests.sh (12 checks), test_fd_sequence.sh,
test_fd3d_gate.sh, test_fd3d_dip.sh: PASS.
Sandbox note: the sandbox was reset on 2026-10-24 and rebuilt from SANDBOX.md (PETSc
3.25.5 from the pip sdist, the sandbox makefile configuration); the graded-row code
written just before the reset was re-implemented.

### Dipping faults in 3-D, step B8: SCEC TPV10, first comparison (2026-10-25)

Reference data (provided by Thor, SCEC site): Kaneko, Lapusta and Ampuero, SPECFEM3D,
100 m elements, on-fault stations st000dp000 (surface, along-strike centre) and
st000dp120 (12 km down dip) plus two off-fault stations; Dunham, WENORUP-2D (a 2-D,
plane-strain version), 18.75 m, on-fault at the surface and 4.5 km down dip, off-fault
at +-0.5 km. fully_dynamic/tpv10_compare.py reads the SCEC HTML files directly.
Code additions (fd3d.c): slip weakening with cohesion (f->coh) and tension treated as
zero normal stress (strength coh + mu max(sigma, 0)); weakening by the accumulated slip
path (f->path, f->dpath; TPV10: backward motion adds to it); a saturated branch of the
local solve for steps that weaken faster than radiation damping can balance (always on,
it only acts where the linearised solve is ill posed) and optionally (f->sw_saturate)
whenever the solution passes d0 within a step (exact for the piecewise-linear
strength; off by default because the prototype references use the linearised solve).
Driver: fd3d_test tpv10 <dx> <strike length> <q2d> <surface margin> <prefix>: the TPV10
geometry and stresses, the full virtual surface (uniform rows, elements of dx: cS dt/hs
= 0.5, stable), the extension row above the trace, normal-stress coupling, window from
the geometry.
Results (centre column; rupture time = slip rate first above 1 mm/s; slip in m):

| station | reference | rupture | peak rate | slip at 5 / 8 / 12 / 15 s |
|---|---|---|---|---|
| surface | Kaneko 3-D, 100 m | 2.320 s | 1.760 m/s | 1.103 / 2.397 / 3.001 / 3.072 |
| | ours 3-D, 1 km | 2.576 | 0.480 | 0.132 / 0.476 / 0.529 / 0.740 |
| | ours 3-D, 500 m | 2.576 | 1.114 | 0.092 / 1.779 / 2.809 / 3.140 |
| 12 km | Kaneko 3-D, 100 m | 0.005 | 6.379 | 3.756 / 4.451 / 4.857 / 4.863 |
| | ours 3-D, 1 km | 0.152 | 6.129 | 3.728 / 3.728 / 3.980 / 4.082 |
| | ours 3-D, 500 m | 0.076 | 5.967 | 3.296 / 3.801 / 4.370 / 4.484 |
| surface | Dunham 2-D, 18.75 m | 2.265 | 6.390 | 3.821 / 5.367 / ... |
| | ours quasi-2-D (nucleation along strike, 30 km), 1 km | 2.273 | 2.918 | 2.767 / 2.999 / ... |
| 4.5 km | Dunham 2-D | 1.901 | 5.184 | 4.666 / 5.955 |
| | ours quasi-2-D, 1 km | 1.818 | 2.712 | 4.993 / 5.160 |

Reading: the machinery runs the full TPV10 problem (60 degree normal fault reaching the
surface, tension at shallow depth, normal-stress changes from the surface) and
converges strongly from 1 km to 500 m: final surface slip 0.74 -> 3.14 m against 3.07,
slip at 12 km 4.08 -> 4.48 against 4.86. At 500 m the surface slip arrives later than
in SPECFEM3D (1.8 m against 2.4 m at 8 s) and the peak slip rate is 37 percent low; the
cohesive zone (about 1 km at mid depth, less near the hypocentre) is resolved by about
two elements at 500 m, against ten at Kaneko's 100 m, so this is the expected
under-resolution, not a sign error. The quasi-2-D runs reproduce Dunham's 2-D rupture
times (2.27 against 2.27 s at the surface, 1.82 against 1.90 s at 4.5 km), but with a
30 km strike the along-strike edges limit the slip after about 5 s, so 2-D final slip
is not comparable; a 16 km strike (q2d16) is worse, as expected.
Cost and limits: 500 m with the full surface is about 2.2 GB of tables and 9 minutes in
the sandbox (7.8 outer and 72 inner iterations per step); 250 m needs a workstation
(about 4 times the elements and twice the steps; tables of order 15 to 20 GB with the
current dense row blocks; single-precision tables or a narrower surface would halve or
quarter that). Next: 250 m (and 100 m if feasible) on a workstation, the off-fault
stations, and the rupture-time contours.

### 3-D, step B9: SCEC TPV101 against Kaneko et al. (2026-10-25)

Reference (provided by Thor, SCEC site): Kaneko, Lapusta and Ampuero, SPECFEM3D, 100 m,
stations st000dp075 (hypocentre), st-090dp075 and st120dp030. (Three further uploads,
meant as Dunham's TPV101 series, were the SCEC login page; no data.)
fully_dynamic/tpv101_compare.py reads the SCEC files; fd3d_test tpv101 now also records
st120dp030 (the first three stations, and the prototype reference test, are unchanged)
and with FD3D_STFILES=1 writes tpv101_<station>.dat. Whole space, FD3D_SMALL geometry,
15 s; strike slip; rupture time at 1 mm/s; slip in m:

| station | run | rupture | peak rate | slip at 3 / 6 / 9 / 15 s |
|---|---|---|---|---|
| st000dp075 | Kaneko 100 m | 0.570 s | 4.347 m/s | 3.674 / 5.504 / 7.312 / 7.419 |
| | ours 250 m, first order | 0.577 | 4.110 | 3.425 / 4.599 / 6.242 / 7.250 |
| | ours 250 m, second order | 0.577 | 4.333 | 3.511 / 4.961 / 6.890 / 7.286 |
| | ours 500 m, second order | 0.650 | 4.558 | 3.227 / 3.774 / 3.989 / 5.248 |
| st090dp075 (vs -090) | Kaneko 100 m | 4.730 | 5.711 | 0 / 2.360 / 4.914 / 4.937 |
| | ours 250 m, first order | 5.810 | 2.720 | 0 / 0.186 / 3.387 / 4.651 |
| | ours 250 m, second order | 5.268 | 3.724 | 0 / 1.427 / 4.106 / 4.747 |
| st120dp030 | Kaneko 100 m | 6.280 | 8.274 | 0 / 0 / 3.712 / 3.736 |
| | ours 250 m, first order | 7.975 | 4.143 | 0 / 0 / 2.076 / 3.534 |
| | ours 250 m, second order | 7.181 | 5.024 | 0 / 0 / 3.017 / 3.605 |

Reading: at 250 m the second-order scheme gives final slip within 1.8 to 4 percent of
SPECFEM3D at all three stations, and arrivals 0.5 s (st090) and 0.9 s (st120) late
against 1.1 and 1.7 s with the first-order scheme; peak slip rates are 35 to 40 percent
low away from the hypocentre. The earlier first-order resolution series (section WP8:
arrival at (9, 7.5) km 6.32, 5.81, 5.51, 5.22, 5.02, 4.94 s from 300 to 75 m,
converging at first order to about 4.74 s) extrapolates to Kaneko's 4.73 s, which
confirms that converged solution independently. At 500 m the rupture does not
propagate properly (the aging-law cohesive zone, a few hundred metres, is not resolved).
Dunham's MDSBI (spectral boundary integral, 100 m) series at the same three stations
agree with SPECFEM3D: rupture 0.570 / 4.715 / 6.275 s against 0.570 / 4.730 / 6.280 s,
final slip 7.400 / 4.939 / 3.733 m against 7.419 / 4.937 / 3.736 m, peak slip rates
4.30 / 5.48 / 7.68 against 4.35 / 5.71 / 8.27 m/s. So the two references agree to 15 ms
in arrival and 0.3 percent in final slip, while peak slip rates still differ by 1 to
7 percent between them at 100 m; our comparison should judge peak rates with that in
mind. MDSBI is the method closest to ours (a boundary integral code), so it is the
natural reference once our resolution approaches 100 m.
Second-order run at 100 m (Thor, 2026-10-25: FD3D_NL=1 FD3D_SMALL=1 FD3D_STFILES=1
bin/fd3d_test tpv101 100 - 15; 390 x 240 elements, 1039 steps):

| station | Kaneko / Dunham (100 m) | ours 100 m | ours 250 m |
|---|---|---|---|
| (0, 7.5): rupture, peak, slip at 15 s | 0.570 / 0.570 s, 4.35 / 4.30 m/s, 7.419 / 7.400 m | 0.577 s, 4.25 m/s, 7.351 m | 0.577, 4.33, 7.286 |
| (9, 7.5) | 4.730 / 4.715 s, 5.71 / 5.48 m/s, 4.937 / 4.939 m | 4.807 s, 5.37 m/s, 4.884 m | 5.268, 3.72, 4.747 |
| (12, 3) | 6.280 / 6.275 s, 8.27 / 7.68 m/s, 3.736 / 3.733 m | 6.409 s, 7.61 m/s, 3.706 m | 7.181, 5.02, 3.605 |

At 100 m: arrivals 0.08 to 0.09 s and 0.13 s late (against 0.55 and 0.90 s at 250 m:
the errors fall by 6 and 7 for a factor 2.5 in element size, i.e. second-order
convergence, which would leave about 0.02 and 0.03 s at 50 m); final slip within 0.3
to 1.1 percent of both references; peak slip rates within 1 to 2 percent of Dunham's
MDSBI and 6 to 8 percent of SPECFEM3D (the references differ by 1 to 7 percent among
themselves). Time series, RMS difference over 0 to 15 s relative to the reference: slip
1.1 / 2.8 / 4.9 percent at the three stations (Dunham against Kaneko: 0.2 percent);
slip rate 3.5 percent at the hypocentre, 57 and 83 percent at the outer stations, where
the sharp rupture-front pulse arrives 0.10 and 0.17 s late (the best-fitting shift), so
that misfit is mostly the arrival delay. The first-order scheme at 100 m (section WP8)
had the (9, 7.5) arrival at 5.02 s; the second-order run's 4.81 s confirms the gain.
The rupture-time grid (tpv101_rt.dat) shows a 43 ms asymmetry between (-9, 7.5) and
(9, 7.5) km, from the half-element offset of the sampling points relative to the
nucleation centre.
Status: TPV101 is reproduced at 100 m to within about 0.1 s in arrival, 1 percent in
slip and a few percent in peak slip rate, converging at second order; a 50 m run would
close the remaining arrival gap. The driver's default is still the first-order scheme so
that the prototype reference test is unchanged; use FD3D_NL=1.

### Speed, step S1: profile of the 3-D step, active-set Jacobi (2026-10-26)

Profile (fd3d_test tpv101 with FD3D_PROF=1, second order, one core): at 500 m over 6 s,
per step 3.1 ms in the known part (history convolution), 56 ms in the local solves,
0.7 ms in finishing; 27.5 Jacobi sweeps per step, each solving every one of the 3744
cells (103,000 local rate-and-state solves of 7.7 Newton iterations per step). So the
local solve, not the convolution, dominated.
Active set: after the first sweep, a cell is recomputed only if its own rate or one of its
3 x 3 stencil neighbours changed by more than 1e-16 relative in the previous sweep (an
unchanged known part gives an unchanged solution; Newton restarted from a nearby guess
changes the last bits even for unchanged input, which must not keep a cell active).
FD3D_NO_ACTIVE=1 switches it off. Results: local solves per step 103,000 -> 16,300 at
500 m (solve 56 -> 7.0 ms per step), 6 times faster per step; at 250 m over 8 s per step
28.7 ms known part, 43 ms solve (28 percent of the cells active on average), 3.9 ms
finishing. The reference tests are unchanged to the printed digits (tpv3, tpv101, mirror
differences identical to before). BP5-FD at 2 km, first event: 21.2 s instead of 70 s of
wall time, 2791 steps and peak 3.814 m/s as before, mean slip 1.9380 against 1.9378 m
(the quasi-steady hand-back test is sensitive at the 1e-4 level).
MPI plan (decided 2026-10-26: MPI only, no OpenMP, to match rsf_solve's
parallelism): the history convolution is distributed by wavenumber slabs (each rank
owns a block of kx rows of the half spectrum: its part of the tables and history rings,
so memory also divides by the number of ranks); the local solves by blocks of fault
columns along strike with a one-column halo exchanged per Jacobi sweep; between the two,
a slab-decomposed 2-D FFT with one all-to-all transpose each way per step. The surface
tables are already per along-strike wavenumber and follow the same slab split. Results
must be bit-identical for any number of ranks (sums over lags per wavenumber in the same
order, local solves identical), which is the test.

### WP7, C implementation in rsf_solve: status (2026-10-18)
Done, each verified in the sandbox before delivery (details in
fd_c_design.md sections 6 to 10, usage in rsf_solve_fd.md):
- M1 source reorganisation (src/rsf_solve/), skeleton, built-in FFT with
  an FFTW-ready interface, bit-identity and timing harness.
- M2 kernels, tables, segment detection, static gate against Is/In
  (1e-16 agreement, including cross-segment normal blocks and the
  half-plane image).
- M3 marcher core, LAPACK banded solve, slip weakening, standalone
  fd_dynamic: TPV205-2D equals the prototype to 1e-10, 25 times faster.
- M4 rate-and-state, images, history seeding, switch as a TS event with
  the phase in the post-event callback, hand-back: BP1-FD three events and
  the single-fault sequence equal the prototype to the printed precision.
- M5 two parallel segments with normal-stress coupling through the phase:
  the prototype's stepover case (1 km step, 2 km overlap, 0.5 km margins)
  reproduced (segment 1 follows segment 0 after 53 min, slips per event
  within 0.3 percent); event log ordering fixed.
Harness: identity PASS on every patch; the timing column is noise at the
run lengths the sandbox allows.

Next, in order:
1. BP1 at 25 m with rsf_solve (run_bp1_fd.sh 25 520 12000, 2026-10-19):
   QD events 196.364, 271.983, 349.804, 428.102, 506.448 yr (prototype
   196.36, 271.98, 349.82, 428.12, 506.48; group 4 at 50 m 196.36, 271.97,
   349.78, 428.07, 506.40). FD phases at 196.3638, 302.1023, 420.8015 yr
   with 11077, 9784, 10554 steps and peak slip rates 13.59, 11.56,
   12.52 m/s (prototype 196.36, 302.11, 420.81 yr; 11074, 9783, 10551
   steps; 13.59, 11.56, 12.51 m/s). The C code reproduces the prototype at
   25 m. Wall time about 1 h per phase on one core, no faster than the
   Python: at 25 m the per-step cost is the history convolution (L M per
   source, here 12000 x 4096 x 2 complex products per step), which is the
   same algorithm in both. This makes M6 (partitioned convolution) the
   next code item. The unequal-segment sweep is done and pinned (WP6
   above).
1c. Coarse history for aligned long substeps (2026-10-21, fd_c_design.md
   section 15): exact; BP1 50 m phase 183 s to 50 s (312 s before M6).
   -fd_back_eps 0.05 at 25 m: 1e-6 relative in event times, 12 percent
   faster; stays opt-in.
1b. M6 delivered (2026-10-20, fd_c_design.md section 14): hybrid partitioned
   convolution (default), banded solves, single-FFT pushes, optional
   quasi-steady hand-back. BP1 50 m phase 312 s to 183 s at unchanged
   results, 110 s with -fd_back_eps 0.05. Extra tests: Harris and Day case D
   (two segments, normal coupling) in C against prototype traces to 5e-11.
   To check on the workstation: BP1 at 25 m with -fd_back_eps 0.05 against
   the default criterion, before considering it as a default.
1a. History alignment fix (2026-10-20, fd_c_design.md section 13): found
   while designing M6; prototype and C corrected together; BP1-FD at 50 m
   moves by 0.1 yr in recurrence (105.33 and 118.79 yr against 105.42 and
   118.74; sem2dpack 106.5 and 118.7).
2. M6 performance: partitioned FFT convolution over the lag axis
   (O(log L) amortised per step, factor 20 to 50 at L of a few thousand),
   then OpenMP over Fourier modes; identity to 1e-10 against M5.
3. M7 scope decisions: non-parallel segment pairs (dense lag tables, needed
   for branches and TPV14/15-2D), VE memory carried through a phase
   (-ve_mode with -fd_mode), MPI distribution of the tables, slip law and
   flash heating in the C friction (prototype has them), in-plane half
   plane, 3-D elements.
4. Small items: pre-phase onset line in rsf_events.dat is done; station
   output of the QD stations (-rsf_stations) during phases; checkpoint
   that includes the FD record; a 3-D case in the harness.

### WP7, C implementation in rsf_solve (original plan)
Design in fd_c_design.md (2026-10-17): monitor-attached FD phase, no
changes to the QD code path, bit-identity and timing harness run at every
step, milestones M1 to M7 with their acceptance tests. M1 (source
reorganisation, skeleton, FFT, harness) delivered 2026-10-17. M2 (kernels,
tables, segment detection, static gate against Is; C tables agree with
the prototype dumps to 1e-12) delivered 2026-10-17.
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
