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

Next: NL with the variable substeps in 2-D (prototype, then C), a
random-perturbation stability check at cS dt/dx = 0.5, the full BP5-FD
first event and sequence with -fd_nl 1, then making it the default.

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
