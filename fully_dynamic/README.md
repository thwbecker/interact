# fdkern: in-plane elastodynamic kernels for the FD prototype (WP1)

Python prototype of the 2-D plane-strain (mode II) space-time BIEM
planned in fd_plan.md. Status: WP1 (kernels), WP2 (windowed tables, history
convolution) and WP3 (fixed-step dynamic stepper, slip-weakening and
rate-and-state friction) done and tested against exact solutions and the
SCEC TPV205-2D benchmark; WP4 (adaptive QD integrator and hand-off) and the
rate-and-state sequence machinery (WP4) in place and tested on a single
fault; the stepover ensembles (WP6) are next.

Files
- gen_inplane.py: derivation and verification. Builds the 2-D full-space
  Green's function by integrating the 3-D one along strike, double
  time-integrates it (ramp response), forms the point slip-patch field via
  the equivalent double couple, and generates inplane_kernel.py. Running it
  performs three checks: the hand antiderivatives (second time derivative
  returns the impulse response), the Navier equation for the 2-D Green's
  function in both arrival regions, and a direct numerical z-integration of
  the 3-D Green's function.
- inplane_kernel.py: generated code (do not edit). sigma_point_ramp(x,y,t,
  alpha,beta,rho) gives (sxx,sxy,syy) for a point slip patch of unit length
  and unit slip rate created at t=0. Carries the opposite sign to the
  physical field; the segment routines account for that.
- inplane.py: segment kernels with two engines. "closed" (default,
  production): the Tada and Madariaga (2001, IJNME 50, eqs A4-A6) primitives,
  identical to Ando et al. (2007, eqs 4-6), transcribed from the typeset
  pages; about 2.7e6 lag entries per second. "quad" (cross-check):
  quadrature of the independently derived point field along the source
  segment (Gauss-Legendre on sub-intervals split at the P and S front
  crossings, with a cubic endpoint map that removes the inverse square root
  singularities at the fronts; on-segment receivers use a complement
  identity, whole-cone integral equals the plane S-wave field, to avoid the
  hypersingular integrand); about 2e4 entries per second. Also: static
  Hirth-Lothe references, rotation and traction projection for arbitrarily
  oriented segments, lag kernels K^(l) and a tabulation routine.
- antiplane.py: antiplane (mode III) point field built with the same
  pipeline (2-D SH Green's function, double time integral, double couple)
  and segment routines reusing inplane.segment_ramp_stress.
- ref_ando2007.py: the closed-form primitives (in-plane: Tada and
  Madariaga 2001 eqs A4-A6 = Ando et al. 2007 eqs 4-6; antiplane: Ando et
  al. 2007 eqs 30-31 = Tada and Madariaga eqs A2-A3). Both reproduce the
  static dislocation limits to 1e-7 and serve as the production engines.
- test_antiplane.py: plane SH-wave identity, static limits, reference
  transcription check, entry-by-entry comparison with Ando et al. (2007):
  median relative difference below 1e-10 and maximum below 1e-6 away from
  the fronts, below 1e-3 within 2 percent of an arrival.
- gen_series.py / inplane_series.py: large-time expansion of the in-plane
  primitives (t^-1, t^-3, t^-5 terms; the t^3 terms are shown to cancel and
  the t term to equal the static primitive). Used for the dynamic correction
  W once beta t / r > 10, where the closed form suffers cancellation.
- shells.py (WP2): Fault geometry, static operators Ks/Kn, per-pair windowed
  tables of W[l] = K[l] - Ks dt stored as FFTs along the Toeplitz offset
  axis (parallel segment pairs), slip-rate history ring buffers, and the
  history convolution (O(L M) complex multiply-adds per pair and step). The
  window per pair ends where |W| drops below tol times the static
  interaction envelope mu dx^2 dt / (2 pi (1-nu) r^2); for 2-D kernels this
  gives L of order (1/sqrt(tol)) times the S travel time across the pair.
- stepper.py (WP3): fixed-step dynamic time marching with end-of-step
  collocation (e_t = 1, beta dt = dx/2). The first-lag operator K[1]
  (tridiagonal for a straight fault at this time step) is solved exactly:
  directly for the prescribed-tip model, by a banded active-set solve for
  linear slip weakening; rate-and-state friction still uses a two-pass
  predictor-corrector for the neighbour terms. QD mode (radiation damping
  only) on the same code.
- spectral_ap.py: independent spectral BIEM (Perrin-Rice-Zheng kernel,
  periodic domain, trapezoid time convolution) for the antiplane
  self-similar crack; shares nothing with the space-time code.
- run_tpv205.py, ref/: TPV205-2D driver and the SCEC reference series.
- run_hd.py, test_hd.py: Harris and Day (1993) stepover experiments (two
  parallel segments, forced nucleation, normal-stress coupling).
- cycle.py (WP4): earthquake sequences. Adaptive QD integrator (scipy
  DOP853 on slip and state, slip rate from the elementwise closure with a
  bracketed log-Newton), backslip loading in shear and normal stress,
  QD-to-FD hand-off at a slip-rate threshold with the fine-grid history
  rebuilt from the QD dense output, variable FD steps of k dt_min (k a
  power of two up to 64, from xi dc / V_max) using cumulative kernels
  K_k = sum_{l<=k} K[l] as the implicit operator, FD-to-QD return after
  the slip rate has stayed below a threshold for one full window.
  A dynamic=False option runs the same sequence with radiation damping
  only.
- run_cycle1.py, test_cycle.py: single in-plane fault with a 16 km
  velocity-weakening patch (h* about 5 km, cohesive zone about 800 m,
  dx = 100 m).
- test_shells.py, test_stepper.py: table consistency, FFT versus dense
  convolution, window behaviour, kinematic self-similar crack (constant
  stress in time, zero normal stress on the plane), dynamic self-similar
  crack with prescribed tip speed against the exact Kostrov (mode III) and
  Broberg (mode II) factors, QD static limit, TPV205-2D at 100 m against the
  Dunham reference (about 20 s).
- test_inplane.py: pytest suite (both engines where applicable): static limits
  (point and segment, off-line, collinear, self), causality and arrival
  regions, plane-wave identity, self element early-time value -mu/(2 beta),
  Gauss quadrature against adaptive quadrature, lag-kernel telescoping,
  rotation invariance of tractions, reciprocity for parallel equal elements,
  convergence under time-step refinement, table-versus-single consistency.

Conventions are listed at the top of inplane.py. Slip positive when the +n
side moves along +s; stresses are stress-from-slip, so interact's Is/In
(stress-drop convention) differ by a global sign.

Observed numerical properties worth keeping in mind
- The point field diverges as (t - t_arr)^(-1/2) at the fronts; the segment
  integral is finite with a square-root onset. The quadrature map handles
  this; accuracy of segment values at front-crossing lags is about 1e-7
  relative, limited by double-precision cancellation in the generated
  expressions within about 1e-8 relative time of an arrival.
- Static limits are approached as (c t / R)^-2; at very large c t / R the
  generated expressions lose precision, so static checks use c t / R ~ 1e3.
- Closed-form tabulation: about 2.7e6 lag entries per second; a 400 by
  4000 Toeplitz table takes well under a second. The closed form loses
  relative accuracy for c t / r >> 1 (cancellation of t^3 terms between the
  P and S parts, Tada and Madariaga section 3.5); for the velocity
  formulation correction W = K - K_s dt this matters only far outside the
  1e-2 window adopted in the plan, where the asymptotic series (Ando et al.
  2007 eqs 24-26 and Appendix A) would be the remedy.

Validation status (two independent routes agree)
- Route 1: closed-form primitives from the literature (Tada and Madariaga
  2001; Ando et al. 2007), transcribed from the typeset PDFs.
- Route 2: derivation from the 3-D Green's function (gen_inplane.py,
  antiplane.py) plus segment quadrature.
- In-plane: median relative difference 4e-14 and maximum 6e-6 over 6157
  point-time combinations away from the fronts, including self and
  collinear receivers; up to a few percent of a tiny value within 0.1
  percent of an arrival (quadrature cancellation).
- Antiplane: median 6e-15, maximum below 1e-6 away from the fronts.
- Physical identities satisfied by both: static edge/screw dislocation
  limits, plane S-wave identity, exact early-time self value -mu/(2 beta),
  zero before the P wave from the segment ends (uniform slip radiates P
  only from the slip-gradient discontinuities at the ends).
- Time marching, exact solutions: the symmetrically expanding crack with
  prescribed tip speed. Mode III: Kostrov (1964) / Broberg (1999, eq.
  6.9.145), slip = (2 tau/mu) sqrt(a^2 - x^2)/E(k), k^2 = 1 - (v/beta)^2;
  exact C(v)/C(0) = 0.8257 and 0.7052 at v/beta = 0.5 and 0.8, against
  0.8265 and 0.7051 from the space-time code at dx = 100 m and 0.8258 and
  0.7042 from the independent spectral solver (spectral_ap.py). Mode II:
  Broberg (1999, eqs 6.9.85-6.9.87, elliptic-integral function g2); exact
  0.9631, 0.8651, 0.7662 at v/beta = 0.2, 0.5, 0.8, against 0.9630, 0.8651,
  0.7649 (kinematic route, converged in dx) and the stepper within 6e-3 at
  dx = 100 m.
- Time marching, community benchmark: TPV205-2D (SCEC; 2-D slice of TPV5,
  linear slip weakening, overstressed nucleation, two stress asperities)
  against E. Dunham's MDSBI solution at 6.25 m (run_tpv205.py, reference
  series in ref/). With the exact coupled solve of the first-lag system
  (see below), at dx = 100, 50, 25 m the main-front arrival (slip rate
  exceeding 0.5 m/s) at the stations -4.5, -7.5 and +12 km differs from
  the reference by +37/+12/-18 ms, +23/+4/+11 ms and +12/+4/+8 ms, i.e.
  one to three time steps; final slip differs by 0.1 to 0.4 percent at
  25 m; rms slip difference over 12 s is 0.030/0.019/0.011 m and rms shear
  stress difference 0.04 to 0.5 MPa. Peak slip rates are within 3 to 14
  percent at 25 m and converge more slowly. At +12 km a precursory creep
  appears where the stress hovers at the yield value ahead of the front
  (40, 14 and 2.5 mm at 100, 50 and 25 m; the reference stress touches
  the yield value there too but accumulates no slip).
- Harris and Day (1993) stepovers (100 m, their Table 1): case D reproduces
  their Fig. 7d (0.5 km steps jumped in both senses with trigger delays
  +0.23 s dilational and +1.73 s compressional against their +0.3 and
  +1.5 s, at the same locations; 1 and 2 km steps not triggered). Case A
  matches delays and trigger locations within 0.2 s and 0.4 km for
  dilational 1, 2, 4 km and compressional 1, 2 km steps, and runs through
  the supershear transition stably. Two differences: our case A rupture
  runs at 5.8 to 6.0 km/s (close to alpha) and crosses 14 km in 2.7 s
  against their 3.35 s; and a compressional 4 km step is triggered in ours
  (+2.3 s) but never in theirs. Both are unchanged at 50 m and are
  insensitive to the forced-nucleation radius (1 to 5 km) and speed (0.3
  to 0.8 beta), because the spontaneous supershear front overtakes the
  forced one within the first kilometre (L_c is about 200 m at S = 0.49).
  Their 250 m grid, coarser than L_c, and artificial damping are the
  likely cause; the published transition lengths (Andrews 1976b; Dunham
  2007) would settle it. Case C has the right jump pattern but triggers
  earlier than theirs; its cohesive zone is about 120 m, so 100 m is
  under-resolved and a 50 m rerun is pending.
- Sequences (run_cycle1): switch threshold 1e-3 versus 1e-2 m/s changes
  the first event's peak slip rate by 0.4 percent and mean slip by 1
  percent. Over four events the FD sequence settles to recurrence 82.9 yr,
  mean slip 1.149 m, peak slip rate 2.66 m/s, extent 22.2 km; the QD
  sequence on the same fault to 68.1 yr, 0.886 m, 1.06 m/s, 20.4 km.
  An event costs about 2000 variable steps and 30 to 40 s at N = 300.
- The first-lag coupling matters. With beta dt = dx/2 the P wave from an
  element's own ends arrives within one step, so the self term K1_ii is
  about 1.46 eta and the neighbours supply the rest (the row sum of K1 is
  exactly -eta, a discrete plane-wave identity that is tested). A two-pass
  predictor-corrector for the neighbour terms recovers only 90 percent of
  the first-step slip rate of a newly activated element, which produced a
  first-order rupture delay (0.27/0.14/0.07 s at 100/50/25 m in the first
  TPV205-2D runs). The slip-weakening solver now solves the coupled banded
  system exactly with an active-set loop (scipy solve_banded, O(N) per
  iteration); the rate-and-state solver still uses the predictor-corrector
  and should get the same treatment (banded Newton) in WP4.
- Observation: imposing the elliptical slip kinematically produces a clean
  uniform stress in antiplane (1e-4) but a grid-scale oscillation of several
  percent in-plane, consistent with Rayleigh waves excited by the discrete
  tip advance (mode II has the Rayleigh pole); the stepper, where the slip
  rate adjusts self-consistently, does not show it in the slip.

Open items for WP1
- Numerically safer evaluation near arrivals (series expansion) if the
  front accuracy turns out to matter in the stepper tests.
