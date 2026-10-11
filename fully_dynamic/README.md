# fdkern: in-plane elastodynamic kernels for the FD prototype (WP1)

Python prototype of the 2-D plane-strain (mode II) space-time BIEM
planned in fd_plan.md. Status: WP1 (kernels), WP2 (windowed tables, history
convolution) and WP3 (fixed-step dynamic stepper) done and tested; WP4
(adaptive QD integrator and hand-off) and the rate-and-state runs are next.

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
  collocation (e_t = 1, beta dt = dx/2), the first-lag self term implicit,
  neighbour coupling through K[1] (P-wave reach within one step) by
  predictor-corrector for rate-and-state friction and exactly for the linear
  prescribed-tip model; QD mode (radiation damping only) on the same code.
- spectral_ap.py: independent spectral BIEM (Perrin-Rice-Zheng kernel,
  periodic domain, trapezoid time convolution) for the antiplane
  self-similar crack; shares nothing with the space-time code.
- test_shells.py, test_stepper.py: table consistency, FFT versus dense
  convolution, window behaviour, kinematic self-similar crack (constant
  stress in time, zero normal stress on the plane), dynamic self-similar
  crack with prescribed tip speed against the spectral solver (antiplane,
  agreement 1e-3) and against the kinematic factor (in-plane), QD static
  limit.
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
- Time marching: the self-similar (Kostrov) crack with prescribed tip speed
  gives dynamic factors C(v)/C(0) = 0.826 (v = 0.5 beta) and 0.705 (0.8
  beta) in antiplane from the space-time code, against 0.826 and 0.704 from
  the independent spectral solver. In-plane, stepper and kinematic routes
  agree to 0.3 percent (shared kernels). Not yet done: TPV205-2D
  (rate-independent slip weakening) if reference data become available.
- Observation: imposing the elliptical slip kinematically produces a clean
  uniform stress in antiplane (1e-4) but a grid-scale oscillation of several
  percent in-plane, consistent with Rayleigh waves excited by the discrete
  tip advance (mode II has the Rayleigh pole); the stepper, where the slip
  rate adjusts self-consistently, does not show it in the slip.

Open items for WP1
- Numerically safer evaluation near arrivals (series expansion) if the
  front accuracy turns out to matter in the stepper tests.
