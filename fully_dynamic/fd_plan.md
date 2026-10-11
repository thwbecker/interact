# Fully dynamic extension of `rsf_solve`: plan

Status: draft plan, 2026-10-10. Target application: rupture jump
probabilities across stepovers, comparing quasi-dynamic (QD) and fully
dynamic (FD) sequences on the same fault system. Scope of this
document: formulation, feasibility of fast (H-matrix) evaluation, the
Python prototype work packages, the stepover experiment design, and
open decisions. Numbers in section 3 come from the scripts
`hcheck.py` and `hcheck2.py` delivered with this note; they are
2-D antiplane results and should be re-checked for the in-plane
kernels once those are implemented.

## 1. Goal and route

Route chosen: space-time boundary integral formulation with
causality-based truncation (option 3a of the discussion). The static
operators (`Is`, `In`), friction laws, state laws, loading, and the
QD integrator of `rsf_solve` are kept; the FD extension adds a
hereditary dynamic-correction term and a fixed-step dynamic stepper
that is active only during events.

First target is 2-D plane strain (map view, in-plane, mode II) rather
than antiplane, because stepover jumping is controlled by dynamic
normal-stress change and by the stopping phase radiated at segment
ends, neither of which exists in antiplane. The antiplane kernel is
cheap to carry along for integrator tests and for comparison with the
existing antiplane cycle tests and the SEAS BP1-FD benchmark.

Known limitation of the 2-D stage: infinitely deep faults
overestimate static stress transfer across a step and remove the
seismogenic-depth saturation of rupture. Jump statistics from the 2-D
stage should be read as QD-versus-FD differences, not as absolute
probabilities. A 2.5-D or 3-D stage is needed before absolute
statements.

## 2. Formulation

### 2.1 Exact decomposition

For a receiver element i, the elastodynamic shear traction change due
to slip history on all elements j is

    tau_i(t) = sum_j sum_T K_T[i,j] V_j(t - T)  -  eta V_i(t)

with eta = mu/(2 c_s) the exact instantaneous (radiation damping)
term, K_T the kernel for unit slip rate on element j over one time
step at lag T (pre-integrated over the step), and V the slip rate.
Normal traction has the same structure with its own kernels. The QD
approximation replaces the double sum by K_s delta(t), where K_s is
the static kernel.

### 2.2 Velocity formulation and retarded-static split

Write K_T = K_s dt + W_T. Then

    sum_T K_T V(t - T) = K_s delta(t) + sum_T W_T[i,j] V_j(t - T).

Before the first arrival T_arr(i,j) = r_ij / c (c = c_s in antiplane;
c_p for the first arrival in-plane) the kernel is zero, so
W_T = -K_s dt there. Summing the pre-arrival part analytically gives

    sum_j K_s[i,j] delta_j(t)  +  sum_{T < T_arr} (-K_s dt) V_j(t - T)
      = sum_j K_s[i,j] delta_j(t - T_arr(i,j)).

Hence the exact FD traction is

    tau_i(t) =  sum_j K_s[i,j] delta_j(t - T_arr(i,j))          (A: retarded static)
              + sum_j sum_{T >= T_arr} W_T[i,j] V_j(t - T)       (B: transient)
              - eta V_i(t).                                       (C: damping)

Term A is the static kernel applied to slip delayed by the pair travel
time. Term B is nonzero only in a window after arrival: in 3-D it ends
exactly when the S front has crossed the element pair (Ando 2016); in
2-D it decays algebraically (approximately as (c t / r)^-2, section
3.1) and must be truncated. This split is the same partition as
Ando's Fast Domain Partitioning. It matters for implementation because
A and B have very different data structures (section 4).

Interseismic phase: with V of order 1e-9 m/s both A minus K_s delta
and B are negligible, so the QD integrator is exact to working
precision and no history needs to be evaluated. The FD terms are
switched on by a slip-rate threshold (section 5, WP4).

### 2.3 Conventions to fix early

- `Is` and `In` in `interact` are assembled in the stress-drop
  convention; the dynamic kernels from the literature are in the
  stress-from-slip convention. One global sign flip, gated by a static
  limit check per pair exactly as the `-ve_prony_k0tol` gate does.
- The time step of the pre-integrated kernels is a fixed
  dt_min = beta dx / c_s with beta <= 1/2. Kernel tables are tied to
  this dt_min; changing beta means re-tabulating.
- Slip rate history is stored on the dt_min grid during events;
  on entry into an event the history over the pre-event coarse steps is
  interpolated onto the fine grid (V is slow there, so this is
  accurate).

## 3. Numerical facts from the antiplane check

Geometry: two parallel straight segments of 10 km, element size
dx = L/n, parallel step width W, zero overlap or a 1 km gap, 2-D
antiplane, full space. Kernel: finite differences of the Ando et al.
(2007) primitive I_32 as reproduced in Ciardo and Romanet (2026),
pre-integrated over dt. Static kernel: analytic finite-segment screw
dislocation. Checks and results:

1. Static limit. The lag-summed dynamic kernel converges to the static
   kernel: ratio 1.0017 at t = 50 s, 1.00002 at 500 s, 1.000000 at
   5000 s for a pair at 1 km. The decay of the deviation is close to
   t^-2. This validates the primitives, the finite-difference
   construction and the sign/normalisation of (A).
2. Front amplitude. The first lag after arrival carries a value 250 to
   530 times larger than K_s dt (sqrt singularity at the S arrival,
   integrated over dt); subsequent lags fall off as 3, 1, 0.5, ... of
   K_s dt. Consequence: the kernel must be pre-integrated exactly over
   the step (the primitives do this); point sampling in time is not an
   option. The large alternating-sign front values are the reason the
   velocity formulation (which removes K_s dt) is preferable to the
   slip formulation numerically.
3. Transient window length (2-D only). After arrival, the maximum of
   |W_T| over a pair block exceeds 10 % of max |K_s dt| for about
   0.5 s (1 km step) to 1.4 s (2 km step, 1 km gap), 1 % for 1.8 to
   5.5 s, 0.1 % for 6 to 19 s. The window scales with pair distance
   (relative decay is a function of c t / r), so in 2-D the memory per
   pair is roughly proportional to r_ij / c_s divided by the square
   root of the tolerance. Truncation sensitivity must be tested
   (Ciardo and Romanet report insensitivity of peak slip rate and
   cumulative slip to windows as short as 1 % beyond the last arrival,
   for intersecting faults, slip formulation). In 3-D this issue
   disappears.
4. Compressibility, lag frame (what HTOOL or HACApK would do if handed
   the per-lag matrix W_T through a kernel callback): with a binary
   cluster tree, leaf 16, admissibility eta_H = 1.5, tolerance 1e-4
   relative to |K_s dt|, the low-rank storage over all lags is 0.41 to
   0.46 of dense, with full-rank blocks whenever the front crosses the
   block. The front is a band in (i, j) at fixed T and is intrinsically
   high rank.
5. Compressibility, arrival-shifted frame (relative lag
   s = T - T_arr(i,j), FDP=H style): storage 0.02 to 0.06 of dense,
   maximum block rank 4 to 5, median rank 1, out to s = 4 times the
   last arrival. In the shifted frame the kernel is a smooth function
   of the pair geometry and compresses as well as the static kernel.

## 4. H-matrix assessment

Conclusions for our stack (dense, HTOOL/MATHTOOL, H2OPUS, HACApK
lattice, hmmvp; all accept a kernel callback):

- The static operators stay as they are: term A uses K_s, term C is
  local. Nothing changes for the interseismic phase.
- Per-lag H-matrices built with the existing backends are not worth
  it: a factor 2 at best (fact 4), with hundreds to thousands of
  lag matrices to assemble and store. This rules out the simple
  "one HTOOL matrix per lag" design.
- The arrival-shifted frame compresses extremely well (fact 5), but
  turning that into a fast matvec requires handling pair-dependent
  delays inside low-rank blocks: y_i = sum_k U_ik sum_j V_jk
  V_j(t - s - r_ij / c), and the delay r_ij / c is not separable in
  (i, j) at the accuracy needed (within an admissible block at
  eta_H = 1.5 the first-order expansion of r_ij misses by of order
  d / (18 c), many dt_min). This is the problem the plane-wave
  time-domain (PWTD) family and Sato and Ando's FDP=H solve with
  directional expansions. It is a custom algorithm, not available in
  any backend we use, and the 3-D lattice H-matrix report (JpGU 2024)
  notes accuracy depends on the admissibility parameter. Treat it as a
  stage-3 optimisation, not a prerequisite.
- Term A (retarded static) does not need low rank: binning pairs by
  arrival lag gives sparse shells whose total storage is one dense
  N x N matrix (each pair belongs to exactly one shell). Cost per
  dynamic step is O(N^2), the same as one dense matvec.
- Term B in 3-D: a few lags per pair (P-to-S interval plus the front
  band), so sparse-shell storage of order N^2 times a small factor and
  O(N^2) work per step. In 2-D the window is longer (fact 3), which is
  why the 2-D prototype should use the Toeplitz structure below.
- 2-D prototype: for straight segments in a homogeneous full space all
  kernels depend only on (x_i - x_j) within a segment pair, so each
  shell is Toeplitz. Storage per segment pair is O(n_lags x N) and the
  convolution is O(n_lags x N log N) per step with FFTs, or O(n_lags x
  N^2) direct. This is sufficient for the stepover ensembles and is the
  same translational-invariance idea already used in the VE assembly.

Recommended staging:

| stage | problem | term A | term B | expected cost per dynamic step |
|---|---|---|---|---|
| 1 | 2-D prototype, N ~ 1e3 | Toeplitz per segment pair | Toeplitz per segment pair, truncated window | ~ 1e7 flop |
| 2 | 3-D in rsf_solve, N ~ 1e4 | sparse lag shells | sparse shells, exact window | ~ N^2 x few |
| 3 | 3-D, N >~ 1e5 | H-matrix with directional delays (FDP=H / PWTD) | same | ~ N log N x window, custom code |

Honest reading: the H-matrix approach is usable for the static part
today and for the dynamic part only through a custom shifted-frame
algorithm. For the jump-probability study in 2-D and for 3-D at BP5
resolution (N of order 1e4), stages 1 and 2 are enough and need no
new compression technology. The dynamic phase at stage 2 scales as
dx^-3 (N^2 times the number of steps per event), so 250 m resolution
in 3-D is where stage 3 becomes necessary.

## 5. Python prototype: work packages

Language and layout: standalone Python (numpy, scipy), one module per
work package, pytest for all tests, following the `bp1_pt.py` pattern.
Friction formulas copied from `rsf_engine.c` (regularized law
`v = 2 v0 exp(-psi/a) sinh(tau/(sigma a))`, aging law in psi form) so
that QD runs of the prototype can be compared against `rsf_solve`
output directly.

WP1 Kernels.
- Antiplane: Ando et al. (2007) primitives (already implemented in the
  check scripts, to be moved into the module with tests).
- In-plane (plane strain): Tada and Yamashita (1997) closed forms for
  piecewise-constant slip rate on a straight segment, shear and
  normal components, P and S arrivals. Receiver traction obtained by
  projecting the stress tensor on the receiver normal and tangent.
- Static kernels: analytic finite-segment edge dislocation (shear and
  normal), compared against `interact` 2-D in-plane `Is`, `In` on a
  test geometry (sign convention gate).
- Tests: (a) lag sum converges to the static kernel for several pair
  distances and orientations (both components); (b) causality: zero
  before c_p arrival, shear component of antiplane zero before c_s;
  (c) telescoping identity of the finite differences; (d) comparison
  with a point-source Green's function at large distance; (e)
  symmetry/reciprocity relations expected for the segment pair;
  (f) time-step consistency: results independent of beta within the
  pre-integration accuracy when the response to a smooth slip-rate
  history is compared at fixed physical time.

WP2 History data structure and convolution.
- Shell representation of A and B, Toeplitz specialisation for
  straight segments, per-pair window from a tolerance on |W_T| /
  |K_s dt| (section 3, fact 3), ring buffer of V on the dt_min grid.
- Tests: dense reference convolution versus Toeplitz/FFT path; window
  truncation error against a long-window reference for a prescribed
  slip-rate history (rupture-like pulse).

WP3 Dynamic stepper.
- Algebraic form at each dt_min: tau_i = tau_0 + loading - K_s delta
  (retarded) + B - eta V, sigma_i likewise, V inverted per element
  with the regularized law (Newton, as in `vel_from_rsf`), second
  order predictor-corrector (Lapusta et al. 2000 structure).
- Tests: single-segment dynamic rupture against TPV205-2D (SCEC
  dynamic rupture benchmark, 2-D in-plane, slip-weakening friction, many
  code solutions available); antiplane two-fault problem of Ciardo and
  Romanet (2026, parameters in their Table 1) for the antiplane path.

WP4 QD/FD hand-off.
- QD phase: adaptive explicit integrator equivalent to the `rsf_solve`
  default (RK with error control) in the tau-ODE form.
- Switch to FD when max V exceeds V_sw (to be scanned, start at
  1e-3 m/s and 1e-4 m/s); on entry, fill the fine-grid history from the
  coarse steps; switch back when max V has dropped below V_sw and the
  norm of the B contribution is below tolerance; the history buffer is
  then released.
- Tests: an FD run with V_sw set below the interseismic slip rate (FD
  always on) must agree with the switched run; a QD run of the
  prototype must reproduce the `rsf_solve` BP1-QD anchor.

WP5 Verification of the complete cycle code.
- Antiplane BP1-FD (image method for the free surface, exact in
  antiplane) against the SEAS archive solutions.
- In-plane: no community SEAS benchmark exists for 2-D in-plane FD
  cycles; verification relies on WP3 single-event benchmarks, the
  static limit, convergence in dx and beta, and QD agreement with
  `rsf_solve`.

WP6 Stepover experiments (section 6).

WP7 Design note for the C implementation in `rsf_solve`: shell
assembly from `interact` element kernels, coupling to the IMEX split
(local damping and state implicit, nonlocal history explicit), restart
of the history buffer, output of dynamic fields.

Order: WP1, WP2, WP3 with TPV205-2D, WP4, WP5 antiplane, then WP6.
WP1 is where most of the risk sits (in-plane kernels are long and
easy to get wrong); the tests in WP1 are the main deliverable of the
first iteration.

## 6. Stepover experiment design

Physics expected from the literature: dynamic jumps are promoted by
dynamic normal-stress reduction (extensional steps jump wider than
compressional ones) and by the stopping phase; over cycles the jump
probability tracks receiver prestress, which is set by prior slip and
afterslip; QD treatment and oversized cells both change the
probability while frequency-magnitude statistics and static stress
drops can still match (Lambert and Lapusta 2021). Design consequences:

- Resolution: cohesive zone (about 9 pi / 32 L_b) resolved with at
  least 5 cells on both segments; convergence check at 2x refinement
  on a subset.
- Geometry parameters: step width W (0.5 to 5 km), overlap (negative
  for gaps), extensional versus compressional sense (sign of slip
  relative to the step), segment lengths relative to L_nuc.
- Frictional parameters: a/b contrast, D_c, possibly a mild VS zone
  at segment ends as arrest control; no enhanced weakening in the
  first round.
- Pairing: QD and FD runs share loading, spin-up and parameters; the
  FD run is the QD run with A and B switched on. Events are compared
  from the spun-up cycle, not from matched initial conditions only.
- Ensemble size: tens of events per configuration; with the 2-D cost
  estimates of section 4 this is feasible on a workstation.
- Classification of a jump on the receiver segment after an event on
  the source segment: dynamic (seismic slip rate reached within the
  wave travel time plus a few seconds of source arrest), postseismic
  (within hours to days, afterslip-driven; QD captures this), or
  independent (later). Delayed dynamic triggering of several seconds
  is expected from the dynamic rupture literature and should be
  resolvable.
- Outputs per event: moment per segment, rupture speed, peak slip
  rate, prestress on the receiver at source nucleation, time of
  receiver nucleation relative to source arrest, and the receiver
  stress time series at the arrival of the stopping phase.

## 7. Risks and open questions

- In-plane kernel correctness (WP1); mitigated by the test set.
- 2-D tail truncation: window tolerance versus cost; test on the
  two-fault problem by varying the window.
- Hand-off thresholds: in a BIEM the switch only decides whether a
  small term is evaluated, so sensitivity should be weak, but it must
  be measured (WP4).
- Self-element first lag: with beta <= 1/2 the self kernel over one
  step is dominated by eta; the residual is part of B and must not be
  double counted with C.
- 2-D depth bias (section 1); candidate mitigations are a 2.5-D
  stiffness term for the static part or going to 3-D for term A only
  while keeping B in 2-D, both approximate; decide after the 2-D
  QD-versus-FD differences are known.
- 3-D free surface for stage 2: half-space static K_s with full-space
  dynamic correction misses surface-reflected phases; the consistent
  alternative (auxiliary traction-free surface) adds degrees of
  freedom. Not needed for the 2-D stage.

## 8. Decisions requested

1. Confirm in-plane (map view) as the first target and antiplane as a
   secondary test path.
2. Window tolerance for the 2-D transient term (proposal: start at
   1e-2 relative, test 1e-3 on the two-fault problem).
3. Whether to copy friction and state formulas verbatim from
   `rsf_engine.c` (recommended, for one-to-one QD comparison) or to
   implement only the aging law in the prototype.
4. Whether the stepover ensemble should include a VS arrest zone at
   segment ends or rely on free segment ends.

## 9. WP1 status (2026-10-10)

Done. In-plane and antiplane segment kernels are implemented and validated
by two independent routes that agree to machine precision away from the
wave fronts: the closed-form primitives of Tada and Madariaga (2001) /
Ando et al. (2007), and a derivation from the 3-D Green's function with
segment quadrature. The closed form is the production engine (about 3e6
lag entries per second); the derivation is kept as the cross-check.
Package: fdkern (README.md there). Changes to the plan:
- WP1 kernel tests now include an external closed-form comparison, so the
  first external check of the kernels no longer waits for TPV205-2D; that
  benchmark remains the first check of the time-marching scheme (WP3).
- Kernel tabulation cost is no longer a consideration for the 2-D stage.
- Implementation detail for the C code (WP7): the Tada-Madariaga primitives
  are the natural element kernel for straight 2-D segments, four primitive
  evaluations per (pair, lag). The 3-D stage will need the triangular
  element kernels of Tada (2006) or numerical space-time quadrature.

## 10. WP2 and WP3 status (2026-10-10, later)

Done and tested (see fdkern/README.md). Findings that touch the plan:
- The closed-form primitives need the large-time series for the dynamic
  correction W once beta t / r exceeds about 10; without it roundoff near
  nodal directions of the radiation pattern prevents the window criterion
  from closing. The series (t^-1, t^-3, t^-5) is generated symbolically and
  verified against the closed form at moderate beta t / r.
- The window criterion is defined relative to the static interaction
  envelope at the pair distance, not the static value itself, because the
  radiation pattern has nodal directions. At tol = 1e-2 the window is about
  5 S-travel times across the pair for the 2-D kernels, growing like
  1/sqrt(tol).
- A stepper benchmark that needs no external data is in place: the
  self-similar crack with prescribed tip speed, checked in antiplane against
  an independent spectral solver (agreement 1e-3 in the dynamic factor).
  TPV205-2D remains desirable for the rate-independent friction path.
- Per-step cost of the history convolution at N = 120, L = 300 is well
  below a millisecond; the stepover ensembles at N of a few hundred per
  segment are feasible as planned.
Next: WP4 (adaptive QD integrator on the same operators, hand-off on a
slip-rate threshold, fine-grid history reconstruction) and the rate-and-state
single-fault cycle against the rsf_solve QD anchor.

## 11. WP3 external validation (2026-10-11)

- Exact self-similar solutions (Broberg 1999, section 6.9): mode III factor
  1/E(k) and mode II factor 2(1-k^2)(k^2-b^2)/g2(b) are reproduced to 1e-3
  and better by the space-time code; the mode II check covers the in-plane
  kernels, the history convolution and the stepper together.
- SCEC TPV205-2D against Dunham's MDSBI 6.25 m solution: first-order
  convergence of arrival times (errors 0.27/0.14/0.07 s at 100/50/25 m,
  Richardson extrapolation within 1 ms of the reference), final slip within
  0.2 to 1.4 percent at 25 m.
- Implications: the stepover ensembles should use dx of 50 m or finer for
  TPV5-like parameters (the cohesive zone is about 600 m static); a 25 m
  run of a 30 km fault over 12 s takes about 10 minutes in the Python
  prototype on one core, dominated by the history convolution (L = 3300
  lags, M = 4096), so the 2-D ensembles are feasible but the C port will
  matter for the finer resolutions.
- Harris and Day (1993) setup noted for WP6: 2-D plane strain, 28 km
  segments, overlap 5 km, step widths -5 to +5 km in 0.5 km increments,
  slip weakening, four friction cases (A: mu_s 0.75, mu_d 0.3, supershear;
  D: mu_d 0.51, subshear, jumps only 0.5 km). Their Table 1 holds the
  initial stresses needed to reproduce it.
