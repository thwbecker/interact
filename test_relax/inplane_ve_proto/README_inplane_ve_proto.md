# Standalone prototype: plane-strain (in-plane) viscoelastic kernels for thrust faults in cross section

Status 2026-08-29.  DELIBERATELY OFF THE MAIN BRANCH: this directory
is a self-contained python prototype plus its validation gates, the
groundwork for implementation-plan step 7b (plane-strain layered
kernels for rsf_solve via -ve_prony_file).  Nothing here touches
tracked code.  Gravity IS included in the standard local-buoyancy
form (see below).

## What it solves

2-D plane-strain quasi-static deformation of an ELASTIC PLATE
(thickness H, free surface) over a MAXWELL HALF-SPACE (Maxwell in
shear, elastic in bulk), driven by dip-slip dislocations, e.g. a
thrust fault in cross section.  Wavenumber-domain state-vector
formulation Y = [ux, uz, sxz, szz](k, z):

- exact analytic basis per uniform region (verified against Navier
  symbolically with sympy, gate T0); for k*width < 2 the pure
  exponential columns degenerate and are replaced by matrix-
  exponential propagator columns (a genuine failure mode found
  during validation: the ill-conditioned matcher amplified rounding
  into a smooth ~10 percent long-wavelength stress error, invisible
  to quadrature refinement);
- point dislocations inserted as state-vector jumps derived from the
  moment-tensor equivalent-force representation ([ux] = m_xz/(2 pi
  mu), [uz] = m_zz/(2 pi (lam+2mu)), [sxz] = (ik/2 pi)(m_xx - lam
  m_zz/(lam+2mu)), [szz] = 0); a fault = a line distribution of
  point moments (segment_sources);
- viscoelasticity by correspondence, mu2(s) = mu2 s/(s + 1/tau),
  lam2(s) = K2 - (2/3) mu2(s), inverted with a fixed Talbot contour;
  on the contour the moduli are complex and the real-axis conjugate
  symmetry of the k-integral FAILS: the negative-k half line is
  computed via the parity relation A(-k) = S A(k) S with
  S = diag(1,-1,1,-1) (a second solve per k).  Forgetting this is
  another failure mode the gates catch (T2 fails by factors).

## Gravity

Advected-boundary (local buoyancy) gravity, the standard treatment
for crustal deformation (the same class as Rundle 1982 and Wang's
PSGRN): continuity of TOTAL stress across the deformed free surface
and interface gives, to first order (z positive down, uz positive
down),
    free surface:  szz(0) = rho1 g uz(0)
    interface:     szz_below - szz_above = (rho2 - rho1) g uz(H),
two extra linear rows in the wavenumber-domain matcher; cost zero.
NOT included: self-gravitation of the perturbation and internal
buoyancy from compressible density gradients (both negligible at
crustal scales), and any lateral pre-stress advection.  What gravity
buys physically: the elastic fields barely change (rho g/(mu k),
sub-percent at seismogenic wavelengths, gate T4a), but the RELAXED
state changes qualitatively: the through-plate neutral mode (a fault
cutting the whole plate over a fluid substrate has zero-energy block
motions without gravity) is cured, with the long-wavelength relaxed
response capped by buoyancy instead of blowing up (gate T4c: at
k = 1e-6 the g = 0 relaxed surface response exceeds the g-on one by
6e5; the relaxed limit with gravity is the buoyant flexure state,
verified against a direct solve at 1e-7, gate T4b).

## Gates (run_inplane_proto_test; all PASS on delivery)

T0  basis ODE residual 6e-10; source conjugate symmetry; free
    surface, source jump, interface continuity exact in the matcher.
T1a k-solver vs interact's OKADA kernels (4000-km 3-D vertical
    dip-slip fault at mid-strike = the plane-strain limit; the
    repository binary): sxx/szz/sxz rel misfit 2e-3 to 7e-3 over a
    window spanning plate and substrate depths.  This is the
    community-grade elastic anchor.
T1b independent FEM (fem2d.py, bilinear quads, split nodes):
    agrees to ~2.5 percent near the source but carries a ~9 percent
    bias at depth: bilinear quads SHEAR-LOCK in this bending-
    dominated free-surface problem (they match full-space fields to
    ~1 percent where nothing flexes, and Okada sides with the
    k-solver, firmly attributing the bias).  Kept as a soft gate and
    as scaffolding for future VE time-stepping checks with better
    elements.
T2  uniform Maxwell medium: every stress component of a held
    dislocation must follow the single exact correspondence factor
    f(t) = L^-1[mu(s)(lam(s)+mu(s))/((lam(s)+2mu(s)) s)]/D0 (the
    full-space field is proportional to mu/(1-nu)); deep-source
    configuration, sympy-exact partial fractions vs the full
    machinery: agreement 3e-9 at t = 0.3, 1, 3 tau.
T4  gravity: elastic fields insensitive at the expected rho g/(mu k)
    level (Okada anchor intact); relaxed-with-gravity limit matches
    the direct buoyant-flexure solve to 1e-7; the through-plate
    zero-energy mode is suppressed by 6e5 at k = 1e-6 (and the
    suppression grows toward k -> 0, the buoyancy signature).
T3  plate over Maxwell substrate: t->0+ reproduces the unrelaxed
    layered elastic solution (2e-4); t = 300 tau matches the
    mu2 -> 0 relaxed elastic solution (3e-3); interface tractions
    sxz, szz continuous through time at machine precision while sxx
    develops the legitimate modulus-contrast jump; plate-point
    stress evolves monotonically between the limits.
T7  the batched matched-pair evaluation path (what kernel generation
    actually runs) against the original per-wavenumber loop, kept as
    fields_pairs_ref: separated bases reproduce basis/basis_cs to
    4e-16; the batched matcher agrees with solve_k to 7e-11 at
    condition numbers up to 2e8, with both solutions leaving a scaled
    residual of 2e-15, i.e. both are backward stable and the
    difference is the conditioning; assembled fields agree to 2e-13,
    and to 4e-14 against the diagonal of the untouched fields_xz.
T6  Maxwell-time scaling: the relaxation part of both traction
    families, sampled at matched t/tau for Maxwell times in integer
    and non-integer ratios, agrees to 1e-13 of the relaxation scale
    with gravity on and off, and the relaxed limit does not move at
    all.  See the section below; the file-level counterpart is
    run_tau_scaling_test.

## Literature anchors

Rani & Singh (GJI 1992): closed forms, long dip-slip fault, uniform
elastic half-space.  Freund & Barnett (BSSA 1976): 2-D surface
deformation for dip-slip.  Singh & Rosenman (PEPI 1974):
viscoelastic half-space dislocation via correspondence (the T2
construction).  Savage (JGR 1998): edge dislocation in a layered
half-space.  Rundle (JGR 1978, 1982): layered viscoelastic(-
gravitational) dislocation solutions and thrust cycles.  Fukahata &
Matsu'ura (GJI 2006): multilayered viscoelastic internal
deformation, equivalence theorem (target for the layered extension).
Sen & Debnath (AJCAM 2012) give the uniform-VE dip-slip creep
variant of regime figures.

## Literature replication: Rundle (JGR 1982) Figure 2

rundle82_compare.py reproduces Rundle's Figure 2 configuration (his
purely viscoelastic, no-gravity case): 30-deg dipping thrust,
surface-breaking, down-dip width W = H, elastic layer over Maxwell
half-space, lambda = const, times 5 tau_a and 45 tau_a with
tau_a = 2 eta/G = 2 tau_M.  Against points digitized from his
figure (+-2 units reading error): the 5 tau_a change profile agrees
at 5.8 percent rms of peak (peak -27.2 vs his -26), including the
zero crossings and both far lobes; the 45 tau_a profile matches in
shape with the 2-D basin ~15 percent deeper than his (-87 vs -75).
That residual is quantitatively consistent with his FINITE fault
length (2L = 6.7H): interact-Okada elastic finite-length factors at
mid-strike are 0.8-1.1 over the basin and collapse (sign flips)
beyond ~2.5H, and late-time relaxation weights exactly the long
wavelengths that finite length suppresses.  Figure:
rundle82_fig2_compare.png.

THE GRAVITY CASE, his Figure 3 (rundle82_fig3_compare.py): same
fault, H = 30 km, rhoL = 3.3, rhoH = 3.8 g/cc, muL = lamL = 30 GPa
(from Rundle 1981, J. Phys. Earth 29, 173).  With the advected-
boundary buoyancy rows on, the 45 tau_a basin depth matches to 1.3
percent (-44.4 vs his -45.0), core rms 6.5 percent of peak; 5 tau_a
at 9 percent rms.  Gravity is exactly what cut his basin from -75
(Fig 2) to -45 (Fig 3), and the same suppression emerges here, so
the buoyancy implementation is validated against the literature, not
just internally.  Note the finite-fault-length discrepancy of the
no-gravity comparison largely disappears with gravity on, because
buoyancy itself removes the long wavelengths that finite length
would otherwise clip: the 2-D and his 3-D solutions converge.
Figure: rundle82_fig3_compare.png.  Thatcher & Rundle (1984) cycle
curves remain as follow-up.

## Modern benchmark: PSGRN/PSCMP (Wang et al. 2006, 2020 version)

psgrn_compare.py cross-validates against the community-standard
layered viscoelastic-gravitational dislocation code (source:
github.com/RongjiangWang/PSGRN-PSCMP_2020, cloned and gfortran-built
in the sandbox at ~/psgrn_pscmp; input templates in psgrn_inputs/).
Case: H = 30 km elastic plate (mu = lam = 30 GPa, rho 3300) over a
Maxwell half-space (30 GPa, rho 3800, tau_M = 1 yr), 30-deg
surface-breaking thrust, W = 30 km, 1 m slip; PSCMP fault 600 km
long, mid-strike profile (2-D limit), gravity ON and OFF via PSGRN's
gravity factor vs the prototype's buoyancy rows.  Postseismic uplift
CHANGE profiles at 10 and 90 tau_M, 101 receivers, curve-on-curve:

    gravity ON : 2.65 percent of peak rms (10 tau_M), 1.56 (90)
    gravity OFF: 1.09 percent of peak rms (10 tau_M), 1.56 (90)

(residual budget: their finite 600-km fault, their fuller gravity
theory vs local buoyancy, GF grid and FFT time sampling, our
quadrature).  The gravity suppression of the relaxed basin (-44 vs
-87 units) agrees between the codes and with the Rundle Fig 2/3
pair.  PSGRN also slots in later as the 3-D kernel generator for
-ve_prony_file.  Figure: psgrn_compare.png.

## BP3-VE: external kernels driving rsf_solve (-ve_mode 3)

bp3_ve_kernels.py closes the loop from this prototype to earthquake
cycles: it samples the RELAXATION part of the shear kernel,
dK_ij(t) = K_ij(t) - K_ij(0), for a dipping-fault discretization
(differences of VE and elastic solves on one truncated k-grid, so
the singular direct parts cancel and dK is regular everywhere), fits
per-pair amplitudes on a fixed shared ladder with the relaxed limit
as a constraint and one held-out sample as the gate (4.6e-4 of the
relaxation scale in the demo), and writes the -ve_prony_file that
rsf_solve's new generator-agnostic -ve_mode 3 consumes.  Runtime
gates in rsf_solve: the file's high-k elastic K0 block must match
the assembled Is at off-diagonal entries (0.65 percent in the demo;
note interact assembles Is in the stress-DROP convention, hence the
generator's default sign = -1), and a zero-amplitude file reproduces
the elastic catalog BIT-IDENTICALLY.  Demonstration (dip-60 BP3
thrust, ds = 667 m, H = 40 km plate, gravity on, tM = 45 yr, sigma_n
fixed): large-event recurrence lengthens from 124.5 yr (elastic
mean) to 145 yr (+17 percent) at similar slip per event, the
in-plane analog of the contained-fault loading-pathway relaxation.
The same file contract is the intended landing pad for
PSGRN/PSCMP-derived 3-D layered viscoelastic-gravitational kernels
(sample step responses with PSGRN per source patch, fit the same
ladders, write the same file).  Both traction families are
implemented (see above); what remains open is 3-D kernels via the
PSGRN route.

## Maxwell-time scaling: one sampling pass per sweep

The substrate enters the wavenumber-domain system only through
mu2(s) = mu2 s/(s + 1/tau) and lam2(s) = K2 - 2 mu2(s)/3, and the two
buoyancy rows carry no s, so substituting s = sigma/tau leaves the
system dependent on s only through s*tau.  The step response is then
a function of t/tau alone,

    dK(t; tau) = dKtilde(t/tau),     dK(t) = K(t) - K(0+),

with K(0+) elastic and Maxwell-time independent in any case.  The
fitted amplitude blocks C_p are therefore the same for every Maxwell
time and only the ladder tau_p scales, so a sweep over Maxwell times
costs ONE sampling pass plus a fit per Maxwell time.  This is exact
for a single Maxwell time with an elastic plate; a second relaxation
time (Maxwell plate, Burgers or bi-viscous substrate,
depth-dependent viscosity) leaves only the ratios free and would
need one pass per case.  Gravity does not spoil it, which is worth
stating explicitly because buoyancy is the one place where an
s-independent term enters the matcher.

What bp3_ve_kernels.py does with it:

- the sampling times and the ladder are fixed multiples of tM, so in
  units of tM they are the same set for every Maxwell time;
- the sampled cache is named after a signature of what it actually
  depends on (geometry, H, gravity, sign, sampling parameters) rather
  than after the output file, so every Maxwell time written into one
  directory shares it, across separate invocations and across the
  serial and parallel paths;
- tM_yr may be a comma-separated list, with %TM% in the output name,
  to emit a whole sweep from one call;
- a cache written before this path existed is adopted rather than
  resampled, and the cache is written through a temporary file so
  that an interrupted write cannot cost the sampling.

Measured (run_tau_scaling_test, dip-60 BP3 geometry, gravity on and
off): the amplitude blocks of a scaled file and of a directly
generated one agree to 1e-11 to 1e-10 of the relaxation scale,
against a held-out fit residual of 1e-3 to 1e-2 in the same units,
i.e. the difference is roundoff and not a modelling choice.  A file
written for the Maxwell time the sampling was done at is
bit-identical to what a single-Maxwell-time run produced before this
path existed; a file for another Maxwell time differs from a direct
generation of it only through the rounding of the scaled ladder.
The array-level statement is gate T6 (1e-13), which is the cleaner
number because it is not limited by the 8-digit kernel file format.

## Cost of kernel generation

The wavenumber quadrature is evaluated in batched form.  Both region
bases are separable in depth,

    B(z) = sum_m f_m(z) C_m ,     C_m independent of z,

with f = (ch, sh, kz ch, kz sh) in the cosh/sinh form and
(ed, eg, kz ed, kz eg) in the exponential one, so the state vector at
a receiver is a sum of four small matrix-vector products rather than a
4x4 basis built per receiver, and the quadrature over wavenumbers
becomes a matrix product.  The factor that carries the geometry,

    g_m[k, i] = w_k e^{i k x_i} f_m[k, z_i],

contains no material parameters, so it is built once per source and
reused for every node of the Laplace contour and every sample time;
the negative-k parity contribution reuses it as its conjugate.  What
remains per contour node is the batched 10x10 matcher solve and one
matrix product per source, region and parity.

Measured on the SEAS BP3 dip-60 geometry at ds = 25 m (1600 columns,
one core of a 2-core container): 11.4 s per column before, 1.16 s
after, i.e. 0.5 h instead of 5.1 h for a full kernel file, and about
15 minutes on two cores.  Combined with the Maxwell-time scaling
above, a three-Maxwell-time sweep for one case costs one such pass
rather than three.  Smaller consequences of the same change: the
gate suites now run in tens of seconds (run_tau_scaling_test 16 s,
run_ve_normal_test 19 s), and the ds = 2 km demo including both
viscoelastic cycle runs finishes in 9 s.

The sampled values move by about 1e-11 relative to the previous
per-wavenumber loop, which is four orders below the fit's own
held-out residual; kernel files are therefore no longer bit-identical
to ones generated before this change, and the cache version was
deliberately NOT bumped so that a part-way sampling run survives.
The per-column cache write was also put on a wall-clock cadence
(SAVE_EVERY_S): at 25 m the cache is 450 MB, so writing it after
every column had become the dominant cost once the sampling itself
was fast.  An interrupted run now loses at most that much sampling.

Not batched: fields_xz, the full source-by-receiver-grid routine used
by the validation gates T1 to T5 and by the PSGRN and Rundle
comparisons.  It is not on the kernel-generation path, and leaving it
alone keeps those gates as an independent reference.  A worthwhile
follow-up for large runs is the memory of the parallel path: each
worker currently allocates the whole (nt, n, n) sample array although
it fills only its own columns, which is 450 MB per worker at 25 m.

## Hereditary NORMAL-stress relaxation (two-family kernels)

On a dipping fault, slip changes the fault-normal traction as well as
the shear traction, and when the substrate relaxes BOTH evolve.
bp3_ve_kernels.py therefore fits a second family of amplitude
matrices, from the same field solves (one extra projection,
n.sigma.n instead of t.sigma.n), on the SAME tau ladder: the
relaxation spectrum is a property of the medium, not of the traction
component.  The file declares the number of families in its header
and carries an N0 gate block plus the normal amplitudes.  rsf_solve
(-ve_mode 3) activates the family exactly when the elastic normal
path is on (-calc_sigma_dot), and refuses that flag with a
shear-only file, since relaxing shear with frozen normal traction is
not a consistent medium.  NOTE the sign conventions differ between
the two families: interact assembles In compression positive (it is
scaled by -1 internally) while the projection here is tension
positive, so the normal blocks carry the opposite global sign; the
N0 gate is what pins this (a pure flip reads exactly 2.0, which is
how it was found).

Gates: run_ve_normal_test (format and guards, N0 vs In agreement,
and a regression that a two-family file with the normal family
inactive is bit-identical to the shear-only file).  Measured on the
demonstration case (dip-60 BP3 thrust, ds = 1 km, H = 40 km,
gravity, tM = 45 yr): held-out fit residual 1.1e-3 for the normal
family (4.6e-4 shear), N0 vs In 1.2 percent, K0 vs Is 0.7 percent.

HOW MUCH IT MATTERS IS NOT SETTLED BY THIS DEMO.  Normal-stress
coupling changes measured recurrence by tens of percent, so it is
clearly a first-order ingredient for dipping faults, but at
demonstration resolution the result is not stable enough to quote a
number, let alone a sign.  Mean recurrence over the same four runs
as a function of the slip threshold that defines a "large event"
(ds = 1 km, stop 1200 yr):

    slip >   elastic            viscoelastic
             frozen  evolving   frozen  evolving
    0.3 m     15.7    11.7       10.8    15.8
    0.5 m     56.3    73.4       60.7    68.1
    1.0 m     80.3   128.3       97.1   136.2
    2.0 m    120.5   171.0      161.9   181.6

The elastic-vs-viscoelastic difference changes SIGN between the
0.5 m and 1.0 m thresholds, so no sign can be claimed from these
runs.  The reason is resolution: ds = 1 km is about 3 x Lb here
(BP3 asks for 25 m), so the catalog is dominated by partial-rupture
chatter and "recurrence" depends on where the event-size cut is
drawn.  Any physical statement about how normal-stress relaxation
shifts recurrence needs the converged mesh; this demo exercises the
machinery, it does not measure the effect.

## Caveats and next steps

- Gravity is local-buoyancy only (see above): adequate for the
  earthquake-cycle band; full viscoelastic-gravitational normal-mode
  machinery (Rundle 1982; Yu, Rundle & Fernandez 1996) only becomes
  necessary for very long wavelengths/timescales or self-consistent
  geoid work.
- The normal-stress hereditary path: for dip slip the fault-normal
  stress relaxes too; rsf_solve integration needs In(t) ladders
  alongside Is(t) (doubles the memory states; -calc_sigma_dot
  machinery exists on the elastic side).
- Next: sample K_ij(t) = stress-at-patch-i from unit-slip-at-patch-j
  step responses on a fault discretization using ve_fields_xz, fit
  the per-pair Prony ladders with the same held-out gate as the
  antiplane mode 2, and exercise the -ve_prony_file
  interface end to end.
- FEM upgrade (quadratic or incompatible-mode elements) if an
  independent time-domain VE check is wanted beyond the T2/T3
  limits.

## Files

inplane2d.py  the solver module (basis, matcher, sources, k
    integration, Maxwell moduli, Talbot, ve_fields_xz; and the
    batched matched-pair path pair_ctx / fields_pairs_ctx, with the
    original loop kept as fields_pairs_ref for the T7 gate)
fem2d.py      independent FEM reference (bilinear quads, split nodes)
t0_basis.py t1_elastic.py t2_maxwell_corr.py t3_plate_maxwell.py
t4_gravity.py t6_tau_scaling.py t7_batched_fields.py
run_inplane_proto_test   gate driver
run_tau_scaling_test     file-level Maxwell-time scaling gate
    (bp3_ve_kernels.py scaled emission vs direct generation)
tau_scale_compare.py     the comparison it uses (also usable on any
    two kernel files that should differ only in the Maxwell time)
