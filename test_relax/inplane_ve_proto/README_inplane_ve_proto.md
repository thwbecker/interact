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
ladders, write the same file).  SHEAR ONLY so far: hereditary
normal-stress relaxation (In(t)) remains open, so -calc_sigma_dot
VE runs are refused.

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
  antiplane mode 2, and exercise the planned -ve_prony_file
  interface end to end.
- FEM upgrade (quadratic or incompatible-mode elements) if an
  independent time-domain VE check is wanted beyond the T2/T3
  limits.

## Files

inplane2d.py  the solver module (basis, matcher, sources, k
    integration, Maxwell moduli, Talbot, ve_fields_xz)
fem2d.py      independent FEM reference (bilinear quads, split nodes)
t0_basis.py t1_elastic.py t2_maxwell_corr.py t3_plate_maxwell.py
run_inplane_proto_test   gate driver
