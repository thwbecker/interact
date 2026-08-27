# Loading conditions for visco-elastic earthquake-cycle models: what changes, what is implemented

TWB / Claude, August 2026.  Companion to README_ve_antiplane.md (the
validation ledger) and noda_mn_tests.md (the Kato 2002 and Miyake &
Noda 2019 replications, where most of the numbers below were
established).  All statements are backed by runs cross-validated
against independent exact-kernel integrators
(kato02_exact_chain.py, mn19_exact_chain.py).

## Why loading matters more once the medium relaxes

In an ELASTIC medium the loading representation is nearly immaterial:
backslip loading (stress from the slip deficit against steady plate
motion) and the equivalent constant stressing rate produce the same
cycle to numerical precision (137.87 vs 137.86 yr recurrence on the
Kato-2002 geometry), because recurrence is stress-budget controlled;
the substrate depth drops out of the elastic kernel entirely.  Once
the medium relaxes, the loading PATHWAY is itself a relaxing object,
and the question "how does viscoelasticity change the cycle" is
ill-posed until the drive is specified.  This resolves the apparent
contradictions in the literature: published results are slices
through the following map.

## The four regimes (2-D antiplane, rate-and-state faults)

1. BACKSLIP, FAULT CONTAINED IN THE ELASTIC PLATE (the well-posed
   configuration; -ve_mode 2 with fault bottom above -ve_plate_h).
   The drive is anchored in the deep creeping fault section inside
   the unrelaxing plate.  Elastic and VE cycles both have unique,
   initial-condition-independent attractors (virgin and spun memory
   converge to identical recurrence).  Substrate relaxation mildly
   weakens the elastic coupling of the deep creep to the seismogenic
   zone: mild recurrence LENGTHENING, +1.3 / +1.0 / +0.6 / +0.1
   percent for (tM = 2.5 yr, substrate 2 km below a 20-km fault),
   (5 yr, 2 km), (5 yr, 5 km), (5 yr, 20 km); elastic reference
   250.9 yr; robust to rtol, Prony ladder size, and initial state.
   (An earlier sweep reported SHORTENING here; that sign was a bias
   of the pre-fix lagged memory sink, see -ve_h_stage.)

2. BACKSLIP, DEGENERATE GEOMETRY (fault cutting the entire plate,
   the Kato 2002 configuration; also any uniform-Maxwell whole-space
   backslip problem).  The relaxed operator annihilates a uniform
   slip deficit, leaving a neutral mode: NO steady cycle attractor.
   Early, protocol-matched cycles are meaningful and show strong
   transient SHORTENING (-22 percent at tr = 1-10 yr for Kato's
   geometry: the afterslip zone abuts the substrate and relaxation
   re-energizes post-event reloading), which heals back toward the
   elastic recurrence over ~100 cycles.  Multi-kyr trajectories
   integrate kernel-representation detail along the neutral mode and
   should not be interpreted.  For physics, contain the fault.

3. CONSTANT FAR-FIELD STRESS (the Miyake & Noda 2019 configuration:
   whole fault rate-and-state, rate-strengthening surroundings creep
   at Vpl, uniform driving stress; realized here as backslip on a
   fully frictional fault, since the uniform mode carries no stress).
   With relaxation the medium restores the driving stress after every
   event and nothing ties fault slip to far-field displacement.
   Recurrence LENGTHENS with stronger relaxation and DIVERGES,
   t_rec = alpha/(t_c/t_load - t_cr) + beta, with the patch
   permanently stuck (their ST class) below t_cr: 0.86 for an
   isolated patch strip, about 1.1 with 4.096R-periodic patches,
   1.52 in their fully dynamic periodic case (the remainder is
   inertia).  Quasi-dynamically there is also a mild t_rec MINIMUM
   (about 0.75x elastic near t_c/t_load = 2.4) before the divergence.

4. CONSTANT STRESSING RATE on a relaxing kernel (-rsf_stress_rate_file
   with the backslip drive off): ILL POSED.  The external pump never
   relaxes while the elastic constraint does; cycles accelerate
   without bound (28 -> 4 yr and falling in the test).  Constant-rate
   terms are meaningful only as PERTURBATIONS on a well-posed backslip
   configuration (stress steps, seasonal or postseismic transients),
   or in purely elastic runs, where they reproduce the backslip cycle
   exactly.

## The literature, placed on this map

Savage-Prescott-lineage imposed-cycle models (Savage & Prescott 1978;
Hetland & Hager 2005) see viscoelasticity only in surface-deformation
timing, by construction.  Kato (EPS 2002) is regime 2, and his
reported null effect on recurrence is not reproducible from his
equations; it is quantitatively explained by a dynamically inert
single-image-family kernel (noda_mn_tests.md).  Miyake & Noda (EPS
2019) is regime 3; their divergence and stuck transition reproduce
here.  Bonafede-style and Allison & Dunham (2018, 2021)-style
shortening lives where relaxing material sits directly beneath the
afterslip zone (regime 2's transient; their thermomechanical and
power-law couplings can push either sign).  Lambert & Barbot (GRL
2016) demonstrate the coupled fault-plus-strain-volume framework
(early-onset postseismic flow) without recurrence systematics; their
Vpl = 0 remark defines a mantle-driven loading class.  Shi, Wei &
Barbot (JGR 2022) show viscoelastic stress transfer synchronizing
asperities, the 3-D analog of the multi-fault experiments in
../anti_plane_cycles.

## What is implemented in rsf_solve

- Backslip hereditary loading: the default; -ve_mode 1 (uniform
  Maxwell, exact single pole) and 2 (plate over Maxwell half-space);
  memory initialization -ve_h_init (1 spun, 0 virgin); sink accuracy
  -ve_h_stage (1 stage-consistent default, 0 cheap lagged).
- Constant far-field stress class: build the fault fully frictional
  with creeping surroundings (gen_mn19.py) under backslip.
- Constant stressing rates and per-cell loading-rate perturbations:
  -rsf_stress_rate_file (tau_dot, sigma_dot per patch), with the
  regime-4 caveat above.
- Stress-step and transient-perturbation experiments: restart from a
  checkpoint with modified initial conditions, or add a transient via
  the stress-rate file in elastic runs.
- NOT implemented: mantle-flow driving through evolving strain-volume
  elements (the Lambert-Barbot Vpl = 0 class); this is the natural
  extension alongside -ve_prony_file for 3-D kernels, and interact's
  eigenstrain machinery is the hook.
