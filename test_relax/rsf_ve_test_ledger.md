# rsf_solve visco-elastic implementation: test ledger and status

TWB / Claude, August 24, 2026. Supersedes the status portions of
rsf_ve_theory_and_status.md; the theory sections there remain valid.
Companion: rsf_ve_layered_fit_design.md (7b fit design),
relax_fault_ve.md (testbed usage), ve-note-v5/v6 (theory note).
Numbers are specific to the tested configurations named below.

## What is validated, with the numbers

Uniform Maxwell half space (the exact case; steps 1 and 2):
- Three-pole Prony representation of all stress kernels (rates t_M,
  1.8 t_M, 1.2 t_M at nu = 1/4), amplitudes from sampled
  effective-moduli assemblies: elastic-sum, held-out-pole, and
  relaxed-limit identities pass at rounding (1e-13 to 1e-15) on
  Okada, plane-strain, and triangular kernels (ve_check), including
  a nu = 0.35 cross-check; displacement kernels carry no t_M pole
  and relax to the nu = 1/2 field.
- Time histories against Gaver-Stehfest numerical inversion at 8e-5
  (ve_laplace_check; the inverter floors near 1e-4).
- Prescribed-slip testbed (relax_fault_ve): closed-form route
  1.5e-14; interior stresses are exact by the same machinery. The
  section-1 effective-modulus scheme (relax_fault) is a
  one-relaxation-time approximation (late e-folding 1.005 t_M vs the
  exact 1.8 t_M) and is retained for comparison only.
- rsf_solve step 3 plumbing (delivered as a patch; adoption status
  to confirm): -maxwell_time sampling, SI scaling, matvec
  verification shared by dense, hmmvp, and HACApK paths at 4e-14 to
  2.5e-7 depending on tolerance; option-absent runs byte-identical.

2-D antiplane layered (elastic plate over Maxwell substrate; exact):
- His dip = -90 screw element verified against the closed-form
  surface profile at 2e-8 (test_twod/anti_plane/run_test, corrected).
- Image machinery through the element against the independent
  Rybicki arctan series at 1e-15 to 2e-15, all anchors and Gamma in
  [-1, 1], general rigidity contrast Gamma(s) = Gamma_0 +
  (1 - Gamma_0) b / (s + b).
- Prony fits against the independent Nur-Mavko time series: np =
  2..6 gives 3e-2 / 4e-3 / 6e-4 / 1.4e-4 / 1.4e-5 of coseismic on
  the [b/15, 1.5b] ladder; ladder floor and image truncation are
  paired (n_img 40 to 80 with that floor). Profile figures
  (run_layered_test): elastic dots exact at print precision,
  displacement 1 to 2e-4 per curve, velocity 2e-4 to 5e-3 per curve
  via the velocity-consistent s-domain fit (absolute 1e-5 level).

3-D layered, what the truth looks like (Kaj propagator, surface):
- True relaxation spectrum: three real branches over kH; the slow
  flexural branch reaches rate 0 at k = 0, so any finite exponential
  set truncates the far field.
- Shared-rate fits of the true responses: over 24 mechanisms/depths
  x 20 stations x 3 components (1140 signals), one optimized
  six-rate ladder (0.038, 0.186, 0.441, 0.909, 1.542, 2.459 in
  units of a = 1/(2 t_M)) fits at median 1.2e-4, worst 3.2e-3; the
  plain geometric [a/30, 2a] is within a factor two; np = 4 median
  1.5e-3. The ladder transfers across H = 10/20/40 unchanged in
  a units. Velocities consistent with finite differences at 1e-3.

## What does not work, deliberately on the record

- The scalar image construction (2-D-exact) transplanted to 3-D:
  recovers 5 to 66 percent of the true transient amplitude, thrust
  horizontals with near-zero shape correlation; relaxed mid-plate
  strike slip 3 percent vs the true 73 percent of coseismic. Not
  usable as a 3-D layered kernel at any tested depth fraction or
  mechanism. ve_pom_compare.c documents this.
- Rate ladders placed on the cubic-root branches: worse than blind
  geometric (amplitude weighting, not root support, governs).
- np = 8 on the fixed 2-D ladder span: degrades (basis crowding).

## Interior stresses (the current work front)

The pom maxwell propagator extension (-zr): interior transient
displacements plus szz, sxz, syz from the state vector, propagation
direction and physical z-up conventions established by derivation
from the generator matrix and finite-difference calibration.
- VALIDATED at the FD floor (1e-6): strike slip on vertical planes
  and tensile sources (azimuthal orders 0 and 2), t = 10/100/500 yr,
  zr = 2/8/14 km in a 20 km plate; free-surface tractions vanish at
  zr -> 0; t -> 0 null.
- OPEN DEFECT: azimuthal order |M| = 1 content (dip slip on any
  plane, strike slip on dipping planes) fails at the 0.1 to 0.9
  level. Two candidate fixes tested negative (removing the M = -1
  synthesis flips; a diag(1,-1,1,-1) parity conjugation, which
  improved syz twentyfold while worsening szz). Next: mask azimuthal
  orders one at a time in the wavenumber loop and read the required
  transform off measured per-order ratios, or derive the convention
  from pom_source_vectors against the m-file. Until fixed, interior
  stresses are trusted only for the validated content classes.
- Not yet built: in-plane components sxx, syy, sxy (ingredients
  already synthesized in code: the horizontal divergence, the uz
  gradient, and the generator row-2 identity for sxx); interface
  continuity check at z = H; the uniform-Maxwell interior anchor
  against the exact three-pole stresses (needs the layer variant
  extended with tR1 = tR2, or an independent route); finite-fault
  quadrature (the layer code's Gauss template).

## What to do next with rsf_solve, in order

1. Close the |M| = 1 defect and add sxx/syy/sxy; then the
   uniform-Maxwell anchor decides whether the interior extension is
   trusted (acceptance: match the exact three-pole interior stresses
   at the fit-accuracy level over the time grid).
2. Confirm/apply the step 3 plumbing patch in rsf_solve (dense and
   H-matrix verification banner; np > 1 H-matrix check on a real
   host is still pending from the container-only test).
3. Implement the -ve_prony_file route: header (np, shared rates),
   per-pair amplitude rows, held-out TIME samples in the
   verification banner. Generator script wraps the extended pom for
   layered 3-D; the uniform case needs no file (exact poles).
4. Step 4 proper: per-pole state vectors h^(p) (np = 6 layered, 3
   uniform), exact-exponential update, checkpoint layout sized by
   VE_MAX_NP = 6; slider and small-array regression against the
   elastic path and against relax_fault_ve limits.
5. First layered cycle applications in 2-D antiplane (exact
   machinery, already validated end to end) while 3-D interior
   kernels mature; Savage-Prescott as the cycle benchmark.
