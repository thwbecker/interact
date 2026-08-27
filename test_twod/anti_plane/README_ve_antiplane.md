# Visco-elastic 2-D antiplane machinery: implementation and validation ledger

Status 2026-08-26.  This file summarizes WHAT exists for visco-elastic
(VE) earthquake-cycle modeling on 2-D antiplane shear zones, and every
validation case conducted, with the drivers that reproduce each.  The
step-by-step design rationale is in src/relax/rsf_ve_implementation_plan.md
and the code-level notes in src/rsf_ve.c.

## The machinery

Hereditary VE stressing in rsf_solve via a Prony/exponential
representation of the interaction kernel: tau_dot = K_el (v - vpl)
- sum_p h_p/tau_p with per-patch memory states h_p forced by
C_p (v - vpl); the relaxed operator C_inf = K_el - sum_p C_p is
implicit, so the coseismic response is exactly the assembled elastic
matrix and everything is bit-identical when VE is off.  Generators:
-ve_mode 1, uniform Maxwell (single exact pole, C_1 = K_el,
tau_1 = -ve_tmaxwell_yr = eta/G); -ve_mode 2, elastic plate over
Maxwell half-space (two-family Bonafede image kernel sampled through
the actual elements, per-pair 6-pole ladder fit with a held-out gate).
Options: -ve_mode, -ve_tmaxwell_yr, -ve_plate_h, -ve_g2fac, -ve_np,
-ve_h_init (1 backslip-spun memory, 0 virgin), -ve_h_stage (1 default:
stage-consistent sink forcing under TS error control; 0: cheaper lagged
forcing).  VE-aware checkpoint/restart.  A future -ve_prony_file slots
3-D kernels into the same architecture.

## Validation ledger (all antiplane cases conducted)

1. Layered kernel correctness (test_relax + test_twod/anti_plane,
   run_layered_test, ve_layered_check): the two-family image
   construction (families at 2nH+-[c1,c2] per reflection order with
   Erlang weights) against a finite-difference solve, a k-space closed
   form U(k) ~ (1-E)(1+Gamma E)/(1-Gamma E^2), and Nur-Mavko; anchors
   at 1e-15.  This work FOUND the missing second family in the
   repository's original construction (order-unity far-field errors)
   and fixed all four kernel functions.

2. Savage-Prescott cycle testbed (run_sp_cycle_test, ve_sp_cycle.c,
   sp_analytic.py): imposed periodic events on a through-plate fault;
   the Prony/h machinery against exact-in-time Erlang superposition
   and an independent numpy analytic; velocities to 4-9e-3 Vpl, stress
   to 6e-3 of the cycle swing, spun cycle-mean velocity 0.497-0.499 of
   the exact 0.5.

3. RSF cycles, BP5-style depth-dependent a-b (run_bp2d_cycle_test,
   gen_bp2d.py): elastic bit-identity when VE off; VE elastic limit
   (t_M = 1e7 yr) reproducing elastic event times to 1e-4 yr; the
   locked-fault hereditary loading gate against the closed-form
   two-family Erlang analytic (1.5-3e-3 of a 22 MPa swing); sustained
   layered cycles vs uniform-Maxwell loading saturation; -imex vs RK
   and restart reproducibility.  Physics, CORRECTED after the
   -ve_h_stage fix (the earlier sweep reported steady-state SHORTENING
   for contained faults; that sign was a bias of the lagged memory
   sink): contained faults over a relaxing substrate LENGTHEN mildly
   and monotonically with relaxation strength and substrate proximity,
   +1.3 / +1.0 / +0.6 / +0.1 percent at (tM = 2.5 yr, H = 22 km),
   (5, 22), (5, 25), (5, 40) for a 20-km fault with elastic recurrence
   250.9 yr; robust to rtol, np = 5/6, and virgin-vs-spun starts
   (identical attractor).  The loading-pathway relaxation (Miyake-Noda
   mechanism) wins for contained faults; the strong afterslip-reloading
   SHORTENING is a through-plate phenomenon (Kato replication).

4. Kato (2002) replication (run_noda_test gates 1-3, gen_kato02.py,
   kato02_exact_chain.py, noda_mn_tests.md): elastic Tcy/us match his
   Table 1 with zero free parameters (digitized Fig. 2 profiles); his
   viscoelastic null result does NOT replicate (true effect -23 to
   -4 percent over tr = 1..300, -31 percent for h = 15/12) and is
   quantitatively explained by a single-image-family (family A only)
   implementation, the same bug class fixed in this repo; confirmed
   against an exact-Erlang-chain integrator of his equations,
   including his truncation, his Eq. (3) misprint, and his Fig. 4
   afterslip signature.

5. Miyake & Noda (2019) replication (run_noda_test gates 4-6,
   gen_mn19.py, mn19_exact_chain.py): uniform-Maxwell rate-weakening
   patch; recurrence divergence t_rec = alpha/(t_c/t_load - t_cr)
   + beta and the permanently stuck (ST) class reproduced
   quasi-dynamically; t_cr = 0.86 (isolated strip) rising to about 1.1
   with 4.096R-periodic patches, vs their 1.52 (fully dynamic,
   periodic); mode 1 agrees with the exact single-state integrator to
   better than 1 percent after the -ve_h_stage fix.

6. SEAS BP1-QD community anchor (test_twod/anti_plane_cycles,
   run_cycles2d_test): the elastic antiplane benchmark of Erickson et
   al. (SRL 2020) reproduces with zero tuning; recurrence 77.88 /
   78.21 / 78.24 / 78.19 yr at 200 / 100 / 50 / 25 m cells against the
   community-converged 78.2 yr, demonstrating RSF resolution
   convergence explicitly (the suggested criterion ds < Lb/3 holds
   with margin).

7. Experiment suite (test_twod/anti_plane_cycles, README_cycles2d.md):
   single- and multi-fault (tested to 32 faults) drivers with slip-rate
   AND shear-stress field frames (-field_stress; step-cadence, dense
   through events), spacetime and stress-cross-section rendering, and
   per-fault catalogs; VE assembly accelerated by a translational-
   invariance sample cache (exact; 8x at 16 faults).  Multi-fault runs
   desynchronize into migrating sequences with partial ruptures;
   substrate relaxation adds long-range interaction (cf. Shi, Wei &
   Barbot, JGR 2022).

8. Off-fault stress cross sections (run_xsection, plot_xsection.py,
   this directory): sigma_xy(x, z, t) and sigma_yz through plate AND
   relaxing substrate, reconstructed from -field_slip frames via the
   validated image series (plate: Erlang-weighted two-family images;
   substrate: the transmitted series, whose stress weights reduce to
   Poisson masses p_n(bt) obtained from the same memory chains).
   Driver: fault reaching a fraction of H (default 0.5), elastic plus
   tM/T_rec = 0.1, 0.5, 1, 2, 10, all spun; snapshots through the last
   full cycle plus the coseismic change.  Built-in checks: on-fault
   kernel match, free-surface traction, and interface-traction
   continuity at every snapshot time (jointly tests image positions
   and both media's time weights).  Note this configuration produces
   period-2 cycles (alternating small/large events); T_rec is the mean
   clustered interval.

9. Loading conditions (documented in ve_loading_conditions.md, this
   directory):
   backslip is well posed for faults contained in the elastic plate;
   through-plate and uniform-medium geometries have a stress-free
   relaxed slip mode (protocol-matched early cycles meaningful,
   multi-kyr drift not); constant far-field stress reproduces the M&N
   stability question; a constant stressing rate reproduces the
   elastic backslip cycle exactly (137.86 vs 137.87 yr) but is ill
   posed on a relaxing kernel (unbounded acceleration) and should only
   perturb a well-posed backslip configuration.

## Document map

- README_ve_antiplane.md (this file): validation ledger, entry point.
- ve_loading_conditions.md: the four loading regimes, elastic vs VE,
  the literature placed on that map, and what is implemented.
- noda_mn_tests.md: the Kato (2002) and Miyake & Noda (2019)
  replications and their diagnosis; driver run_noda_test.
- run_xsection + plot_xsection.py: off-fault stress cross-section
  evolution (plate and substrate), elastic vs Maxwell-time series.
- ../anti_plane_cycles/README_cycles2d.md: the experiment suite
  (single/multi-fault, snapshots, figures); driver run_cycles2d,
  gates run_cycles2d_test.
- src/relax/rsf_ve_theory_and_status.md and
  rsf_ve_implementation_plan.md: theory lineage and build order, with
  dated status appendices.
- test_relax/rsf_ve_test_ledger.md: the original step-1/2 and 3-D
  ledger, with a dated update superseding its pre-two-family-fix rows.
- rsf_solve.md (repo root): solver usage, H-matrix backends, and the
  VE options/cost section.

## Known limitations

- Quasi-dynamic only (radiation damping); near-critical relaxation
  problems are inertia-sensitive (M&N t_cr).
- The 6-pole ladder spans [t_M/1.5, 45 t_M]; multi-century Erlang
  tails of the layered kernel are folded into the implicit relaxed
  operator.  Irrelevant for contained faults; through-plate NEUTRAL
  MODE trajectories beyond a few relaxation times are
  representation-sensitive (see noda_mn_tests.md caveats).
- -ve_h_stage 0 (lagged sink) is only safe when recurrence is far
  from relaxation-critical; the default (1) costs np matvecs per RHS.
- interact's classic one-step I-matrix path still lacks the antiplane
  stress-BC pairing (rsf_solve and the compressed-operator paths are
  wired).
