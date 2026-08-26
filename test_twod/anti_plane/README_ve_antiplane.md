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
   and restart reproducibility.  Physics: recurrence shortens
   monotonically with stronger relaxation and shallower substrate
   (-0.5 to -1.2 percent for a substrate 20 km below the fault, to
   -6.4 percent at 2 km below), with virgin and spun starts converging
   to the same attractor for contained (non-through-plate) faults.

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

6. Loading conditions (documented in kato_mn_replication_notes.md):
   backslip is well posed for faults contained in the elastic plate;
   through-plate and uniform-medium geometries have a stress-free
   relaxed slip mode (protocol-matched early cycles meaningful,
   multi-kyr drift not); constant far-field stress reproduces the M&N
   stability question; a constant stressing rate reproduces the
   elastic backslip cycle exactly (137.86 vs 137.87 yr) but is ill
   posed on a relaxing kernel (unbounded acceleration) and should only
   perturb a well-posed backslip configuration.

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
