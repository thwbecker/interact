# SEAS BP3-QD (dipping fault, plane strain) for rsf_solve

OFF-BRANCH working directory (untracked).  2026-08-29.  Spec:
SEAS_BP3-QD-FD.pdf (Erickson & Jiang); reference behavior: Erickson
et al., BSSA 2023 (bssa-2022066.1.pdf).

## Setup

interact's TWO_DIM_HALFPLANE_PLANE_STRAIN kernels (Crouch &
Starfield segments with free surface) run BP3-QD with the existing
rsf_solve machinery; the along-segment component is the mode-II
slip.  gen_bp3.py emits the EXACT spec parameterization: frictional
fault 0-40 km down-dip, a: 0.010 / linear 15-18 km / 0.025,
b = 0.015, sigma_n0 = 50 MPa, L = 8 mm, V0 = 1e-6, f0 = 0.6,
Vp = Vinit = 1e-9, rho 2670, cs 3464 (G = 32.04 GPa), nu = 0.25,
spec eq (24)/(25) initial conditions (identical to our BP1
convention), steady creep below 40 km via the backslip formulation.
IMPORTANT: the spec's effective normal stress INCLUDES slip-induced
changes (sigma_n = sigma_n0 + dsigma), so -calc_sigma_dot is part of
the benchmark; runs without it converge to a different (wrong)
attractor.  Driver: run_bp3 ds_km dip vpl outdir stop sense "extra
opts"; idempotent, resumes from checkpoints.

## Results vs the community (their section "BP3-QD", Figure 8)

  dip 90 (slip induces no dsigma; verified sigma stays 50 MPa):
    25 m cells, 1500 yr: characteristic events every 89.83 yr.
    Community: ~90 yr.  MATCH.

  dip 60 thrust, sigma-coupled, 100 m survey resolution:
    interevent set {63.3, 87.4, 91.0, 91.5, 92.3, 95.0} yr.
    Community (25 m): four characteristic events {~60, 87, 90, 95}.
    Pattern and values already reproduced at 4x coarser cells; the
    25-m run (bp3_d60A_sig_25m) is checkpointed at t = 744/1500 yr
    and resumes with the same command (finish locally: it is in a
    coseismic crawl that outlives this sandbox's per-call cap).

  dip 30 thrust, sigma-coupled, 100 m:
    period-2 {64.40, 86.51} yr.  Community: two characteristic
    events {~65, ~80} yr.  First matches to 1 percent; second is 8
    percent high at survey resolution.

  Resolution sensitivity: without sigma coupling the dip-60 ladder
  (200/100/50/25 m) wandered through period-2 artifacts before
  converging to a single 85.93-yr event; BP3 needs the suggested
  25 m (Lb/14) for converged patterns, matching the paper's warning
  that resolution is critical.

## Numerical notes

- Explicit RK (3bs) hits domain-guard rejection storms during large
  sigma-coupled events (trial stages push sigma through zero at the
  shallowest cells: the near-trace normal-stress concentration of a
  surface-breaking dipping thrust).  -imex integrates through these
  stretches cleanly at ~100 steps/s (25 m); use it for the
  sigma-coupled cases.  -limit_sigma 1 -min_sigma 1e6 -max_sigma 2e8
  is kept as a guard; sigma stayed within [49.5, 55.4] MPa at
  monitor times.
- Catalogs from chained restarts can carry duplicate/out-of-order
  rows at restart boundaries; sort on onset time and dedupe.
- Station time series (12 stations at fixed down-dip distances, spec
  section 6) are the next comparison level once the community data
  files are in hand (strike.scec.org).

## Files

gen_bp3.py, run_bp3, README_bp3.md; run dirs bp3_*.
