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

## Sense of faulting

branch +1 (vpl = +1e-9, the default of run_bp3 and gen_bp3.py) is
THRUST: positive along-segment slip against the down-dip tangent
lifts the hanging wall (checked with interact on the 2 km dip-60
geometry: +0.51 m at 5 km on the hanging-wall side, -0.29 m on the
footwall), and the fault-normal stress change of that slip agrees
with an independent Okada (1992) evaluation to 0.1 percent at every
cell for dip 30 and 60 (positive slip unclamps the shallow fault).
This is also the spec's convention (thrust = positive slip and
positive Vp).  branch -1 is NORMAL.  Between 2026-09-01 and
2026-09-06 the scripts and notes carried these labels the other way
round; directories stamped in that period with branch -1 and named
*_thrust or *_t_* hold NORMAL-sense results and vice versa.  The
physics of rsf_solve was never affected, only the labels.

## Results vs the community (their section "BP3-QD", Figure 8)

  The paper's Figure 8 gives interevent times for the THRUST cases
  only (normal cases are in its supplement).  25 m, 1500 yr, 3bs.

  dip 90 (slip induces no dsigma; verified sigma stays 50 MPa):
    both branches identical, characteristic events every 89.83 yr.
    Community: ~90 yr (first event ~185 yr).  MATCH.

  dip 60 thrust (branch +1), rtol 1e-5 and 1e-6 agree to 0.06 %:
    first event 179.3 yr, then 90.9, 59.0, 85.6, 88.7, 66.4, 85.8,
    88.5, settling on a 3-cycle {66.5, 85.8, 88.5} yr.
    Community: first event ~178 yr, then ~98, 57, 87, 92, 97, 62, 87,
    92, 97, ... i.e. a 4-cycle {~60, 87, 92, 97}; the codes disagree
    among themselves at this level (sbplib settles on a 1-cycle near
    87, TriBIE on a 2-cycle 62/95).  The first five events match to
    a few years; the long-term attractor is within the community
    spread but not the majority one.

  dip 30 thrust (branch +1), rtol 1e-5:
    first event 182.3 yr, then 65.5, 68.7, then a 1-cycle at 68.5.
    Community: first event ~172 yr, then a 2-cycle {~63, ~82}.
    NOT reproduced at 25 m.  The 100 m run gives first event 177.2
    and a 2-cycle {64.5, 86.7}, closer to the community; 50 m gives a
    deep partial first event at 164.7 yr (7-20 km down-dip, not
    surface breaking, 0.63 m) and then a 1-cycle at 69.4; a 12.5 m
    run reached its first event at 169.6 yr.  The first event is
    thus not converged with cell size at either 25 or 12.5 m, while
    the eight community codes agree on ~172 yr.  These differences
    are tolerance-independent (100 and 50 m identical at rtol 1e-5,
    1e-6, 1e-7).  Open.

  dip 60 normal (branch -1): 1-cycle 86.75 yr at 25 m (rtol 1e-5 and
    1e-6), against 94.0 at 100 m and 94.5 at 50 m; first event 174.9
    (25 m) vs 195.0 (100 m), 180.1 (50 m), 187.6 (12.5 m).
  dip 30 normal (branch -1): 1-cycle 80.63 yr at 25 m; 100 m gives a
    2-cycle {74.7, 106.0}.
    No community numbers in the main text for either.

  Resolution: without sigma coupling the dip-60 ladder (200/100/50/
  25 m) converged to a single 85.93-yr event.  With sigma coupling the
  first-event time and the attractor type still move between 50, 25
  and 12.5 m on both branches, so the spec's 25 m is not sufficient
  for convergence of the sequence here, although the interval scale
  (65-95 yr) is stable from 200 m down.

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
