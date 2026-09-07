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
  only (normal cases are in its supplement).  25 m, 3bs, rtol 1e-5,
  -calc_sigma_dot.  Community station files can be downloaded per
  station from the SEAS platform without a login; the HBI (Ozawa)
  fltst_dp125 file for dip 30 thrust is the reference used below.

  dip 90 (slip induces no dsigma; verified sigma stays 50 MPa):
    both branches identical, characteristic events every 89.83 yr.
    Community: ~90 yr (first event ~185 yr).  MATCH.

  dip 30 thrust (branch +1): onsets 176.18, 263.23, 327.90, 414.48,
    479.13, 565.71 yr against HBI 176.18, 263.23, 327.93, 414.51,
    479.19, 565.77: within 0.06 yr over six events; interevent
    2-cycle {87.05, 64.67} against HBI {87.06, 64.68}.  At dp125 the
    interseismic tau and sigma agree with HBI to 1e-4 MPa.  MATCH.
    (50 m gives the same onsets to 0.05 yr; 100 m to 1 yr.)

  dip 60 thrust (branch +1): onsets 177.1, 277.2, 334.2, 421.2, 513.6
    yr, intervals 100.1, 57.0, 86.9, 92.4, against the community's
    ~178, then ~98, 57, 87, 92, 97 (a 4-cycle {60, 87, 92, 97}; sbplib
    and TriBIE differ among themselves at this level).  MATCH through
    the events run so far; the long-term attractor is to be confirmed
    with a full 1500 yr run.

  dip 60 and 30 normal (branch -1): no community numbers in the main
    text; the runs made before 2026-09-07 (86.75 and 80.63 yr 1-cycles
    at 25 m) are affected by the geometry-precision artifact below and
    need repeating.

  Geometry precision (2026-09-07).  Until this date gen_bp3.py wrote
  the segment centres with %.6e, i.e. to about 1e-2 m at 12 km.  The
  centres of a dipping fault are then collinear only to that level,
  and a glide segment, which produces exactly zero normal traction on
  its own plane, produces a spurious normal traction on its slightly
  misaligned neighbours that grows as 1/ds: at 25 m it was ten times
  the physical free-surface term at the adjacent cells.  Dip 90 is
  immune (x = 0 exactly), which is why it matched while the dipping
  cases drifted with resolution: dip 30 thrust first events 177 / 165
  / 182 / 170 yr at 100 / 50 / 25 / 12.5 m, a deep partial first
  event at 50 m, a 1-cycle instead of the 2-cycle at 25 m.  With
  full-precision centres the sequence converges and matches HBI.  All
  dipping-fault results produced before this date, including the
  viscoelastic demo's elastic references at 100 m, carry the artifact
  (small at 100 m: first event 177.2 against 176.2 yr).  Any
  generator that writes centres of a non-vertical fault should write
  them at full precision; make_thrust.py was changed as well.

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
