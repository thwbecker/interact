# Becker and Schmeling (1998) redone with rate-and-state friction

Becker, T. W. and Schmeling, H. (1998), Earthquake recurrence-time
variations with and without fault-zone interactions, Geophys. J. Int.,
135, 165-176, used 2-D plane-strain finite elements with split-node
faults and static-kinetic friction under progressive simple shear.
This folder contains rsf_solve runs set up as analogues of their
single-fault models (their Figs 4 to 6) and of their leftstep and
rightstep systems A to D (Figs 7 to 10), the scripts that ran and
plotted them, and the resulting figures.

All figures here were made from rsf_solve output, not taken or
digitized from the paper. The comparison with the paper below is by
eye against the published figures and text.

## Files

    run_bs98          runs all eight models with ../run_model and calls
                      plot_bs98.py; usage: run_bs98 [rundir] [stride]
    plot_bs98.py      makes the four figures from the run directories;
                      usage: plot_bs98.py [rundir]
    bs98_fig4_5.png   mean tau - f0 sigma and event moments of the single
                      faults at beta = 45 and 29.5 degrees, first 550 yr
    bs98_fig6.png     event time against event number for the single
                      faults at beta = 29.5, 45, 54.6, 60.5 degrees, with
                      the paper's eq. (11)
    bs98_fig7_8.png   maps and per-fault moments of models A and B
    bs98_fig9_10.png  maps and space-time pattern of slipping cells
                      (|v| > 1 cm/s at the series records) of models C, D

## How the figures were made

The runs were made on 2026-11-08 in a sandbox (Ubuntu 24.04, one core)
with `bin/rsf_solve` built, as described in ../COMPILE.md (PETSc
3.19.6 from the Ubuntu package), from interact c4e2c67 plus the
changes later committed as 756f9f8. The only rsf_solve source
difference to 756f9f8 is the fully dynamic code changed in 5f0e42b,
which these runs do not use; a build of 756f9f8 reproduced the
beta45 series byte for byte. The runs were made with:

    becker_and_schmeling98/run_bs98 rundir 8

which calls `../run_model` once per model (rtol 1e-4, Runge-Kutta 3bs,
series output at every 8th accepted step during events, final step
recorded) and then `plot_bs98.py rundir`. Model A was first started
with a negative strike list that gen_system.py did not accept; it was
rerun on its own afterwards with the command now in run_bs98
(`--strikes=-15.5,0,-15.5`). The full set takes about 25 minutes of
wall time on one core and about 830 MB of series output.

| run | model time | cells (cell size) | wall time | events (all faults) |
|---|---|---|---|---|
| beta45 | 2000 yr | 226 (442 m) | 16 s | 32 |
| beta29.5 | 900 yr | 226 (442 m) | 9 s | 27 |
| beta54.6 | 2000 yr | 526 (190 m) | 34 s | 14 |
| beta60.5 | 2000 yr | 713 (140 m) | 45 s | 9 |
| modelA | 500 yr | 678 (442 m) | 120 s | 43 |
| modelB | 1500 yr | 1128 (266 m) | 698 s | 55 |
| modelC | 300 yr | 1200 (333 m) | 276 s | 36 |
| modelD | 300 yr | 1200 (333 m) | 241 s | 26 |

Common gen_system.py parameters (run_bs98 has the full commands):

    --G 20e9 --sigma 10e6 --taudot 6.4e-4 --dc 0.01 --a 0.010 --b 0.015
    --L 100e3 --W 15e3 --sigma_min 0.5e6

beta is 45 degrees plus the relative strike in `--strikes` (beta =
29.5, 54.6, 60.5 for -15.5, 9.6, 15.5). Models A and B are chains of
three 100 km faults with 5 km gaps (`--layout chain --gap 5e3`), the
side faults at strikes -15.5 (A) and +15.5 (B). Models C and D are
chains of five faults of 100, 50, 100, 50 and 100 km, the 50 km
segments at -15.5 (C) and +15.5 (D).

## Correspondence with the paper, and differences

- Loading: the pure-shear far-field stressing rate of gen_system.py
  gives the same shear and normal stress histories on a fault as their
  simple shear (constant mean stress, maximum shear stress growing at
  G times the strain rate, 3.2e-14/s, i.e. 20.2 kPa/yr). Background
  normal stress 10 MPa, G 20 GPa, fault length 100 km and W 15 km as
  in the paper; f0 = 0.6 takes the place of their mu_s = 0.6.
- Friction: rate and state (aging law, a 0.010, b 0.015, dc 0.01 m;
  h* about 10 km at 10 MPa) instead of static-kinetic friction with
  mu_s 0.6, mu_d 0.5. The stress drop is an outcome; the recurrence of
  the beta = 45 run corresponds to an effective mu_s - mu_d of 0.126,
  against 0.1 in the paper.
- Geometry: one-cell-deep strips in an elastic half space (2.5-D)
  rather than plane strain, so the moments come out smaller (about
  0.5e20 N m for the 45 degree fault, against 0.9e20 N m).
- Start: all faults start at steady state at the plate rate, so all
  rupture within the first few years. The paper loaded from zero shear
  stress and skipped ahead to the first failure, so first-event times
  are not comparable.
- The tilt of the side faults in models A to D is not given in the
  paper; 29.5 and 60.5 degrees to S1 were used here, the same values
  as for the single faults. The 5 km gaps are also a choice.
- Normal stress interaction between faults is included
  (`--sigma_interact 1`), as in the paper's elastic models.
- Cells are sized for the largest normal stress expected from the far
  field, except in model B, where the clamped side faults reach 33 MPa
  and the cells there end at about L_b/1.5. The normal stress limiter
  (0.5 MPa) was not reached in any run.

## Numerical and initial-condition checks

`run_checks [rundir]` reruns the models with finer cells (res), with
rtol 1e-5 instead of 1e-4 (tol), and with another random initial
slip-rate perturbation (seed2, `--seed 2`); for model C also with
point-symmetric gaps (`--gap_dir bisect`, runs Cs_*) without
(`--v_pert 0`) and with the perturbation. `compare_checks.py
bs98_rundir checks_rundir` prints per-fault statistics and writes the
event rasters checks_modelA.png, checks_modelB.png and
checks_modelC.png (one row per variant, tick size scaled with moment).
The runs were made on 2026-11-08 as described above, with the series
written at every 16th dense step; model B at the finer resolution was
interrupted twice by sandbox restarts and continued from its
checkpoints, as was the model B run with the other perturbation; the
duplicated catalog rows of the interrupted parts were removed before
the comparison (run_checks does not include these restarts). The model
B tolerance run was not made.

- Single faults: at half the cell size and at rtol 1e-5 the event
  sequences of beta = 60.5, 54.6 and 29.5 agree with the base runs to
  1 to 2 yr; the largest difference is in the last three intervals of
  beta = 54.6 at the finer cells (174, 184, 194 against 175, 187,
  199 yr). The slower than predicted recurrence growth on the clamped
  faults and the end of events on the 29.5 degree fault at 786 yr are
  therefore not discretization effects.
- Model A: finer cells (1014 instead of 678) and rtol 1e-5 give the
  same event counts per fault and mean intervals within 0.4 yr. Another
  random initial perturbation leaves the main fault and side fault 2
  nearly unchanged, but side fault 3 then has 19 instead of 16 events,
  more scatter in the intervals (CV 0.40 against 0.20) and some small
  events.
- Model C: finer cells (1805 instead of 1200), rtol 1e-5 and another
  perturbation give the same event counts and mean intervals within
  0.3 yr, except the connecting fault 3 at the finer cells (34.8
  against 32.8 yr). The offset of up to 20 yr between the outer faults
  5 and 4 in all of these comes from the geometry: with
  `--gap_dir previous` (used in run_bs98) the connecting segments are
  about 1.3 km out of point symmetry. With symmetric gaps and no
  perturbation, faults 5 and 4 stay within 0.13 yr and faults 2 and 3
  within 0.3 yr of each other over 300 yr; with the perturbation, 5 and
  4 stay within 0.4 yr while 2 and 3 drift up to 3.7 yr apart. In none
  of the variants does the middle fault 1 rupture within 1 yr of faults
  5 and 4; it tends to fail with the connecting faults, 3 to 13 yr
  apart from the outer faults. The loss of the whole-zone ruptures of
  the paper is thus a property of these rate-and-state models, not of
  the numerics or the initial perturbation.
- Model B: at 1.8 times finer cells (2025 instead of 1128; L_b/2.7
  instead of L_b/1.5 at the highest normal stress on the clamped side
  faults) the side faults have the same event counts, mean intervals
  within 1 yr and the same moment ranges, with event times within about
  6 yr after 1300 yr. The main fault has 33 instead of 34 events. The
  first two small main-fault events after the side-fault event at about
  1145 yr occur in both, with the same moments; after that the
  sequences differ (four more small events between 1194 and 1250 yr in
  the base run, a small pair at 1349 yr in the finer run). The onset of
  this irregular phase is resolved; its details are not.
  With another random initial perturbation (base resolution; this run
  was also continued once from a checkpoint) the model B sequence is
  nearly the same as the base run, including the six small main-fault
  events between 1157 and 1250 yr (times within 0.5 yr, moments within
  0.02e20 N m), so the late irregular phase depends on the resolution
  in its details but not on the initial perturbation.

## Results, with the comparison to the paper (one run per model)

- Figs 4 and 5: the 45 degree fault gives periodic events every 61 to
  62 yr of constant size; the 29.5 degree fault gives events that
  become more frequent and smaller (intervals from 52 to about 35 yr,
  moment down by more than half over 550 yr). Both as in the paper.
- Fig. 6: with the effective friction drop of the 45 degree run, the
  paper's eq. (11) matches the 29.5 and 45 degree runs closely. For the
  clamped faults (54.6, 60.5 degrees) the recurrence also grows with
  time, but less than eq. (11) predicts: the first intervals are 10 to
  20 percent longer than predicted, the later ones shorter. The 29.5
  degree fault stops producing events after about 790 yr, when the
  normal stress has fallen to about 2 MPa and the nucleation length
  exceeds the fault length; static-kinetic friction has no such
  transition.
- Figs 7 and 8: in model A (leftstep) the side faults fail more and
  more often with smaller events (intervals from 33 to 36 yr down to 21
to 27 yr), and the
  main fault recurs every 51 to 52 yr, 16 percent faster than the
  single 45 degree fault; the paper found 15 percent. In model B
  (rightstep) the side faults fail less and less often (intervals from
  about 90 to 110 yr up to 200 to 250 yr) with moments growing from 0.5 to 1.2e20 N m, as
  in the paper. The main fault of model B recurs every 49 to 54 yr,
  faster than the single fault, where the paper found no significant
  change, and it has a few smaller events late in the run.
- Figs 9 and 10: in model C the short connecting segments fail about
  every 30 yr with small events; in model D they are mostly locked and
  fail only a few times, both as in the paper. Unlike the paper, the
  three long faults of model C do not stay synchronized after the first
  events.

In short, the trends with fault orientation and with left- versus
rightstep interaction come out as in the paper, and two numbers (the
45 degree recurrence, by construction of the friction drop, and the
15 to 16 percent shortening of the main-fault recurrence in model A)
agree closely; the clamped-fault recurrence growth, model B's main
fault and the synchronization in model C differ. These are single
runs with one choice of friction parameters and side-fault angles.
