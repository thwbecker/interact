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
