# 2-D antiplane checks against Wu & Chen (2014) and Cattania (2019)

The strike_slip_networks runs are 2.5-D (one-cell-deep strips, in-plane
mode II slip along strike). The two published parameter studies of
single-fault complexity used here as references are 1-D faults in 2-D
antiplane (mode III) geometry, so this folder sets up rsf_solve's 2-D
antiplane elements (as in ../../test_twod/anti_plane_cycles) to
reproduce them directly before comparing with the 2.5-D results
(../single_fault_regimes).

- `gen_antiplane.py model X [--law aging|slip] [--cpl n] [--vs_W f]
  [--depth_W d] [--out dir]` writes geom.in, rsf.dat, ic.in, run.opts
  and params.txt:
  - `wuchen`: Wu & Chen (2014) model A, uniform fault of length W, a
    0.008, b 0.016, dc 0.1 mm, sigma 50 MPa, G 30 GPa, cs 3464 m/s, f0
    0.6, v0 = vpl = 0.035 m/yr, loaded by the steady slip around the
    fault (backslip at vpl on the fault cells); their whole space is
    approximated by burying the fault 50 W below the free surface.
    X = W/hc, hc = (2/pi) G dc/((b-a) sigma) = 4.77 m; their regime
    boundaries W/h* = 1.2, 5.5, 33, 105 are W/hc = 2.4, 11, 66, 210 at
    b/a = 2.
  - `cattania`: Cattania (2019), vertical fault from the free surface,
    velocity weakening to depth W (a - b = -0.005), velocity
    strengthening below (a - b = +0.005) to depth 2 W, b 0.02, sigma 50
    MPa, dc 0.1 mm, loaded by creep from below (backslip at vpl).
    X = W/L_inf, L_inf = b G dc/(pi (b-a)^2 sigma). The paper does not
    state G, cs, vpl or the extent of the creeping region; here G 30
    GPa, cs 3464 m/s, vpl = v0 = 1e-9 m/s and a creeping region of W.
  - `--law slip` uses the slip law with otherwise identical parameters.
  - initial state: steady state at the plate rate with 0.1 percent random
    perturbations of the slip rate.
- `run_antiplane model X law ncycles [rundir] [extra rsf_solve options]`
  runs one case (rtol 1e-5, 3bs, catalog at 1e-3 m/s, all-cell series at every 16th dense step; the runs reported below used every 4th, which made the W/hc 63 series 3 GB)
  for ncycles times a rough recurrence estimate; `CPL=n` sets the cells
  per L_b = G dc/(b sigma).
- `run_antiplane_sweep [rundir] [ncycles] [laws]`: Wu & Chen W/hc 3,
  8.6, 20, 31.5, 63 (cells L_b/8, i.e. hc/10 as in their study) and 120
  (L_b/4); Cattania W/L_inf 3.2, 6.5, 13 and 26 (L_b/4); completed runs
  are skipped, so the sweep can be restarted.
- `plot_antiplane.py [rundir]` writes antiplane_bifurcation.png
  (interevent times / mean and rupture lengths / W against X),
  antiplane_nucleation.png (nucleation position along the fault, which
  separates Wu & Chen's centre- and side-origin regimes),
  antiplane_summary.png, antiplane_spacetime.png and
  antiplane_summary.txt. Events starting within one hour of the end of
  the previous one are merged (60 s in the first version, which left
  brief threshold crossings about a minute before some slip-law events
  as separate events), and the first quarter of each run is
  discarded.

## Results so far (2026-11-09, 40 recurrence estimates per run, single runs)

Completed: Wu & Chen aging at W/hc 3 to 63, slip law at W/hc 3 to 31.5;
Cattania aging at W/L_inf 3.2 to 13, slip law at 3.2 and 6.5. Not
finished: Wu & Chen slip law at 63, W/hc 120, Cattania slip law at 13,
W/L_inf 26.

Wu & Chen (2014) model A, aging law:

| W/hc | this model | Wu & Chen |
|---|---|---|
| 3, 8.6 | periodic full ruptures nucleating at the centre | periodic, centre origin (W/hc 2.4 to 11) |
| 20 | periodic full ruptures nucleating near the ends | periodic, side origin (W/hc 11 to 66) |
| 31.5 | full ruptures every 5.2 events with partial ruptures of 0.2 to 0.5 W between them, irregular (CV 0.78) | side origin; partial ruptures (narrow distribution) only above W/hc of about 66 |
| 63 | a periodic three-event cycle (full rupture plus two of about 0.12 W) | narrow event distribution |

The centre-to-side transition is reproduced. Partial ruptures appear
already at W/hc 31.5 here, below the boundary of about 66 in Wu & Chen;
whether this comes from the buried half-space standing in for their
whole space, the catalog threshold, the initial state or the run
length has not been tested.

Wu & Chen model A, slip law: periodic full ruptures at W/hc 3, 8.6 and
20. Nucleation is central at 3 and 8.6 and at 0.27 W at 20, between the
centre and the end position of the aging law at the same W/hc (0.15 W),
i.e. the move of the nucleation towards the ends is delayed with the
slip law (antiplane_nucleation.png; the f_cent column of the table counts
positions within 0.25 W of the centre and so lists 0.27 W as central).
At W/hc 31.5 the slip law gives periodic full ruptures nucleating near
the end (0.86 W), where the aging law already has irregular partial
ruptures. Both transitions (centre to side origin, and onset of partial
ruptures) thus occur at larger W/hc with the slip law, as in the 2.5-D
runs (../single_fault_regimes).

Cattania setup, slip law: periodic full ruptures at W/L_inf 3.2 and 6.5,
as with the aging law. The slip-law catalogs again contain many short
bursts right after large events (247 and 688 catalog entries, merged
into 14 and 22 events).

Cattania (2019), aging law:

| W/L_inf | this model | Cattania |
|---|---|---|
| 3.2 | periodic full ruptures | periodic full ruptures (her Fig. 1a) |
| 6.5 | periodic full ruptures | onset of partial ruptures at W/L_inf 5 +- 2 |
| 13 | two ruptures per cycle: a full one and one of 0.38 W, alternating intervals (0.62 and 1.38 times the mean) | two ruptures per cycle (her Fig. 1b, d) |

Her estimate is 0.45 (W/L_inf)^1/2 = 1.6 events per cycle at 13, and
she found two in her simulations, as here. Positions along this fault
run from the free surface (0) to the bottom of the velocity-weakening
part (1); at W/L_inf 3.2 the nucleation zone is a large fraction of the
fault, so the position is only loosely defined.

## Parameter-space maps (`run_parameter_maps`, `plot_parameter_maps.py`)

`run_parameter_maps [rundir] [laws] [ncycles] [maxcells]` scans the Wu &
Chen setup over W/hc = 2, 5, 11, 22, 45, 90 and b/a = 1.33, 2, 3 (a
varied at fixed b), and the Cattania setup over W/L_inf = 2, 4, 8, 16,
32 and a/b = 0.5, 0.625, 0.75, 0.875, for each state evolution law
(aging, slip, PRZ selectable), cells L_b/4, 30 recurrence estimates per
run; cases needing more than maxcells cells are skipped and completed
runs are skipped on restart. `plot_parameter_maps.py` classifies each
run as stable, slow, periodic with central or side nucleation, or with
partial ruptures of a narrow or wide size distribution, and draws
parameter_map_wuchen.png in the layout of Wu & Chen's regime figure
(with their boundaries W/h* = 0.6, 1.2, 5.5, 33, 105, h* = hc/(1 - a/b))
and parameter_map_cattania.png in the layout of Cattania's Figure 2a
(events per full-rupture cycle by colour); the numbers are in
parameter_maps.txt.

`compare_papers.py [rundir]` makes the comparison quantitative:
for the Wu & Chen setup it predicts each run's regime from their
boundaries and counts agreement (runs within 10 percent of a boundary
match either neighbouring regime), and gives per b/a the W/h* interval
of each transition between neighbouring grid points; for the Cattania
setup it gives the ratio of events per full-rupture cycle to her
estimate 0.45 (W/L_inf)^1/2 and the W/L_inf interval of the onset of
partial ruptures. It writes comparison.txt, comparison_wuchen.png and
comparison_cattania.png. plot_antiplane.py caches the events of each run
in events.npz, so that the large series files can be deleted once a run
has been analysed.

### Results (2026-11-09/10, aging and slip law; single runs)

Runs: Wu & Chen setup W/hc 2 to 45 at b/a 1.33, 2, 3 and W/hc 90 at
b/a 2 (aging; W/hc 90 at b/a 3 and with the slip law were not finished);
Cattania setup W/L_inf 2 to 32 at a/b 0.5 and 0.625, 2 to 16 at 0.75,
2 and 4 at 0.875 (the slip-law run at a/b 0.75, W/L_inf 16 was stopped
at 13.6 of 38.7 yr because of cost). PRZ runs are not included yet.

Wu & Chen (2014) setup, regime predicted by their boundaries:

| law | runs in the predicted regime |
|---|---|
| aging | 10 of 16 |
| slip | 13 of 15 |

Transition intervals in W/h* (their values 1.2, 5.5, 33):

| b/a | law | events become seismic | nucleation moves to the ends | partial ruptures begin |
|---|---|---|---|---|
| 1.33 | aging | 0.5 to 1.2 | 2.8 to 5.5 | 5.5 to 11 |
| 1.33 | slip | 0.5 to 1.2 | 2.8 to 5.5 | not up to 11 |
| 2 | aging | 1 to 2.5 | 5.5 to 11 | 11 to 22 |
| 2 | slip | at or below 1 | 5.5 to 11 | not up to 22 |
| 3 | aging | at or below 1.3 | 3.3 to 7.3 | 7.3 to 15 |
| 3 | slip | at or below 1.3 | 3.3 to 7.3 | not up to 30 |

With the aging law the first two transitions agree with Wu & Chen within
the grid resolution at all three b/a, and partial ruptures begin a
factor of 2 to 5 below their W/h* = 33; a wide size distribution
appears at W/h* 30 (b/a 3) and 45 (b/a 2), below their 105. In W/hc
the onset of partial ruptures lies between about 20 and 45 at all three
b/a, i.e. it does not follow a constant W/h* as their boundary does
(the grid is coarse). Possible reasons for the earlier onset, none
tested: the catalog counts all events above 1 mm/s while their
classification may consider only larger ones, the buried half space
standing in for their whole space, and the run length (30 recurrence
estimates); cell size is unlikely, since the b/a 2 scan at their
resolution (hc/10) gave the same.

With the slip law, events become seismic at a smaller W/h* (b/a 2:
events at W/h* 1, where the aging law has slow events, consistent with
its smaller nucleation size), the nucleation moves to the ends at the
same W/h*, and no partial ruptures occurred up to the largest W/h* run
(30 at b/a 3), i.e. at least a factor 2 to 4 later than with the aging
law.

Cattania (2019) setup:

| a/b | law | events per cycle / her estimate (W/L_inf 8, 16, 32) | onset of partial ruptures |
|---|---|---|---|
| 0.5 | aging | 0.79, 0.90, 0.79 | W/L_inf 8 to 16 |
| 0.625 | aging | 0.79, 1.11, 0.81 | 8 to 16 |
| 0.75 | aging | 1.03, 1.16 | 4 to 8 |
| 0.5 | slip | 0.79, 0.56, 0.39 | not up to 32 |
| 0.625 | slip | 0.79, 0.56, 0.39 | not up to 32 |
| 0.75 | slip | 0.79, 0.60 | not up to 16 |

(her onset: 5 +- 2; a ratio of 0.79 at W/L_inf 8 means one event per
cycle where her estimate is 1.27.) With the aging law her estimate
matches the number of events per cycle to about 20 percent where
partial ruptures occur, and her onset at a/b 0.75; at a/b 0.5 and 0.625
the onset lies between W/L_inf 8 and 16, later than her 5 +- 2. With
the slip law no partial ruptures occurred up to W/L_inf 32.

In both setups the slip law therefore leaves the onset of seismic
behaviour and the move of nucleation to the ends nearly where the
aging law has them (or earlier), but delays partial ruptures by at least
a factor of 2 in fault size, consistent with the 2.5-D strips
(../single_fault_regimes). The slip-law runs at a/b >= 0.75 are
expensive (up to 20 min each on one core).

Wu & Chen setup, aging law, against their map:

- the slow-event region (no catalog event, slip rate well above the
  plate rate) at W/hc 2 for b/a 1.33 and 2 and periodic events at b/a 3
  fall where their boundaries put them (at b/a 1.33 their stable region
  begins; the run there had maximum slip rates of about 10 vpl);
- the change from central to side nucleation follows their W/h* = 5.5
  curve at all three b/a (centre at W/hc 11 for b/a 1.33 and 2, side at
  11 for b/a 3, side at 22 for all);
- partial ruptures (narrow distribution) appear at W/hc 22 for b/a 3,
  i.e. W/h* of about 15, well before their W/h* = 33 curve, as in the
  single-parameter scan at b/a 2.

Slip law, same setup: seismic events already at W/hc 2 for b/a 2, where
the aging law gives only slow events (consistent with the smaller
nucleation size of the slip law, Ampuero & Rubin 2008); the centre-side
change at the same places as with the aging law; and periodic full
ruptures at W/hc 22 for b/a 3, where the aging law has partial ruptures.
The slip law thus makes small faults seismic earlier but delays partial
ruptures.

Cattania setup: at W/L_inf 2 and 4 all runs with both laws give periodic
full ruptures. With the slip law at a/b 0.75 the catalog listed a
one-cell threshold crossing about 60 s before each large event; with
the one-hour merging these are part of the large events.

## Cattania setup for the aging, slip and two-state slip laws

`gen_antiplane.py --law two` uses two state variables with the slip law:
b1 + b2 = b, b2/b1 = 0.84, dc2 = dc1/0.048 (the ratios of Gu et al. 1984
and Becker 2000), and dc1 such that the critical stiffness of the
two-state spring slider equals the one-state (b - a) sigma/dc, so that
W/L_inf refers to the same linear-stability length as for one state
variable (a choice, not an established equivalence). Cells are
G dc2/(b2 sigma)/32, which does not resolve the first variable's
G dc1/(b1 sigma); in the 2.5-D check this kept the character of the
sequences but changed recurrence by about 5 percent.

`run_cattania_space.py [rundir] [maxcells] [ncycles] [laws]` runs the
Cattania setup at W/L_inf = 2, 3, 4, 6, 8, 11, 16, 23, 32, 45, 64 and
a/b = 0.5, 0.625, 0.75 for the aging, slip and two-state laws in order of
estimated cost, skipping runs above maxcells cells (default 800) and
completed runs. `plot_cattania_space.py [rundir]` draws
cattania_space.png in the layout of Cattania's figures (W/L_inf-a/b map
coloured by events per full-rupture cycle, events per cycle against
0.45 (W/L_inf)^1/2, CV of interevent times, rupture-length survival
functions) and writes cattania_space.txt.

Results (2026-11-10; single runs, 30 recurrence estimates each; runs
completed: aging to W/L_inf 45, 32, 16 for a/b 0.5, 0.625, 0.75; slip to
32, 32, 16; two-state to 16, 8, 4, limited by cost):

| a/b | law | first W/L_inf with partial ruptures | events per cycle at the largest W/L_inf run | CV there |
|---|---|---|---|---|
| 0.5 | aging | 16 (11: irregular timing, CV 0.18, all full) | 1.96 (45) | 0.43 |
| 0.625 | aging | 11 | 2.06 (32) | 0.45 |
| 0.75 | aging | 8 | 2.08 (16) | 0.45 |
| 0.5 | slip | none up to 32 | 1.00 | 0.00 |
| 0.625 | slip | none up to 32 (32: CV 0.15, all full) | 1.00 | 0.15 |
| 0.75 | slip | 16 (1.08 events per cycle) | 1.08 (16) | 0.16 |
| 0.5 | two-state | none up to 16 | 1.00 | 0.00 |
| 0.625 | two-state | none up to 8 | 1.00 | 0.00 |
| 0.75 | two-state | none up to 4 (CV 0.10 at 4) | 1.00 | 0.10 |

- Aging law: the onset of partial ruptures moves from W/L_inf 6 to 8 at
  a/b 0.75 to 11 to 16 at a/b 0.5, against a single band of 5 +- 2 in
  Cattania's figure. Above the onset the number of events per cycle
  stays near 2 from W/L_inf 16 to 45 while her estimate grows to 3.0,
  so the ratio drops from about 1 to 0.65 at 45. CV grows to 0.4 to 0.5
  (she reports up to 1.5 at W/L_inf 415, outside this range).
- Slip law: periodic full ruptures over almost the whole range; the
  first partial ruptures at a/b 0.75 and W/L_inf 16, i.e. a factor of
  about 2 later than with the aging law at the same a/b, and none up to
  W/L_inf 32 at a/b 0.5 and 0.625. Cattania noted that the slip law's
  smaller nucleation length favours partial ruptures while its smaller
  fracture energy favours full ones; here the second effect dominates
  throughout.
- Two state variables (slip law): periodic full ruptures wherever runs
  were affordable, like the one-state slip law; the largest W/L_inf
  could not be reached because of the finer cells.
