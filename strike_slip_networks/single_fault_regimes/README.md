# Single-fault regimes against W/h*: aging, slip and two-state slip law

`run_regimes [rundir] [tmax_yr] [variants] [Xlist] [extra rsf_solve options]`
runs one strike-slip strip (60 km long, 15 km deep, 10 percent
velocity-strengthening ends with a 0.03, sigma 50 MPa, optimal
orientation, constant normal stress) loaded by uniform plate-rate
backslip, so that the creeping ends load the locked velocity-weakening
patch from its edges, as in the 1-D antiplane models of Wu & Chen (2014,
loaded by steady slip around the fault) and Cattania (2019, loaded by
down-dip creep). W = 48 km is the velocity-weakening length and
h* = (2/pi) G/(1-nu) dc b/((b-a)^2 sigma), which is twice the L_inf of
Rubin & Ampuero (2005) and corresponds to the h* of Wu & Chen at b/a = 2.
dc is set from W/h* = 1, 2, 4, ..., 64; b/a = 2 throughout.

- one: one state variable, aging law, a 0.01, b 0.02; cells L_b/4,
  i.e. L_inf/6.8
- slip: the same with the slip law
- two: two state variables, slip law, a = b1 = 0.01, b2 = 0.84 a,
  dc2 = dc1/0.048 (Becker 2000), dc1 = 0.0803 dc so that the
  spring-slider critical stiffness equals the one-state one (a choice,
  see ../single_fault_sweep/README.md); cells G dc2/(b2 sigma)/32, about
  half the one-state size

`plot_regimes.py [rundir]` merges events that start within one day of
the end of the previous one (this absorbs the bursts of one or two
cells at the strip tips seconds to minutes after large slip-law events,
which the catalog lists separately), discards the first 20 percent of
each run, and writes:

- regimes_summary.png: events per full-rupture cycle (full rupture:
  length at least 0.9 W), CV of the intervals of all events and of full
  ruptures, and the fraction of full ruptures against W/h*, with the
  estimate alpha = 0.45 (W/L_inf)^1/2 of Cattania (2019) and the regime
  boundaries of Wu & Chen (2014) at W/h* = 1.2, 5.5, 33 and 105;
- regimes_bifurcation.png: every normalized interevent time and every
  rupture length / W against W/h*, one column per law;
- regimes_intervals.png: histograms of the normalized interevent times;
- regimes_lengths.png: survival functions of rupture length / W;
- regimes_spacetime.png: log10 slip rate along strike for selected runs;
- regimes_summary.txt: the numbers.

## Runs and results (2026-11-08/09, single runs)

Runs: 8000 yr for W/h* 1 to 16; 4000 yr for W/h* 32 (aging, slip law).
The two-state run at W/h* 16 was stopped at about 4250 yr because of its
cost, and two state variables at W/h* 32 and 64 and the one-state runs
at W/h* 64 are not included here (cost of 30 minutes or more per run on
one core).

| W/h* | aging law | slip law | two states, slip law |
|---|---|---|---|
| 1 | no events (stable) | no events | no events |
| 2, 4, 8 | periodic full ruptures (CV 0.000) | periodic full ruptures | periodic full ruptures (CV up to 0.003) |
| 16 | partial ruptures: 3.2 events per full rupture, CV 0.68, intervals of full ruptures irregular (CV 0.57) | periodic full ruptures | periodic full ruptures (4250 yr) |
| 32 | 3.1 events per full rupture, CV 0.78, but full ruptures periodic (CV 0.000): a repeating three-event pattern with a small rupture (0.11 W) shortly after one of them | periodic full ruptures (43 of 43, interval 75.3 yr) | not run |

- With the aging law, partial ruptures appear between W/h* 8 and 16
  (W/L_inf 16 to 32), and the number of events per cycle (3.1 to 3.2)
  is of the order of Cattania's estimate (2.5 to 3.6 for W/L_inf 32 to
  64). The transition lies below the W/h* of about 33 at which Wu &
  Chen (2014) found partial ruptures; their 1-D antiplane fault is
  loaded kinematically on both sides, ours is a strip of finite depth
  in a half space loaded through creeping ends, and which of these
  differences matters has not been tested. Between W/h* 16 and 32 the
  sequence goes from irregular to a periodic three-event pattern, i.e.
  the route is not monotonic, as also found by Sun & Becker (2026) for
  dip.
- With the slip law, one or two state variables, all W/h* up to 32 (one
  state) and 16 (two states) give periodic full ruptures. Cattania
  (2019) noted that the slip law's smaller nucleation length should
  favour partial ruptures but its smaller fracture energy full ruptures;
  in these runs the second effect dominates. The slip-law runs instead
  have many short bursts at the strip tips right after each large event
  (for example 429 catalog entries merged into 43 events at W/h* 32);
  their number depends on resolution (../single_fault_sweep/README.md).
- Resolution: L_inf is resolved by 6.8 cells (one state) and the slip-law
  sequence at W/h* 20 did not change at L_b/8 in the earlier sweep;
  resolution checks at the W/h* where the regime changes have not been
  made here.
