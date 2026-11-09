# Single-fault sweep: friction assumptions and orientation in the far-field loading

`run_single_sweep [rundir] [tmax_yr] [variants] [dx_div_two]` runs one
strike-slip strip (L 60
km, W 15 km, 10 percent velocity-strengthening ends with a 0.03, sigma_0
50 MPa, far-field stressing rate 21.5 kPa/yr, normal stress evolving
with the far-field rate, limiter at 5 MPa) for 3000 yr in 30
configurations:

- friction (`variants`): one state variable, aging law, a 0.01, b 0.02
  (`one`); the same with the slip law (`slip`); or two state variables,
  slip law, a = b1 = 0.01, b2 = 0.84 a, dc2 = dc/0.048 as in Becker (2000)
  (`two`);
- L/h* of 10, 20, 30 (dc varied). For two state variables, whose
  weakening acts over dc2 because b1 = a, dc is chosen so that the
  spring-slider critical stiffness (0.0803 a sigma/dc) equals the
  one-state (b - a) sigma/dc at the same L/h*; this matching is a
  choice, not an established equivalence;
- strike relative to the optimal orientation -20, -10, 0, 10, 20 degrees
  (angle to S1 of 25 to 65 degrees; negative unclamping, positive
  clamping, cf. Becker & Schmeling 1998).

Cells: L_b/4 at sigma_0 for one state variable (L_b/2 or less at the end
of the +20 degree runs, where sigma reaches about 90 MPa); for two state
variables a cell size of G dc2/(b2 sigma) divided by `dx_div_two`
(16 by default), which leaves the first variable's L_b = G dc/(b sigma)
under-resolved (0.6 to 0.7 cells at 16, 1.3 at 32).

`plot_sweep.py [rundir]` makes sweep_timeseries.png (mean tau/sigma of
the velocity-weakening cells per run), sweep_recurrence.png (interval to
the next event against time, by orientation), sweep_summary.png and
sweep_summary.txt. Recurrence statistics use events with Mw >= 6: the
two-state runs have, after nearly every large event, one to four
re-ruptures of one or two cells within seconds (Mw about 4.3 to 4.7),
listed as separate events in the catalog, which would otherwise
dominate the interval statistics.

## Runs shown

The figures and sweep_summary.txt in this directory combine the aging-law
runs (`run_single_sweep dir 3000 "one two" 16`), the slip-law runs
(`run_single_sweep dir 3000 "slip two" 32`) and two-state runs at
`dx_div_two` 32 for L/h* 10 and 20 but 16 for L/h* 30 (the finer runs
took about 14 min each there and were not completed).

## Results (2026-11-08, single runs, 3000 yr each)

- Orientation dominates within 3000 yr. Unclamped faults (-10, -20)
  have recurrence intervals that shorten steadily and smaller events,
  clamped faults (+10, +20) lengthening intervals; at -20 degrees the
  normal stress approaches the limiter near the end of the run and
  events become very frequent (one state) or give way to steady
  sliding (two states, L/h* 10). These are the trends of Becker &
  Schmeling (1998), now with rate-and-state friction.
- After removing the linear trend, the Mw >= 6 intervals are nearly
  periodic in every run (detrended CV 0.000 to 0.114). The largest values,
  0.07 to 0.11 at L/h* 20 to 30 with one state variable, come from
  alternating event sizes: at L/h* 20 and optimal orientation the
  sequence alternates intervals of about 338 and 404 yr and Mw 7.37 and
  7.42 (period 2), as in the 10,000 yr pilot run.
- No partial ruptures (Mw >= 6 with less than half the velocity-weakening
  area) occurred under this far-field loading at L/h* up to 30, for
  either friction law. The single-fault pilot (../single_fault_pilot)
  shows that uniform backslip loading produces partial ruptures and
  irregular sequences at L/h* 20 and 30.
- Two state variables (with the stiffness matching above) give shorter
  recurrence (about 306 against 404 yr at L/h* 20, optimal orientation)
  and the small re-ruptures, but not a different regime in these runs.
- The one-state slip law gives nearly the same large-event sequences as
  the aging law (7, 8 and 10 Mw >= 6 events at the optimal orientation
  for L/h* 10, 20, 30, against 7, 8 and 11), the same orientation trends,
  detrended CV of at most 0.073, and no period 2 at L/h* 20 (CV 0.000).
- Both slip-law variants add many small events (Mw about 3.8 to 5.0):
  9 to 98 per run for one state variable, rising with L/h*, and 9 to 49
  for two. In the runs checked (L/h* 20, optimal orientation) they occur
  within seconds to about 40 min after a large event, at the very ends
  of the strip (outermost 1 to 3 km, inside the velocity-strengthening
  ends), i.e. bursts of post-seismic slip at the strip tips that exceed
  the 1 mm/s catalog threshold. Their number and size depend on the
  resolution: for the one-state slip law at L_b/8 there were 26 of
  Mw 4.1 to 4.6 against 36 of Mw 4.3 to 5.0 at L_b/4, while the eight
  large events kept their onsets within 14 yr over 2600 yr. They are
  not used in the statistics above (Mw >= 6 only).
- The permutation entropy values in the table rest on 2 to 24 large
  events and on peak stresses that drift with the normal stress; they
  should not be read as measures of complexity here.

Resolution check (`run_resolution_check`, L/h* 20, optimal
orientation): the one-state run at L_b/8 reproduces the period-2
sequence (same onsets within 10 yr over 2600 yr, same alternating Mw);
the two-state run at half the cell size has the same character but a 5
percent shorter recurrence (291 against 306 yr) and more small
re-ruptures, so the two-state runs are not converged in recurrence.

The 3000 yr runs hold only 2 to 24 large events, too few for regime
statistics of the stationary (optimal) cases; the drifting cases are
not stationary by construction.
