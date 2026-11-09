# Single-fault pilot (PLAN.md 10.5)

All results below are from rsf_solve runs made on 2026-11-08 in the
sandbox (one core, PETSc 3.19.6), single runs unless noted; they are
first indications, not converged statistics.

## P0, P1: two-state spring slider

See ../../test_slider/two_state (README.md, slider2sv.py,
bifurcation_scan.py). rsf_solve with `-rsf_two_state` on a single patch
reproduces the nondimensional spring slider (Gu et al. 1984; Becker
2000) where the latter is bounded. The kappa' scan (1e5 dc/v0, analysis
of the last 5e4, start at 1.001 v0, bifurcation_scan.py, figure
p1_bifurcation.png there) gives 2 distinct stress extrema down to
kappa' = 0.87, 4 from 0.865 to 0.8574, 8 at 0.857 and 0.856, then 15 to
28 near 0.8552 to 0.8545, a periodic window with 12 at 0.854, and many
(47 to 89, PE of the maxima 5.4 to 6.0 bits) from 0.8536 to 0.852. At
kappa' = 0.8555, 0.851 and 0.85 the run settled into large stick-slip
events (stress drops of about 17 a sigma, slip rates near 1 m/s, made
finite by radiation damping) instead of the small-amplitude attractor.
At 0.8555 this depends on the start: 1.0001 v0 stayed on the small
attractor, 1.01 and 0.99 v0 had large events early and returned to it,
1.001 v0 stayed in large events. At 0.85 all four starts tried ended in
large events. The scan points are too coarse to compare the bifurcation
values with Becker (2000) and Sun et al. in detail.

## P2: single strip, L/h* scan

p2.sh: one fault, L 60 km, W 15 km, a 0.01, b 0.02, 10 percent VS ends
(a 0.03), sigma 50 MPa, aging law, L_b/4, L/h* set through dc, 10,000 yr,
series at every 8th dense step. Two loadings (gen_system.py `--load`):
farfield (uniform far-field stressing rate, as in the other examples)
and backslip (uniform plate-rate backslip, so that the creeping ends load
the locked part from its edges). Measures from ../complexity.py over the
last three quarters of each run (first fault group):

| run | events | mean T [yr] | CV | cluster-N | distinct extrema | PE [bits] | D2 |
|---|---|---|---|---|---|---|---|
| L/h* 5, farfield | 15 | 482 | 0.00 | 1 | 2 | 0.0 | 1.01 |
| L/h* 5, backslip | 28 | 269 | 0.00 | 1 | 2 | 0.0 | 1.00 |
| L/h* 10, farfield | 17 | 446 | 0.00 | 1 | 2 | 0.0 | 1.00 |
| L/h* 10, backslip | 32 | 233 | 0.00 | 1 | 2 | 0.0 | 1.00 |
| L/h* 20, farfield | 20 | 373 | 0.09 | 2 | 4 | 1.0 | 1.27 |
| L/h* 20, backslip | 90 | 84 | 0.69 | 8 | 88 | 5.1 | 1.49 |
| L/h* 30, farfield | 28 | 274 | 0.00 | 1 | 2 | 0.0 | 1.00 |
| L/h* 30, backslip | 111 | 68 | 0.74 | 3 | 4 | 3.8 | 1.04 |

p2_overview.png shows tau/sigma over the last 4000 yr and the return
maps. Under far-field loading the strip is periodic at L/h* 5, 10 and 30
and alternates two event sizes (period 2) at L/h* 20, i.e. the sequence
is not monotonic in L/h*. Under backslip loading partial ruptures appear:
irregular at L/h* 20, and a repeating three-event pattern at L/h* 30
whose CV is high because it contains near-coincident event pairs
(intervals near zero in the return map; whether these are separate
ruptures or one rupture split by the catalog's velocity threshold has
not been checked). The analysis windows hold 15 to 110 events, fewer
than the 100 cycles planned; the far-field runs need about four times
longer runs.

Not done in this session: L/h* 40 (started), the twin run of L/h* 20
backslip with a 1 Pa kick (`--ic_kick 1`) for divergence, 40,000 yr
far-field runs at L/h* 20 and 30 (p2c.sh), resolution checks, P3, P4.

The L/h* used here is h* = (2/pi) G/(1-nu) dc b/((b-a)^2 sigma) with L
the full strip length; it differs by constant factors from the
W' = W/L_inf of Sun & Becker (2026), so the threshold values are not
directly comparable.
