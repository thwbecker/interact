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
  antiplane_summary.txt. Events starting within 60 s of the end of the
  previous one are merged, and the first quarter of each run is
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
