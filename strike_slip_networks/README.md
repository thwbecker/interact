# strike_slip_networks: quasi 2-D (2.5-D) strike-slip earthquake cycle ensembles with rsf_solve

Driver scripts for generating many rate-and-state earthquake cycle
models of a strike-slip fault system with `rsf_solve` from
https://github.com/thwbecker/interact, intended as training data for
machine learning (transformer or graph neural network) on the dynamics
of such systems.  Status as of 2026-11-08: first version, tested in a
sandbox on the configurations listed below only.  The implementation
plan is in PLAN.md.

## Model

Each fault is a strip of rectangular Okada patches along strike, one
patch deep from the free surface to the locking depth W.  Slip is thus
a function of the along-strike coordinate only, while stress transfer
is the 3-D half-space solution for finite depth extent (this is the
"2.5-D" part; a true 2-D plane strain model would have no depth
scale).  All faults are vertical and slip in mode II (strike slip,
`-rsf_slip_mode 0`) in their own strike direction (`--strikes`, one
value per fault), with the regularized aging law by default.

Loading is a uniform far-field horizontal stressing-rate tensor,
resolved on each fault into an along-strike shear rate tau_dot
(positive = left-lateral drive, interact's slip sign) and a normal
rate sigma_dot (compression positive), written per cell to
`-rsf_stress_rate_file` and used with `-vpl 0`, i.e. no backslip.  The
tensor is given by principal rates `--sr_s1`, `--sr_s3` (compression
positive, vertical rate 0) and the S1 azimuth `--sr_az`; the default
is pure shear with magnitude tau_ref = `--taudot` or G v_eq/(pi W)
and S1 at azimuth 135, so that a strike-0 fault gets tau_dot = tau_ref
and no normal rate.  With alpha the strike measured from the S1 axis,

    tau_dot   = (s1 - s3)/2 sin(2 alpha)
    sigma_dot = (s1 + s3)/2 - (s1 - s3)/2 cos(2 alpha)

A fault 45 degrees clockwise of S1 is driven left-lateral, 45 degrees
counterclockwise right-lateral (as the Garlock and San Andreas faults
under roughly north-south S1).  The handedness was derived from the
strain field and interact's documented convention (positive strike
slip = left-lateral); the resolution is checked to be consistent with
interact's patch frames by the rotation test below, but no independent
check of the absolute left/right-lateral sense has been made.

Normal stress evolves by default (`--sigma_evolve 1` puts
`-calc_sigma_dot -limit_sigma` into run.opts): the far-field sigma_dot
(switch off with `--ff_normal 0`) plus the normal stress changes from
slip on the other faults through the normal interaction matrix (zero
on a planar vertical strike-slip fault itself).  `--sigma_interact 0`
(rsf_solve `-sigma_dot_interact 0`) keeps the far-field rate and drops
the changes from slip on the faults, which removes the normal stress
drift at stepovers described under Observations.  The limiter keeps
sigma within `--sigma_min` and `--sigma_max` (5 and 200 MPa).  Two
consequences of a constant far-field normal rate:

- A fault oblique to the optimal orientation has no steady state: its
  normal stress drifts linearly (13.8 kPa/yr for 20 degrees off, 18.6
  for 30 degrees in the reference case), so over a 3000 yr run sigma
  changes by tens of MPa, until the limiter clips it.  This is
  intended (recurrence variations from evolving normal stress, cf.
  Becker and Schmeling, 1998, GJI); the secular drift means oblique
  faults are non-stationary, so for stationary statistics keep runs
  short relative to sigma/sigma_dot, or use `--ff_normal 0`.  In a
  single-fault test, 10 degrees of obliquity changed the recurrence
  from 219 yr (horizontal) to 166-187 yr (unclamping) and 305-378 yr
  (clamping) within 1200 yr, with Mw drifting by about 0.1.
- Under pure shear and f0 = 0.6 the Coulomb stressing rate
  tau_dot - f sigma_dot becomes negative for strikes more than about
  27 degrees clockwise of the optimal orientation (sigma increasing);
  such faults never fail (a strike-30 test produced no event in 300 yr
  while sigma rose by 5.6 MPa, as expected).  Faults rotated the other
  way unclamp and fail more often.

`--rot` rotates the whole system (geometry, strikes and S1 azimuth);
a pure rotation leaves the results unchanged.  The default is
`--rot 90`, which lays the faults out left to right (strike 90 for
`--strikes 0`, S1 at azimuth 225, left-lateral); the default right
step then goes toward -y.

Further consequences of stress-rate loading worth keeping in mind:

- The long-term slip rate is an outcome, not an input.  For a finite
  stress-driven fault it is below the "equivalent" deep slip rate
  v_eq used to set tau_dot = G v_eq/(pi W) (the default; `--taudot`
  sets the rate directly).  In the 60 x 15 km reference case the
  long-term rate was about half of v_eq.
- Velocity-strengthening cells have no imposed creep rate; they creep
  at the rate at which their own stress transfer balances tau_dot, and
  are otherwise driven by the neighbouring events.
- The whole fault is loaded uniformly, which appears to favour
  characteristic full-fault ruptures (see Observations).

For nfault > 1 the faults are arranged en echelon: fault k is shifted
by k times `--step_x` perpendicular to strike and starts at
y = k (L - overlap) (before the rotation), so `--overlap` > 0 gives
overlapping tips and < 0 a gap.  `--stepover right` (default) places
each next fault to the right when looking along strike, `--stepover
left` to the left.  For left-lateral slip a right step is restraining
and a left step releasing, and the reverse for right-lateral slip;
gen_system.py prints the type of each step and writes it to
params.json (`step_types`).  Faults given a non-zero entry in
`--strikes` are rotated about their centre by default, or with
`--pivot inner` about their end nearest the middle of the system, so
that the geometry at the inner tips stays that of the parallel case.
`--load backslip` replaces the far-field rate by uniform plate-rate
backslip (rsf_solve -vpl), and `--ic_kick` adds a shear stress to one
cell's initial value for twin runs.

`--layout chain` instead places the faults tip to tip along a
segmented trace: fault k starts `--gap` beyond the end of fault k-1
and runs along its own strike from `--strikes`; `--lengths` gives one
length per fault (chain only).  By default the gap points along the
previous fault; `--gap_dir bisect` points it along the bisector of the
two strikes, which keeps point-symmetric chains symmetric.  Lists with a negative first value
have to be written as `--strikes=-15.5,0,-15.5`.

Geometry resolution: the cell size is set from the cohesive zone
Lb = G dc/(b sigma) as Lb/`--cells_per_lb` (default 3), consistent
with the checks in `rsf_solve`.  The nucleation length reported is
h* = 2/pi G/(1-nu) dc b/((b-a)^2 sigma) with the in-plane modulus.

## Files

    gen_system.py   writes geom.in, rsf.dat (a b), ic.in (steady state
                    at v_init with log-normal perturbation), srate.in
                    (tau_dot), params.json, run.opts; all model
                    parameters are its command line options (--help)
    run_model       name tmax_yr ncore [gen args] [-- rsf_solve args]
                    creates the run directory, generates the inputs
                    and runs rsf_solve with catalog, per-fault
                    monitors and the all-cell series output
    run_sweep       prefix nrun tmax_yr seed0 [nparallel]: random
                    parameter draws (edit draw_params) and one
                    run_model per member, table of parameters in
                    prefix.table
    rsf_series.py   reader for series.bin, catalogs and geometry;
                    per-fault mean/std, uniform resampling, event
                    windows, export_training() to a compressed npz
    plot_faults.py  per-fault figure: map view with the sense of shear
                    (arrows from the sign of the resolved shear rate)
                    and the S1 axis on top, then one column per fault
                    with rows log10 V (mean, std band, max), tau with
                    sigma on the right axis, log10 theta, and
                    cumulative slip with event Mw; --spacetime adds a
                    log10 V along-strike image per fault
    plot_run.py     older overview figure, all faults in one panel per
                    quantity, plus space-time image of fault 0 and map
    PLAN.md         implementation plan and status
    complexity.py   cycle complexity measures (CV, cluster-N, distinct
                    extrema, permutation entropy, D2, TwoNN, PCA)
    single_fault_pilot/
                    L/h* scan of a single strip, far-field against
                    backslip loading (README.md there)
    COMPILE.md      build notes for rsf_solve (PETSc versions, flags,
                    HACApK and hmmvp, sandbox specifics)
    examples/       make_examples (runs and plots the example systems)
                    and its figures
    becker_and_schmeling98/
                    rsf_solve runs that redo the models of Becker and
                    Schmeling (1998, GJI) with rate-and-state friction,
                    the scripts that made them, figures, and a
                    README.md comparing them with the paper

`run_model` sources an optional `config.sh` next to it (not part of
the repository) that may set the variables `rsf_solve` and `mpirun`,
e.g.

    rsf_solve=$HOME/progs/src/interact/bin/rsf_solve
    mpirun="mpirun --oversubscribe"

## Changes to rsf_solve made for this work

All are new options; with their defaults the output of earlier runs
is reproduced byte for byte.

    -rsf_two_state <bool>, -rsf_b2 <b2>, -rsf_dc2 <m>
                             second state variable (Ruina 1983; Gu et
                             al. 1984; Becker 2000), both following
                             -state_law; the series then has a sixth
                             field, log10 theta2.  gen_system.py: --b2,
                             --dc2, --state_law.  Checked against the
                             nondimensional spring slider in
                             ../test_slider/two_state.  Not with -imex
                             or -fd_mode.

    -sigma_dot_interact <bool>  with -calc_sigma_dot: 1 (default)
                             normal stress also changes from slip on
                             the faults; 0 only from the loading rates
                             (backslip, -rsf_stress_rate_file).  The In
                             matrix is still assembled, since the
                             backslip loading uses it, so memory is not
                             saved, but one matrix-vector product per
                             right-hand side evaluation is.  Not
                             implemented with -fd_mode (stops with an
                             error).

With `-sigma_dot_interact 0` the normal stress of each fault in the
angled example followed sigma_0 + sigma_dot t to float precision (in
both the explicit and the `-imex` integration), stayed at exactly
50 MPa on parallel faults, and 2 MPI ranks gave the same series as 1.

### Series output

`rsf_solve` wrote per-cell fields only as one small binary per frame
and group, and the state variable was only available through
per-station ASCII files.  A single compact binary stream, written on
rank 0, was added for this work:

    -rsf_series_file <f>     header int32 magic 0x52534631, int32 ncell,
                             int32 nfield (5); then per record a float64
                             time [s] and five float32[ncell] blocks:
                             log10|v| [m/s], tau [MPa], log10 theta [s],
                             slip [m], sigma [MPa] (compression positive)
    -rsf_series_dt_yr <dt>   cadence while quiescent (default 1 yr)
    -rsf_series_vdense <v>   max|v| above which steps count as dense
                             (default 1e-6 m/s)
    -rsf_series_stride <n>   every n-th accepted step in a dense phase,
                             counted from its first step (default 1)
    -rsf_series_mode <m>     0 (default): the three options above;
                             1: change based, see below
    -rsf_series_dlogv <d>    mode 1: log10|v| change (default 0.1)
    -rsf_series_dstress <d>  mode 1: tau or sigma change [MPa] (default 0.1)
    -rsf_series_dtmax_yr <dt> mode 1: maximum interval (default 5 yr)
    -rsf_series_final <bool> also record the final accepted step

All new options default to the earlier behaviour; with the defaults the
series file was byte-identical to the one written before they were
added.  In mode 1 a record is written at the first accepted step at
which some cell's log10|v| has changed by dlogv, or its tau or sigma by
dstress, since the last record, or dtmax has passed.

The gather is collective, so the option works under MPI; the write
reads the state only.  theta is computed from psi with the uniform dc
(as in the station output), so it is approximate if `-rsf_dc_file` is
used.  On a restart (`-rsf_restart`) the existing series file is cut
back to its complete records up to the checkpoint time, which removes
what the interrupted run wrote after the checkpoint and a partial last
record (as left by a crash), and the restarted run appends from there
without writing the checkpoint record twice.  A restart test (stop
after 2500 steps, a partial record appended by hand, restart to the
end) gave a series identical to that of the uninterrupted run.  The
catalogs still keep the rows the interrupted run wrote after the
checkpoint (see `-rsf_restart` in `rsf_solve -h`).  Records are written at accepted
solver steps only, so while quiescent the spacing is the larger of the
requested cadence and the step length, which is capped by
`-dt_monitor_yr` (5 yr in `run_model`).  The final accepted step can
lie after `-stop_time_yr` (0.18 yr after it in the single-fault case),
since rsf_solve steps over the stop time.

### Choice of sampling

Since the step size is error controlled, the change per accepted step
is already roughly uniform: in the uniform single-fault case (every
step recorded, 14,602 records over 3000 yr) the largest change of
log10 V over all cells per step had a median of 0.1 to 0.3 in all
phases, and the interseismic tau change per step was 0.1 MPa, set by
the 5 yr step cap.  88 percent of the steps fell in the coseismic
phases (max V above 0.1 m/s, steps of about 0.01 s, about 1400 steps
per event).  Thresholds below the per-step change therefore select
nearly every step, and the two ways of thinning the series behave
similarly.  Interpolating each thinned series linearly in time back to
all steps gave (largest error over cells per step, 99th percentile and
maximum over steps):

| sampling | records | log10 V p99 / max | tau [MPa] p99 / max |
|---|---|---|---|
| every step | 14,602 | 0 / 0 | 0 / 0 |
| dense stride 2 | 7908 | 0.06 / 0.10 | 0.11 / 0.13 |
| dense stride 4 | 4561 | 0.19 / 0.28 | 0.38 / 0.45 |
| dense stride 8 | 2886 | 0.54 / 1.08 | 1.29 / 1.54 |
| mode 1, dlogv 1, dstress 1 | 8645 | 0.21 / 0.34 | 0.22 / 0.93 |
| mode 1, dlogv 1, dstress 2 | 5190 | 0.25 / 0.34 | 0.38 / 1.56 |
| mode 1, dlogv 2, dstress 5 | 2734 | 0.77 / 1.00 | 1.18 / 3.76 |

(the stride rows keep every interseismic step.)  For a similar record
count the stride has comparable or smaller typical errors, and mode 1
bounds the log10 V error more tightly at strong thinning.  Stride
sampling depends on the solver settings (rtol, Runge-Kutta scheme),
mode 1 thresholds do not.  This is one case with one fault; the
choice for training data is open (PLAN.md 5.3).  `run_model` uses
stride 1 unless `-- -rsf_series_stride n` is given; the examples use
stride 4.

The output was tested by the runs below and by reproducing the first
BP5 2 km event.  On 2 MPI ranks the series of the uniform single-fault
case was bit-identical to the 1-rank series, in mode 0 with stride 1
and 4 and in mode 1.  The record selection of the stride and of mode 1
was compared with an emulation in Python on the every-step series and
selected the same steps.  Larger rank counts, multi-fault systems
under MPI, and restarts with mode 1 were not tested.

## Test runs (sandbox, 1 core, rtol 1e-4, 3bs)

Reference parameters: G 32.04 GPa, sigma 50 MPa, a 0.010, b 0.015,
dc 0.02 m (Lb 854 m, 3 cells of 284 m, h* 6.5 km), L 60 km, W 15 km,
v_eq 1e-9 m/s (tau_dot 21.5 kPa/yr).

| run | cells | tmax | wall | events | behaviour |
|---|---|---|---|---|---|
| single uniform | 211 | 3000 | 4.6 s | 9 | periodic, full ruptures, Mw 7.4, T 327 yr |
| 3 faults en echelon, uniform, step 5 km, no overlap | 633 | 3000 | 53 s | 30 | each fault periodic, sequential triggering within years |
| single, 15 % white noise on a, 10 % VS margins | 211 | 3000 | 5 s | 14 | periodic |
| single, 30 % noise on a, 5 km correlation length | 211 | 3000 | 5 s | 15 | periodic |
| single, L 120 km, L/h* 11 | 703 | 3000 | 40 s | 14 | periodic |
| sweep member, L/h* 33, VS margins | 498 | 600 | 9 s | 3 | periodic |
| single, 4 asperities, 20 % barriers | 211 | 3000 | 10 s | 144 | irregular, partial and multi-asperity ruptures, Mw 5.1 to 6.6 |
| 3 faults, 3 asperities each, random lengths, 10 km overlap | 633 | 1000 | 46 s | 75 | irregular on all faults |
| single uniform, sigma evolving (In matrix on) | 211 | 1500 | 2.3 s | 5 | identical catalog to the first row |
| single uniform, `--strikes 90` (right-lateral for either `--rot`) | 211 | 800 | 2.2 s | 3 | mirror image of the first row, negative slip |
| 3 faults, strikes 0, 20, -20, 15 km overlap, 3 km step, sigma evolving | 633 | 500 | 28 s | 6 | fault at +20 (clamping) stops after its first event, fault at -20 (unclamping) fails every 170 yr |
| same, whole system rotated by 37 and by 90 degrees | 633 | 500 | 28 s | 6 | same catalog to 3 digits (37) and 4 digits (90) |
| single, strike 30 | 211 | 300 | | 0 | Coulomb rate negative, no events, sigma + 5.6 MPa as expected |
| sweep member, 3 faults, 6 asperities, 1392 cells | 1392 | 600 | | | reached 139 yr in about 170 s before a timeout; not completed |
| examples/single: one asperity, 10 % VS ends | 211 | 3000 | 5.7 s | 14 | periodic, Mw 7.2 every 218 yr |
| examples/three_right: 3 parallel faults, right steps 4 km, overlap 10 km (restraining) | 633 | 3000 | 104 s | 40 | near-periodic; sigma at the overlapping tips rises steadily, to 83 MPa |
| examples/three_left: as above with left steps (releasing) | 633 | 3000 | 114 s | 45 | near-periodic, 209 yr on all faults; sigma at the tips falls steadily, to 6 MPa |
| examples/three_angled: right steps, outer faults at +15 and -15 degrees (inner pivot), L_b/5 | 1056 | 3000 | 393 s | 57 | clamped fault slows to 476 yr; unclamped fault speeds up and turns irregular (Mw 5.5 to 7.2) after about 2000 yr, with cells at the 5 MPa limiter |
| examples/three_right_ff: as three_right, `--sigma_interact 0` | 633 | 3000 | 69 s | 41 | near-periodic; sigma stays at 50 MPa |
| examples/three_angled_ff: as three_angled, `--sigma_interact 0` | 1056 | 3000 | 234 s | 48 | sigma linear in time per fault; clamped fault slows, unclamped fault speeds up with Mw decreasing gradually, without the irregular phase of three_angled |

Observations, with the caveat that this is a handful of runs with the
aging law:

- Uniform or smoothly heterogeneous VW faults produced characteristic,
  periodic cycles over the whole range tried (L/h* 5 to 33), including
  with VS margins and with 30 percent correlated heterogeneity in a.
  Under uniform far-field loading there is no stressing-rate gradient
  to break the symmetry, and the one-cell-deep discretization removes
  the depth structure that produces partial ruptures in SEAS-type
  models.  Whether larger L/h* or the slip law changes this was not
  tested (`-- -state_law 2` passes the slip law through run_model).
- Segmentation into asperities separated by VS barriers produced
  irregular sequences readily, as did asperity faults in an en echelon
  arrangement.  For an ensemble with a range of complexity the asperity
  count, barrier fraction and asperity length randomness are the
  levers that worked; a/b and dc set event size and recurrence.
- Cost scales roughly as N^2 times the number of accepted steps; the
  step count grows with the number of events.  Up to about 600 cells a
  3000 yr run takes seconds to a minute on one core.  The 1392 cell
  sweep member was not run to completion (139 yr in about 170 s, a
  naive extrapolation gives about an hour per 3000 yr, but the event
  rate varies along the run); see PLAN.md section 4 for the cost model
  and the larger-network estimates.  Note that cells_per_lb below 3 is
  not recommended (L_b/3 is already the resolution floor used here).

- In the stepover examples, with parallel faults and therefore no
  far-field normal rate, the normal stress at the overlapping fault
  tips changes by roughly the same amount at each event of the
  neighbouring fault, and since slip accumulates without bound so does
  this change: over 3000 yr the cells of the middle fault in the
  overlap gained 21 MPa in the restraining case and lost 27 MPa in the
  releasing case, and the extremes reached 83 and 6 MPa.  The model has
  no mechanism (off-fault yielding, fault growth) that would relieve
  this, so stepover systems drift like oblique faults do, and long runs
  reach the limiter.  The sign (compression in restraining, tension in
  releasing steps) is the one expected for the slip sense labelled in
  gen_system.py, which supports, but does not by itself prove, the
  left/right-lateral convention.  The tip values are grid dependent
  (see PLAN.md).
- Creeping tips (`--tip_creep`, velocity strengthening graded from
  `--a_vs` to `--tip_a` toward each end) did not change this drift.
  In the restraining case over 1500 yr, the overlap cells gained 10.5
  to 11.1 MPa, about 0.5 MPa per metre of mean slip on the fault, with
  10 percent VS margins, with 10 and 20 km graded tips (a 0.02 to 0.04),
  with 10 km weakly strengthening tips (a 0.016 to 0.02), and with a
  third of each fault weakly VS.  The VS tips creep and in the long
  run accumulate slip like the rest of the fault, and the slip taper
  toward the ends is set by elasticity, so the tip stress grows with
  the cumulative slip in all variants.  `--tip_creep` is kept as an
  option, but it is not a remedy for the drift; `--sigma_interact 0`
  is.
- Under rising normal stress the resolution degrades as described in
  PLAN.md 2.1: at the restraining tips the cells reached L_b/1.8 at
  83 MPa, while the angled example, run at L_b/5 for sigma_0, stayed
  at L_b/2.3 or better.  In the angled example the unclamped fault
  reached the 5 MPa limiter in some cells, so its late irregular
  phase is affected by the limiter.
- The per-fault event counts above come from the per-group catalogs,
  which share the global event bracket; events on different faults
  that overlap in time are merged (PLAN.md 6).

Checks made:

- Resolution: the 4-asperity case at Lb/3 and Lb/6 gave the same first
  cycle (onset 38.8 vs 38.4 yr, same rupture extents) and the same
  event size distribution afterwards, with trajectories diverging after
  a few cycles as expected for this kind of system.  Conclusions from
  individual trajectories at Lb/3 should therefore be limited to
  statistics.
- Tolerance: rtol 1e-5 tracked the 1e-4 run to about 400 yr (ten
  events) before diverging, with the same event sizes; cost was about
  1.7 times higher.
- Loading: the interseismic mean stress rate on the uniform fault
  (22.3 kPa/yr from the series) matches tau_dot plus a small afterslip
  contribution; the normal stress drift of the strike-30 fault matches
  sigma_dot times time; rotating a three-fault system with mixed
  strikes reproduces its catalog (consistency of the resolved loading
  with interact's patch frames).
- Normal stress evolution: the extra matvec did not change the single
  fault result and roughly doubles the per-step cost for multi-fault
  systems (28 s for 500 yr of the three-strike case).
- Series output: the BP5 2 km first event and record count were
  consistent with the monitor and catalog; the reader round-trips the
  header and record layout.

## Usage

The examples above are produced by `examples/make_examples [tmax_yr]
[stride] [rundir]`, which runs them and writes the figures
single.png, three_right.png, three_left.png, three_right_left.png,
three_angled.png, three_right_interact_ff.png and
three_angled_interact_ff.png (one column per fault: log10 V, tau with sigma,
log10 theta, slip with event Mw, and log10 V along strike against
time).


    # single fault, 3000 yr, default parameters
    ./run_model s1 3000 1
    # four asperities, slip law
    ./run_model asp4 3000 1 --n_asp 4 --barrier_frac 0.2 --vs_margin 0.05 -- -state_law 2
    # three en echelon faults, left steps
    ./run_model e3 3000 1 --nfault 3 --step_x 4e3 --overlap 10e3 --n_asp 3 --asp_rand 0.3
    # three faults with different strikes, S1 axis north-south, no far-field normal rate
    ./run_model e3s 2000 1 --nfault 3 --strikes 300,320,335 --sr_az 0 --ff_normal 0
    python3 plot_faults.py e3 e3.png --spacetime
    python3 plot_run.py asp4 asp4.png --tmin 1000 --tmax 2000
    # ensemble of 50 members, 2000 yr each, four at a time
    ./run_sweep sweep 50 2000 100 4
    # training export (series interpolated to 0.5 yr, plus event windows)
    python3 -c "from rsf_series import Run, export_training; export_training(Run('asp4'), 'asp4/train.npz')"

`run_model` expects `rsf_solve` on the PATH or in `config.sh`
(variable `rsf_solve`).  The test runs used two builds: PETSc 3.23.4
built from source on Ubuntu 22.04, and the Ubuntu 24.04 PETSc 3.19.6
package with HACApK and hmmvp.  The second build reproduced the
reference catalog of the first to the digits given in the table, with
a slightly different count of RHS evaluations, so trajectories from
different builds should not be expected to agree exactly.  With the
Debian/Ubuntu PETSc package, the distribution flags in its pkg-config
file (including -flto) made the link fail; overriding them worked:

    make bin/rsf_solve PETSC_DIR=/usr/lib/petscdir/petsc3.19/x86_64-linux-gnu-real \
        PETSC_ARCH= PETSC_HTOOL_USED= SCARGS="-Wno-discarded-qualifiers" \
        SFARGS="-ffixed-line-length-132" CXXFLAGS=

## Open points

- The training export stores events as ragged object arrays; a fixed
  length, resampled event representation may suit a transformer better.
- Per-cell dc and sigma heterogeneity, dipping faults, and normal
  stress evolution are supported by rsf_solve but not exposed in
  gen_system.py.
- The default series cadence produces large files for runs with many
  creeping barriers; `-rsf_series_vdense` and `-rsf_series_dt_yr` are
  the knobs.
- No check against a 3-D (multi-row) model of the same fault has been
  made, so the effect of the one-cell-deep approximation on event
  statistics is not quantified.
