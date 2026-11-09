# Implementation plan: 2.5-D fault system cycle ensembles for learned embeddings

Revision 7, 2026-11-08. Status markers: [done], [todo], [open].

Changes from revision 6: regime mapping and learning program (section
10), with what Sun & Becker (2026, GJI 247) and Sun, Gilpin & Becker
(submitted, JGR) imply for it; two-state-variable rate-and-state law in
rsf_solve, checked against the nondimensional spring slider; new code
requirements (10.4) and the single-fault pilot (10.5).

Changes from revision 5 to 6: options for far-field-only normal stress
(rsf_solve -sigma_dot_interact 0, gen_system --sigma_interact 0) and
creeping tips (--tip_creep); tests show that creeping tips do not
change the stepover drift and the far-field-only option removes it
(2.3); compile notes in COMPILE.md.

Changes from revision 4 to 5: output sampling options implemented in
rsf_solve (stride, change based, final record) with defaults that
leave the output unchanged, and compared on the single fault (5.3);
stepover and pivot options in gen_system.py; example runs (single
fault, three faults with right and left steps, three faults with
angled outer faults); normal stress drift at stepover tips and
limiter engagement observed in them (2.1, 2.3, 8).

Changes from revision 3 to 4: code moved to interact/strike_slip_networks;
theta reconstruction stated to need b and dc as well as a, with the
visibility of a, b and dc chosen per experiment and both variants run
(5.1); depth check limited to uniform rows down to W (8); output
sampling restated as an open design question with a proposal (5.3,
6); H-matrix backends HACApK and hmmvp now build in the sandbox (4.3).

Changes from revision 2 to 3: channel sets restructured around observability
(theta is redundant given V, tau, sigma, a), window length added as an
experiment axis; cell size set from the maximum expected normal stress
and limiter engagement recorded; cost model corrected against the logs
and unmeasured numbers marked as such; per-cell output kept for all
tiers; order of work changed so that the depth check and a tier 0 to 2
learning baseline come before any tier 5 machinery. Matrix truncation
is deliberately not pursued: full dense first, then H-matrix
approximations.

## 1. Goal

Generate a large, physically consistent ensemble of quasi-dynamic
rate-and-state earthquake cycle simulations of strike-slip fault
systems, from one fault to networks of up to hundreds of faults,
spanning periodic to irregular behaviour, and learn low-dimensional
embeddings of the observable fields in order to study the system
dynamics in embedding space. Models are trained on the fields, not on
the interaction matrices; whatever coupling structure the embedding
recovers has to come from the data.

## 2. Model family [done, README.md]

- Vertical strike-slip faults, one Okada cell deep (locking depth W),
  strips along strike, half space, mode II, aging law, quasi-dynamic.
- Far-field horizontal stressing-rate tensor resolved per fault into
  shear and normal rates; no backslip. Normal stress evolves from the
  far-field rate and from slip on the other faults. Oblique faults
  therefore have secular recurrence and size drift (intended, Becker
  and Schmeling 1998).
- Fixed across the ensemble: G, c_s, sigma_0 = 50 MPa, b = 0.015, f0,
  v0, aging law.
- Varied: number of faults, L, W, step, overlap, strike relative to S1,
  a/b, dc, asperity count, barrier fraction and length randomness, VS
  margins, tau_ref, initial-state seed.

### 2.1 Resolution under evolving normal stress [todo]

L_b = G dc/(b sigma) and h* scale with 1/sigma. Cells sized from
sigma_0 lose resolution on clamping faults: at 20 degrees from optimal
(13.8 kPa/yr) sigma reaches about 91 MPa after 3000 yr and L_b/3
becomes about L_b/1.65. Rule: size cells from
sigma_max_expected = min(sigma_0 + max(sigma_dot, 0) t_max, sigma_max)
per fault, with the global cell size from the most demanding fault,
and record sigma_max_expected and the resulting cells per L_b in
params.json.

This rule covers the far-field rate only. At stepover tips the normal
stress also drifts from the accumulated slip on the neighbouring fault
(2.3), which reached 83 MPa in a restraining example with parallel
faults (cells at L_b/1.8 there). An estimate for the tips, e.g. the
per-event tip change times the expected number of events, or a
measured value from a pilot, is needed in addition.

### 2.2 Limiter and stability transitions [todo]

The normal stress limiter (5 to 200 MPa) is not physical. Choose run
length per member so that sigma_0 + sigma_dot t_max stays inside the
limits with margin, and record in the series whether the limiter is
engaged anywhere (count of clipped cells per record, or a flag in the
per-group output). Unclamping faults approach low sigma, where h* grows
toward the fault length (fivefold at 10 MPa) and the fault moves
toward stable sliding; this transition is a legitimate regime and is
labelled (h*/L_VW over time per fault) rather than avoided.

The limiter was engaged in the angled three-fault example (the
unclamped fault reached 5 MPa in some cells after about 2000 yr), which
makes the limiter flag in the series (6) more urgent than planned.

### 2.3 Normal stress at stepover tips [options in place, choice open]

With parallel faults and no far-field normal rate, the normal stress
at overlapping tips changes at each event of the neighbouring fault,
and since slip accumulates without bound, so does the change: over
3000 yr the overlap cells of the middle fault gained 21 MPa in a
restraining case and lost 27 MPa in a releasing one (extremes 83 and
6 MPa). Nothing in the model (off-fault yielding, fault growth)
relieves this, so stepover systems are non-stationary like oblique
faults, and the drift is concentrated at grid-dependent tips. Options
to decide between: accept it as part of the model family and limit
run length per member as in 2.2; or limit the tip effect, e.g. by
velocity-strengthening tips that creep and spread the interaction, or
by keeping the normal stress interaction only from the far-field rate.
The sign of the drift (compression in restraining steps) is the
expected one.

Both options were implemented and tested:

- Creeping tips (`--tip_creep`, graded velocity strengthening toward
  the ends) did not change the drift: 10.5 to 11.1 MPa over 1500 yr in
  the overlap of the restraining case, about 0.5 MPa per metre of mean
  slip, for five tip variants. The VS tips accumulate slip like the
  rest of the fault in the long run, and the slip taper toward the
  ends is elastic, so the tip stress grows with cumulative slip
  regardless.
- Far-field-only normal stress (`--sigma_interact 0`, rsf_solve
  `-sigma_dot_interact 0`) removes the drift by construction; sigma
  then follows sigma_0 + sigma_dot t per fault.

What remains to decide is which normal stress model the ensemble
uses: full interaction (with the drift, run length limited per member
as in 2.2), far-field only (stationary stepovers, but no normal stress
coupling between faults, so clamping and unclamping by neighbouring
events is lost), or both as an ensemble axis.

## 3. Complexity ladder

| tier | configuration | expected behaviour (tests so far) | cells |
|---|---|---|---|
| 0 | single fault, uniform VW, optimal orientation | periodic characteristic events | 100 to 200 |
| 1 | single fault, VS margins, oblique | periodic with recurrence and size drift | 100 to 300 |
| 2 | single fault, 2 to 6 asperities | irregular, partial and multi-asperity ruptures | 200 to 300 |
| 3 | 3 en echelon faults, simple a-b, parallel | each near-periodic, triggering, phase locking | 600 |
| 4 | 3 faults, mixed strikes, asperities, drift | irregular on all faults plus interaction | 600 to 900 |
| 5 | networks of 10 to 100 faults, mixed strikes and lengths | collective behaviour, unknown | 1000 to 10000 |

Tiers 0 and 1 are controls: an embedding that does not recover cycle
phase and sigma drift there is not worth training on tiers 4 and 5.
Tier 1 runs are non-autonomous unless sigma is among the inputs (set
D below), which is part of what the channel comparison tests.

## 4. Computational cost

### 4.1 Measured (sandbox, one core, dense, 3bs, rtol 1e-4)

| cells | faults | sigma | model time | events | RHS evals | RHS per event | wall |
|---|---|---|---|---|---|---|---|
| 211 | 1 uniform | fixed | 3000 yr | 9 | 65,885 | 7300 | 4.6 s |
| 211 | 1, 15 % noise | fixed | 3000 yr | 14 | 99,181 | 7100 | 5 s |
| 703 | 1, L 120 km | fixed | 3000 yr | 14 | 194,967 | 13,900 | 40 s |
| 211 | 1, 4 asperities | fixed | 3000 yr | 144 | 227,385 | 1600 | 10 s |
| 633 | 3 uniform | fixed | 3000 yr | 30 | 207,444 | 6900 | 53 s |
| 633 | 3, 3 asperities | evolving | 1200 yr | 95 | | | 93 s |
| 633 | 3, simple | evolving | 1200 yr | 17 | | | 71 s |

Not measured, and to be treated as such: a 1392-cell three-fault run
reached 139 yr in about 170 s before a timeout (the "11 min per 600 yr"
in revision 2 was an extrapolation); a 2110-cell ten-fault run spent
4.8 s and 3.8 s assembling Is and In and then filled the sandbox disk
with per-cell output before completing, so no wall time exists.

Wall time per RHS evaluation: 70 microseconds at 211 cells, 255 at 633
cells (one matrix), i.e. about 4 GFlop/s effective dense matvec,
memory-bandwidth bound at larger N. RHS evaluations per event range
from about 1600 (small asperity events) to 7000 (full-fault Mw 7.4
events) and 14,000 (120 km fault), i.e. they grow with event size and
rupture duration in steps.

### 4.2 Cost model and tier 5 estimate

wall = n_RHS x t_matvec x (1 or 2 matrices).
n_RHS = (events per model time) x (RHS per event) + quiet-time steps;
with many faults there is little quiet time, since some fault is
nearly always nucleating or relaxing and the time step is global.
Events per time scale with the number of faults and strongly with
segmentation (a 4-asperity fault produced 16 times as many events as a
uniform fault of the same length). t_matvec scales as N^2 dense.

100 faults, 50 to 80 cells each (5000 to 8000 cells), 1000 yr:
- events: about 50 per 1000 yr per segmented fault (144 per 3000 yr
  measured), 3 to 5 for simple faults; 300 to 5000 per member.
- RHS evaluations: simple faults 300 to 500 events x 7000, segmented
  faults 5000 events x 1600; 2 to 10 million per member.
- dense, two matrices (Is, In), 5000 to 8000 cells: 0.4 to 1 GB of
  matrix read per RHS, about 0.05 to 0.13 s on one core; 30 to 300
  core hours per member.
- with 16 MPI ranks dense (speed-up limited by memory bandwidth,
  assume about 10): 3 to 30 h per member.
These are extrapolations by one to two orders of magnitude in N from
the measured runs and should be replaced by the 10 and 30 fault pilot
numbers (section 8). The upper end is impractical for an ensemble of
100 members with the dense operator; this is where H-matrices are
needed.

### 4.3 Measures, in the order they will be tried

1. Fewer cells per fault through larger dc and shorter faults for
   tier 5 (dc 0.03 to 0.05 m, L 20 to 40 km, 30 to 80 cells per fault,
   events Mw 6 to 7), within the resolution rule of 2.1.
2. Full dense operator with MPI (16 to 24 ranks) as the reference for
   the tier 5 pilot.
3. H-matrix operator (`-use_hmatrix 3`, HACApK, or 1, HTOOL), already
   production in interact for BP5; validated against the dense
   reference on the pilot networks for event statistics before use.
   Neither the speed-up nor the accuracy for strip networks has been
   measured. HACApK and hmmvp build in the current sandbox and
   reproduced the single-fault catalog (onsets within 0.005 yr over
   1000 yr); multi-fault systems and the normal-stress matrix are
   untested.
4. Shorter model time for tier 5 (1000 yr); long-period statistics
   come from tiers 0 to 4.
5. Output volume control (section 5.3).

Truncating far-field interactions is not pursued.

## 5. Data products

### 5.1 Channel sets and the observability question

The per-cell dynamical state is (psi, tau, sigma); V follows from it
through the friction law V = 2 v0 exp(-psi/a) sinh(tau/(sigma a)), and
theta follows from psi. Given a, sigma and V, tau, theta is therefore
determined exactly; and tau is itself determined by the V history
(tau = tau_0 + tau_dot t - K slip). Channels are thus not independent
pieces of information, and the experiment is one of observability
from finite windows:

| set | channels | question |
|---|---|---|
| A | log10 V | how long a slip-rate history recovers the embedding of B or D (delay embedding) |
| B | log10 V, tau | complete local state at fixed sigma, given a |
| D | log10 V, tau, sigma | complete local state with evolving normal stress |

The second experiment axis is the input window length (number of past
samples per training example), from single snapshots to several
cycles. The expected outcome is that A converges to B and D as the
window grows, and the interesting quantities are how fast, and whether
the latent dimension needed is the same.

The state variable is handled as an explicit two-way comparison, and
every experiment states which variant it uses:

- theta reconstructed (default): theta is not an input; it is a target
  the model must recover from the other channels. With V, tau, sigma,
  a, b and dc available (f0 and v0 are fixed across the ensemble), the
  reconstruction should be exact to float precision (the conversion
  uses the run's uniform dc), so failure indicates a pipeline problem;
  with any of a, b, dc hidden, or from set A alone, how well and from
  how long a window theta is recovered is a result.
- theta known: theta is an input channel (sets C = V, tau, theta and
  D' = V, tau, theta, sigma). Since theta is redundant given V, tau,
  sigma and a, any difference from the reconstructed variant in latent
  dimension or forecast skill measures what the model gains from being
  handed the state rather than inferring it, not new physics.

Whether a, b, dc and sigma are visible to the model (as static or
dynamic inputs) is likewise stated per experiment, and both variants
(frictional parameters given, frictional parameters hidden) are run.
dc varies across the ensemble, so with dc hidden theta carries
information about dc, and with a hidden about a.

### 5.2 Representations [todo]

- Per-cell tensors for all tiers, including all tier 5 members:
  cells ordered along strike within faults, faults in system order,
  with static attributes (fault id, along-strike coordinate, x, y,
  strike, a - b, W, cell size, resolved tau_dot and sigma_dot) as side
  information.
- Per-fault aggregates (mean, std, min, max per channel, fraction of
  cells above 1e-6 m/s) computed from the per-cell data in the export
  step, used as a comparison representation only: they are themselves
  a hand-made low-dimensional embedding and cannot be the only form,
  or the embedding result would be predetermined by the aggregation.
- Time sampling, both exported, with the time increment as an input:
  (a) uniform grid, 0.5 to 1 yr; (b) event windows at the solver's
  step density, resampled onto a fixed number of samples in log time
  from onset, one tensor shape for all events.
- Labels for evaluation only: catalog (onset, Mw, extent, fault), cycle
  phase per fault, sigma and its drift rate, limiter flag, h*/L_VW,
  tier, parameter vector.
- Normalization statistics per channel over the ensemble.

Graph structure, if a GNN is used: along-strike chain plus tip-to-tip
links across steps, purely geometric. A transformer with along-strike
positional encoding is the simpler baseline.

### 5.3 Sampling and volume [open]

The dynamics span time scales of centuries (interseismic loading),
years down to seconds (nucleation), tens of seconds (coseismic) and
hours to years (afterslip), so no fixed output cadence fits. As written now, the series is recorded at every accepted step
while max V > 1e-6 m/s and otherwise at the larger of 0.5 yr and the
step length (capped at 5 yr by -dt_monitor_yr). In the uniform
single-fault case this gave about 1500 records per event, 95 percent of
all records, and interseismic spacings of 5 yr.

Status: rsf_solve now has both a stride over accepted steps in the
dense phase (`-rsf_series_stride`) and change-based recording
(`-rsf_series_mode 1`), as new options that leave the default output
unchanged. On the uniform single fault the error-controlled steps
already carry a roughly uniform change (per-step max |d log10 V| of
0.1 to 0.3 in all phases), so the two behave similarly at equal
record counts (README.md, Choice of sampling). The coseismic phases
hold 88 percent of the steps (about 1400 per event at 0.01 s), which
is where any thinning acts. Which sampling, and how much coseismic
detail, the learning needs is to be decided from the tier 0 to 2
baseline; the examples use stride 4.

Remaining proposals, not yet agreed:
- Replay rather than storage: keep checkpoints at fixed model-time
  intervals (numbered files, about 3N doubles each) together with the
  inputs, so that any window can be regenerated at any output density.
  Requires a model-time checkpoint option and a test of whether a
  restart reproduces the continuous run.
- Per member by default: catalog, per-fault aggregates at every
  accepted step, the change-based per-cell series, checkpoints.
  Full-density per-cell output for a validation subset only.

For training, the stored sequence can be presented as an irregular
sequence with continuous time encoding, as a slow uniform stream plus
fixed-shape event windows in log time (5.2), or on the change-based
clock with dt as a channel. These are compared in the tier 0 to 2
baseline. Note that with change-based sampling the sample spacing
itself carries information about the dynamics (small dt during fast
change), as do event windows aligned on catalog onsets; for the
observability experiments of 5.1 the time representation is therefore
stated per experiment like the channel set.

Volume estimates as written now:

Per-cell records cost 4 bytes x 5 channels x N. At 200 cells: 50 to 150
MB per 3000 yr member as written now, 20 to 60 MB compressed. At 10^4
cells a record is 200 kB; with a yearly interseismic cadence and a
strided event phase (every 10th to 50th accepted step) a 1000 yr tier
5 member is of order 1 to 10 GB uncompressed. Feasible for 100
members on a cluster file system, not on a laptop; keep full-cadence
series only for a validation subset.

## 6. Code changes

rsf_solve (interact), output layer:

- [done] all-cell series with sigma (`-rsf_series_file`).
- [done] output sampling options (5.3): `-rsf_series_stride`,
  `-rsf_series_mode 1` with `-rsf_series_dlogv`, `-rsf_series_dstress`,
  `-rsf_series_dtmax_yr`; defaults unchanged. The choice for the
  ensemble is open.
- [done] `-rsf_series_final`: record the final accepted step.
- [done] series file cut back to the checkpoint on `-rsf_restart`
  (previously a crash could leave a partial record and misalign all
  appended records). Catalog rows written after the checkpoint by the
  interrupted run are still kept.
- [done] `--gap_dir bisect` for point-symmetric chains (gen_system.py).
- [done] `-sigma_dot_interact 0`: normal stress from the loading rates
  only. The In matrix is still assembled (the backslip loading uses
  it); skipping it when -vpl is 0 would halve the matrix memory for
  tier 5 and is a possible later change.
- [open] model-time checkpoints kept as numbered files, for replay
  (5.3).
- [todo, raised in priority by the examples] limiter flag: number of
  clipped cells per series record (one extra int32, or a sixth float
  block with the clip state per cell).
- [todo, before tier 5] per-group event brackets: the catalog split by
  group reuses the global onset/arrest bracket (max |v| over all
  cells), so simultaneous events on different faults merge. Needed for
  tier 5 labels.
- [later, optional] per-group aggregated series in the solver, only if
  per-cell output proves too large in practice.
- [optional] `-rsf_dump_matrices` as a diagnostic.

gen_system.py:

- [todo] cell size from sigma_max_expected (2.1); run length check
  against the limiter (2.2); write both to params.json.
- [todo] per-fault lists for a, dc, n_asp, W (lengths: done for
  `--layout chain`, `--lengths`).
- [todo] `--strike_rel`: strike relative to the S1 axis.
- [done] `--stepover right|left` and `--pivot center|inner`; step types
  (restraining, releasing) in params.json.
- [done] `--tip_creep`, `--tip_a` (graded VS tips) and
  `--sigma_interact` (far-field-only normal stress).
- [done] `--layout chain`, `--gap`, `--lengths`: segmented traces with
  faults placed tip to tip (used for the Becker and Schmeling 1998
  analogues in becker_and_schmeling98/).
- [todo] `--max_cells` guard.
- [todo] tier and drawn parameter vector into params.json.
- [todo, before tier 5] network layout mode (`--network <file>` or
  random layout with fault lengths log-normal, strikes about the
  optimal orientation, minimum spacing from the existing closeness
  check).
- [todo] depth option `--nrow` (cells in depth, default 1) for the
  2.5-D versus multi-row check (section 8); rows of equal height
  between the surface and W, same friction per column.

run_sweep:

- [todo] tiered draws, per-tier member counts and model time, N cap,
  backend options per tier (dense for tiers 0 to 4; dense MPI then
  H-matrix for tier 5).
- [todo] validation subset per tier at L_b/6 and rtol 1e-5.

rsf_series.py / export:

- [todo] HDF5 per member: all channels stored once, channel set and
  window length selected at load; both time representations; per-fault
  aggregates; labels; static attributes; ensemble index.
- [todo] event-window resampling in log time.
- [todo] ensemble statistics: recurrence mean and CV per fault, Mw
  distribution, partial-rupture fraction, sigma drift, limiter
  engagement, event rate against number of faults.

Plotting: `plot_faults.py` (per-fault columns, several runs side by
side) for tiers 0 to 4; a fault-by-time raster (`--raster`) for tier 5.

## 7. Learning side (outline)

- Baselines per channel set and window length: PCA on per-cell and
  per-fault snapshots and on delay-embedded windows; intrinsic
  dimension estimates.
- Autoencoder plus a latent dynamics model (neural ODE, Koopman or
  SINDy on latent coordinates, or a latent transformer with dt as
  input). Diagnostics: smoothness of latent trajectories across
  events, latent dimension needed, whether the latent flow is
  autonomous (expected only with sigma visible for oblique faults).
- Transformer over time with per-cell tokens; GNN on the chain graph
  plus a temporal model; objectives: next state, masked
  reconstruction, time to next event per fault.
- Evaluation: cycle phase and sigma drift on tiers 0 and 1 without
  labels; theta reconstruction as a pipeline check; clustering by
  tier; forecast skill on held-out seeds and geometries; per-cell
  versus per-fault-aggregate embeddings.

## 8. Checks

- [done] resolution L_b/3 vs L_b/6 and rtol 1e-4 vs 1e-5 (asperity
  case): same statistics, diverging trajectories.
- [done] rotation invariance of the resolved loading.
- [todo, first] one-cell-deep approximation: tier 0 and tier 2
  configurations against 3 to 4 rows of equal height down to W, with
  the same friction per column, compared in recurrence and Mw. If the
  statistics differ qualitatively, the model family (section 2) is
  revisited before anything else. This tests the depth discretization
  of the locked strip; the effect of a velocity-strengthening region
  below W is not part of the check.
- [todo] left/right-lateral sense: single-cell slip, sign of the
  off-fault stress change compared with a hand calculation from the
  Okada solution. Partial support so far: the normal stress drift in
  the stepover examples has the sign expected for the labelled slip
  sense (compression in restraining steps).
- [todo] stepover tips: one tier 3 stepover case at L_b/6, since the
  normal stress at fault tips is grid dependent and drives the wide
  sigma spread in the three-fault runs; compare the tip drift rate of
  2.3 between L_b/3 and L_b/6.
- [todo] resolution rule of 2.1 on a 20 degree clamping fault, L_b/3 at
  sigma_0 against L_b/3 at sigma_max_expected.
- [todo] seed decorrelation time per tier.
- [todo, before tier 5] 10 and 30 fault pilots, dense, for scaling and
  per-group event detection; then H-matrix against dense on the same
  pilots.

## 9. Order of work

1. [done] Sandbox reset and rebuild (HANDOFF.md). Decide the output
   sampling (5.3), since the depth runs multiply the cell count. Then
   the checks that can change the model family: depth (`--nrow`),
   handedness, stepover tips, resolution under evolving sigma.
2. rsf_solve: final record, output sampling, limiter flag. gen_system.py:
   sigma-based cell size, limiter check, per-fault lists, strike_rel,
   max_cells, metadata.
3. Tiers 0 to 2 pilot (20 members each), HDF5 export, ensemble
   statistics.
4. First learning baseline on tiers 0 to 2: PCA and autoencoder over
   channel sets A, B, D and window lengths, each with theta
   reconstructed and theta known.
   Revise the data specification (cadence, window lengths, channels)
   from what this shows.
5. Tiers 3 and 4 pilot and full tiers 0 to 4 ensemble.
6. Tier 5: per-group event brackets, network generator, 10 and 30
   fault dense pilots, H-matrix validation against dense, then the
   tier 5 ensemble on a cluster.
7. Latent dynamics models, transformer, GNN.

## 10. Regime mapping, complexity measures and the learning program

### 10.1 Question

Can this model family span periodic, characteristic, period-doubled and
chaotic sequences, and can learned low-dimensional embeddings of the
field dynamics recover the structure of those regimes (dimension,
bifurcations, synchronization) and approach their intrinsic
predictability? The regime structure has to be established with
standard nonlinear-dynamics measures first, as ground truth for the
learned embeddings.

### 10.2 What the runs so far show, and what is missing

- Uniform velocity-weakening strips under uniform far-field loading were
  periodic and characteristic for L/h* from 5 to 33. Irregularity came
  from segmentation, fault interaction, and the normal stress drift of
  oblique faults and stepovers. No clean period-two sequence has been
  identified, and no run was long enough to characterize chaos.
- The 2.5-D strips have no depth structure, and the VS margins are not
  driven at the plate rate, so the edge loading that drives partial
  ruptures and period doubling in finite-fault models (e.g. Kato 2014;
  Barbot 2019; Cattania 2019) is absent.

### 10.3 Lessons from Sun & Becker (2026) and Sun, Gilpin & Becker (submitted)

Both use HBI (Ozawa et al. 2023) on a single planar fault, 50 x 20 km,
triangular elements of 80 m, aging law, backslip loading, constant
normal stress.

- Controls. W' = W/L_inf (seismogenic width over nucleation length) is
  the primary control: semiperiodic to irregular sequences for W'
  above about 20, strongly irregular above about 45. Fault dip (the
  free-surface asymmetry), the gradient of the down-dip a - b
  transition and lateral velocity-strengthening bounds also change the
  regime, and their effect depends on the regime of the reference model
  (the same change can stabilize or destabilize). Dip acts as a
  stiffness-like control parameter with period-2, 4 and higher periods
  and irregular windows, but not a clean Feigenbaum cascade; periodic
  and irregular windows lie close together in dip.
- Numerics. Under-resolution changes the regime: L_b/dx of 2.2 gave a
  spurious cluster-10 sequence, 1.8 a more irregular one, while 6.9 and
  5.5 agreed statistically. Mesh topology changed recurrence by 1 to 3
  percent. Irregular runs diverge from round-off (different thread
  counts) while their statistics hold. Our runs use L_b/3 by default;
  for regime mapping L_b/6 and a resolution check per regime are
  needed.
- Measures. CV of recurrence with classes (periodic < 0.1,
  semiperiodic < 0.3, irregular < 0.6, highly irregular), cluster-N
  counting of repeated intervals (5 percent of the median), recurrence
  return maps (Lorenz maps), distinct stress extrema (prominence 5
  percent of the stress range, small swings removed, tolerance 3e-5),
  order-6 permutation entropy of the peak-stress sequence, correlation
  dimension D2 (Grassberger & Procaccia) of (max log V, mean tau/sigma,
  mean state) on randomly drawn times with monotone cubic interpolation
  (removes the over-weighting of coseismic steps), Theiler window of two
  mean recurrence intervals, checks with TwoNN and DANCo, and PCA of
  the double-centred slip-rate field with the participation ratio.
  Lyapunov exponents for the spring slider (Benettin et al.), basins of
  attraction in the unstable eigenplane.
- Findings to build on. D2 stays near 1 for periodic orbits of any
  period and rises to about 1.7 to 1.9 only in irregular regimes; the
  extrema count and PE resolve the periods. Spatial complexity (PCA
  modes, PR) and temporal complexity are only partly coupled. The
  two-state spring slider shows final-state sensitivity (coexisting
  period-4 and period-6 or chaotic attractors in narrow stiffness
  ranges, interwoven basins).
- Implications here. The one-cell-deep strips can only represent the
  W' control along strike (L/h*) and cannot carry down-dip transitions
  or dip effects; a dip option for the strips and the multi-row option
  are therefore needed to connect to these results (10.4). The measures
  above become the regime labels of the ensemble, and D2 alone is not a
  sufficient label.

### 10.4 Code requirements

rsf_solve:

- [done] two state variables: `-rsf_two_state 1 -rsf_b2 -rsf_dc2`;
  friction psi = f0 + b ln(v0 theta1/dc) + b2 ln(v0 theta2/dc2), both
  variables follow `-state_law`; state vector of five entries per cell
  (psi total, tau, sigma, slip, psi2 = f0 + phi2); series output gains
  log10 theta2 as a sixth field. Not with `-imex` or `-fd_mode` (stops
  with an error). Defaults unchanged; outputs of runs without the
  option byte-identical to before; 2 MPI ranks identical to 1. Checked
  against the nondimensional spring slider of Gu et al. (1984) and
  Becker (2000) on the single-patch test (test_slider/two_state): same
  number of stress extrema and extrema values to about 1e-3 for
  kappa' = 0.95, 0.9, 0.86 and 0.8552.
- [todo] per-group event brackets (multi-fault event labels).
- [todo] catalog rows of an interrupted run removed on restart (the
  series is already cut back).
- [todo] limiter flag in the series.
- [later] per-cell b2 and dc2; tangent-linear Lyapunov estimates are
  not planned (twin runs instead).

gen_system.py:

- [done] `--b2`, `--dc2`, `--state_law` (steady-state initial stress
  with b + b2).
- [done] `--ic_kick`, `--ic_kick_cell`: perturb one cell's initial
  stress, for twin runs with otherwise identical inputs.
- [done] `--load farfield|backslip|both`: uniform plate-rate backslip
  (rsf_solve -vpl) instead of or with the far-field rate. Backslip on
  the VS cells alone does not load the locked part (the creeping cells
  then carry no slip deficit), so the backslip option applies to all
  cells.
- [todo] `--dip` for the strips (dipping one-cell-deep faults, strike
  slip), as the stiffness-like control of Sun et al.
- [todo] `--nrow` with a depth profile of a - b (VW zone, transition
  gradient, VS below), for the depth check and for comparison with
  down-dip transition effects.
- [todo] frictional bands (VS band of given width and position across a
  fault).

Analysis (new `complexity.py`):

- [done] complexity.py: recurrence intervals, CV and classes, cluster-N;
  distinct stress extrema (prominence, small-swing removal, tolerance);
  order-6 permutation entropy; D2 on randomly drawn, interpolated times
  with a Theiler window (Lorenz 2.04 against the published 2.05, closed
  orbit 1.00, `complexity.py --test`); TwoNN; PCA of the double-centred
  field with participation ratio. Return maps in
  single_fault_pilot/plot_p2.py.
- [todo] projection on reference modes; per-fault measures for
  multi-fault runs.
- [todo] twin-run divergence: separation in (log V, tau) and in event
  times against time, largest Lyapunov exponent estimate and
  predictability horizon in cycles.

Run management:

- [todo] scan driver for one-parameter scans with long runs, resumable
  from checkpoints (sandbox restarts between sessions kill running
  jobs), recording wall time per run.

### 10.5 Single-fault pilot

Goal: find out whether the 2.5-D strip family reaches period doubling
and chaos, and which control parameter does it, before any ensemble.

Status (single_fault_pilot/README.md): P0 done; P1 first pass done
(cascade 2, 4, 8, then irregular, with a periodic window, and large
stick-slip events reached from some starts near the accumulation point);
P2 first pass done at L/h* 5 to 30 (far-field loading: periodic,
period 2 at L/h* 20, periodic again at 30; backslip loading: partial
ruptures, irregular at 20, three-event pattern at 30), runs too short
for 100 cycles.

Single-fault sweep (single_fault_sweep/README.md): one against two
state variables, L/h* 10 to 30, strike -20 to 20 degrees from optimal,
3000 yr, far-field loading: orientation dominates (secular recurrence
trends), detrended sequences nearly periodic, period 2 at L/h* 20 with
one state variable (resolved at L_b/8), no partial ruptures under
far-field loading; two-state runs not yet converged in recurrence. The
one-state slip law gives nearly the same large-event sequences as the
aging law; both slip-law variants add small post-seismic bursts at the
strip tips that are counted as events at the 1 mm/s threshold and
depend on resolution.

Regime sweep against W/h* with backslip loading
(single_fault_regimes/README.md): aging law periodic up to W/h* 8,
partial ruptures from W/h* 16 (about 3 events per cycle, of the order
of Cattania 2019); slip law, one and two state variables, periodic full
ruptures up to W/h* 32 and 16. W/h* 64 and two-state runs above 16 to
do.

- P0 [done]: two-state law against the nondimensional spring slider.
- P1: two-state spring slider in rsf_solve over kappa' from 0.95 to
  0.85 with long runs (5e4 dc/v0 settling, 5e4 analysis): bifurcation
  diagram of stress extrema and comparison with Becker (2000) and Sun
  et al. Tests the analysis tools on a known cascade. Note that at
  kappa' = 0.85 runs started at 1.001 to 1.05 v0 left the small-
  amplitude attractor and settled into large stick-slip events (with
  radiation damping), which needs a look.
- P2: single strip, aging law, uniform VW with VS ends, L/h* from about
  5 to 60 (dc or sigma varied at fixed L), L_b/6, at least 100 cycles
  after spin-up, loading (a) uniform far-field rate, (b) loaded margins.
  Measures of 10.3 per run; twin runs for a few cases.
- P3: as P2 with the two-state law (b2, dc2 varied), to see whether the
  second state variable produces period doubling on an extended fault.
- P4: dip scan of a single strip (90 to about 10 degrees) at fixed L/h*,
  the analogue of Sun et al. (needs `--dip`).
- Resolution: one L_b/3 against L_b/6 comparison per regime found.
- Cost: single strips of 200 to 700 cells take seconds to minutes per
  2000 yr at L_b/3 (section 4); L_b/6 roughly quadruples that.
