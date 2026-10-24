# Viscoelastic earthquake-cycle demos: literature systematics with rsf_solve

Two self-contained bash + python demos showing how to drive
rsf_solve's viscoelastic machinery for the two 2-D configurations,
each reproducing published cycle systematics.  Both are DEMOS:
survey resolution, minutes-to-tens-of-minutes runtimes; the full
replications with convergence studies live in
test_twod/anti_plane/noda_mn_tests.md (antiplane) and the ledger
README_ve_antiplane.md.

## 1. Antiplane, uniform Maxwell medium: Miyake & Noda (EPS 2019)

    ./run_mn19_demo                # ~15 min, writes mn19_demo/mn19_demo.png

A rate-weakening patch in a uniform Maxwell medium (-ve_mode 1, the
single EXACT relaxation pole) under constant far-field loading.
Sweeps t_c/t_load and plots the recurrence against the divergence
law t_rec = alpha/(x - x_cr) + beta, with Miyake & Noda's
fully-dynamic periodic-fault reference curve (alpha 2.63, beta 7.28,
x_cr 1.52) overlaid.  What to expect: at large x a mild QD
recurrence MINIMUM (~0.75x elastic near x ~ 2.4, visible already in
the first three points), then the upturn and divergence near
x_cr ~ 0.86 for these quasi-dynamic isolated strips (the offset from
their 1.52 is inertia plus fault periodicity: regime and scaling
reproduce, the critical point is physics-sensitive).  Below x_cr the
patch is permanently stuck (no catalog: their ST class).

## 2. In-plane dipping fault, plate over Maxwell half-space (+gravity)

    ./run_bp3ve_demo               # one case directory, kernels cached per tM
    ./run_many_bp3ve_demo          # the three coupling cases, both branches, both dips

The BP3 dipping fault (dip 60, surface-breaking to 34.6 km, or dip 30
to 20 km with the tenth argument of run_bp3ve_demo) contained in a
40-km elastic plate over a Maxwell half-space, run through -ve_mode 3
with external kernels from test_relax/inplane_ve_proto/
bp3_ve_kernels.py.  Three coupling cases: ve_g1 (shear and normal
tractions relax, gravity as interface buoyancy), ve_g0 (same without
gravity), ve_shear (shear only, sigma_n frozen, also in the elastic
reference, so its reference cycle differs).  Branch +1 (vpl > 0) is
THRUST, -1 is NORMAL (_normal suffix); dip 30 adds _dip30.  Between
2026-09-01 and 09-06 the branch labels were the other way round.  The
driver runs with -ve_h_init 0 (zero initial memory); the earlier
default produced a secular slip deficit over the tau ladder and the
"+14 percent" quoted in older versions of this file came from it and
is withdrawn.  The hereditary integration is checked against the
analytic response of the kernel file (both traction families, 4
digits); the kernel physics (plate over Maxwell half-space, buoyancy
rows) is checked against PSGRN/PSCMP to 0.5 percent without gravity
(test_relax/psgrn_rundle/README.md), with the known caveat that the
buoyancy-row gravity overstates the early-time gravity effect by
about a factor two relative to PSGRN's full treatment.

Analysis: bp3ve_analyze.py reads the catalogs of every case
directory, finds the settled state of each run as the longest
periodic tail of the interevent series (tolerance 0.5 yr, periods up
to 40, means over whole periods), flags cycles whose members still
creep, and writes a summary table and figures plus a small tarball
of the catalogs.  bp3ve_overview.py draws event slip vs time for
elastic and selected Maxwell times (bp3ve_overview.png).  Supporting
tests: test_restart_ve.sh (continuous vs checkpoint-extended run,
bitwise identical at 100 m) and run_res_test.sh with bp3ve_compare.py
(50 m vs 25 m: settled states and transition times agree to within
one cycle).

### What the 25 m, 15000 yr sweep shows (rtol 1e-5, 3bs; read 2026-10-24)

Times in Maxwell times t_M = eta/mu; T is the elastic mean recurrence
(85 yr thrust, 95 yr normal at dip 60; 76 and 92 yr at dip 30).  Every
run was taken to 15000 yr and only the settled cycles are quoted;
earlier readings of 5000 yr runs were window averages across slow
transitions and are superseded.

Dip 30 (elastic cycles are simple 2-cycles, {64.7, 86.6} thrust and
{74.1, 109.6} normal, matching README_bp3.md):

- Thrust: no change of recurrence within 0.3 percent at any t_M from
  25 to 1000 yr, with or without the normal-stress channel or
  gravity; per-event slip changes by less than 1 percent.
- Normal: the long interval shortens while the short one stays at
  74 yr; mean change -4.9, -2.8, -1.6, -0.6, -0.1 percent at t_M =
  25, 50, 100, 250, 1000 yr without gravity (0.27 to 11 T), about half
  that with gravity, zero with the normal-stress channel frozen.  So
  for a simple cycle the effect is a few percent at most, largest
  near t_M ~ 0.3 T, roughly 1/t_M beyond, gone by 10 T, and carried
  by the normal-stress channel.

Dip 60 (the elastic thrust cycle is a 4-cycle {63.3, 87.1, 91.2,
97.5}, the community-majority BP3 state, and sits near a bifurcation):

- Thrust, sigma coupled (gravity on or off): for every t_M from 50
  (g0) / 75 (g1) to 1000 yr the settled state is the 3-cycle {65,
  86.5, 91}, mean recurrence -4.45 to -5.1 percent, with 30 to 51
  settled cycles per run.  The effect does not decay up to t_M = 12 T.
  What grows with t_M is the time to reach it, t_settle ~ 2300 +
  4.7 t_M yr (g0) and 2050 + 5.7 t_M (g1): the elastic 4-cycle is
  unstable at all these t_M with a growth rate ~ 1/t_M, and the end
  state is the same.  Fast relaxation (t_M = 5, 15) gives exact
  1-cycles (-0.5, +3.5 to +4.9 percent).  At t_M = 25 to 50 the
  system is at the edge of the 4-cycle's stability: a creeping
  5-cycle (g1, 25), type-I intermittency (g0, 35), or the 4-cycle
  retained (g1, 35 and 50), and which of these occurs depends on
  perturbations at the 1e-15 level; no single number should be
  quoted there.
- Thrust, shear only: -7.2 percent (4-cycle {69, 87, 91.5, 103}) for
  t_M = 50, 75, 100 with transition times 2400, 3550, 11060 yr;
  t_M = 150 not settled by 15000 yr; t_M >= 250 on the elastic
  2-cycle within 0.4 percent.  The transition time diverges between
  t_M = 100 and 250.
- Normal, sigma coupled: 1-cycle 95.1 -> 2-cycle {75, 99.5}, -8.8 to
  -7.6 percent for t_M = 35 to 500 yr; at t_M = 500 the transition
  comes at 10200 yr.
- Per-event slip and stress drop scale with the interval; the
  ruptured extent changes by less than 1 percent.  None of this
  identifies a mechanism; the station-series reading (HANDOFF 6)
  excluded the interseismic sigma path.

Reading.  For a simple cycle, viscoelastic relaxation of the
substrate changes recurrence by at most a few percent with the
expected t_M dependence and only through the normal-stress channel.
For a multi-periodic cycle near a bifurcation it acts as a selector
between coexisting attractors: a step of 5 to 8 percent that persists
to arbitrarily long t_M while the time to switch grows linearly with
t_M.  The dip-60 BP3 thrust is in the second regime.  A recurrence
versus t_M curve is therefore the wrong summary for it; the attractor
map is.

## Conventions

No environment variables; positional arguments with defaults
throughout.  Both drivers skip completed catalogs on re-invocation
and extend a case from its own checkpoint when called with a longer
stop time (delete cat_*.dat to force reruns).  Event extraction uses the
rsf_catalog files with a slip threshold for the in-plane case
(under-resolved meshes produce small-event chatter; the demos
compare LARGE-event recurrence at fixed mesh).
