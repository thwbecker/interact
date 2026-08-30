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

## 2. In-plane dipping thrust, plate over Maxwell half-space (+gravity)

    ./run_bp3ve_demo               # kernels ~10 min per tM (cached), runs minutes

The BP3 dip-60 geometry (surface-breaking to 34.6 km depth)
contained in a 40-km elastic plate over a Maxwell half-space with
gravity, run through -ve_mode 3 with external kernels from
test_relax/inplane_ve_proto/bp3_ve_kernels.py (the PSGRN-validated
plane-strain machinery; see tools/psgrn/README_psgrn.md for the 3-D
PSGRN route with the same file contract).  Plots large-event
recurrence vs relaxation strength: contained-fault backslip loading
lengthens recurrence (+14 percent at t_M = 45 yr ~ 0.36 T_rec in the
demo), the loading-pathway mechanism of the antiplane
ve_loading_conditions.md map in its in-plane form, i.e. the
Rundle (1982) / Thatcher & Rundle (1984) regime with resolved
rate-and-state events instead of imposed ones.  sigma_n held fixed
(hereditary normal-stress relaxation not yet implemented).

## Conventions

No environment variables; positional arguments with defaults
throughout.  Both drivers skip completed catalogs on re-invocation
(delete cat_*.dat to force reruns).  Event extraction uses the
rsf_catalog files with a slip threshold for the in-plane case
(under-resolved meshes produce small-event chatter; the demos
compare LARGE-event recurrence at fixed mesh).
