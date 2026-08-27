# 2-D antiplane cycle experiments: single and multiple faults, elastic and visco-elastic

This directory is the experimentation suite on top of the VE machinery
validated in ../anti_plane (see README_ve_antiplane.md there for the
full validation ledger).  One driver runs everything and renders the
figures:

    ./run_cycles2d profile ds_km D_km nfault spacing_km tM_yr H_km stop_yr outdir [extra rsf_solve options]

profile bp1 is the SEAS BP1-QD benchmark friction (community-anchored:
our recurrence 78.2 yr at fine resolution matches Erickson et al., SRL
2020); profile bp5 is the BP5-flavored depth profile (VS cap 0-4 km,
VW core, VS below; nucleation-seeded).  tM_yr = 0 runs elastic; tM_yr
> 0 adds a Maxwell half-space below an elastic plate of thickness H_km
(-ve_mode 2; keep the fault bottom D_km comfortably above H_km, see
CAVEATS).  Every run writes per-fault event catalogs
(rsf_catalog.gNNN.dat), monitors, and FIELD FRAMES of both slip rate
and shear stress (tmp_rsf/rsf_{vel,tau}.gNNN.FFFFFF.bin, float32
(x1,x2,value) triples on the solver's own step cadence, dense through
events; frame times in rsf_vel.times; -field_stress is on by default
here).  plot_cycles2d.py renders per run:

    spacetime_v.png / spacetime_tau.png         depth-time evolution (calendar time)
    spacetime_v_frames.png / _tau_frames.png    same on the frame axis (coseismic resolved)
    stress_profiles.png                         tau(z) cross sections through the last cycle
    event_raster.png                            (multi-fault) onsets per fault vs time

Gating test: ./run_cycles2d_test (2-3 min: BP1-QD anchor, resolution
convergence, VE regression, two-fault smoke test).

## Recipes for the planned experiments

Maxwell time vs recurrence time (single fault; elastic reference first):

    ./run_cycles2d bp5 0.5 20 1 10 0  22 2500 el_ref   -dc 0.08
    for tm in 1 2.5 5 12 25 50;do
       ./run_cycles2d bp5 0.5 20 1 10 $tm 22 2500 ve_tm$tm -dc 0.08
    done

Fault depth / plate thickness ratio (D/H):

    for H in 21 22 25 30 40;do
       ./run_cycles2d bp5 0.5 20 1 10 5 $H 2500 ve_H$H -dc 0.08
    done

Multiple faults, elastic vs viscoelastic (2, 8, 16, 32):

    for nf in 2 8 16 32;do
       ./run_cycles2d bp5 0.5 20 $nf 10 0 22 1500 el_${nf}f -dc 0.08
       ./run_cycles2d bp5 0.5 20 $nf 10 5 22 1500 ve_${nf}f -dc 0.08
    done

Reference numbers from this setup (steady, spun memory, elastic
recurrence 250.90 yr): VE lengthens mildly, +1.3 / +1.0 / +0.6 / +0.1
percent at (tM=2.5, H=22), (5, 22), (5, 25), (5, 40); robust to rtol,
np = 5/6, and virgin-vs-spun initialization (identical attractor).
Multi-fault systems desynchronize into migrating sequences and partial
ruptures; the VE runs additionally show long-range interaction through
the substrate (cf. Shi, Wei & Barbot, JGR 2022, for the 3-D analog).

## Practical notes

- Resolution: the generator prints Lb/3 and h*; stay below Lb/3 (bp5 at
  ds = 0.5 km has a 4x margin; bp1 needs ds <= 0.1 km, and even 0.2 km
  is only 0.5 percent off, see run_cycles2d_test).  If in doubt run a
  ds/2 pair; the BP1 series here converged 77.88 / 78.21 / 78.24 /
  78.19 yr at 200 / 100 / 50 / 25 m.
- Memory initialization: the driver defaults to -ve_h_init 1 (spun
  backslip memory), the right start for steady-cycle studies; virgin
  (VE_HINIT=0) needs ~90 tM of spin-up but converges to the same
  attractor for contained faults.
- Cost: elastic runs are fast (a 32-fault, 640-cell, 600-yr run is
  minutes).  VE runs cost more: assembly (seconds to ~1 min thanks to
  the translational-invariance sample cache) plus np = 6 matvecs per
  RHS stage for the stage-consistent memory sink.  -ve_h_stage 0 is
  about 4x cheaper but biases recurrence where relaxation is dynamically
  important (percent-level here, tens of percent in relaxation-critical
  problems); use it only for exploration and re-run keepers with the
  default.  Multi-fault VE runs spend most time in event cascades.
- Long runs: chain calls with the same outdir; the driver auto-restarts
  from rsf_checkpoint.bin, and checkpoints every 10000 accepted steps
  or 120 s of wallclock, whichever comes first (CKPT_INT / CKPT_WALL
  override; -rsf_checkpoint_wall is the underlying rsf_solve option),
  so interrupted runs resume with at most ~2 minutes lost.  The driver reports "reached t = ... of ...
  requested" at the end; an INCOMPLETE run was killed externally
  (queue/time limit), not by any internal step limit (max steps is
  effectively unbounded).  Note the VE runs are 3-5x slower than
  elastic at the same fault count; at 32-64 faults budget accordingly
  or chain restarts.

## CAVEATS

- Keep the fault contained (D_km < H_km): a fault cutting the entire
  plate has a stress-free relaxed slip mode, no steady cycle attractor,
  and its multi-kyr behavior depends on kernel-representation detail
  (README_ve_antiplane.md; the Kato 2002 replication in ../anti_plane).
- Contained faults over a shallow substrate LENGTHEN recurrence mildly
  (loading-pathway relaxation); through-plate faults shorten strongly
  but transiently.  Earlier notes reporting steady-state SHORTENING for
  contained faults were biased by the pre-fix lagged memory sink.
- rsf_solve is quasi-dynamic (radiation damping): near relaxation-
  critical regimes (recurrence diverging), inertia shifts thresholds
  (see the Miyake & Noda replication).
