# psgrn_rundle: PSGRN/PSCMP reference runs for the 2-D plate-over-Maxwell thrust

3-D reference (Wang et al. PSGRN/PSCMP, 2020 version, built by
tools/psgrn/install_psgrn) for the surface-breaking 30 deg thrust of
Rundle (JGR 1982) Figs. 2 and 3: H = 30 km elastic plate (mu = lam =
30 GPa, rho 3300) over a Maxwell half-space (30 GPa, rho 3800, tau_M =
eta/mu = 1.0 yr), fault 30 km down-dip (bottom at 15 km), 1 m slip,
mid-strike profile from 3 H on the footwall to 5 H on the hanging
wall.  Times are Maxwell times tM = eta/mu; Rundle's tau_a = 2 eta/mu,
so his 5 and 45 tau_a are 10 and 90 tM.  Displacements x 100 per metre
of slip, uplift positive, horizontal positive towards the hanging wall.

Scripts (positional arguments with defaults, no environment variables):

- run_psgrn_rundle.sh: Green's functions with gravity factor 0 and 1,
  PSCMP faults of length 200 (Rundle's 2L = 6.7 H), 600 and 2000 km
  (2-D limit); psgrn_rundle_compare.py plots and tabulates against
  Rundle's digitised points and against the 2-D code
  (../ve_surface/ve_thrust_relax.py).  Output psgrn_rundle.png,
  psgrn_rundle_horizontal.png, comp.log.
- run_psgrn_convergence.sh: numerical-parameter variants of the 2-D
  limit case (integration accuracy, distance sampling and range, source
  depth grid, PSCMP patch size, gravity); psgrn_convergence_compare.py
  tabulates PSCMP minus 2-D.  Output psgrn_convergence.png.

Results (2026-10-05/06, logs in this directory):

- Rundle Fig. 2 (no gravity): with L = 200 km PSCMP gives basins of
  -25.1 and -73.6 at 10 and 90 tM against his -26 and -75 (rms 1.2
  and 7.5 of peaks 26 and 75); the finite-length depth factors
  relative to L = 2000 km are 0.983 and 0.860.  The 2-D codes
  (infinite fault) give -26.1 and -87.1; the 15 percent depth
  discrepancy with Rundle at 45 tau_a is his fault length.  His
  hanging-wall flank at 45 tau_a is steeper than PSCMP's (zero crossing
  3.0 vs 3.8 H); unexplained, within what a digitised 1982 figure can
  carry.  Fig. 3 (gravity): L = 200 km gives -42.6 against his -45.
- 2-D limit vs the 2-D codes: the PSGRN source-depth grid of the
  inplane_ve_proto templates (8 depths from 1 km) extrapolated the
  shallowest PSCMP patch row (0.5 km) and made PSCMP 1.7 to 2.5 percent
  shallow; with 29 depths on the patch centres (0.5 to 14.5 km) the
  no-gravity agreement is -0.2 percent at 10 tM and +0.5 percent at
  90 tM (rms 0.05 and 0.27 x 1e-2 slip), unchanged by integration
  accuracy 0.05 to 0.005, distance sampling 126 to 605, range 500 to
  1200 km, or patch size 5 x 2 to 2.5 x 1 km.  The templates were
  corrected (2026-10-06).
- Gravity formulation: with the depth grid fixed, PSCMP with gravity
  is +4.5 percent deeper at 10 tM and -1.3 percent at 90 tM than the
  interface-buoyancy 2-D codes (-24.7 vs -23.6, -43.7 vs -44.3).
  Relative to no gravity, PSGRN's gravity reduces the 10 tM basin by
  1.3 units, the buoyancy rows by 2.5: the buoyancy-row formulation
  (used by inplane2d.py, bp3_ve_kernels.py and hence rsf_solve
  -ve_mode 3, and by ve_thrust_relax.py --gravity 1) overstates the
  early-time gravity effect by about a factor two and agrees late.
  The dropped term is the body force of the density perturbation of
  compression, rho g div u (relative size rho g H / mu = 0.03).
  See ../ve_surface/README_ve_thrust_relax.md for the partial
  (ill-posed without self-gravitation) pre-stress-advection test.
