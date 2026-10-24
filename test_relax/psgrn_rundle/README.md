# psgrn_rundle: comparison of the 2-D plate-over-Maxwell codes with PSGRN/PSCMP

This directory holds the scripts, inputs and logs of the comparison
between interact's 2-D in-plane viscoelastic machinery (the physics
behind rsf_solve `-ve_mode 3` kernels, `test_relax/inplane_ve_proto/
inplane2d.py`, and the independent surface-displacement code
`test_relax/ve_surface/ve_thrust_relax.py`) and the 3-D layered
viscoelastic-gravitational dislocation code PSGRN/PSCMP (Wang et al.
2006; 2020 version, source github.com/RongjiangWang/PSGRN-PSCMP_2020,
built by `tools/psgrn/install_psgrn`).  It also contains the test of
the finite-fault-length explanation for the one literature anchor the
2-D codes did not reproduce, Rundle (JGR 1982) Figure 2 at 45 tau_a.

Carried out 2026-10-05 and 06.  Numbers below are for the stated
configuration and should not be read as general error bounds.

## 1. Problem and parameters

Elastic plate of thickness H = 30 km (vp 5.2223, vs 3.0151 km/s, rho
3300 kg/m3, i.e. mu = lambda = 30 GPa) over a half-space that is
elastic in bulk and Maxwell in shear (vp 4.8676, vs 2.8103, rho 3800;
eta = 9.468e17 Pa s, so tau_M = eta/mu = 1.0003 yr).  The relaxation
time is the Maxwell time eta/mu throughout interact; Rundle's tau_a
and Johnson's tR are 2 eta/mu, so Rundle's 5 and 45 tau_a are 10 and
90 tau_M here.  PSGRN relaxes the shear modulus at constant bulk
modulus, as inplane2d.py and ve_thrust_relax.py do by default; the
alternative lambda = const convention changes the 90 tau_M basin by
0.1 percent and the 10 tau_M basin by 4 percent of itself
(ve_thrust_relax.py `--lam-const`), so it is not a factor at the level
of this comparison once stated.

Fault: 30 degree thrust, surface breaking, 30 km down-dip (bottom at
15 km = 0.5 H), uniform slip of 1 m, applied as a step at t = 0.  This
is the geometry of Rundle (1982) Figures 2 and 3 and of
`inplane_ve_proto/rundle82_compare.py`.

Gravity: PSGRN's gravity factor 0 or 1 (0 switches all gravitational
terms off; 1 is its full treatment with pre-stress advection and the
gravitational potential perturbation).  The 2-D codes have no gravity
or the interface-buoyancy formulation: sigma_zz(0) = rho_1 g u_z at
the free surface and a jump (rho_2 - rho_1) g u_z in sigma_zz at the
layer base, with the Eulerian incremental stress in the bulk equations
and no bulk gravitational terms.

## 2. How the 3-D runs are set up so that they have a 2-D limit

PSCMP frame: strike 90 (east), dip 30, so the fault dips south;
reference point of the fault plane at the trace start (O_depth 0), at
latitude 0 and a longitude such that the fault is centred on
longitude 0; length L along strike, 30 km width, patches 5 x 2 km
(2.5 x 1 km in one test), slip_ddip = -1 m (PSCMP's down-dip slip is
positive for normal faulting).  Receivers: 101 points along longitude
0 from 90 km north (footwall) to 160 km south (hanging wall), 2.5 km
spacing; one of them (latitude 0) sits on the trace.

Mapping to the 2-D frame of the codes: x (positive towards the
hanging wall, in units of H) = -latitude x 111.195 km / H; u_x
towards the hanging wall = -U_north; uplift = -U_down.

Fault lengths: L = 200 km is Rundle's fault (half-length 10 H / 3,
2L = 6.7 H); L = 600 km is the length used in the older
`inplane_ve_proto/psgrn_compare.py`; L = 2000 km (67 H) is treated as
the 2-D limit.  The L = 600 and 2000 km profiles differ by 0.5
percent at 90 tau_M and 0.1 percent at 10 tau_M, so 2000 km is 2-D to
better than that.

Time: PSGRN's Green's functions are sampled uniformly over 90 yr
(128 or 256 samples; the two give identical basins to 0.1); PSCMP
writes the receiver time series, whose first row (t = 0) is the
coseismic field.  Values at 10 and 90 tau_M are linear interpolations
in time.  The compared quantity is the postseismic CHANGE, u(t) -
u(0), as in Rundle's figures, times 100 per metre of slip.

Comparison metrics (psgrn_rundle_compare.py, psgrn_convergence_compare
.py): basin depth (minimum of the uplift change) and its position,
zero crossings of the profile, rms difference over the 101 receivers,
depth ratios between fault lengths, and for the coseismic field the
maximum and rms difference with the receiver on the trace excluded
(|x| > 0.1 H), since that point sits on the coseismic step of 0.45
(vertical) and 0.85 (horizontal) and the two codes resolve it
differently.  Rundle's values are points digitised from his Figures 2
and 3 (reading error about 2 units), listed in the script.

## 3. Numerical parameters of PSGRN/PSCMP and what they changed

The templates inherited from `inplane_ve_proto/psgrn_inputs` had 126
distances from 0 to 500 km, integration accuracy 0.05, and 8 source
depths from 1 to 15 km.  Variants were run on the L = 2000 km,
no-gravity case (run_psgrn_convergence.sh; basin depths x 1e-2 slip,
2-D code values -26.1 at 10 tau_M and -87.1 at 90 tau_M):

    variant                      10 tau_M   90 tau_M   note
    base (templates)              -25.6      -85.6     2.0 / 1.7 percent shallow
    accuracy 0.01                 -25.6      -85.6     no change
    accuracy 0.005                -25.6      -85.6     no change
    252 distances                 -25.6      -85.6     no change
    504 distances                 -25.6      -85.5     no change
    distances to 1200 km          -25.6      -85.5     no change
    29 depths 0.5 to 14.5 km      -26.0      -87.6     0.3 / 0.6 percent
    all of the above              -26.0      -87.6     same as the depth grid alone
    all, 2.5 x 1 km patches       -26.0      -87.6     same

The depth grid was the only parameter that mattered: the shallowest
PSCMP patch row has its centres at 0.5 km depth, below the first
Green's function depth of 1 km, and was extrapolated; since the
shallow rows carry the largest surface response, a 0.5 km geometric
error produced a 2 percent amplitude error.  With 29 depths placed on
the patch centres the no-gravity agreement is -0.2 percent at
10 tau_M and +0.5 percent at 90 tau_M (rms 0.05 and 0.27 x 1e-2 slip,
zero crossings within 0.01 H), and the sign of the residual has
changed, so what remains is a different and smaller term that none of
the tested parameters touch.  The psgrn_inputs templates were
corrected to the 29-depth grid on 2026-10-06; the older
`inplane_ve_proto/psgrn_compare.png` was made with the 8-depth grid
and its 1 to 3 percent residuals are mostly this.  The default in
run_psgrn_rundle.sh is now the 29-depth grid as well, so a fresh run
of it gives the corrected numbers; comp.log in this directory is from
the first run with the template grid.

The 2-D codes' own numerics were converged separately
(ve_surface/README_ve_thrust_relax.md: Talbot terms, source counts,
k grids and window at 1e-5 to 1e-4 slip; inplane2d.py and
ve_thrust_relax.py agree with each other to 1e-3 slip).

## 4. Results

Without gravity.  PSCMP in its 2-D limit agrees with the 2-D codes to
0.2 percent of the basin at 10 tau_M and 0.5 percent at 90 tau_M.  At
Rundle's fault length PSCMP gives -25.1 and -73.6 against his -26 and
-75 (rms 1.2 and 7.5 of peaks 26 and 75); the depth factors relative
to the 2-D limit are 0.983 and 0.860, which is exactly the factor by
which the 2-D infinite-fault basin (-87) exceeded his.  The 15
percent discrepancy of the 2-D codes with Rundle at 45 tau_a is
therefore his finite fault length.  His hanging-wall flank at 45
tau_a is steeper than PSCMP's at the same length (zero crossing 3.0
vs 3.8 H), which remains unexplained and is within what a digitised
1982 figure may carry.  A finite-length factor of 0.868 at 45 tau_a
obtained independently from the Johnson plate-over-Maxwell code
(test_relax/plate_over_maxwell_c, with the caveats in its README)
agrees with PSCMP's 0.860.

With gravity.  With the corrected depth grid, PSCMP (full gravity) is
4.5 percent deeper than the interface-buoyancy 2-D codes at 10 tau_M
and 1.3 percent shallower at 90 tau_M (-24.7 vs -23.6 and -43.7 vs
-44.3).  Relative to the no-gravity values, PSGRN's gravity reduces
the 10 tau_M basin by 1.3 units and the buoyancy rows by 2.5, so the
buoyancy-row formulation overstates the early-time gravity effect by
about a factor two while the two agree on the late-time (isostatic)
effect; the same signature appears in `inplane_ve_proto/
rundle82_fig3_compare.py` (9 percent rms at 5 tau_a, 1.3 percent at
45 tau_a).  At Rundle's fault length PSCMP gives -42.6 against his
Figure 3 value of -45.  The gravity suppression of the late-time
basin (-87 to -44) is reproduced by PSGRN, both 2-D codes and Rundle.

Coseismic.  PSCMP's t = 0 field against the 2-D code, trace receiver
excluded: uplift rms 0.23 and maximum 0.43, horizontal rms 0.08 and
maximum 0.34 (x 1e-2 slip, peaks 38 and 51) without gravity; with
gravity the uplift maximum is 1.09, confined to the near field, where
PSGRN's gravity and the buoyancy rows differ already at t = 0.

## 5. Interpretation of the gravity difference

The interface-buoyancy formulation keeps the Eulerian incremental
stress and the two interface conditions but drops the body force of
the density perturbation of compression, -rho g (div u) z, of
relative size rho g H / mu = 0.03 for these parameters.  That term
acts while the plate is being compressed, i.e. early, and not in the
isostatic end state, which is the pattern observed.  Writing the
problem instead with the Hookean (incremental Lagrangian) stress,
which is then continuous at interfaces and zero at the free surface,
puts the gravity into bulk terms (+rho g d_x u_z in the x equation,
-rho g d_x u_x in the z equation); ve_thrust_relax.py `--gravity 2`
does this and reproduces PSGRN's 10 tau_M basin and zero crossings,
but without self-gravitation that system has imaginary eigenvalues
below k H ~ rho g H / mu and is ill-behaved for a few times that, so
its late-time result depends on a wavenumber cutoff and it is marked
experimental.  A well-posed version needs the gravitational potential
perturbation (a 6 x 6 system, as in PSGRN); it has not been written.
For the rsf_solve kernels, which use the buoyancy rows, the
consequence is an error of a few percent in the gravity part of the
relaxation, which in the BP3 cycle sweeps is itself a one-to-two
point modifier of recurrence.

## 6. Reproducing

Build PSGRN/PSCMP once: `tools/psgrn/install_psgrn <dir>` (gfortran).
Then, from an empty scratch directory with a few GB free (Green's
functions are large and reusable):

    <repo>/test_relax/psgrn_rundle/run_psgrn_rundle.sh <dir>/bin
    python3 <repo>/test_relax/psgrn_rundle/psgrn_rundle_compare.py . "200 600 2000" "0 1" \
        <repo>/test_relax/ve_surface/ve_thrust_relax.py
    <repo>/test_relax/psgrn_rundle/run_psgrn_convergence.sh <dir>/bin "base nz29 all all_g1"
    python3 <repo>/test_relax/psgrn_rundle/psgrn_convergence_compare.py . "base nz29 all all_g1" \
        <repo>/test_relax/ve_surface/ve_thrust_relax.py

All arguments are positional with defaults (see the script headers);
PSGRN runs take tens of minutes to an hour per Green's function set,
PSCMP minutes.  The compare scripts run ve_thrust_relax.py themselves
for the 2-D reference (about 20 s per gravity setting).

Files: run_psgrn_rundle.sh, psgrn_rundle_compare.py,
run_psgrn_convergence.sh, psgrn_convergence_compare.py; psgrn_*.log
and pscmp_*.log of the first run; comp.log (first run, template depth
grid); psgrn_rundle.png and psgrn_rundle_horizontal.png (uplift and
horizontal change, four panels: gravity off and on at 10 and 90
tau_M, every fault length, the 2-D curve, Rundle's points).

## 7. Open items

- The 0.5 percent residual at 90 tau_M without gravity, smooth and
  independent of every PSGRN/PSCMP parameter tested; not pursued.
- The hanging-wall flank shape of Rundle's 45 tau_a curve.
- A gravity formulation with the potential perturbation for the 2-D
  codes, if the early-time gravity effect is ever to be quoted to
  better than a factor two.
- PSCMP at L = 100 km and the coseismic epoch could also serve to
  check the Johnson-code dip-slip orientation question recorded in
  test_relax/plate_over_maxwell_c/README.md.
