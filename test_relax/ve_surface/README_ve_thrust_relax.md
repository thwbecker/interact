# ve_thrust_relax.py: postseismic relaxation of a 2-D thrust in an elastic layer over a Maxwell half-space

Standalone Python (numpy, matplotlib).  Written 2026-10-05.

Note on provenance: the repository's test_relax/inplane_ve_proto
(inplane2d.py) solves the same problem with the same approach
(propagator in k, buoyancy rows, correspondence principle, Talbot)
and is validated against Rundle (1982) and PSGRN; it was written as
the kernel generator for rsf_solve and evaluates fields on a receiver
grid.  This script was written independently and later compared with
it (see Verification); the two agree to 1e-3 slip.  It differs in
having the elastic half-space part in closed form in x (so the trace
is resolved), a graded k grid with a Filon rule for the far field, a
command-line interface for dip, extent, taper and times, and in
running in about 10 s per case instead of a minute.

## Problem

Plane strain.  Elastic layer of thickness H (z down, free surface at
z = 0) over a half-space that is elastic in bulk and Maxwell in shear
with relaxation time tM = eta / mu2.  A fault dipping at `dip` towards
+x from the surface trace at x = 0 (or from depth `top` H) to depth
`extent` H carries a prescribed slip distribution (uniform, cosine
tapered at the bottom, or half sine), applied as a step in time.
Thrust sense: the hanging wall (the +x side at the surface) moves
up-dip.  Surface displacements are computed at t = 0 (coseismic), at
the requested multiples of tM, and for the relaxed limit mu2 -> 0.

Gravity enters as buoyancy at the two density contrasts,
sigma_zz(0) = rho1 g u_z(0) at the free surface and
[sigma_zz] = (rho2 - rho1) g u_z at the layer base.  This is the usual
interface-buoyancy approximation (as in Rundle-type plate-over-
asthenosphere models); there is no self-gravitation and no advection
of pre-stress inside the layers.  With `--gravity 0` both terms are
off.

## Method

Fourier transform in x.  In z the exact exponential basis of the
plane-strain system is used (eigenvalues +-|k| with Jordan chains;
the four vectors are written out in `basis()` and checked against the
system matrix), so no growing exponential is ever formed and the
linear system is well conditioned at every wavenumber.  The fault is
a line of point double couples; each enters as a jump of the
displacement-stress vector at its depth (derived from the body-force
equivalent, `source_jump()`).  The Maxwell half-space enters through
mu2(s) = mu2 s tM / (1 + s tM), lambda2(s) = K2 - 2 mu2(s)/3, and the
time domain is recovered by the fixed Talbot inversion (Abate and
Valko 2004; applied separately to the real and imaginary parts of the
complex k-domain response, which needs the transform at the conjugate
nodes as well).

The response is assembled as the homogeneous elastic half-space
response plus the layered viscoelastic difference.  The half-space
part is evaluated in closed form in x: the k-domain response of a
point double couple at depth d is (m/mu) e^{-|k| d}(p0 + p1 |k| d)
with coefficients that are rational in lambda and mu (read off the
6x6 layer solution, which it matches to 1e-15), and its inverse
transform is a sum of Lorentzians and their Hilbert partners
(`halfspace_surface_x`).  The layered difference decays like
exp(-k H (2 - extent)) and is transformed on a graded k grid (dense
below k H = 0.2 for the slow long-wavelength relaxation, then
pi/(4 xmax)) with a Filon-trapezoid rule, which integrates e^{ikx}
exactly against a piecewise-linear U(k) so that the accuracy does not
degrade towards the edge of the x window (a plain trapezoid rule gave
5e-4 wiggles there, which is what the first version of the far-field
figures showed).  The far-field constant offset and the odd 1/x tail
of the horizontal field are removed in k and added back analytically.

## Verification

- Basis vectors satisfy A v1 = -k v1, (A + k) v2 = v1 and the
  growing counterparts to 1e-16, also for complex moduli.
- Talbot: 1/(s(s+1)) -> 1 - e^{-t} to 3e-12; complex-valued test
  function to 1e-12.
- Closed-form half-space path vs the general layer solver with
  mu2 = mu1: 1e-15; the general solver is independent of where a
  fictitious interface is placed in a homogeneous medium (1e-15); the
  closed form in x vs the k-integration of the same response: 1.5e-6
  (the k-integration's error).
- Coseismic homogeneous half-space vs cutde (Nikkhoo and Walter
  triangular dislocations, half-space, 6000 km long rectangle): dip
  30, bottom 0.5 H, uniform slip, 40 points from x = -100 H to 100 H
  including +-0.2 H: max difference 4.6e-5 slip in u_x and 3.9e-5 in
  u_z.  Jumps at the trace: u_z 0.45 (slip sin 30 = 0.5 smoothed over
  the source spacing), u_x 0.85 (cos 30 = 0.866).
- Literature and the interact generator (rundle82_compare_both.png,
  made with `--H 30 --mu1 30 --rho1 3300 --rho2 3800 --g 9.8
  --lam-const --times 10 90`): Rundle (JGR 1982) Figures 2 and 3 are
  this configuration (30 deg surface-breaking thrust to 0.5 H, lambda
  held constant), at 5 and 45 tau_a = 10 and 90 tM.  Against the
  repository's test_relax/inplane_ve_proto (inplane2d.py), the change
  in uplift agrees to 0.1 to 0.14 x 1e-2 slip rms, peaks -27.1 vs
  -27.2 and -87.1 vs -87.3 (no gravity), -23.8 vs -23.9 and -44.2 vs
  -44.4 (gravity).  Against Rundle's digitised points the 2-D codes
  match at 5 tau_a (1.5 rms of a 26 peak) but give a basin of -87 at
  45 tau_a where he shows -75, with wider flanks.
- PSGRN/PSCMP, 3-D, at Rundle's fault length (test_relax/psgrn_rundle
  in the repository, run 2026-10-05 with run_psgrn_rundle.sh): with
  L = 200 km (his 2L = 6.7 H) PSCMP gives -25.1 at 5 tau_a and -73.6
  at 45 tau_a without gravity, against his -26 and -75; the
  finite-length depth factor relative to L = 2000 km is 0.983 and
  0.860, i.e. exactly the 0.86 the 2-D codes needed, so the depth
  discrepancy is his finite fault length.  What remains is the shape
  on the hanging-wall flank at 45 tau_a: PSCMP's zero crossing is at
  3.8 H against his 3.0 H (rms 7.5 of 75 there, 1.2 of 26 at 5 tau_a);
  the footwall flank and the far lobes agree.  With gravity PSCMP
  L = 200 gives -42.6 against his Fig. 3 -45.  In the 2-D limit
  (L = 2000 km) PSCMP is 1.7 to 2.5 percent shallower than both 2-D
  codes in all four cases (-85.6 vs -87.1, -43.1 vs -44.2 at 90 tM;
  -25.6 vs -26.1, -24.3 vs -24.6 at 10 tM, constant bulk modulus),
  the same offset the earlier psgrn_compare showed.  Resolved with
  run_psgrn_convergence.sh (test_relax/psgrn_rundle): the PSGRN
  source-depth grid of the templates (8 depths from 1 to 15 km)
  extrapolates the shallowest PSCMP patch row (centres at 0.5 km);
  with 29 depths on the patch centres (0.5 to 14.5 km) the offset is
  -0.3 percent at 10 tM and +0.6 percent at 90 tM (rms 0.05 and 0.31
  x 1e-2 slip, zero crossings within 0.01 H), while integration
  accuracy (0.05 to 0.01), distance sampling (126 to 252) and distance
  range (500 to 1200 km) each change the result by 0.1 percent or
  less.  The inplane_ve_proto psgrn_compare templates carry the same
  depth grid, so their 1 to 3 percent residuals are mostly this.
- Gravity formulation (run_psgrn_convergence.sh, all_g1): with the
  depth grid fixed, PSCMP with gravity is 4.5 percent deeper than the
  buoyancy-row 2-D codes at 10 tM and 1.3 percent shallower at 90 tM
  (-24.7 vs -23.6 and -43.7 vs -44.3).  Against the no-gravity values,
  gravity reduces the 10 tM basin by 1.3 units in PSGRN and by 2.5 in
  the buoyancy-row formulation, which therefore overstates the
  early-time gravity effect by about a factor two while agreeing on
  the late-time (isostatic) effect; the same pattern appears in the
  repository's Rundle Fig. 3 comparison (9 percent rms at 5 tau_a, 1.3
  percent at 45).  The term the buoyancy rows drop is the body force
  of the density perturbation of compression, rho g div u, of relative
  size rho g H / mu = 0.03.  `--gravity 2` adds the pre-stress
  advection terms (Hookean stress continuous, traction-free surface,
  bulk terms +rho g d_x u_z and -rho g d_x u_x) and reproduces the
  PSGRN 10 tM basin and zero crossings (-24.5, -0.88/2.15) when used
  above k H = 0.25, but without self-gravitation the system has
  imaginary eigenvalues below k H ~ rho g H / mu and is ill-behaved
  for a few times that, so the late-time result depends on the cutoff
  (-45.5 vs -50.6 for kcut 0.25 vs 0.11); it is EXPERIMENTAL.  A
  well-posed version needs the gravitational potential perturbation
  (6x6 system, as in PSGRN).  For the cycle kernels this matters at
  the level of half the early-time gravity effect, i.e. a few percent
  of the postseismic signal; the no-gravity and late-time physics are
  verified to 0.5 percent.  The Johnson-code finite-length factor (0.868, with
  its mirrored-mechanism caveat) agrees with PSCMP's 0.860.
- Limits: t = 0.001 tM reproduces the coseismic layered solution to
  4e-5; t = 500 tM reproduces the relaxed solution to 1e-5 (gravity
  on).
- Relaxed far-field horizontal offset between the two sides equals
  the Saint-Venant value (integral of the horizontal Burgers vector
  over depth) / H in every case tried: 0.4327 (dip 30, extent 0.5:
  0.866 x 0.5), 0.8654 (extent 1), 0.2498 (dip 60), 0.2757 (half
  sine: 0.866 x 2/pi x 0.5).
- Convergence at the defaults (gravity on and off, x = +-1, 5, 20,
  90 H, t = 2 and 50 tM): Talbot 24/32/48 terms 1e-9; layer sources 24
  vs 48 8e-5 (near the trace); dense small-k points 400 vs 1600 2e-6;
  k oversampling 4 vs 8 and window 100 vs 200 H both 1.3e-5 (far-field
  u_x at 50 tM); k cutoff 40 vs 60 /H 1e-6.  Overall accuracy about
  1e-4 slip, set by the source discretisation near the trace.

## Driver

`run_thrust_relax.sh [dip] [extent] [times] [H] [xmax] [prefix]` runs
one fault (default: the SEAS BP3 geometry of run_bp3ve_demo, dip 60 to
0.866 H in the 40 km plate; the Rundle 1982 geometry, dip 30 to 0.5 H,
is given as a commented example) without and with interface-buoyancy
gravity at the given Maxwell times, writes GMT-ready profile tables
and calls `plot_ve_thrust_gmt.sh` for a near-field and a full-window
PDF with u_x and u_z against gravity off and on.  Run it from a
scratch directory; it finds the Python tool through its own location.

## Usage

    python3 ve_thrust_relax.py -h
    python3 ve_thrust_relax.py --gravity 1 --xmax 100 --out g1      # default case
    python3 ve_thrust_relax.py --gravity 0 --xmax 100 --out g0
    python3 ve_thrust_relax.py --extent 1.0 --out e1                 # fault cuts the layer
    python3 ve_thrust_relax.py --taper sin --out sin1
    python3 ve_thrust_relax.py --dip 60 --out d60

Defaults: H 40 km, dip 30, bottom at 0.5 H, surface breaking, uniform
slip, mu1 = mu2 = 32.04 GPa, nu 0.25, rho 2800 / 3300, times 0.5 1 2
5 10 50 tM, x window +-50 H.  About 25 s per case (4001 x points).  Output <out>.npz
(x in H, times, ux, uz in units of slip, relaxed and half-space
fields), <out>.txt (far-field table), <out>_near.png (+-3 H),
<out>_far.png (full window, symlog).  Conventions: x positive towards
the hanging wall, u_x positive towards +x, u_z positive up.

## Results for the default case (dip 30, bottom at H/2, uniform slip)

Coseismic: hanging-wall uplift 0.41 slip at the trace falling to
-0.09 at x = 1 H, footwall subsidence -0.04; horizontal field
antisymmetric about the trace at large distance (footwall towards
+x, hanging wall towards -x), decaying as 1/x: 0.013 at 20 H.

With gravity, the surface evolves monotonically towards the relaxed
state, which is reached to within 10 percent by 50 tM at distances
below 5 H (horizontal) and 2 H (vertical) and more slowly further
out.  Horizontal: the two sides converge by a total of slip x cos(dip)
x (fault depth / H) = 0.433 (0.216 each side), constant out to the
edge of the window, i.e. the shortening of the cut upper half of the
plate is distributed over the whole thickness once the substrate no
longer transmits shear.  Vertical: the coseismic hanging-wall uplift
reverses; the fault zone sags to -0.48 slip at x = 1 H and -0.2 at
x = -0.5 H, with flanking bulges of +0.07 (x = -5 H) and +0.12
(x = +5 H) that decay to zero by 20 H.  The reason is that the thrust
shortens the upper half of the plate and leaves the lower half
intact, so the plate carries a bending moment that the relaxed
substrate cannot resist; the flexural wavelength (about 1.7 H for
these parameters) sets the width of the sag and bulges.  The
isostatic uplift of the thickened zone (of order slip x sin(dip) x
(rho2 - rho1)/rho2 = 0.08) is small against this.

Without gravity there is no bounded relaxed state: the plate kinks
about the fault and the far field grows without limit as the
substrate relaxes.  At finite time the fields are finite and differ
from the gravity case only beyond 5 tM: at 50 tM the sag at x = 1 H
is -0.71 instead of -0.47, the bulges at +-5 H are 0.12 to 0.14 and
still growing, and the horizontal convergence has spread to 100 H
(0.1 slip at x = -100 H) instead of levelling at 0.216.  The code
suppresses the mu2 -> 0 curve in this case.

Fault cutting the whole layer (extent 1, with gravity): the two plate
halves are decoupled, the full slip cos(dip) = 0.866 appears as
rigid convergence in the far field (0.433 each side), there is no
bending moment, and the hanging-wall uplift persists in the relaxed
state (0.22 at x = 0.5 H, 0.17 at 1 H, footwall -0.17 at -0.5 H).

## Limits

Linear, small displacements; buoyancy only at the two interfaces;
Maxwell in shear and elastic in bulk for the half-space; the layer
is elastic; slip is prescribed and does not evolve.  The source
discretisation smooths the trace step over the fault spacing (H/192).
For `extent` 1 the deepest source sits at the interface; the case
runs but the intended comparison is a fault that just reaches the
base of the layer.  Not included: a viscoelastic layer, a Burgers or
power-law substrate, afterslip, 3-D finite fault length.  The
plate-over-fluid bending result depends on the lower half of the
layer being uncut; in a setting where the whole lithosphere
shortens, extent 1 is the closer analogue.
