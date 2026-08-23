# relax_fault_ve: prescribed-slip visco-elastic relaxation testbed

TWB / Claude, August 2026. Documents src/relax/relax_fault_ve.c and
the supporting machinery (src/relax/prony_kernel.c,
src/includes/prony_kernel.h, src/testing/ve_check.c,
src/testing/ve_laplace_check.c, test_relax/run_ve_tests,
test_relax/run_fault_relax_test) as of commit 99b3a13. Background
theory is in rsf_ve_design.md and the staged development plan in
rsf_ve_implementation_plan.md; this program is steps 1 and 2 of that
plan, the machinery test in front of the rsf_solve integration. The
starting point for the visco-elastic effort was DAM's note
(ve-note-v4.pdf, copy in src/relax/): its section 1 effective-modulus
scheme is implemented in the companion program relax_fault, and the
exact Prony representation used here grew out of working through the
note's sections 2 to 4, as documented in rsf_ve_design.md.

## What the program computes

A slip distribution is prescribed on the fault (elliptical taper over
a fraction range of the patch index), applied either instantaneously
at t = 0 or as a linear ramp of duration t_ramp, and then held. The
medium is homogeneous, Maxwell visco-elastic in shear with elastic
bulk response, with Maxwell time t_M. The program outputs the
time-dependent resolved shear stress on every patch and the
time-dependent displacement at a profile of free-surface observation
points.

The stress and displacement histories come from the correspondence
principle, evaluated exactly rather than through a time-stepping
approximation of the constitutive law: in the Laplace domain the
visco-elastic solution is the elastic solution at effective moduli
G_bar(s) = G s/(s + 1/t_M) and nu_bar(s) from the fixed bulk modulus,
and for this material the resulting step-response kernels are finite
sums of exponentials. For the kernel families tested here the pole
set consists of at most three simple poles at material-determined
times, t_M, 3(1-nu)/(1+nu) t_M, and 3/(2(1+nu)) t_M, plus a constant
(fully relaxed) term. Amplitudes are obtained by evaluating the
existing elastic Green's functions at a small number of sampled
effective moduli and solving a small weight system; two extra sample
points are held out and used to verify per entry that the pole set is
complete, so an incomplete representation would be detected rather
than silently fitted. The constant term is carried throughout even
where it vanishes (homogeneous stress kernels) because displacement
kernels and the planned layered extension both need it.

Two properties of the homogeneous configuration are worth knowing
when interpreting output. Stress on the fault relaxes completely,
with the late-time decay governed by the slowest pole,
3(1-nu)/(1+nu) t_M, not by t_M itself. Displacements do not relax to
zero but to the elastic field of the effectively incompressible
(nu = 1/2) limit, the permanent deformation, and their transient
carries no t_M pole because displacements from dislocations depend on
the moduli only through nu. A consequence that surprised us in
testing: for a vertical (dip 90) fault the surface displacement field
appears to be independent of nu entirely (verified numerically to all
printed digits for both slip senses on the geometry below), so such a
configuration shows no postseismic surface deformation at all in a
uniform Maxwell medium; a dipping fault, or in the future a layered
medium, is needed for a nontrivial surface transient. Strike slip on
a vertical but finite fault does produce a transient through the
finite-length end effects.

## Program structure

The computation has two phases. The assembly phase evaluates, once,
the per-term amplitude arrays: for every slipping source patch, the
stress amplitudes resolved in the slip direction on every receiver
patch, and the displacement amplitudes at every observation point.
This phase contains all Green's function evaluations. The time loop
that follows is kernel-free: per step it combines the precomputed
amplitudes with the per-term time bases (computed once per step), so
its cost is a small multiple of the output size. On a 500-patch case
the two phases cost about half a second each; the loop share is
dominated by formatted output.

The fault stress is advanced by three routes that are compared
against each other every step: route A, the closed form; route B, the
per-pole state variables h_p integrated with a fixed-step RK4, a
stand-in for the PETSc TS integration that rsf_solve will use; and
route C, the exact-exponential state update that will become the
alternative h mode there. For a step or a ramp resolved by the time
step, route C agrees with the closed form at rounding level and route
B at its truncation order (fourth order in Dt was confirmed by
halving the step). The maximum relative deviations of B and C from A
are written per step, so a regression in the state machinery would be
visible directly in the output.

## Usage

    relax_fault_ve [Dt] [mode] [lfrac] [rfrac] [t_M] [t_ramp] [t_max_fac]

with defaults 0.1, 0 (strike), 0.4, 0.6, 1, 0 (step), 50. The
geometry is read from geom.in. Output files: rf_ve_stress.dat (per
step: time, step index, resolved stress on each patch, the layout of
relax_fault output), rf_ve_disp.dat (per step: time, then per
observation point its x coordinate and ux uy uz), and rf_ve_route.dat
(per step: time and the route B and C deviations). The observation
profile is currently 101 points spanning x from -5 to 5 at the free
surface. A startup banner reports the worst held-out amplitude
residual; values above about 1e-8 would warrant attention, though on
the configurations tested so far they sit at 2e-11 or below, with one
Okada dip-slip pairing at 1.3e-8 that appears to be sample-placement
conditioning rather than a missing pole.

test_relax/run_fault_relax_test runs this program (method 2) next to
relax_fault (method 1, the effective-modulus implicit-Euler scheme of
section 1 of DAM's note)
on a buried vertical fault of 500 patches and plots the patch-time
stress field and surface displacement profiles. The two methods agree
exactly at t = 0 and differ afterwards at leading order, not by
discretization: the implicit-Euler scheme's continuum limit decays
all deviatoric stress as exp(-t/t_M), whereas the exact solution's
late decay time is 3(1-nu)/(1+nu) t_M, a factor 1.8 slower at
nu = 1/4. In the tested configuration the schemes stay within about
15 percent of each other through ten Maxwell times and then diverge;
which level of disagreement matters will depend on the application.

## Verification

Three test programs cover the machinery at increasing depth, and all
of them read geom.in from the working directory.

ve_check (test_relax/run_ve_tests drives it over 3-D Okada
rectangles, genuine w = 0 plane-strain segments, and triangular
elements) asserts at rounding level that the amplitude terms sum to
the elastic kernel at t = 0, that the held-out Laplace samples are
reproduced, that the relaxed stress amplitude vanishes, and that the
displacement constant term matches a Richardson extrapolation of
direct small-s evaluations while the t_M displacement amplitude
vanishes. Observed residuals on the tested geometries are 1e-12 or
smaller.

ve_laplace_check verifies the full time histories, not just
endpoints, against Gaver-Stehfest numerical inversion of the
Laplace-domain solution built from direct real-s evaluations, a route
that makes no assumption about the exponential representation.
Agreement on the tested strike-slip configuration is at the 1e-4
level of the coseismic amplitude, which appears to be the Stehfest
method's own double-precision floor (the deviation is non-monotone in
the Stehfest order, the signature of inverter roundoff rather than of
an error in the representation); the Stehfest-free endpoint anchors
agree at 1e-9 or better. The pass threshold is set accordingly at
3e-4.

The route comparison inside relax_fault_ve itself, described above,
runs on every invocation.

One repair to shared code came out of this work: the GC_DISP_ONLY
branch of eval_triangle_nw.c destroyed the just-computed triangle
displacements by calling set_stress_and_disp_nan, whose mode argument
selects the requested outputs. Nothing else in the code exercised
that path; it now marks only the unused stress output.

## Known limitations and next steps

The exactness of the three-pole representation is established for the
homogeneous shear-Maxwell, elastic-bulk material and the kernel
families tested (Okada rectangles in half and full space, plane
strain segments, Nikkhoo triangles); other rheologies or a layered
medium are outside its scope, which is why the machinery carries the
constant term and a general term count. A uniform Maxwell half space
fully relaxes, so cycle loading by backslip saturates there;
sustained loading with a nonzero relaxed kernel is what the planned
plate-over-substrate extension (step 7 of the implementation plan)
provides, and the antiplane image-series case there has a closed-form
reference. The immediate next step is the rsf_solve integration
(step 3 of the plan): option parsing and sampled operator assembly,
dense and H-matrix, with the option-absent path verified
bit-identical.
