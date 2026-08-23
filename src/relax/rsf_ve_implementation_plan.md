# Visco-elastic machinery for interact/rsf_solve: step-by-step implementation plan

TWB / Claude, August 2026. Against master at 5be436e. Revised to split step 7 into 7a/7b as agreed. Companion to
rsf_ve_design.md, which holds the theory; this document is the build
order. Design constraints, agreed: (a) the machinery must extend to an
elastic plate over a viscous half space (layered case) without
rework, so the Prony structure carries a relaxed (constant) kernel
term and externally supplied amplitudes from the start; (b) both
treatments of the internal variables h are provided, as additional
integration variables under the TS error controller and as an
approximate exact-exponential update between accepted steps; (c) the
machinery is exercised first in a relax_fault-style prescribed-slip
program, with stress relaxation on the fault and surface deformation
tracked, before any cycle run; (d) rsf_solve without the Maxwell
option stays bit-identical at every stage; (e) H-matrix operator
compression works for all new operators.

One theoretical addition relative to rsf_ve_design.md, needed for the
surface-deformation requirement and convenient for the layered case:
displacement kernels do not fully relax. Displacements from a
dislocation in a homogeneous medium depend on the moduli only through
nu, so under the correspondence substitution their step response is

    U_ij(t) = U^relaxed_ij + sum_{p in {2,3}} D^p_ij exp(-t/tau_p),

with U^relaxed the displacement kernel evaluated at the effective
incompressible limit nu_bar(0) = 1/2 (permanent deformation) and no
tau_m term, since G cancels. This is exactly the constant-plus-Prony
form the layered stress kernels will need, so the fitting and state
machinery is written once with a constant term that happens to vanish
for homogeneous stress kernels. The stress-kernel pole census (three
simple poles for half space, two for full space, no repeated poles)
stands as measured; the same census procedure should be rerun for the
displacement kernels as part of step 2 rather than assumed.

The steps. Each ends with stated acceptance criteria; the bit-identity
and H-matrix checks recur rather than being saved for the end.

Step 0, baseline and harness (no source changes). Freeze reference
catalogs at 5be436e: the existing test_hmatrix two-fault case run
dense and H-matrix, np 1 and np 4, with rsf_monitor.dat and
rsf_events.dat archived; the test_twod analytic case; a relax_fault
run. Add a small comparison script under a new test_ve/ directory that
byte-compares monitor and event files, since this check will run after
every subsequent step. Commit the pole-census sampler (the
kernel_scan program from the analysis phase) into test_ve/ as the
standing verification tool. Acceptance: catalogs reproduce themselves.

Step 1, the Prony core library, self-contained. New files, suggested
src/ve/prony_kernel.c plus a header, defining: a prony_spec structure
holding n_p, the rates 1/tau_p, an optional constant-term flag, and
the sampling weights; a routine that, given (nu, tau_m), fills the
homogeneous rates from the closed formulas and solves the small
Vandermonde-type system for the weights W_pk over n_p + 2 sample
points s_k (three plus two held out for stress, two plus two for
displacement, constant term included as the s -> 0 basis column); and
the central entrywise evaluator: given source patch, receiver point or
patch, component, and the spec, evaluate the existing elastic Green's
function at the n_s sampled effective moduli (calc_medium_elastic_
parameters at G_bar(s_k), nu_bar(s_k)) and return the amplitude set
C^p_ij (and the held-out residual on request). This evaluator is
deliberately shaped like a kernel function of (i, j), because that is
what both dense assembly and the H-matrix backends consume. The
layered hook lives here too: the spec can alternatively be filled from
a file (rates, and either precomputed amplitude matrices or a request
to Prony-fit externally supplied step-response samples), with the
homogeneous path being just one generator among two. A standalone
unit-test program (test_ve/ve_check.c) builds dense C^p for a small
mixed geometry and asserts, to rounding, that sum_p C^p plus the
constant term equals the elastic kernel, that the held-out residual is
at rounding level, and that the reconstructed R(t) matches direct
effective-moduli evaluations at arbitrary t via the inverse-transform
identity. Acceptance: ve_check passes for 2-D segments, Okada
rectangles, and triangles, half and full space; no other binary
changes; step 0 catalogs unchanged (nothing links the new code yet
except ve_check).

Step 2, the prescribed-slip testbed with bulk relaxation and surface
deformation. Extend the relax_fault family, preferably as a new
src/relax/relax_fault_ve.c beside the existing tool so the section 1
scheme remains available for comparison. Inputs as now (slip from
fault[].u, held after t = 0), plus tau_m and an observation
specification, either a surface grid in the style of interact's
bc.in output plane or a file of observation points. The program
computes, using the step 1 library: fault tractions
tau_i(t) = sum_p exp(-t/tau_p) [C^p u]_i via three routes that must
agree, the closed form, the h-as-ODE stepping (a small fixed-step
integrator standing in for the TS), and the exact-exponential
between-step update, which exercises both h modes outside the cycle
solver as requested; and surface displacements
u_obs(t) = [U^relaxed u]_obs + sum_p exp(-t/tau_p) [D^p u]_obs on the
observation set, written as time series. Validation within this step:
at t = 0 stress and displacement match a one-step elastic interact
evaluation on the same geometry; as t grows, deviatoric stress
relaxes to zero and surface displacement converges to the elastic
evaluation at nu = 1/2, both to rounding; the three h routes agree to
their respective truncation orders; and a side-by-side with the
existing relax_fault documents where the section 1 scheme's
volumetric shortcut departs. This step is where the machinery is
declared trusted or not, before rsf_solve is touched. Acceptance: the
above checks pass on the 2-D 500-patch case and one 3-D Okada case;
run_fault_relax_test-style scripts and expected outputs committed
under test_ve/.

Step 3, rsf_solve plumbing without dynamics. Add -maxwell_time (and
the parallel -ve_prony_file for the layered path) to the option
parsing in rsf_init.c; when present, build the prony_spec and
assemble the sampled operators next to the existing
calc_petsc_Isn_matrices calls in rsf_solve.c, for the shear operator
and, when -calc_sigma_dot is active, the normal operator, in both the
dense and use_hmatrix branches; solve for the amplitude operators;
run the held-out check on a random probe vector (dense: entrywise on
a sample of entries as well) and print a banner with the rates and
residuals, aborting on a structured residual. For the H-matrix
branch, the amplitude operators are compressed by handing the
backends the step 1 entrywise evaluator wrapped per pole, so each C^p
is compressed exactly like the elastic operator is today; the
matrix-free alternative, applying C^p as the W-weighted sum of the
n_s compressed sampled operators, should also be wired since it
reuses the existing compression path unchanged and gives a cross-check
on the compressed-C^p route. No RHS change yet; with the option given,
the code may run the elastic dynamics after printing the VE banner, or
stop under an explicit -ve_setup_only. Acceptance: option-absent runs
byte-identical to the step 0 catalogs, dense and H-matrix, np 1 and 4;
with the option, dense and H-matrix amplitude operators agree in
matvec to the compression tolerance on random vectors.

Step 4, rsf_solve dynamics, in-state h. Grow the per-patch block from
its current layout by n_p entries (plus n_p for sigma when coupled),
gated on the option so the option-absent layout is untouched; extend
rsf_ODE_RHSFunction in rsf_engine.c with the sink terms
-sum_p h^p_i/tau_p in the traction rate and the h equations
dh^p/dt = -h^p/tau_p + [C^p (V - Vpl)]_i, reusing the already-computed
elastic matvec where sum_p C^p = K^elastic allows; place the linear
decay of h into the implicit side of the rsf_imex.c split, which
keeps the block-diagonal stage solver structure since the decay is
diagonal; extend the checkpoint header with a version and the block
size so restarts refuse a mismatched layout, and carry h through
restart; optionally add per-pole norms to the monitor output.
Acceptance, in order: option-absent bit-identity re-run; slider
(single patch) with -maxwell_time against an independent integration
of the corresponding small ODE system to integrator tolerance;
-maxwell_time set to an enormous value against the elastic catalog at
tolerance level; a locked-fault (velocity pinned near zero) rsf_solve
run against the step 2 closed form; and a checkpoint-split run
reproducing an uninterrupted one, as in the earlier restart tests.

Step 5, the approximate h mode. Add -ve_h_mode {state, step} with
state the default from step 4; step keeps h outside the TS vector at
the option-absent block layout and applies the exact-exponential
update with the accumulated slip increment between accepted steps,
which is the memory-light, checkpoint-light variant and the natural
mode if the enlarged state ever bothers the error controller.
Acceptance: slider and the two-fault case agree between the two modes
within a tolerance consistent with the step mode's
piecewise-constant-velocity assumption, with the difference shrinking
under tightened rtol; documentation of when each mode is preferable.

Step 6, H-matrix coverage as a first-class test. Extend the
test_hmatrix sweep to include -maxwell_time runs: dense versus
H-matrix cycle comparisons for both h modes, at the tolerances the
existing sweep uses, plus operator-level matvec comparisons for each
C^p. Any backend that cannot compress the Prony-wrapped kernel
cleanly is identified here rather than in production. Acceptance:
dense/H-matrix cycle differences for VE runs are commensurate with the
elastic dense/H-matrix differences at the same compression tolerance.

Step 7, the layered path, in two sub-steps ordered by how much of
the answer is exact.

Step 7a, the antiplane elastic plate over a Maxwell half space. This
configuration has a closed-form time-dependent solution as an image
series (the Nur-Mavko / Savage-Prescott construction): the image at
depth 2 n H turns on with a degree-(n-1) polynomial in t/tau_m times
exp(-t/tau_m), so the step-response kernel is K_relaxed plus a
geometrically convergent series of t^k exp(-t/tau_m) terms,
truncatable at known accuracy for a given source-receiver distance
over the plate thickness H. Implement this as the first per-entry
generator behind the -ve_prony_file interface, validate it internally
against the elastic limit at t = 0, the relaxed limit, and the known
screw-dislocation-in-a-plate fields, and Prony-fit it per entry with
the fit residual reported; the required Prony term count versus
target accuracy is measured here rather than assumed, since this is
the one layered configuration with an exact reference. The
Savage-Prescott cycle solution for the same geometry then provides a
literature-grade end-to-end benchmark for a layered rsf_solve run
(or its spring-slider reduction), a stronger acceptance test than
anything the homogeneous case offers. Testbed-first as elsewhere: a
layered relax_fault_ve run demonstrating surface velocity migration
and stress re-loading, both absent in uniform Maxwell, precedes any
layered cycle run. Acceptance: image-series generator matches the
closed forms at rounding level, Prony residual versus term count
documented, Savage-Prescott comparison at a stated tolerance,
option-absent bit-identity re-run.

Step 7b, plane strain and 3-D layered kernels, gated on 7a. Here
exactness is not promised; amplitudes come from Prony-fitting step
responses produced by two independent routes that must agree: the
Plate-over-Maxwell (Kaj) port as the primary generator, and numerical
inverse Laplace transformation (a Talbot-type contour) of the layered
elastic solution assembled with s-dependent substrate moduli as the
cross-check, with care near the branch points the layering
introduces. Agreement between the routes is the validation;
disagreement localizes the defect to one side. The file format from
step 3 (-ve_prony_file: n_p, rates, relaxed-kernel flag, amplitudes
as precomputed matrices for small problems or, preferably, through
the per-entry generator interface so the H-matrix backends compress
layered amplitude operators the same way) is exercised end to end.
Validate surface deformation for a strike-slip source against the
MATLAB reference and a locked-fault loading curve against direct
evaluation of the layered hereditary integral; the nonzero relaxed
kernel is what keeps backslip loading alive, so a long two-fault
cycle here is the first physically complete layered cycle run. The
Prony term count for target accuracy is an empirical outcome and may
exceed the antiplane case, particularly for shallow sources near the
interface. Acceptance: generator/ILT agreement at a stated
tolerance, fit residuals at a stated target,
surface-deformation match to the reference, and, once more,
option-absent bit-identity.

Standing rules across all steps: every commit re-runs the step 0
byte-comparison in dense and H-matrix modes; new code lives in
src/ve/ and test_ve/ so the elastic paths are touched only at the
explicitly listed integration points (option parsing, operator
assembly call sites, RHS, IMEX split, checkpoint header); and no step
proceeds while the previous step's acceptance list has an open item.
The natural first deliverable is steps 1 and 2 together, since step 2
is the requested prescribed-slip demonstration with surface
deformation and settles trust in the machinery before rsf_solve is
touched at all.
