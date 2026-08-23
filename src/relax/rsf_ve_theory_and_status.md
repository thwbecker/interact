# Visco-elastic relaxation for interact/rsf_solve: theory, plan, and status

TWB / Claude, August 2026. This note documents how the visco-elastic
machinery grew out of DAM's note (ve-note-v4.pdf, copy in src/relax/),
what was adopted from it and what was changed, the staged
implementation plan, and where the work stands at the completion of
the prescribed-slip testbed, i.e. before the rsf_solve integration.
Companion documents: rsf_ve_implementation_plan.md (the full step
list with acceptance criteria) and relax_fault_ve.md (the testbed
program itself).

## What was adopted from DAM's note

The note poses the right problem: a boundary-element cycle code whose
elastic interaction kernels should acquire Maxwell visco-elastic
relaxation without abandoning the precomputed-operator structure, and
it proposes two complementary routes.

Section 1 gives an incremental effective-modulus scheme: over a time
step Dt, the shear modulus is replaced by G' = G alpha with
alpha = t_M/(t_M + Dt), the bulk modulus held fixed, and the stress
updated from the change of the slip field evaluated at the effective
moduli plus an alpha decay of the accumulated deviatoric stress. We
implemented this essentially as written; it is src/relax/relax_fault.c
and remains in the code as the comparison method. Two properties were
established in testing and documented there: the deviatoric step
ratio is exactly alpha on all kernel families, i.e. the scheme does
what its derivation says, and its volumetric part inherits a
Dt dependence from the Green's-function shortcut (the note's K(E),
G(E) helpers with nu confirmed by DAM to mean the Poisson ratio).

Sections 2 to 4 propose the direction the final implementation took:
represent the time-dependent interaction kernel by exponentials
(there called a Prony series) so that the hereditary convolution
integral collapses into a small number of internal variables with
local decay equations, which is the only structure compatible with a
cycle solver. That architectural idea, exponential kernels carried by
per-patch state variables h_p with dh_p/dt = -h_p/tau_p + forcing, is
adopted unchanged and is what everything below implements.

## What was changed relative to the note

Working through sections 2 to 4 surfaced two substantive corrections,
recorded here because the final equations differ from the note's.

First, the object that relaxes is the step response to slip, not to
slip velocity: the stressing rate is the elastic kernel times the
velocity plus the relaxation of accumulated past slip, so the state
forcing must be C_p (V - V_pl), the kernel convolved with slip
history, rather than a kernel applied to velocity as such. Second,
the single-exponential form K_el exp(-t/t_M) assumed in the note is
exact only in the antiplane limit. For the general (plane strain and
3-D) kernels the correspondence principle gives more structure, which
turned out to be small and exact rather than approximate, as follows.

## The exact representation

In the Laplace domain the visco-elastic solution is the elastic
solution at effective moduli, G_bar(s) = G s/(s + 1/t_M) for shear
Maxwell behavior with the bulk modulus K_b held elastic, and
nu_bar(s) = (3 K_b - 2 G_bar)/(6 K_b + 2 G_bar). Because the elastic
kernels are rational in the moduli, the step response is a finite sum
of exponentials whose rates are material properties, independent of
geometry. Inspection of the modulus combinations that occur in the
kernels admits at most three simple poles, at t_M,
3(1-nu)/(1+nu) t_M, and 3/(2(1+nu)) t_M (1.8 t_M and 1.2 t_M at
nu = 1/4), plus a constant (fully relaxed) term. A numerical pole
census against the actual interact kernels (Okada rectangles in half
and full space, plane-strain segments, Nikkhoo triangles) confirmed
this set is complete at rounding level, with the full-space kernels
needing only two of the poles and no repeated-pole terms appearing
anywhere; a cross-check at nu = 0.35 confirmed the rate formulas.
For the homogeneous medium the relaxed stress term vanishes (a
uniform Maxwell half space relaxes completely), but the constant term
is carried throughout because displacement kernels need it (they
relax to the nu = 1/2 field, the permanent deformation, and carry no
t_M pole since G cancels in displacements from dislocations) and the
layered extension needs it (a plate over a substrate does not fully
relax, which is what keeps cycle loading alive there).

Amplitudes are obtained without any new Green's functions: the
existing elastic kernels are evaluated at a small number of sampled
effective moduli and a small weight system maps the samples to the
per-pole amplitudes; two extra samples are held out and must be
reproduced, so an incomplete pole set is detected rather than fitted
over. For the cycle solver this becomes: a handful of startup
assemblies, three amplitude operators C_p plus the elastic operator,
per-patch states h_p with dh_p/dt = -h_p/tau_p + [C_p (V - V_pl)],
and the stressing rate gaining the sink -sum_p h_p/tau_p. The same
machinery with externally supplied rates and amplitudes (and a
nonzero relaxed operator) is the layered case, which is why the
representation carries a general term count and the constant term
from the start.

## The implementation plan, in brief

The staged plan (rsf_ve_implementation_plan.md, agreed August 2026)
is: step 0, frozen reference catalogs and a byte-comparison harness,
since the standing rule is that rsf_solve without the Maxwell option
stays bit-identical at every stage; step 1, the self-contained Prony
core (spec, weights, entrywise amplitude evaluators shaped as kernel
functions so H-matrix backends can consume them later); step 2, a
prescribed-slip testbed exercising the full machinery, including both
planned treatments of the internal variables and surface deformation,
before any cycle run; step 3, rsf_solve plumbing only (-maxwell_time,
sampled assemblies, verification banner, a setup-only mode), with no
dynamics change; step 4, the dynamics with h as additional
integration variables; step 5, the alternative exact-exponential h
update between accepted steps; step 6, H-matrix coverage as a
first-class test; step 7a, the layered antiplane plate over a Maxwell
half space via the Nur-Mavko / Savage-Prescott image series, which
has closed forms to validate against, including a cycle benchmark;
step 7b, plane-strain and 3-D layered kernels from the
Plate-over-Maxwell (Kaj) port, cross-checked by numerical inverse
Laplace transformation.

## Status at the completion of the testbed

Steps 1 and 2 are implemented, adopted into the repository
(src/relax/prony_kernel.c, src/includes/prony_kernel.h,
src/relax/relax_fault_ve.c, src/testing/ve_check.c,
src/testing/ve_laplace_check.c, driven by test_relax/run_ve_tests and
test_relax/run_fault_relax_test), and validated on the configurations
tested so far as follows. The unit checks hold at rounding level on
all three kernel families: amplitude terms sum to the elastic kernel,
held-out samples are reproduced, the relaxed stress amplitude
vanishes, and the displacement constant term matches small-s
extrapolation. The testbed's three routes for the internal variables
(closed form, stepped ODEs, exact-exponential update) agree at their
respective orders, for step and ramp loading. The full time histories
were verified against Gaver-Stehfest numerical inversion of the
Laplace-domain solution, an independent route that assumes nothing
about the exponential representation, to the inverter's own
double-precision floor. The comparison against the section 1 scheme
shows exact agreement at t = 0 and a leading-order difference in the
late-time decay (t_M versus 1.8 t_M e-folding in the tested
configuration), which quantifies what the exact representation adds.
One shared-code repair came out of the work (the GC_DISP_ONLY branch
of the triangle evaluator destroyed its own displacements; nothing
else exercised that path), and one physical null result is worth
remembering when designing tests: a vertical fault in a uniform
Maxwell medium produces no postseismic surface deformation, because
its surface displacement field appears to be nu-independent, so
dipping faults or layering are needed for a surface transient.

The next stage is the rsf_solve integration, beginning with the
step 3 plumbing (operator assembly and verification behind
-maxwell_time, option-absent path bit-identical), followed by the
step 4 dynamics.
