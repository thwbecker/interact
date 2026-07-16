# slider: single-patch benchmark for state evolution laws and ODE solvers

A single 10 x 10 km Okada patch (`geom_slider10.in`, half-lengths in meters)
with BP5 velocity-weakening friction (`rsf_slider.in`: a = 0.004, b = 0.03)
loaded by backslip is a spring slider with the code's own self-stiffness
(K11 about -3e6 Pa/m here, k/k_c about 0.65) and radiation damping.  It
cycles with about 140 yr recurrence, every run takes seconds to minutes, and
solver behavior can be charted without any of the cost or ambiguity of the
full BP5 problem.  `run_slider` executes the matrix of law variants x solvers
x tolerances (resumable: completed runs are skipped; an optional argument
gives a wall-clock budget in seconds); `analyze_slider.py` produces a physics
table from the tight-tolerance 5dp reference runs and work-precision data
(cost versus the RMS event-time error against the reference).  A 20 km patch
(`geom_slider.in`) is also provided (k/k_c about 0.3, larger drops, about
290 yr recurrence).

All numbers below are from one configuration (dense, serial, PETSc 3.19.6,
tmax 450 yr, three events, phase error over matched events); they should be
regenerated rather than quoted for other setups.

## Physics of the laws (reference runs)

| variant | recurrence [yr] | drop [MPa] | peak v [m/s] | duration [s] |
|---------|-----------------|------------|--------------|--------------|
| aging   | 142.3           | 12.7       | 0.95         | 65           |
| slip    | 135.9           | 12.4       | 2.6          | 16           |
| prz     | 136.5           | 14.2       | 3.3          | 6.5          |
| sato    | 142.1           | 14.7       | 3.1          | 9.6          |
| kt      | 143.3           | 14.7       | 3.1          | 15           |

Observations, specific to this configuration: on a single degree of freedom
the laws differ strongly in coseismic character (PRZ events are ten times
shorter and three times faster than aging, the gated laws sit in between with
aging-like healing) but only weakly in recurrence (spread under about 5
percent).  The much larger recurrence differences seen on BP5 (about 235
versus 172 yr at 2 km) are therefore not a zero-dimensional property of the
laws; they emerge from the spatially extended dynamics (partial ruptures,
arrest levels, front behavior).  An alternative PRZ normalization
(d theta/dt = 1 - Omega^2, healing matched to aging at the cost of a doubled
relaxation rate) was also tested here and lengthened the slider recurrence by
2.3 yr, under 2 percent, consistent with the estimate that the healing
normalization is a minor contributor to recurrence differences; the runtime
option for it was subsequently removed as not useful.

## Solvers (work-precision, phase error of matched events)

Explicit 3bs converges cleanly for every law: for aging, RMS event-time error
falls from 1.0e-1 yr at rtol 1e-3 to 2.0e-4 yr at 1e-6 for 364 to 2103 steps;
PRZ costs only about 40 percent more steps than aging at equal tolerance on
this problem, so the single-cell PRZ stiffness is unremarkable; the BP5-scale
PRZ cost explosion is a rupture-front (many-cell) effect.  rtol 1e-3 is
outside the stable envelope for the slip, PRZ, Sato, and KT laws (PETSc
aborts with DIVERGED_STEP_REJECTED at the first event).

The current plain IMEX (`-imex`, Newton line-search stage solves) is
dominated by explicit 3bs on this problem everywhere it was tested: it costs
more per step (implicit function evaluations two to three times the RHS
count), its error stagnates near 1e-2 yr for aging instead of improving from
rtol 1e-5 to 1e-6, it loses events or wedges at loose tolerances (the aging
runs at rtol 1e-3 and 1e-4 stall at the first nucleation in a domain-check
rejection storm, about 2e5 rejected attempts per 1e2 accepted steps, which
appears not to count toward any PETSc rejection limit; the step ceiling in
run_slider now terminates such runs within seconds and the analysis flags
them), and the PRZ IMEX runs timed out mid-series at all tolerances tried.
This reproduces, on the simplest possible system and for the AGING law, the
stage-problem pathology diagnosed on BP5: without the cell's elastic
self-stiffness in the implicit part, the per-cell stage problem loses its
root at large dt (a saddle-node) and the nonlinear solve fails, capping the
step.  The slider therefore makes a suitable development target for the
stage-solver work: a correct treatment should first make this matrix clean
(IMEX at least matching explicit accuracy law by law) before returning to
BP5.

## Files

    geom_slider10.in   the 10 km patch (default)
    geom_slider.in     a 20 km variant (lower k/k_c, longer cycles)
    rsf_slider.in      a, b
    run_slider         the matrix driver (parameters as plain assignments
                       at the top; resumable; optional wall-clock budget
                       argument)
    analyze_slider.py  physics table and work-precision extraction; writes
                       slider_physics.dat and slider_wp.dat
