# Two-state-variable rate-and-state law in rsf_solve: spring-slider check

`rsf_solve -rsf_two_state 1 -rsf_b2 <b2> -rsf_dc2 <dc2>` adds a second
state variable (Ruina 1983; Gu et al. 1984; Becker 2000). The friction
is psi = f0 + b ln(v0 theta1/dc) + b2 ln(v0 theta2/dc2), and both state
variables follow `-state_law`.

`slider2sv.py` checks the implementation on a single Okada patch (the
10 km patch of `../geom_slider10.in`), which is a spring slider with the
code's self-stiffness. The parameters are mapped onto the
nondimensional system of Gu et al. (1984) as used by Becker (2000) and
in `../../ode_2sv_rsf` (beta1 = 1, beta2 = 0.84, rho = dc1/dc2 = 0.048):
a = b = 0.01, b2 = 0.84 a, sigma 50 MPa, v0 = vpl = 1e-9 m/s, slip law,
dc chosen so that K dc/(a sigma) = kappa' kappa_c. The rsf_solve run
uses the default radiation damping, 5dp and rtol 1e-8. The comparison
is against the nondimensional ODE integrated with scipy from the same
initial velocity (1.05 v0), over 2e4 units of dc/v0, analysis in the
second half.

Results of the check (2026-11-08, PETSc 3.19.6):

| kappa' | rsf_solve distinct extrema, y range | ODE distinct extrema, y range |
|---|---|---|
| 0.95 | 2, -0.653 to 0.162 | 2, -0.652 to 0.161 |
| 0.90 | 2, -1.035 to 0.103 | 2, -1.035 to 0.102 |
| 0.86 | 4, -1.584 to 0.117 | 4, -1.584 to 0.118 |
| 0.8552 | 15, -1.813 to 0.140 | 16, -1.813 to 0.139 |

(y = (tau - tau0)/(a sigma); the 15 against 16 is two extrema closer
than the 2e-3 tolerance used to count distinct values.) At kappa' = 0.85
both leave the small-amplitude attractor from this initial perturbation:
the ODE overflows (no radiation damping), and rsf_solve settles into
large stick-slip events (y down to -17, peak slip rates of about 1 m/s),
also from a 1.001 v0 start. At kappa' = 0.8535, with a 1.001 v0 start
and 6e3 units of settling, rsf_solve stayed bounded with many distinct
extrema; longer runs would be needed to separate period 16 or 32 from
chaos there.

`bifurcation_scan.py rsf_solve outfile [vpert] [ttotal] [tsettle]`
scans kappa' from 0.95 to 0.85 (1e5 dc/v0, analysis of the last 5e4,
start at 1.001 v0); p1_extrema.dat and p1_bifurcation.png are the result
of that scan, summarized in ../../strike_slip_networks/single_fault_pilot/README.md.
