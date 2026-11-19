# Benchmarks of the fully dynamic code

Comparisons with community benchmarks and published results, as opposed to the regression tests
of `src/testing/rsf_solve/run_tests.sh` (which compare the C code with the Python prototype and
with exact solutions). Status as recorded in `../docs/fd_plan.md` and `../proto/README.md`; the
audit of 2026-11-20 re-runs them with the current code and the fastest settings that keep the
results within the tolerances below.

| benchmark | reference | code and driver | recorded result | tolerance used |
| --- | --- | --- | --- | --- |
| SCEC TPV205-2D | Dunham's MDSBI 6.25 m on-fault series (untracked, `proto/ref/`) | prototype `proto/run_tpv205.py`, `proto/tests/test_stepper.py`; C `fd_dynamic tpv205` against the prototype | prototype within 0.05 m slip at three stations (100 m) | slip 0.05 m, rupture time 1 step |
| Harris and Day (1993) case D | their Fig. 7d | prototype `proto/run_hd.py`, `proto/tests/test_hd.py`; C `fd_dynamic hd` against the prototype | the 0.5 km compressional step jumps, the 1 km step does not; fault-1 end 5.28 s (theirs 5.6 s), see fd_plan.md | qualitative; timings within 0.03 s of the prototype |
| Andrews (1976b), Dunham (2007) supershear transition | their transition lengths | prototype `proto/run_transition.py`, `proto/run_dunham.py` | reproduced at 38 to 40 elements per cohesive zone | within the reading of their figures |
| SEAS BP1-QD and BP1-FD | SEAS submissions (untracked, `proto/ref/`) | prototype `proto/run_bp1.py`, `proto/plot_bp1.py`; rsf_solve `benchmarks/run_bp1_fd.sh` | QD event times within 0.05 yr of group 4, FD within 1 yr of sem2dpack (50 m); rsf_solve equals the prototype at 25 m | event times 0.1 yr (QD), 1 yr (FD); slip per event 0.05 m |
| SEAS BP4-QD (3-D) | SEAS submissions (not in the repository) | rsf_solve `benchmarks/run_bp4_qd.sh` | QD only | as in the SEAS comparison |
| SCEC TPV3 (3-D, slip weakening) | SCEC reference | prototype `proto/run_tpv3.py`; C `fd3d_test tpv3` against the prototype | C equals the prototype | 1e-8 relative (C against prototype) |
| SCEC TPV101 (3-D, rate and state) | SCEC reference | prototype `proto/run_tpv101.py`; C `fd3d_test tpv101` against the prototype | C equals the prototype | 1e-8 relative (C against prototype) |
| SCEC TPV14/15 (3-D branch, free surface) | SCEC FaultMod 100 m stations (`tpv14_stcmp.py` reads the saved HTML) | C `fd3d_test tpv14`, `benchmarks/tpv14_stcmp.py` | main fault converges (rupture times 0.13 to 0.22 s late, slip within 5 to 10 percent at 250 m); the branch outcome does not converge (fd_plan.md H13) | rupture time 0.25 s, slip 10 percent (main fault) |
| 2-D half-plane statics | Crouch and Starfield (interact's tdd_coeff) | `run_tests.sh` (unit, rsf tiers) | 1.7e-9 (Cagniard-de Hoop), 7.8e-8 (DWN) | 1e-5 |

Untracked reference data (TPV205-2D and BP1 series) live in `fully_dynamic/proto/ref/` on the
machines that have them; checks that need them are skipped elsewhere.
