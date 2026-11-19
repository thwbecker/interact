# fully_dynamic: the fully dynamic extension of rsf_solve

The C code is in `src/rsf_solve/` (`fd*.c`, `rsf_fd.c`, option `-fd_mode` of `rsf_solve`); its
user documentation is `../rsf_solve_fd.md` and its tests are in `../src/testing/rsf_solve/`
(`run_tests.sh`). This directory holds everything else.

| folder | contents |
| --- | --- |
| `docs/` | `fd_plan.md` (plan, decisions and results, the running record), `HANDOFF.md` (how to resume), `fd_c_design.md` (design of the C port) |
| `proto/` | the Python prototypes (2-D and 3-D kernels, time marching, half-plane DWN and Cagniard-de Hoop, virtual surfaces), the generators of C code and reference data, and the drivers of the benchmark runs; `proto/README.md` describes them |
| `proto/tests/` | tests of the prototypes; `run_proto_tests.sh [quick\|all]` |
| `subduction2d/` | 2-D subduction sections with rsf_solve: `make_2d_subduction.py` (idealized sections), `coarsen_2d.py` (resampling), `check_2d_resolution.py` (Lambda0 and h* per element), `run_2d_cascadia.sh` (QD and FD runs, summary), `compare_runs_2d.py` (events and phases, QD against FD) |
| `stepover/` | two-segment stepover sequences with rsf_solve: `gen_step.py`, `run_step_fd.sh`, `run_many_stepover`, `summarize_rsf_step.py` |
| `benchmarks/` | SEAS BP1 and BP4 with rsf_solve (`run_bp1_fd.sh`, `run_bp4_qd.sh`), SCEC TPV14 station comparison (`tpv14_stcmp.py`) |
| `tools/` | `regress_rsf.sh` (bit identity and timing of rsf_solve against a reference binary), `SANDBOX.md` (build recipe without H-matrix libraries) |

Scripts are run from the repository root unless their header says otherwise. Paths changed on
2026-11-20 (see `docs/fd_plan.md`, "Layout of 2026-11-20"); older entries of the plan use the
former flat layout.
