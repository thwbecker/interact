"""Harris and Day (1993) stepover regression (slow: about 3 minutes).

Case D with 50 m elements, values of the current prototype (2026-11-20), which the C code
(fd_dynamic hd) reproduces. Resolution matters for the trigger: the trigger delay and position
converge from 50 m on (fault 1 end, trigger delay, trigger position):
  100 m: 5.283 s, 0.505 s, 15.45 km (too early: the cohesive zone is under-resolved)
   50 m: 5.261 s, 1.790 s, 14.88 km
   25 m: 5.258 s, 1.786 s, 14.81 km (fd_dynamic only)
Harris and Day (their Fig. 7d): fault 1 ends at 5.6 s, fault 2 triggered 1.5 s later at about
14.6 km. The converged trigger delay and position are within 0.3 s and 0.2 km of theirs; the end of
fault 1 is 0.34 s earlier (not resolved here: the nucleation procedures differ, ours is forced over
2.5 km at 0.5 c_s). At 100 m the prototype matched their 5.6 s only before the fix of the history
indexing of 2026-10-18 (commit a4950d9). See fully_dynamic/docs/fd_plan.md (audit of 2026-11-20)."""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pytest
import run_hd


@pytest.mark.slow
def test_case_D_half_km_compressional_jumps_and_one_km_does_not():
    out = run_hd.run(case="D", W=500.0, dx=50.0, tend=10.0, verbose=False)
    # prototype of 2026-11-20 at 50 m: fault 1 ends at 5.261 s, fault 2 triggered at 7.051 s at
    # 14.875 km and ruptures 18.65 km (time step 0.0072 s)
    assert abs(out["t_end1"] - 5.261) < 0.02
    assert out["t_trig2"] is not None and abs(out["t_trig2"] - 7.051) < 0.02
    assert abs(out["x_trig2"] - 14875.0) < 100.0
    assert out["len2"] > 15.0
    out = run_hd.run(case="D", W=1000.0, dx=50.0, tend=10.0, verbose=False)
    assert out["t_trig2"] is None and out["len2"] == 0.0
