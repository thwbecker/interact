"""Harris and Day (1993) stepover regression (slow: about 20 s to 2 minutes).

Case D, 100 m elements: the values are those of the current prototype (2026-11-20). Harris and Day
(their Fig. 7d) have fault 1 ending at 5.6 s and fault 2 triggered 1.5 s later at x ~ 14.6 km;
the prototype gave 5.6 s before the fix of the history indexing of 2026-10-18 (commit a4950d9,
shells.System.history_stress: lag l refers to the slice pushed at k + lag0 - 1 - l; one step was
missing before) and gives 5.28 s, a trigger 0.51 s later at 15.45 km now. The qualitative result
of Harris and Day holds: the 0.5 km compressional step jumps, the 1 km step does not.
See fully_dynamic/docs/fd_plan.md (M11 notes of 2026-11-20)."""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pytest
import run_hd


@pytest.mark.slow
def test_case_D_half_km_compressional_jumps_and_one_km_does_not():
    out = run_hd.run(case="D", W=500.0, dx=100.0, tend=10.0, verbose=False)
    # prototype of 2026-11-20: fault 1 ends at 5.283 s, fault 2 triggered at 5.788 s at 15.45 km
    # and ruptures 19 km (the time step is 0.0144 s, the element 100 m)
    assert abs(out["t_end1"] - 5.283) < 0.03
    assert out["t_trig2"] is not None and abs(out["t_trig2"] - 5.788) < 0.03
    assert abs(out["x_trig2"] - 15450.0) < 150.0
    assert out["len2"] > 15.0
    out = run_hd.run(case="D", W=1000.0, dx=100.0, tend=10.0, verbose=False)
    assert out["t_trig2"] is None and out["len2"] == 0.0
