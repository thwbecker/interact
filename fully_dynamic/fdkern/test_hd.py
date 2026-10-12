"""Harris and Day (1993) stepover regression (slow: about 2 minutes)."""
import numpy as np, pytest
import run_hd


@pytest.mark.slow
def test_case_D_half_km_compressional_jumps_and_one_km_does_not():
    out = run_hd.run(case="D", W=500.0, dx=100.0, tend=10.0, verbose=False)
    # Harris and Day Fig. 7d: fault 1 end at 5.6 s, trigger 1.5 s later at x ~ 14.6 km, rupture propagates
    assert abs(out["t_end1"] - 5.6) < 0.3
    assert out["t_trig2"] is not None and 1.0 < out["t_trig2"] - out["t_end1"] < 2.2
    assert 14.3e3 < out["x_trig2"] < 15.3e3
    assert out["len2"] > 10.0
    out = run_hd.run(case="D", W=1000.0, dx=100.0, tend=10.0, verbose=False)
    assert out["t_trig2"] is None and out["len2"] == 0.0
