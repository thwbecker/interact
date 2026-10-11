"""Sequence tests (about 1.5 minutes): QD/FD hand-off, switch insensitivity,
FD versus QD on the same fault."""
import numpy as np, pytest
import run_cycle1

YR = 3.15576e7


def first_event(V_sw, dynamic=True):
    model, y0, x = run_cycle1.setup(V_sw=V_sw, dynamic=dynamic)
    model.run(y0, t_end=5 * YR, max_events=1)
    e = model.events[0]
    return e, x


def test_switch_threshold_insensitivity():
    e3, x = first_event(1e-3)
    e2, _ = first_event(1e-2)
    assert abs(e3["Vmax"] - e2["Vmax"]) < 0.05 * e3["Vmax"]
    assert abs(e3["slip"].mean() - e2["slip"].mean()) < 0.02 * e3["slip"].mean()
    assert abs(e3["t_start"] - e2["t_start"]) < 60.0          # same nucleation time to a minute
    assert 2.0 < e3["Vmax"] < 3.0 and 0.9 < e3["slip"].mean() < 1.1


def test_fd_exceeds_qd_in_peak_slip_rate_and_slip():
    efd, x = first_event(1e-3, dynamic=True)
    eqd, _ = first_event(1e-3, dynamic=False)
    assert efd["Vmax"] > 2.0 * eqd["Vmax"]
    assert efd["slip"].mean() > 1.2 * eqd["slip"].mean()
    assert np.sum(efd["slip"] > 1e-2) >= np.sum(eqd["slip"] > 1e-2)
