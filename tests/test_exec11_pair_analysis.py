import importlib.util
import pathlib

import numpy as np

MODULE_PATH = pathlib.Path(__file__).resolve().parents[1] / "analysis" / "exec11_pair_analysis.py"
SPEC = importlib.util.spec_from_file_location("exec11_pair_analysis", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
fourth_hit_and_counts = MODULE.fourth_hit_and_counts


def sample_rows():
    event = np.array([1, 0, 0, 0, 0, 2, 2, 2, 3, 3, 3, 3, 3, 3])
    gid = np.array([28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 29])
    time = np.array([9, 4, 1, 3, 2, 1, 2, 3, 8, 2, 4, 6, 10, 11], dtype=float)
    return event, gid, time


def test_fourth_hit_and_fixed_denominator():
    frame = fourth_hit_and_counts(*sample_rows(), n_events=6)
    assert len(frame) == 6
    assert frame.loc[0, "npe_A"] == 4
    assert frame.loc[0, "t_A_ns"] == 4
    assert frame.loc[1, "npe_A"] == 1
    assert frame.loc[2, "npe_A"] == 3
    assert frame.loc[4, "npe_A"] == 0
    assert np.isnan(frame.loc[4, "t_A_ns"])
    assert frame.loc[5, "npe_A"] == 0


def test_signs_and_pair_threshold():
    event = np.array([0] * 8)
    gid = np.array([28] * 4 + [29] * 4)
    time = np.array([1, 2, 3, 4, 2, 3, 4, 6], dtype=float)
    frame = fourth_hit_and_counts(event, gid, time, n_events=2)
    assert frame.loc[0, "passed_pair_4pe"]
    assert frame.loc[0, "delta_t_ps"] == -2000
    assert frame.loc[0, "R_log_ratio"] == 0
    assert not frame.loc[1, "passed_pair_4pe"]


def test_log_ratio_sign():
    event = np.zeros(12, dtype=int)
    gid = np.array([28] * 8 + [29] * 4)
    time = np.arange(12, dtype=float)
    frame = fourth_hit_and_counts(event, gid, time, n_events=1)
    assert np.isclose(frame.loc[0, "R_log_ratio"], np.log(2))
    assert np.isclose(frame.loc[0, "npe_asymmetry"], 1 / 3)


def test_input_order_independence():
    event, gid, time = sample_rows()
    expected = fourth_hit_and_counts(event, gid, time, n_events=6)
    order = np.array([8, 3, 12, 1, 10, 5, 0, 13, 7, 2, 11, 4, 9, 6])
    actual = fourth_hit_and_counts(event[order], gid[order], time[order], n_events=6)
    for column in expected:
        np.testing.assert_equal(expected[column].to_numpy(), actual[column].to_numpy())


def test_rejects_out_of_range_event_ids():
    with np.testing.assert_raises(ValueError):
        fourth_hit_and_counts(np.array([3000]), np.array([28]), np.array([1.0]))
