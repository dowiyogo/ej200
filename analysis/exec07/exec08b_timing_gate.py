#!/usr/bin/env python3
"""Compare intrinsic same-end SUM4 timing in EndTop and End-only campaigns."""

from __future__ import annotations

import argparse
import math
import pathlib
import sys

import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from exec07.common import END_CLUSTERS  # noqa: E402
from exec07_photon_budget import leading_edge_time  # noqa: E402

POSITIONS = (-400, 0, 400)
ENDTOP_FILES = {
    x: f"photon_hits_x{x}mm.root" for x in POSITIONS
}
END_ONLY_FILES = {
    -400: "photon_hits_run007.root",
    0: "photon_hits_run015.root",
    400: "photon_hits_run023.root",
}


def cluster_times(
    event_id: np.ndarray, global_id: np.ndarray, time_ns: np.ndarray,
    ids: tuple[int, ...], n_events: int,
) -> np.ndarray:
    mask = np.isin(global_id, ids)
    events = event_id[mask]
    times = time_ns[mask]
    order = np.lexsort((times, events))
    events = events[order]
    times = times[order]
    starts = np.searchsorted(events, np.arange(n_events), side="left")
    stops = np.searchsorted(events, np.arange(n_events), side="right")
    result = np.full(n_events, np.nan)
    for event in range(n_events):
        result[event] = leading_edge_time(times[starts[event]:stops[event]])
    return result


def analyze(path: pathlib.Path, campaign: str, expected_x: int) -> list[dict]:
    with uproot.open(path) as root_file:
        arrays = root_file["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns", "gun_x_mm"], library="np",
        )
    event_id = arrays["event_id"].astype(int)
    global_id = arrays["global_id"].astype(int)
    time_ns = arrays["time_ns"].astype(float)
    unique_events = np.unique(event_id)
    unique_x = np.unique(arrays["gun_x_mm"])
    if unique_x.tolist() != [float(expected_x)]:
        raise RuntimeError(f"{path}: gun_x={unique_x}, expected {expected_x}")
    n_events = int(unique_events.max()) + 1
    end_mask = global_id < 16
    counts = np.bincount(
        event_id[end_mask].astype(np.int64) * 16 + global_id[end_mask],
        minlength=n_events * 16,
    ).reshape(n_events, 16)
    triggers = {
        name: cluster_times(event_id, global_id, time_ns, ids, n_events)
        for name, ids in END_CLUSTERS.items()
    }
    rows = []
    for side, ids in (("left", tuple(range(8))), ("right", tuple(range(8, 16)))):
        a = triggers[f"end_{side}_A_SUM4"]
        b = triggers[f"end_{side}_B_SUM4"]
        valid = np.isfinite(a) & np.isfinite(b)
        delta = a[valid] - b[valid]
        event_npe = counts[:, ids].sum(axis=1)
        rows.append({
            "campaign": campaign,
            "x_beam_mm": expected_x,
            "side": side,
            "n_events": n_events,
            "n_triggered": int(np.count_nonzero(valid)),
            "npe_mean": float(np.mean(event_npe)),
            "npe_sem": float(np.std(event_npe) / math.sqrt(n_events)),
            "sigma_group_ps": float(np.std(delta) * 1000.0 / math.sqrt(2.0)),
            "sigma_group_error_ps": float(
                np.std(delta) * 1000.0 / math.sqrt(2.0)
                / math.sqrt(2.0 * (len(delta) - 1))
            ),
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endtop-dir", required=True, type=pathlib.Path)
    parser.add_argument("--end-only-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("analysis/exec07/exec08b_timing_gate.csv"))
    args = parser.parse_args()
    rows = []
    for x_mm in POSITIONS:
        print(f"EndTop x={x_mm}", flush=True)
        rows.extend(analyze(args.endtop_dir / ENDTOP_FILES[x_mm], "EndTop", x_mm))
        print(f"End-only x={x_mm}", flush=True)
        rows.extend(analyze(args.end_only_dir / END_ONLY_FILES[x_mm], "End-only", x_mm))
    raw = pd.DataFrame(rows)
    comparisons = []
    for x_mm in POSITIONS:
        for side in ("left", "right"):
            subset = raw[(raw.x_beam_mm == x_mm) & (raw.side == side)].set_index("campaign")
            endtop = subset.loc["EndTop"]
            end_only = subset.loc["End-only"]
            observed = endtop.sigma_group_ps / end_only.sigma_group_ps
            observed_error = observed * math.hypot(
                endtop.sigma_group_error_ps / endtop.sigma_group_ps,
                end_only.sigma_group_error_ps / end_only.sigma_group_ps,
            )
            predicted = math.sqrt(end_only.npe_mean / endtop.npe_mean)
            predicted_error = 0.5 * predicted * math.hypot(
                end_only.npe_sem / end_only.npe_mean,
                endtop.npe_sem / endtop.npe_mean,
            )
            difference_error = math.hypot(observed_error, predicted_error)
            comparisons.append({
                "x_beam_mm": x_mm,
                "side": side,
                "npe_endtop": endtop.npe_mean,
                "npe_endtop_sem": endtop.npe_sem,
                "npe_end_only": end_only.npe_mean,
                "npe_end_only_sem": end_only.npe_sem,
                "sigma_endtop_ps": endtop.sigma_group_ps,
                "sigma_endtop_error_ps": endtop.sigma_group_error_ps,
                "sigma_end_only_ps": end_only.sigma_group_ps,
                "sigma_end_only_error_ps": end_only.sigma_group_error_ps,
                "observed_sigma_ratio": observed,
                "observed_sigma_ratio_error": observed_error,
                "predicted_sqrt_npe_ratio": predicted,
                "predicted_sqrt_npe_ratio_error": predicted_error,
                "ratio_difference_sigma": (observed - predicted) / difference_error,
                "endtop_not_better": endtop.sigma_group_ps >= end_only.sigma_group_ps,
            })
    output = pd.DataFrame(comparisons)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    raw.to_csv(args.output.with_name("exec08b_timing_gate_raw.csv"), index=False)
    print(output.to_string(index=False))
    if not output["endtop_not_better"].all():
        raise SystemExit("TIMING GATE FAILED: EndTop is better than End-only in at least one comparison")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
