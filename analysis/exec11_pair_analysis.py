#!/usr/bin/env python3
"""EXEC_11 pair-scan analysis from the authoritative per-hit sipm_hits tree."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import subprocess
from collections import Counter

import numpy as np
import pandas as pd
import uproot

HOOK_DATA_DIR = pathlib.Path("/home/reriosto/SHiP/t0minidaq/pairscan_2026-06-11")
HOOK_PAIR_IDS = (28, 29)
HOOK_N_EVENTS = 3000
HOOK_POS_REF_1 = "AUTO"
HOOK_POS_REF_2 = "AUTO"
HOOK_FIT_RANGE_NSIGMA = 2.5
HOOK_CALIB_EXCLUDE_TEST = True
HOOK_SPTR_PS = 0
HOOK_ELEC_JITTER_PS = 0
HOOK_WALK = "none"

N_PE_THRESH = 4
EXPECTED_POSITIONS = tuple(float(x) for x in range(-462, -421))
REQUIRED_BRANCHES = {"event_id", "global_id", "time_ns", "gun_x_mm"}
COUNT_DEFINITION = (
    "number of sipm_hits rows per event and SiPM, using the existing EXEC_08 "
    "PE-equivalent convention"
)
FILE_RE = re.compile(r"pairscan_x([+-]?\d+(?:\.\d+)?)mm\.root$")


def analysis_commit() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], text=True
    ).strip()


def extract_position(path: pathlib.Path) -> float:
    match = FILE_RE.fullmatch(path.name)
    if not match:
        raise ValueError(f"Unexpected pair-scan filename: {path.name}")
    return float(match.group(1))


def fourth_hit_and_counts(
    event_ids: np.ndarray,
    global_ids: np.ndarray,
    times_ns: np.ndarray,
    n_events: int = HOOK_N_EVENTS,
    pair_ids: tuple[int, int] = HOOK_PAIR_IDS,
) -> pd.DataFrame:
    """Build fixed-denominator event observables independent of input row order."""
    event_ids = np.asarray(event_ids, dtype=np.int64)
    global_ids = np.asarray(global_ids, dtype=np.int64)
    times_ns = np.asarray(times_ns, dtype=float)
    if not (len(event_ids) == len(global_ids) == len(times_ns)):
        raise ValueError("Input arrays must have equal length")
    if np.any((event_ids < 0) | (event_ids >= n_events)):
        bad = np.unique(event_ids[(event_ids < 0) | (event_ids >= n_events)])
        raise ValueError(f"event_id outside [0,{n_events - 1}]: {bad[:10]}")

    out: dict[str, np.ndarray] = {"event_id": np.arange(n_events, dtype=np.int32)}
    for label, channel in zip(("A", "B"), pair_ids):
        selected = (global_ids == channel) & np.isfinite(times_ns)
        ev = event_ids[selected]
        ts = times_ns[selected]
        counts = np.bincount(ev, minlength=n_events).astype(np.int32)
        fourth = np.full(n_events, np.nan)
        if len(ev):
            order = np.lexsort((ts, ev))
            sorted_ev = ev[order]
            sorted_ts = ts[order]
            starts = np.r_[0, np.cumsum(counts[:-1])]
            eligible = counts >= N_PE_THRESH
            fourth[eligible] = sorted_ts[starts[eligible] + N_PE_THRESH - 1]
        out[f"npe_{label}"] = counts
        out[f"t_{label}_ns"] = fourth
        out[f"passed_{label}_4pe"] = counts >= N_PE_THRESH

    frame = pd.DataFrame(out)
    frame["passed_pair_4pe"] = frame["passed_A_4pe"] & frame["passed_B_4pe"]
    pair = frame["passed_pair_4pe"].to_numpy()
    frame["delta_t_ps"] = np.where(
        pair, 1000.0 * (frame["t_A_ns"] - frame["t_B_ns"]), np.nan
    )
    positive = (frame["npe_A"] > 0) & (frame["npe_B"] > 0)
    frame["R_log_ratio"] = np.where(
        positive, np.log(frame["npe_A"] / frame["npe_B"]), np.nan
    )
    denominator = frame["npe_A"] + frame["npe_B"]
    frame["npe_asymmetry"] = np.where(
        denominator > 0, (frame["npe_A"] - frame["npe_B"]) / denominator, np.nan
    )
    return frame


def inventory_and_derive(path: pathlib.Path, derived_dir: pathlib.Path) -> dict[str, object]:
    x_true = extract_position(path)
    present = np.zeros(HOOK_N_EVENTS, dtype=bool)
    pair_event: list[np.ndarray] = []
    pair_gid: list[np.ndarray] = []
    pair_time: list[np.ndarray] = []
    gun_counter: Counter[float] = Counter()
    n_entries = 0
    n_nonfinite = 0
    hit_counts = {HOOK_PAIR_IDS[0]: 0, HOOK_PAIR_IDS[1]: 0}

    with uproot.open(path) as root_file:
        if "sipm_hits" not in root_file:
            raise RuntimeError(f"{path}: missing sipm_hits")
        tree = root_file["sipm_hits"]
        missing = REQUIRED_BRANCHES - set(tree.keys())
        if missing:
            raise RuntimeError(f"{path}: missing branches {sorted(missing)}")
        for arrays in tree.iterate(
            ["event_id", "global_id", "time_ns", "gun_x_mm"],
            step_size="100 MB",
            library="np",
        ):
            ev = arrays["event_id"].astype(np.int64, copy=False)
            gid = arrays["global_id"].astype(np.int64, copy=False)
            times = arrays["time_ns"].astype(float, copy=False)
            if np.any((ev < 0) | (ev >= HOOK_N_EVENTS)):
                bad = np.unique(ev[(ev < 0) | (ev >= HOOK_N_EVENTS)])
                raise RuntimeError(f"{path}: event_id outside expected range: {bad[:10]}")
            present[np.unique(ev)] = True
            n_entries += len(ev)
            n_nonfinite += int(np.count_nonzero(~np.isfinite(times)))
            values, counts = np.unique(np.round(arrays["gun_x_mm"], 6), return_counts=True)
            gun_counter.update(dict(zip(values.tolist(), counts.tolist())))
            selected = np.isin(gid, HOOK_PAIR_IDS)
            if np.any(selected):
                pair_event.append(ev[selected].copy())
                pair_gid.append(gid[selected].copy())
                pair_time.append(times[selected].copy())
                for channel in HOOK_PAIR_IDS:
                    hit_counts[channel] += int(np.count_nonzero(gid[selected] == channel))

    event_ids = np.concatenate(pair_event)
    global_ids = np.concatenate(pair_gid)
    times_ns = np.concatenate(pair_time)
    frame = fourth_hit_and_counts(event_ids, global_ids, times_ns)
    frame.insert(0, "x_true_mm", x_true)
    derived_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        derived_dir / f"pair_events_x{x_true:+.1f}mm.npz",
        **{column: frame[column].to_numpy() for column in frame.columns},
    )
    dominant_gun, dominant_count = gun_counter.most_common(1)[0]
    return {
        "file": path.name,
        "x_true_mm": x_true,
        "size_bytes": path.stat().st_size,
        "mtime_ns": path.stat().st_mtime_ns,
        "n_entries": n_entries,
        "n_event_ids_present_anywhere": int(np.count_nonzero(present)),
        "event_id_min": int(np.flatnonzero(present)[0]),
        "event_id_max": int(np.flatnonzero(present)[-1]),
        "n_events_missing_from_entire_tree": int(HOOK_N_EVENTS - np.count_nonzero(present)),
        "gun_x_mm_dominant": dominant_gun,
        "gun_x_mm_dominant_fraction": dominant_count / n_entries,
        "hits_A": hit_counts[HOOK_PAIR_IDS[0]],
        "hits_B": hit_counts[HOOK_PAIR_IDS[1]],
        "nonfinite_time_fraction": n_nonfinite / n_entries,
        "n_events_with_A": int(np.count_nonzero(frame["npe_A"])),
        "n_events_with_B": int(np.count_nonzero(frame["npe_B"])),
        "n_events_with_A_and_B": int(np.count_nonzero((frame["npe_A"] > 0) & (frame["npe_B"] > 0))),
        "n_events_passing_4PE_both": int(frame["passed_pair_4pe"].sum()),
        "efficiency": float(frame["passed_pair_4pe"].sum() / HOOK_N_EVENTS),
    }


def validate_file_set(data_dir: pathlib.Path) -> list[pathlib.Path]:
    files = sorted(data_dir.glob("pairscan_x*.root"), key=extract_position)
    positions = [extract_position(path) for path in files]
    if positions != list(EXPECTED_POSITIONS):
        missing = sorted(set(EXPECTED_POSITIONS) - set(positions))
        duplicates = sorted({x for x in positions if positions.count(x) > 1})
        raise RuntimeError(
            f"Expected 41 integer positions -462..-422; missing={missing}, duplicates={duplicates}"
        )
    return files


def write_metadata(output_dir: pathlib.Path) -> None:
    metadata = {
        "data_dir": str(HOOK_DATA_DIR),
        "pair_ids": HOOK_PAIR_IDS,
        "n_events_expected": HOOK_N_EVENTS,
        "n_pe_threshold": N_PE_THRESH,
        "count_definition": COUNT_DEFINITION,
        "data_commit": "f431c01",
        "analysis_commit_at_generation": analysis_commit(),
        "root_phase_0b": "DEFERRED — ROOT runtime unavailable on MSI",
    }
    (output_dir / "analysis" / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


def derive_all(data_dir: pathlib.Path, output_dir: pathlib.Path) -> None:
    files = validate_file_set(data_dir)
    for subdir in ("analysis", "figures", "tables", "derived", "logs"):
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)
    rows = []
    for index, path in enumerate(files, 1):
        print(f"[{index:02d}/{len(files)}] inventory + derive {path.name}", flush=True)
        rows.append(inventory_and_derive(path, output_dir / "derived"))
    inventory = pd.DataFrame(rows).sort_values("x_true_mm")
    inventory.to_csv(output_dir / "analysis" / "data_inventory.csv", index=False)
    if (
        (inventory["n_event_ids_present_anywhere"] != HOOK_N_EVENTS).any()
        or (inventory["nonfinite_time_fraction"] != 0).any()
        or not np.allclose(inventory["gun_x_mm_dominant"], inventory["x_true_mm"])
    ):
        raise RuntimeError("Inventory gate failed; inspect data_inventory.csv")
    write_metadata(output_dir)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("derive",))
    parser.add_argument("--data-dir", type=pathlib.Path, default=HOOK_DATA_DIR)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if args.stage == "derive":
        derive_all(args.data_dir, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
