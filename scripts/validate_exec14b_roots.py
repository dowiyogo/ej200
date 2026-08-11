#!/usr/bin/env python3
"""Validate the 31 main EJ-230 ROOT files and any completed special controls."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path

import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[1]
POSITIONS = (-690, -670, -650, -600, -550, -500, -450, -400, -350, -300,
             -250, -200, -150, -100, -50, 0, 50, 100, 150, 200, 250, 300,
             350, 400, 450, 500, 550, 600, 650, 670, 690)
WINDOW_META = {
    "photon_hits_run_A_x-652mm.root": (-652, 808201, 808211),
    "photon_hits_run_B_x-642mm.root": (-642, 808301, 808307),
    "photon_hits_run_C1_x-648mm.root": (-648, 808401, 808421),
    "photon_hits_run_C2_x-654mm.root": (-654, 808501, 808519),
}
END_ONLY_META = {
    "photon_hits_run007.root": (-400, 317064, 717094),
    "photon_hits_run015.root": (0, 325136, 725198),
    "photon_hits_run023.root": (400, 333208, 733302),
}


def inspect(
    path: Path,
    expected_x: float,
    expected_events: int,
    expected_max_channel: int,
    configuration: str,
    threads: int,
    seed1: int | str = "",
    seed2: int | str = "",
    machine: str = "MSI",
    session: str = "exec14b_special_controls_20260613",
) -> dict[str, object]:
    event_ids: set[int] = set()
    gun_x: set[float] = set()
    channels: set[int] = set()
    with uproot.open(path) as root_file:
        tree = root_file["sipm_hits"]
        entries = int(tree.num_entries)
        for chunk in tree.iterate(
            ["event_id", "global_id", "gun_x_mm"], step_size="250 MB", library="np"
        ):
            event_ids.update(np.unique(chunk["event_id"]).astype(int).tolist())
            gun_x.update(np.unique(chunk["gun_x_mm"]).astype(float).tolist())
            channels.update(np.unique(chunk["global_id"]).astype(int).tolist())
    valid = (
        entries > 0
        and event_ids == set(range(expected_events))
        and gun_x == {float(expected_x)}
        and channels
        and channels == set(range(expected_max_channel + 1))
    )
    return {
        "path": str(path.resolve()),
        "configuration": configuration,
        "material": "OPSC-106",
        "readout": "End" if configuration == "end_only_End" else "EndTop",
        "jitter_ns": 0,
        "machine": machine,
        "session": session,
        "generation_date": datetime.fromtimestamp(path.stat().st_mtime).astimezone().date().isoformat(),
        "threads": threads,
        "seed1": seed1,
        "seed2": seed2,
        "expected_x_mm": expected_x,
        "expected_events": expected_events,
        "expected_channel_max": expected_max_channel,
        "entries": entries,
        "events": len(event_ids),
        "channel_min": min(channels) if channels else "",
        "channel_max": max(channels) if channels else "",
        "gun_x_values": ";".join(str(value) for value in sorted(gun_x)),
        "valid": valid,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("/home/reriosto/SHiP/t0minidaq/results_ej230/data"),
    )
    parser.add_argument(
        "--special-dir",
        type=Path,
        default=REPO / "results_ej230_analysis/sim_missing",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "results_ej230_analysis/root_validation_exec14b.csv",
    )
    args = parser.parse_args()
    rows = []
    for x_mm in POSITIONS:
        path = args.data_dir / f"photon_hits_x{x_mm}mm.root"
        if not path.is_file():
            rows.append({"path": str(path), "expected_x_mm": x_mm, "valid": False})
            continue
        rows.append(inspect(
            path,
            x_mm,
            2000,
            85,
            "main_scan_EndTop",
            24,
            machine="t0minidaq",
            session="session_20260612_171407",
        ))

    for name, (x_mm, seed1, seed2) in WINDOW_META.items():
        path = args.special_dir / "window_dip" / name
        if path.is_file():
            rows.append(inspect(path, x_mm, 2000, 85, "window_dip_EndTop", 16, seed1, seed2))
        else:
            rows.append({
                "path": str(path.resolve()),
                "configuration": "window_dip_EndTop",
                "expected_x_mm": x_mm,
                "expected_events": 2000,
                "valid": False,
            })
    for name, (x_mm, seed1, seed2) in END_ONLY_META.items():
        path = args.special_dir / "end_only" / name
        if path.is_file():
            rows.append(inspect(path, x_mm, 10000, 15, "end_only_End", 16, seed1, seed2))
        else:
            rows.append({
                "path": str(path.resolve()),
                "configuration": "end_only_End",
                "expected_x_mm": x_mm,
                "expected_events": 10000,
                "valid": False,
            })

    fields = sorted({key for row in rows for key in row})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    main_rows = rows[:len(POSITIONS)]
    print(
        f"main_valid={sum(row.get('valid') is True for row in main_rows)}/{len(POSITIONS)} "
        f"special_valid={sum(row.get('valid') is True for row in rows[len(POSITIONS):])}/"
        f"{len(rows) - len(POSITIONS)}"
    )
    return 0 if all(row.get("valid") is True for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
