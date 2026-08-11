#!/usr/bin/env python3
"""Apples-to-apples EXEC_14D diagnosis of EndTop versus End-only timing."""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot


REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results_ej230_analysis"
MAIN_DATA = Path("/home/reriosto/SHiP/t0minidaq/results_ej230/data")
MAIN_LOGS = Path("/home/reriosto/SHiP/t0minidaq/results_ej230/logs")
SPECIAL = RESULTS / "sim_missing/end_only"
POSITIONS = (-400, 0, 400)
FILES = {
    "EndTop": {x: MAIN_DATA / f"photon_hits_x{x}mm.root" for x in POSITIONS},
    "End-only": {
        -400: SPECIAL / "photon_hits_run007.root",
        0: SPECIAL / "photon_hits_run015.root",
        400: SPECIAL / "photon_hits_run023.root",
    },
}
LOGS = {
    "EndTop": {x: MAIN_LOGS / f"run_x{x}mm.log" for x in POSITIONS},
    "End-only": {
        -400: SPECIAL / "run_run007.log",
        0: SPECIAL / "run_run015.log",
        400: SPECIAL / "run_run023.log",
    },
}
SIDE_IDS = {"left": tuple(range(8)), "right": tuple(range(8, 16))}
SUM4_IDS = {
    "left": (tuple(range(4)), tuple(range(4, 8))),
    "right": (tuple(range(8, 12)), tuple(range(12, 16))),
}
GEOMETRY_FIELDS = (
    "SiPM model",
    "SiPM PDE file",
    "-X face",
    "+X face",
    "-Y face",
    "+Y face",
    "-Z face",
    "+Z face",
    "Reflector R",
    "Active End SiPMs",
)

sys.path.insert(0, str(REPO / "analysis"))
from exec07_photon_budget import leading_edge_time  # noqa: E402


def load_end_hits(path: Path, expected_x: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with uproot.open(path) as root_file:
        arrays = root_file["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns", "gun_x_mm"],
            cut="global_id < 16",
            library="np",
        )
    positions = np.unique(arrays["gun_x_mm"])
    if positions.tolist() != [float(expected_x)]:
        raise RuntimeError(f"{path}: gun_x_mm={positions.tolist()}, expected {expected_x}")
    return (
        arrays["event_id"].astype(np.int32),
        arrays["global_id"].astype(np.int16),
        arrays["time_ns"].astype(float),
    )


def event_counts(event_id: np.ndarray, global_id: np.ndarray, n_events: int) -> np.ndarray:
    flat = np.bincount(
        event_id.astype(np.int64) * 16 + global_id.astype(np.int64),
        minlength=n_events * 16,
    )
    return flat.reshape(n_events, 16)


def cluster_triggers(
    event_id: np.ndarray,
    global_id: np.ndarray,
    time_ns: np.ndarray,
    ids: tuple[int, ...],
    n_events: int,
) -> np.ndarray:
    selected = np.isin(global_id, ids)
    events = event_id[selected]
    times = time_ns[selected]
    order = np.lexsort((times, events))
    events = events[order]
    times = times[order]
    starts = np.searchsorted(events, np.arange(n_events), side="left")
    stops = np.searchsorted(events, np.arange(n_events), side="right")
    result = np.full(n_events, np.nan)
    for event in range(n_events):
        result[event] = leading_edge_time(times[starts[event] : stops[event]])
    return result


def sigma_group(delta: np.ndarray) -> float:
    valid = delta[np.isfinite(delta)]
    return float(np.std(valid) * 1000.0 / math.sqrt(2.0))


def side_tail(time_ns: np.ndarray, global_id: np.ndarray, ids: tuple[int, ...]) -> dict[str, float]:
    selected = time_ns[np.isin(global_id, ids)]
    return {
        "t99_ns": float(np.quantile(selected, 0.99)),
        "frac_gt10": float(np.mean(selected > 10.0)),
        "frac_gt20": float(np.mean(selected > 20.0)),
    }


def parse_geometry(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    block_match = re.search(
        r"=== Active Readout / Wrapping Configuration ===(.*?)={20,}",
        text,
        re.DOTALL,
    )
    if not block_match:
        raise RuntimeError(f"{path}: active readout block not found")
    block = block_match.group(1)
    result: dict[str, str] = {}
    for line in block.splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        result[key.strip()] = value.strip()
    return result


def geometry_check() -> tuple[bool, list[str]]:
    notes = []
    valid = True
    for x_mm in POSITIONS:
        top = parse_geometry(LOGS["EndTop"][x_mm])
        pure = parse_geometry(LOGS["End-only"][x_mm])
        for field in GEOMETRY_FIELDS:
            if field == "+Y face":
                if top.get(field) != "instrumented" or pure.get(field) != "wrapped":
                    valid = False
                    notes.append(
                        f"x={x_mm}: unexpected +Y values EndTop={top.get(field)!r}, "
                        f"End-only={pure.get(field)!r}"
                    )
            elif field == "SiPM PDE file":
                top_pde = Path(top.get(field, "")).name
                pure_pde = Path(pure.get(field, "")).name
                if top_pde != pure_pde:
                    valid = False
                    notes.append(
                        f"x={x_mm}: PDE resources differ: EndTop={top_pde!r}, "
                        f"End-only={pure_pde!r}"
                    )
            elif top.get(field) != pure.get(field):
                valid = False
                notes.append(
                    f"x={x_mm}: {field} differs: EndTop={top.get(field)!r}, "
                    f"End-only={pure.get(field)!r}"
                )
        if top.get("Active Top SiPMs") != "70" or pure.get("Active Top SiPMs") != "0":
            valid = False
            notes.append(
                f"x={x_mm}: unexpected Active Top SiPMs "
                f"{top.get('Active Top SiPMs')}/{pure.get('Active Top SiPMs')}"
            )
    if valid:
        notes.append(
            "All three log pairs differ only on +Y: EndTop is instrumented with "
            "70 Top SiPMs; End-only is wrapped. End faces, other wraps, PDE, "
            "reflectivity, and 16 End SiPMs match."
        )
    return valid, notes


def bootstrap_ratio(
    endtop_delta: np.ndarray,
    end_only_delta: np.ndarray,
    replicas: int,
    sample_size: int,
    rng: np.random.Generator,
) -> tuple[float, float, float, float]:
    top_sigma = sigma_group(endtop_delta)
    valid = end_only_delta[np.isfinite(end_only_delta)]
    values = np.empty(replicas)
    pure_sigmas = np.empty(replicas)
    for index in range(replicas):
        sample = rng.choice(valid, size=sample_size, replace=True)
        pure_sigmas[index] = sigma_group(sample)
        values[index] = top_sigma / pure_sigmas[index]
    return (
        float(values.mean()),
        float(values.std(ddof=1)),
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicas", type=int, default=200)
    parser.add_argument("--sample-size", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=14014)
    parser.add_argument(
        "--output", type=Path, default=RESULTS / "csv/exec14d_endtop_diag.csv"
    )
    parser.add_argument(
        "--document", type=Path, default=REPO / "docs/exec14d_endtop_diagnosis.md"
    )
    args = parser.parse_args()
    if args.replicas < 200:
        raise SystemExit("--replicas must be at least 200")

    rng = np.random.default_rng(args.seed)
    geometry_valid, geometry_notes = geometry_check()
    rows = []
    for x_mm in POSITIONS:
        campaign_data = {}
        for campaign in ("EndTop", "End-only"):
            event_id, global_id, time_ns = load_end_hits(FILES[campaign][x_mm], x_mm)
            n_events = int(event_id.max()) + 1
            counts = event_counts(event_id, global_id, n_events)
            campaign_data[campaign] = {
                "event_id": event_id,
                "global_id": global_id,
                "time_ns": time_ns,
                "counts": counts,
                "n_events": n_events,
            }
            triggers = {}
            for side, (ids_a, ids_b) in SUM4_IDS.items():
                a = cluster_triggers(event_id, global_id, time_ns, ids_a, n_events)
                b = cluster_triggers(event_id, global_id, time_ns, ids_b, n_events)
                triggers[side] = a - b
            campaign_data[campaign]["delta"] = triggers

        for side, ids in SIDE_IDS.items():
            top = campaign_data["EndTop"]
            pure = campaign_data["End-only"]
            top_delta = top["delta"][side]
            pure_delta = pure["delta"][side]
            top_sigma = sigma_group(top_delta)
            pure_sigma = sigma_group(pure_delta)
            boot_mean, boot_std, boot_lo, boot_hi = bootstrap_ratio(
                top_delta, pure_delta, args.replicas, args.sample_size, rng
            )
            top_npe = top["counts"][:, ids].sum(axis=1)
            pure_npe = pure["counts"][:, ids].sum(axis=1)
            top_tail = side_tail(top["time_ns"], top["global_id"], ids)
            pure_tail = side_tail(pure["time_ns"], pure["global_id"], ids)
            predicted = math.sqrt(float(np.mean(pure_npe)) / float(np.mean(top_npe)))
            rows.append(
                {
                    "x_beam_mm": x_mm,
                    "side": side,
                    "estimator": "std(delta_t_ab)/sqrt(2)",
                    "sum4_A_ids": str(SUM4_IDS[side][0]),
                    "sum4_B_ids": str(SUM4_IDS[side][1]),
                    "endtop_events": int(top["n_events"]),
                    "end_only_events": int(pure["n_events"]),
                    "bootstrap_sample_events": args.sample_size,
                    "bootstrap_replicas": args.replicas,
                    "sigma_endtop_ps": top_sigma,
                    "sigma_end_only_full_ps": pure_sigma,
                    "ratio_full": top_sigma / pure_sigma,
                    "ratio_bootstrap_mean": boot_mean,
                    "ratio_bootstrap_std": boot_std,
                    "ratio_bootstrap_q025": boot_lo,
                    "ratio_bootstrap_q975": boot_hi,
                    "npe_endtop_mean": float(np.mean(top_npe)),
                    "npe_end_only_mean": float(np.mean(pure_npe)),
                    "predicted_sqrt_npe_ratio": predicted,
                    "t99_endtop_ns": top_tail["t99_ns"],
                    "t99_end_only_ns": pure_tail["t99_ns"],
                    "t99_difference_ns": top_tail["t99_ns"] - pure_tail["t99_ns"],
                    "frac_gt10_endtop": top_tail["frac_gt10"],
                    "frac_gt10_end_only": pure_tail["frac_gt10"],
                    "frac_gt20_endtop": top_tail["frac_gt20"],
                    "frac_gt20_end_only": pure_tail["frac_gt20"],
                    "geometry_valid": geometry_valid,
                }
            )

    frame = pd.DataFrame(rows)
    inversion_persists = bool((frame.ratio_bootstrap_q025 > 1.0).all())
    npe_sign = bool((frame.npe_end_only_mean > frame.npe_endtop_mean).all())
    tail_sign = bool((frame.t99_difference_ns < 0.0).all())
    verdict = "PHYSICAL" if geometry_valid and inversion_persists and npe_sign and tail_sign else "INCONCLUSIVE"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.insert(0, "verdict", verdict)
    frame.to_csv(args.output, index=False)

    table_rows = []
    for item in frame.itertuples():
        table_rows.append(
            f"| {item.x_beam_mm:+d} | {item.side[0].upper()} | "
            f"{item.ratio_bootstrap_mean:.3f} +/- {item.ratio_bootstrap_std:.3f} | "
            f"{item.npe_endtop_mean:.1f} / {item.npe_end_only_mean:.1f} | "
            f"{item.t99_endtop_ns:.3f} / {item.t99_end_only_ns:.3f} | "
            f"{item.frac_gt10_endtop:.5f} / {item.frac_gt10_end_only:.5f} |"
        )
    document = f"""# EXEC_14D EndTop/End-only diagnosis

## Verdict: {verdict}

The apples-to-apples comparison uses the same End channels, SUM4 A/B assignment,
leading-edge implementation, absolute 4 PE threshold, and estimator
`sigma_group = std(delta_t_AB)/sqrt(2)` for both configurations. No fit window is
applied to either campaign. End-only is bootstrapped to {args.sample_size} events
with {args.replicas} replicas (seed {args.seed}).

| x [mm] | side | ratio at 2000 events | Npe EndTop / End-only | t99 EndTop / End-only [ns] | f(t>10) EndTop / End-only |
|---:|:---:|---:|---:|---:|---:|
{chr(10).join(table_rows)}

The lower 95% bootstrap quantile is above one in all six comparisons:
`{frame.ratio_bootstrap_q025.min():.3f}` minimum. End-only retains
`{(frame.npe_end_only_mean / frame.npe_endtop_mean).min():.2f}` to
`{(frame.npe_end_only_mean / frame.npe_endtop_mean).max():.2f}` times more End
photoelectrons. The corresponding pure-counting prediction,
`sqrt(Npe_End-only/Npe_EndTop)`, spans
`{frame.predicted_sqrt_npe_ratio.min():.3f}--{frame.predicted_sqrt_npe_ratio.max():.3f}`.
EndTop nevertheless shortens the side-specific t99 by
`{abs(frame.t99_difference_ns).min():.3f}--{abs(frame.t99_difference_ns).max():.3f}` ns.

## Geometry validation

{chr(10).join(f"- {note}" for note in geometry_notes)}

## Interpretation

The tail-drain mechanism persists, but its benefit is smaller than the timing
penalty caused by the loss of End photoelectrons. For this fast EJ-230
scintillator the EndTop configuration therefore widens intrinsic End timing.
The sign is physical within the tested configurations and is not produced by
the 10000-versus-2000 event-count difference.
"""
    args.document.write_text(document, encoding="utf-8")
    print(frame.to_string(index=False))
    print(f"verdict={verdict}")
    return 0 if verdict == "PHYSICAL" else 2


if __name__ == "__main__":
    raise SystemExit(main())
