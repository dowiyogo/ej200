#!/usr/bin/env python3
"""Quick EXEC_22b summaries for the backpainted surface test."""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
from dataclasses import dataclass

import numpy as np
import uproot
from scipy.optimize import curve_fit

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SPR_RISE_NS = 0.5
SPR_FALL_NS = 5.0
LEADING_EDGE_THRESHOLD_PE = 4.0
BOOTSTRAP_REPLICAS = 400


@dataclass
class PhysicsPoint:
    x_mm: float
    n_events: int
    npe_end: float
    npe_end_l: float
    npe_end_r: float
    n_triggered_lr: int
    sigma_end_ps: float
    sigma_boot_ps: float
    fit_used: bool


def pulse(slow_state: float, fast_state: float, delta_ns: float) -> float:
    peak = (
        SPR_RISE_NS
        * SPR_FALL_NS
        / (SPR_FALL_NS - SPR_RISE_NS)
        * math.log(SPR_FALL_NS / SPR_RISE_NS)
    )
    norm = 1.0 / (
        math.exp(-peak / SPR_FALL_NS) - math.exp(-peak / SPR_RISE_NS)
    )
    return norm * (
        slow_state * math.exp(-delta_ns / SPR_FALL_NS)
        - fast_state * math.exp(-delta_ns / SPR_RISE_NS)
    )


def leading_edge_time(arrivals: np.ndarray) -> float:
    if arrivals.size == 0:
        return math.nan
    arrivals = np.sort(arrivals)
    slow_state = 0.0
    fast_state = 0.0
    index = 0
    while index < arrivals.size:
        current = float(arrivals[index])
        next_index = index
        while next_index < arrivals.size and arrivals[next_index] == current:
            slow_state += 1.0
            fast_state += 1.0
            next_index += 1
        interval = (
            float(arrivals[next_index] - current)
            if next_index < arrivals.size
            else math.inf
        )
        derivative0 = fast_state / SPR_RISE_NS - slow_state / SPR_FALL_NS
        if derivative0 > 0.0:
            peak_delta = math.log(
                fast_state * SPR_FALL_NS / (slow_state * SPR_RISE_NS)
            ) / (1.0 / SPR_RISE_NS - 1.0 / SPR_FALL_NS)
            reach = min(peak_delta, interval)
            if reach >= 0.0 and pulse(slow_state, fast_state, reach) >= LEADING_EDGE_THRESHOLD_PE:
                low, high = 0.0, reach
                for _ in range(60):
                    middle = 0.5 * (low + high)
                    if pulse(slow_state, fast_state, middle) >= LEADING_EDGE_THRESHOLD_PE:
                        high = middle
                    else:
                        low = middle
                return current + high
        if next_index >= arrivals.size:
            break
        slow_state *= math.exp(-interval / SPR_FALL_NS)
        fast_state *= math.exp(-interval / SPR_RISE_NS)
        index = next_index
    return math.nan


def earliest(first: float, second: float) -> float:
    if not math.isfinite(first):
        return second
    if not math.isfinite(second):
        return first
    return min(first, second)


def gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def mad_sigma(values: np.ndarray) -> float:
    median = float(np.median(values))
    return 1.4826 * float(np.median(np.abs(values - median)))


def fit_sigma_ps(values_ns: np.ndarray) -> tuple[float, bool]:
    values_ns = values_ns[np.isfinite(values_ns)]
    if values_ns.size < 20:
        return (1000.0 * float(np.std(values_ns, ddof=1)), False) if values_ns.size > 1 else (math.nan, False)

    median = float(np.median(values_ns))
    sigma = max(mad_sigma(values_ns), 1.0e-6)
    core = values_ns[np.abs(values_ns - median) <= 2.0 * sigma]
    if core.size < 20:
        return 1000.0 * sigma, False

    counts, edges = np.histogram(core, bins=min(80, max(20, int(math.sqrt(core.size) * 3))))
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak_center = float(centers[int(np.argmax(counts))])
    p0 = (float(np.max(counts)), peak_center, max(mad_sigma(core), 1.0e-6))
    try:
        params, _ = curve_fit(
            gaussian,
            centers,
            counts,
            p0=p0,
            bounds=((0.0, float(np.min(core)), 1.0e-8), (np.inf, float(np.max(core)), np.inf)),
            maxfev=20000,
        )
        return 1000.0 * abs(float(params[2])), True
    except (RuntimeError, ValueError):
        return 1000.0 * float(np.std(core, ddof=1)), False


def bootstrap_sigma_ps(values_ns: np.ndarray, rng: np.random.Generator) -> float:
    values_ns = values_ns[np.isfinite(values_ns)]
    if values_ns.size < 20:
        return math.nan
    estimates = []
    for _ in range(BOOTSTRAP_REPLICAS):
        sample = rng.choice(values_ns, size=values_ns.size, replace=True)
        estimates.append(fit_sigma_ps(sample)[0])
    return float(np.std(estimates, ddof=1))


def metadata_events(input_dir: pathlib.Path) -> int | None:
    path = input_dir / "run_metadata.txt"
    if not path.is_file():
        return None
    for line in path.read_text().splitlines():
        if line.startswith("events_per_position="):
            return int(line.split("=", 1)[1])
    return None


def analyze_root(path: pathlib.Path, configured_events: int | None, rng: np.random.Generator) -> PhysicsPoint:
    with uproot.open(path) as root_file:
        arrays = root_file["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns", "gun_x_mm"], library="np"
        )

    event_id = arrays["event_id"].astype(np.int64)
    global_id = arrays["global_id"].astype(np.int64)
    time_ns = arrays["time_ns"].astype(float)
    positions = np.unique(arrays["gun_x_mm"].astype(float))
    if positions.size != 1:
        raise RuntimeError(f"{path}: expected one gun_x_mm, found {positions.tolist()}")
    x_mm = float(positions[0])
    inferred_events = int(np.max(event_id)) + 1 if event_id.size else 0
    n_events = configured_events if configured_events is not None else inferred_events
    if n_events <= 0 or n_events < inferred_events:
        raise RuntimeError(f"{path}: invalid event count {n_events}")

    end_mask = (global_id >= 0) & (global_id < 16)
    end_event = event_id[end_mask]
    end_id = global_id[end_mask]
    end_time = time_ns[end_mask]
    left_counts = np.bincount(end_event[end_id < 8], minlength=n_events)
    right_counts = np.bincount(end_event[end_id >= 8], minlength=n_events)

    triggers = np.full((n_events, 4), np.nan)
    combined = end_event * 4 + end_id // 4
    order = np.argsort(combined, kind="stable")
    combined = combined[order]
    ordered_time = end_time[order]
    unique, starts = np.unique(combined, return_index=True)
    stops = np.r_[starts[1:], combined.size]
    for key, start, stop in zip(unique, starts, stops):
        event = int(key // 4)
        group = int(key % 4)
        if 0 <= event < n_events:
            triggers[event, group] = leading_edge_time(ordered_time[start:stop])

    left = np.array([earliest(a, b) for a, b in triggers[:, :2]])
    right = np.array([earliest(a, b) for a, b in triggers[:, 2:]])
    selected = np.isfinite(left) & np.isfinite(right)
    t_avg = 0.5 * (left[selected] + right[selected])
    sigma_ps, fit_used = fit_sigma_ps(t_avg)
    sigma_boot_ps = bootstrap_sigma_ps(t_avg, rng)
    return PhysicsPoint(
        x_mm=x_mm,
        n_events=n_events,
        npe_end=float(np.mean(left_counts + right_counts)),
        npe_end_l=float(np.mean(left_counts)),
        npe_end_r=float(np.mean(right_counts)),
        n_triggered_lr=int(t_avg.size),
        sigma_end_ps=sigma_ps,
        sigma_boot_ps=sigma_boot_ps,
        fit_used=fit_used,
    )


def write_csv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def json_clean(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_clean(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_clean(item) for item in value]
    return value


def summarize_speed(speed_csv: pathlib.Path, out_dir: pathlib.Path) -> dict[str, object]:
    total = 0
    scint_total = 0
    scint_super = 0
    max_ratio_scint = 0.0
    max_ratio_air = 0.0
    missing: dict[str, int] = {}
    positions: set[float] = set()
    scint_ratios = []

    with speed_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("x_mm") == "x_mm":
                continue
            total += 1
            x = row.get("x_mm", "")
            if x:
                positions.add(float(x))
            has_rindex = row["has_rindex"] == "1"
            material = row["material"]
            volume = row["volume"]
            ratio = float(row["ratio"]) if has_rindex else math.nan
            if not has_rindex:
                missing[f"{volume}:{material}"] = missing.get(f"{volume}:{material}", 0) + 1
                continue
            material_key = material.lower()
            is_scint = (
                row["is_scintillator"] == "1"
                or material_key in {"opsc-101", "ej-204", "ej-200", "ej-230"}
            )
            if is_scint:
                scint_total += 1
                max_ratio_scint = max(max_ratio_scint, ratio)
                if ratio > 1.001:
                    scint_super += 1
                scint_ratios.append(ratio)
            elif material == "G4_AIR":
                max_ratio_air = max(max_ratio_air, ratio)

    fig_path = out_dir / "speed_ratio_hist.png"
    plt.figure(figsize=(7, 4.5))
    if scint_ratios:
        plt.hist(scint_ratios, bins=80, range=(0.98, 1.02), histtype="step", color="black")
    plt.axvline(1.001, color="red", linestyle="--", linewidth=1)
    plt.xlabel("v / (c / n)")
    plt.ylabel("optical-photon steps")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=160)
    plt.close()

    summary = {
        "positions_mm": sorted(positions),
        "n_steps_total": total,
        "n_steps_scintillator": scint_total,
        "max_ratio_scintillator": max_ratio_scint if scint_total else None,
        "max_ratio_air": max_ratio_air if max_ratio_air > 0 else None,
        "superluminal_fraction_scintillator": (
            scint_super / scint_total if scint_total else None
        ),
        "missing_rindex": missing,
        "pass": (
            scint_total > 0
            and not missing
            and (scint_super / scint_total) < 0.001
            and max_ratio_scint <= 1.001
            and max_ratio_air <= 1.001
        ),
        "histogram": str(fig_path),
    }
    (out_dir / "speed_audit_summary.json").write_text(
        json.dumps(json_clean(summary), indent=2, allow_nan=False) + "\n"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    parser.add_argument("--speed-csv", required=True, type=pathlib.Path)
    parser.add_argument("--events", type=int)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    speed_summary = summarize_speed(args.speed_csv, args.output_dir)

    configured_events = args.events or metadata_events(args.input_dir)
    rng = np.random.default_rng(22022)
    points = sorted(
        (analyze_root(path, configured_events, rng) for path in args.input_dir.glob("photon_hits_x*mm.root")),
        key=lambda point: point.x_mm,
    )
    if not points:
        raise SystemExit(f"no photon_hits_x*mm.root files in {args.input_dir}")

    npe_rows = [
        {
            "x_mm": p.x_mm,
            "n_events": p.n_events,
            "npe_end": p.npe_end,
            "npe_end_l": p.npe_end_l,
            "npe_end_r": p.npe_end_r,
            "n_triggered_lr": p.n_triggered_lr,
        }
        for p in points
    ]
    sigma_rows = [
        {
            "x_mm": p.x_mm,
            "n_events": p.n_events,
            "n_triggered_lr": p.n_triggered_lr,
            "sigma_end_ps": p.sigma_end_ps,
            "sigma_boot_ps": p.sigma_boot_ps,
            "fit_used": p.fit_used,
        }
        for p in points
    ]
    write_csv(args.output_dir / "end_npe_summary.csv", npe_rows)
    write_csv(args.output_dir / "end_sigma_summary.csv", sigma_rows)

    center = min(points, key=lambda point: abs(point.x_mm))
    if not speed_summary["pass"]:
        verdict = "ABORTADO"
    elif center.npe_end >= 20.0 and center.sigma_end_ps < 250.0:
        verdict = "EXITO"
    elif 10.0 <= center.npe_end < 20.0 and center.sigma_end_ps < 250.0:
        verdict = "EXITO PARCIAL"
    else:
        verdict = "NO FUNCIONA"

    result = {
        "exec_id": "EXEC_22b",
        "finish_tested": "polishedbackpainted",
        "reflectivity": 0.95,
        "surface_rindex": 1.0,
        "speed_gate": {
            "pass": speed_summary["pass"],
            "max_ratio_scintillator": speed_summary["max_ratio_scintillator"],
            "superluminal_fraction_scintillator": speed_summary["superluminal_fraction_scintillator"],
        },
        "end_center": {
            "npe_before_exec22": 0.38,
            "sigma_before_exec22_ps": 985,
            "npe_after": center.npe_end,
            "sigma_after_ps": center.sigma_end_ps,
        },
        "points": [p.__dict__ for p in points],
        "verdict": verdict,
    }
    cleaned = json_clean(result)
    (args.output_dir / "results_exec22b.json").write_text(
        json.dumps(cleaned, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(cleaned, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
