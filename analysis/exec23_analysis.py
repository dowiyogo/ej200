#!/usr/bin/env python3
"""EXEC_23 mini-scan analysis for explicit air-gap reflector tests."""

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
BOOTSTRAP_REPLICAS = 300


@dataclass
class Point:
    x_mm: float
    n_events: int
    npe_end: float
    npe_l: float
    npe_r: float
    paired_fraction: float
    n_paired: int
    sigma_end_ps: float | None
    sigma_boot_ps: float | None
    fit_used: bool


def pulse(slow_state: float, fast_state: float, delta_ns: float) -> float:
    peak = (
        SPR_RISE_NS
        * SPR_FALL_NS
        / (SPR_FALL_NS - SPR_RISE_NS)
        * math.log(SPR_FALL_NS / SPR_RISE_NS)
    )
    norm = 1.0 / (math.exp(-peak / SPR_FALL_NS) - math.exp(-peak / SPR_RISE_NS))
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


def earliest(a: float, b: float) -> float:
    if not math.isfinite(a):
        return b
    if not math.isfinite(b):
        return a
    return min(a, b)


def gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def mad_sigma(values: np.ndarray) -> float:
    median = float(np.median(values))
    return 1.4826 * float(np.median(np.abs(values - median)))


def fit_sigma_ps(values_ns: np.ndarray) -> tuple[float | None, bool]:
    values_ns = values_ns[np.isfinite(values_ns)]
    if values_ns.size < 20:
        return None, False
    median = float(np.median(values_ns))
    sigma = max(mad_sigma(values_ns), 1.0e-6)
    core = values_ns[np.abs(values_ns - median) <= 2.0 * sigma]
    if core.size < 20:
        return 1000.0 * sigma, False
    counts, edges = np.histogram(core, bins=min(80, max(20, int(math.sqrt(core.size) * 3))))
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak_center = float(centers[int(np.argmax(counts))])
    try:
        params, _ = curve_fit(
            gaussian,
            centers,
            counts,
            p0=(float(np.max(counts)), peak_center, max(mad_sigma(core), 1.0e-6)),
            bounds=((0.0, float(np.min(core)), 1.0e-8), (np.inf, float(np.max(core)), np.inf)),
            maxfev=20000,
        )
        return 1000.0 * abs(float(params[2])), True
    except (RuntimeError, ValueError):
        return 1000.0 * float(np.std(core, ddof=1)), False


def bootstrap_sigma_ps(values_ns: np.ndarray, rng: np.random.Generator) -> float | None:
    values_ns = values_ns[np.isfinite(values_ns)]
    if values_ns.size < 20:
        return None
    estimates = []
    for _ in range(BOOTSTRAP_REPLICAS):
        sigma, _ = fit_sigma_ps(rng.choice(values_ns, size=values_ns.size, replace=True))
        if sigma is not None:
            estimates.append(sigma)
    return float(np.std(estimates, ddof=1)) if len(estimates) > 2 else None


def analyze_root(path: pathlib.Path, n_events: int, rng: np.random.Generator) -> tuple[Point, np.ndarray]:
    with uproot.open(path) as rf:
        arrays = rf["sipm_hits"].arrays(["event_id", "global_id", "time_ns", "gun_x_mm"], library="np")
    event_id = arrays["event_id"].astype(np.int64)
    global_id = arrays["global_id"].astype(np.int64)
    time_ns = arrays["time_ns"].astype(float)
    positions = np.unique(arrays["gun_x_mm"].astype(float))
    if positions.size != 1:
        raise RuntimeError(f"{path}: expected one gun_x_mm, got {positions.tolist()}")
    x_mm = float(positions[0])
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
    paired = np.isfinite(left) & np.isfinite(right)
    t_avg = 0.5 * (left[paired] + right[paired])
    sigma, fit_used = fit_sigma_ps(t_avg)
    boot = bootstrap_sigma_ps(t_avg, rng)
    return (
        Point(
            x_mm=x_mm,
            n_events=n_events,
            npe_end=float(np.mean(left_counts + right_counts)),
            npe_l=float(np.mean(left_counts)),
            npe_r=float(np.mean(right_counts)),
            paired_fraction=float(np.count_nonzero(paired) / n_events),
            n_paired=int(np.count_nonzero(paired)),
            sigma_end_ps=sigma,
            sigma_boot_ps=boot,
            fit_used=fit_used,
        ),
        t_avg,
    )


def write_csv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv_rows(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def summarize_speed(rows: list[dict[str, str]]) -> dict[str, object]:
    by_material: dict[str, dict[str, float]] = {}
    missing = []
    for row in rows:
        mat = row["material"]
        steps = int(row["n_steps"])
        max_ratio = float(row["max_ratio"])
        super_frac = float(row["superluminal_fraction"])
        missing_count = int(row["missing_rindex"])
        item = by_material.setdefault(mat, {"n_steps": 0, "max_ratio": 0.0, "weighted_super": 0.0})
        item["n_steps"] += steps
        item["max_ratio"] = max(item["max_ratio"], max_ratio)
        item["weighted_super"] += steps * super_frac
        if missing_count:
            missing.append({"volume": row["volume"], "material": mat, "count": missing_count})
    for item in by_material.values():
        item["superluminal_fraction"] = (
            item["weighted_super"] / item["n_steps"] if item["n_steps"] else 0.0
        )
        del item["weighted_super"]
    scint = by_material.get("opsc-101") or by_material.get("OPSC-101")
    speed_pass = bool(scint) and scint["max_ratio"] <= 1.001 and scint["superluminal_fraction"] < 0.001 and not missing
    return {
        "pass": speed_pass,
        "by_material": by_material,
        "max_ratio_scintillator": scint["max_ratio"] if scint else None,
        "superluminal_fraction_scintillator": scint["superluminal_fraction"] if scint else None,
        "volumes_without_rindex": missing,
    }


def summarize_history(rows: list[dict[str, str]]) -> dict[str, object]:
    if not rows:
        return {"n_detected": 0}
    ref_scint = np.array([int(r["n_reflections_scint_air"]) for r in rows])
    ref_refl = np.array([int(r["n_reflections_air_reflector"]) for r in rows])
    tir = np.array([int(r["n_total_internal_reflections"]) for r in rows])
    return {
        "n_detected": len(rows),
        "fraction_visited_air_gap": float(np.mean([int(r["visited_air_gap"]) for r in rows])),
        "fraction_visited_reflector_boundary": float(np.mean([int(r["visited_reflector_boundary"]) for r in rows])),
        "mean_scint_air_reflections": float(np.mean(ref_scint)),
        "mean_air_reflector_reflections": float(np.mean(ref_refl)),
        "mean_total_internal_reflections": float(np.mean(tir)),
    }


def make_plots(points: list[Point], center_t_avg: np.ndarray, history_rows: list[dict[str, str]], out_dir: pathlib.Path) -> None:
    png = out_dir / "png"
    png.mkdir(parents=True, exist_ok=True)
    xs = np.array([p.x_mm for p in points])
    order = np.argsort(xs)
    xs = xs[order]
    npe = np.array([points[i].npe_end for i in order])
    sigma = np.array([points[i].sigma_end_ps if points[i].sigma_end_ps is not None else np.nan for i in order])

    plt.figure(figsize=(6, 4))
    plt.plot(xs, npe, marker="o")
    plt.xlabel("x [mm]")
    plt.ylabel("<Npe> END")
    plt.tight_layout()
    plt.savefig(png / "npe_end_vs_x_quick.png", dpi=160)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.plot(xs, sigma, marker="o")
    plt.xlabel("x [mm]")
    plt.ylabel("sigma END [ps]")
    plt.tight_layout()
    plt.savefig(png / "sigma_end_vs_x_quick.png", dpi=160)
    plt.close()

    plt.figure(figsize=(6, 4))
    if center_t_avg.size:
        plt.hist(1000.0 * center_t_avg, bins=60, histtype="step")
    plt.xlabel("t_avg [ps]")
    plt.ylabel("events")
    plt.tight_layout()
    plt.savefig(png / "time_pdf_center.png", dpi=160)
    plt.close()

    plt.figure(figsize=(6, 4))
    center_rows = [r for r in history_rows if abs(float(r["gun_x_mm"])) < 1.0e-6]
    if center_rows:
        vals = [int(r["n_total_internal_reflections"]) + int(r["n_reflections_air_reflector"]) for r in center_rows]
        plt.hist(vals, bins=50, histtype="step")
    plt.xlabel("detected photon reflection count")
    plt.ylabel("photons")
    plt.tight_layout()
    plt.savefig(png / "detected_reflections_hist_center.png", dpi=160)
    plt.close()


def clean(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    parser.add_argument("--events", required=True, type=int)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(23023)
    points: list[Point] = []
    center_t_avg = np.array([])
    for root_path in sorted(args.input_dir.glob("photon_hits_x*mm.root")):
        point, t_avg = analyze_root(root_path, args.events, rng)
        points.append(point)
        if abs(point.x_mm) < 1.0e-6:
            center_t_avg = t_avg
    points.sort(key=lambda p: p.x_mm)

    write_csv(args.output_dir / "end_npe_summary.csv", [
        {
            "x_mm": p.x_mm,
            "events": p.n_events,
            "npe_end": p.npe_end,
            "npe_l": p.npe_l,
            "npe_r": p.npe_r,
        }
        for p in points
    ])
    write_csv(args.output_dir / "end_pairing_summary.csv", [
        {
            "x_mm": p.x_mm,
            "events": p.n_events,
            "n_paired": p.n_paired,
            "paired_fraction": p.paired_fraction,
        }
        for p in points
    ])
    write_csv(args.output_dir / "end_sigma_summary.csv", [
        {
            "x_mm": p.x_mm,
            "events": p.n_events,
            "paired_fraction": p.paired_fraction,
            "sigma_end_ps": p.sigma_end_ps,
            "bootstrap_ps": p.sigma_boot_ps,
            "fit_used": p.fit_used,
        }
        for p in points
    ])

    history_rows = read_csv_rows(args.output_dir / "detected_photon_history.csv")
    speed_summary = summarize_speed(read_csv_rows(args.output_dir / "speed_audit_steps.csv"))
    history_summary = summarize_history(history_rows)
    boundary_rows = read_csv_rows(args.output_dir / "optical_boundary_counts.csv")
    write_csv(args.output_dir / "end_time_structure_summary.csv", [
        {"metric": key, "value": value} for key, value in history_summary.items()
    ])
    make_plots(points, center_t_avg, history_rows, args.output_dir)

    center = min(points, key=lambda p: abs(p.x_mm)) if points else None
    center_pass_npe = center is not None and center.npe_end >= 10.0
    center_pass_sigma = center is not None and center.sigma_end_ps is not None and center.sigma_end_ps < 250.0
    center_pass_pair = center is not None and center.n_paired >= 20
    verdict = "PASS" if speed_summary["pass"] and center_pass_npe and center_pass_sigma and center_pass_pair else "FAIL"
    if center is not None and center.npe_end >= 20.0 and center_pass_sigma:
        strength = "strong"
    elif center is not None and center.npe_end >= 10.0:
        strength = "partial"
    else:
        strength = "fail"

    result = {
        "exec_id": "EXEC_23",
        "verdict": verdict,
        "mini_validation_strength": strength,
        "geometry": {
            "air_gap_mm": 0.10,
            "reflector_thickness_mm": 0.05,
            "scintillator_air_surface": "dielectric_dielectric/unified/polished/no_REFLECTIVITY",
            "air_reflector_surface": "dielectric_metal/unified/polished/REFLECTIVITY=0.95",
        },
        "speed_gate": speed_summary,
        "center_result": {
            "npe_exec21": 0.37,
            "npe_exec22b": 0.763,
            "npe_exec23": center.npe_end if center else None,
            "sigma_exec23_ps": center.sigma_end_ps if center else None,
            "paired_fraction_exec23": center.paired_fraction if center else None,
        },
        "points": [p.__dict__ for p in points],
        "history_summary": history_summary,
        "boundary_categories": len(boundary_rows),
        "scan11": {
            "ran": False,
            "n_positions": 0,
            "sigma_end_ps_max": None,
            "sigma_end_ps_center": None,
            "v_eff_cm_per_ns": None,
            "sigma_x_cm_center": None,
        },
        "endtop_mini": {
            "ran": False,
            "pass": None,
            "npe_end_center": None,
            "speed_gate_pass": None,
        },
    }
    (args.output_dir / "speed_audit_summary.json").write_text(
        json.dumps(clean(speed_summary), indent=2, allow_nan=False) + "\n"
    )
    (args.output_dir / "optical_history_summary.json").write_text(
        json.dumps(clean(history_summary), indent=2, allow_nan=False) + "\n"
    )
    (args.output_dir / "results_exec23_nominal.json").write_text(
        json.dumps(clean(result), indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(clean(result), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
