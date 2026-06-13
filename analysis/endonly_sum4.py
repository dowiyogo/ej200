#!/usr/bin/env python3
"""END-only attenuation and congruent SUM4 timing analysis."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
import uproot
from scipy.optimize import curve_fit

BAR_HALF_LENGTH_MM = 700.0
SPR_RISE_NS = 0.5
SPR_FALL_NS = 5.0
LEADING_EDGE_THRESHOLD_PE = 4.0
DEFAULT_READOUT_JITTER_PS = 20.0
SQRT2 = math.sqrt(2.0)


@dataclass
class Point:
    x_mm: float
    n_events: int
    n_hits: int
    npe_left_mean: float
    npe_left_sem: float
    npe_right_mean: float
    npe_right_sem: float
    n_triggered: int
    trigger_efficiency: float
    mean_delta_ns: float
    sigma_lr_core_ps: float
    sigma_lr_core_err_ps: float
    sigma_lr_rms_ps: float
    sigma_single_ps: float
    sigma_single_err_ps: float
    fit_used: bool


def pulse(slow_state: float, fast_state: float, delta_ns: float) -> float:
    peak = (
        SPR_RISE_NS * SPR_FALL_NS / (SPR_FALL_NS - SPR_RISE_NS)
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
    """Exact Python port of congruent_sum4_timing.C LeadingEdgeTime."""
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


def modeled_total_sigma_ps(intrinsic_sigma_ps: float) -> float:
    """Quadrature total using the SiPMSD default readout jitter."""
    return math.hypot(intrinsic_sigma_ps, DEFAULT_READOUT_JITTER_PS)


def require_physical_ordering(intrinsic_sigma_ps: float, total_sigma_ps: float) -> None:
    if not total_sigma_ps > intrinsic_sigma_ps:
        raise ValueError(
            f"non-physical resolution ordering: intrinsic={intrinsic_sigma_ps}, total={total_sigma_ps}"
        )


def robust_sigma(values: np.ndarray) -> float:
    median = float(np.median(values))
    return 1.4826 * float(np.median(np.abs(values - median)))


def gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def fit_core(values: np.ndarray) -> tuple[float, float, float, float, bool]:
    """Mirror congruent_sum4_timing.C: 100 bins and four +/-2 sigma core fits."""
    if values.size < 20:
        rms = float(np.std(values, ddof=1)) if values.size > 1 else math.nan
        return float(np.mean(values)) if values.size else math.nan, rms, math.nan, rms, False

    center = float(np.median(values))
    sigma = max(robust_sigma(values), 1.0e-4)
    lo, hi = center - 8.0 * sigma, center + 8.0 * sigma
    counts, edges = np.histogram(values, bins=100, range=(lo, hi))
    centers = 0.5 * (edges[:-1] + edges[1:])
    amplitude = float(np.max(counts))
    covariance = None
    fit_used = False
    for _ in range(4):
        selected = (centers >= center - 2.0 * sigma) & (centers <= center + 2.0 * sigma)
        if np.count_nonzero(selected) < 4:
            break
        try:
            parameters, covariance = curve_fit(
                gaussian,
                centers[selected],
                counts[selected],
                p0=(amplitude, center, sigma),
                bounds=((0.0, lo, 1.0e-6), (np.inf, hi, hi - lo)),
                maxfev=20000,
            )
        except (RuntimeError, ValueError):
            break
        amplitude, center, sigma = map(float, parameters)
        sigma = abs(sigma)
        fit_used = True

    rms = float(np.std(values, ddof=1))
    if not fit_used:
        return float(np.mean(values)), rms, rms / math.sqrt(2.0 * (values.size - 1)), rms, False
    sigma_error = (
        math.sqrt(float(covariance[2, 2]))
        if covariance is not None and covariance[2, 2] >= 0.0
        else math.nan
    )
    return center, sigma, sigma_error, rms, True


def metadata_events(input_dir: pathlib.Path) -> int | None:
    path = input_dir / "run_metadata.txt"
    if not path.is_file():
        return None
    for line in path.read_text().splitlines():
        if line.startswith("events_per_position="):
            return int(line.split("=", 1)[1])
    return None


def analyze_file(path: pathlib.Path, configured_events: int | None) -> Point:
    with uproot.open(path) as root_file:
        tree = root_file["sipm_hits"]
        arrays = tree.arrays(["event_id", "global_id", "time_ns", "gun_x_mm"], library="np")

    event_id = arrays["event_id"].astype(np.int64)
    global_id = arrays["global_id"].astype(np.int64)
    time_ns = arrays["time_ns"].astype(float)
    positions = np.unique(arrays["gun_x_mm"].astype(float))
    if positions.size != 1:
        raise RuntimeError(f"{path}: expected one gun_x_mm, found {positions.tolist()}")
    x_mm = float(positions[0])
    inferred_events = int(np.max(event_id)) + 1 if event_id.size else 0
    n_events = configured_events if configured_events is not None else inferred_events
    if n_events < inferred_events or n_events <= 0:
        raise RuntimeError(f"{path}: invalid event count {n_events}, observed event_id up to {inferred_events - 1}")

    end_mask = (global_id >= 0) & (global_id < 16)
    if np.count_nonzero(~end_mask):
        raise RuntimeError(f"{path}: ntuple contains non-END global IDs")
    end_event = event_id[end_mask]
    end_id = global_id[end_mask]
    end_time = time_ns[end_mask]

    left_counts = np.bincount(end_event[end_id < 8], minlength=n_events)
    right_counts = np.bincount(end_event[end_id >= 8], minlength=n_events)
    npe_left_mean = float(np.mean(left_counts))
    npe_right_mean = float(np.mean(right_counts))
    npe_left_sem = float(np.std(left_counts, ddof=1) / math.sqrt(n_events)) if n_events > 1 else 0.0
    npe_right_sem = float(np.std(right_counts, ddof=1) / math.sqrt(n_events)) if n_events > 1 else 0.0

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
    delta = left[selected] - right[selected]
    mean, sigma, sigma_error, rms, fit_used = fit_core(delta)
    return Point(
        x_mm=x_mm,
        n_events=n_events,
        n_hits=int(end_event.size),
        npe_left_mean=npe_left_mean,
        npe_left_sem=npe_left_sem,
        npe_right_mean=npe_right_mean,
        npe_right_sem=npe_right_sem,
        n_triggered=int(delta.size),
        trigger_efficiency=float(delta.size / n_events),
        mean_delta_ns=mean,
        sigma_lr_core_ps=1000.0 * sigma,
        sigma_lr_core_err_ps=1000.0 * sigma_error,
        sigma_lr_rms_ps=1000.0 * rms,
        sigma_single_ps=1000.0 * sigma / SQRT2,
        sigma_single_err_ps=1000.0 * sigma_error / SQRT2,
        fit_used=fit_used,
    )


def attenuation_model(distance_mm: np.ndarray, n0: float, lambda_mm: float) -> np.ndarray:
    return n0 * np.exp(-distance_mm / lambda_mm)


def attenuation_fit(distance: np.ndarray, npe: np.ndarray, error: np.ndarray) -> dict[str, object]:
    valid = np.isfinite(distance) & np.isfinite(npe) & np.isfinite(error) & (npe > 0.0)
    if np.count_nonzero(valid) < 3 or np.unique(distance[valid]).size < 2:
        return {
            "n0_pe": math.nan,
            "n0_error_pe": math.nan,
            "lambda_eff_cm": math.nan,
            "lambda_eff_error_cm": math.nan,
            "n_points": int(np.count_nonzero(valid)),
            "status": "insufficient_points",
        }
    distance = distance[valid]
    npe = npe[valid]
    error = np.maximum(error[valid], 1.0e-9)
    try:
        parameters, covariance = curve_fit(
            attenuation_model,
            distance,
            npe,
            sigma=error,
            absolute_sigma=True,
            p0=(float(np.max(npe)), 1000.0),
            bounds=((0.0, 1.0), (np.inf, np.inf)),
            maxfev=20000,
        )
    except (RuntimeError, ValueError):
        return {
            "n0_pe": math.nan,
            "n0_error_pe": math.nan,
            "lambda_eff_cm": math.nan,
            "lambda_eff_error_cm": math.nan,
            "n_points": int(distance.size),
            "status": "fit_failed",
        }
    errors = np.sqrt(np.diag(covariance))
    return {
        "n0_pe": float(parameters[0]),
        "n0_error_pe": float(errors[0]),
        "lambda_eff_cm": float(parameters[1] / 10.0),
        "lambda_eff_error_cm": float(errors[1] / 10.0),
        "n_points": int(distance.size),
        "status": "ok",
    }


def write_csv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path)
    parser.add_argument("--events", type=int)
    args = parser.parse_args()

    files = sorted(args.input_dir.glob("photon_hits_x*mm.root"))
    if not files:
        raise SystemExit(f"no photon_hits_x*mm.root files in {args.input_dir}")
    output_dir = args.output_dir or args.input_dir / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    configured_events = args.events or metadata_events(args.input_dir)

    points = sorted(
        (analyze_file(path, configured_events) for path in files),
        key=lambda point: point.x_mm,
    )
    curve_rows = []
    timing_rows = []
    for point in points:
        curve_rows.append(
            {
                "x_mm": point.x_mm,
                "n_events": point.n_events,
                "distance_left_mm": BAR_HALF_LENGTH_MM + point.x_mm,
                "distance_right_mm": BAR_HALF_LENGTH_MM - point.x_mm,
                "npe_left_mean": point.npe_left_mean,
                "npe_left_sem": point.npe_left_sem,
                "npe_right_mean": point.npe_right_mean,
                "npe_right_sem": point.npe_right_sem,
                "n_end_hits": point.n_hits,
            }
        )
        timing_rows.append(
            {
                "x_mm": point.x_mm,
                "n_events": point.n_events,
                "n_triggered_lr": point.n_triggered,
                "trigger_efficiency": point.trigger_efficiency,
                "mean_delta_lr_ns": point.mean_delta_ns,
                "sigma_lr_core_ps": point.sigma_lr_core_ps,
                "sigma_lr_core_error_ps": point.sigma_lr_core_err_ps,
                "sigma_lr_rms_ps": point.sigma_lr_rms_ps,
                "sigma_single_ps": point.sigma_single_ps,
                "sigma_single_error_ps": point.sigma_single_err_ps,
                "sigma_total_default_readout_ps": modeled_total_sigma_ps(point.sigma_single_ps),
                "physical_ordering_pass": (
                    modeled_total_sigma_ps(point.sigma_single_ps) > point.sigma_single_ps
                ),
                "fit_used": point.fit_used,
            }
        )
    write_csv(output_dir / "attenuation_curve.csv", curve_rows)
    write_csv(output_dir / "sigma_t_sum4.csv", timing_rows)

    x = np.asarray([point.x_mm for point in points])
    left = np.asarray([point.npe_left_mean for point in points])
    right = np.asarray([point.npe_right_mean for point in points])
    left_error = np.asarray([point.npe_left_sem for point in points])
    right_error = np.asarray([point.npe_right_sem for point in points])

    figure, axis = plt.subplots()
    axis.errorbar(x, left, yerr=left_error, marker="o", label="END -X")
    axis.errorbar(x, right, yerr=right_error, marker="s", label="END +X")
    axis.set(xlabel="x [mm]", ylabel="mean Npe / event", title="EJ-230 end-only attenuation")
    axis.grid(True)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_dir / "attenuation_curve.png", dpi=160)
    plt.close(figure)

    figure, axis = plt.subplots()
    axis.errorbar(
        x,
        [point.sigma_single_ps for point in points],
        yerr=[point.sigma_single_err_ps for point in points],
        marker="o",
    )
    axis.set(
        xlabel="x [mm]",
        ylabel="sigma(t_L - t_R) / sqrt(2) [ps]",
        title="EJ-230 END SUM4 timing",
    )
    axis.grid(True)
    figure.tight_layout()
    figure.savefig(output_dir / "sigma_t_sum4.png", dpi=160)
    plt.close(figure)
    fit_rows = []
    for side, distance, npe, error in (
        ("left", BAR_HALF_LENGTH_MM + x, left, left_error),
        ("right", BAR_HALF_LENGTH_MM - x, right, right_error),
        (
            "combined",
            np.concatenate((BAR_HALF_LENGTH_MM + x, BAR_HALF_LENGTH_MM - x)),
            np.concatenate((left, right)),
            np.concatenate((left_error, right_error)),
        ),
    ):
        fit_rows.append({"side": side, **attenuation_fit(distance, npe, error)})
    write_csv(output_dir / "attenuation_fits.csv", fit_rows)

    (output_dir / "analysis_metadata.txt").write_text(
        "\n".join(
            (
                "readout=END-only global IDs 0..15",
                "sum4={0..3},{4..7},{8..11},{12..15}",
                f"spr_rise_ns={SPR_RISE_NS}",
                f"spr_fall_ns={SPR_FALL_NS}",
                f"leading_edge_threshold_pe={LEADING_EDGE_THRESHOLD_PE}",
                "end_time=earliest SUM4 threshold crossing per end",
                "estimator=sigma(deltaT_LR)/sqrt(2)",
                "time_origin=Geant4 event start; input time_ns",
            )
        )
        + "\n"
    )
    print(f"END-only analysis complete: {len(points)} positions -> {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
