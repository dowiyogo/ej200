#!/usr/bin/env python3
"""EXEC_08 photon-budget analysis of the EXEC_07 EndTop scan."""

from __future__ import annotations

import argparse
import math
import pathlib
import sys
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from scipy.stats import norm, poisson

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from exec07.common import (  # noqa: E402
    BAR_HALF_LENGTH_MM,
    END_CLUSTERS,
    END_LEFT_IDS,
    END_RIGHT_IDS,
    EXPECTED_POSITIONS_MM,
    LEADING_EDGE_THRESHOLD_PE,
    N_CHANNELS,
    N_EVENTS,
    SIGMA_NUMERATOR_NS,
    SPR_FALL_NS,
    SPR_RISE_NS,
    TAU_D_NS,
    TAU_R_NS,
    TOP_IDS,
    TOP_LEFT_IDS,
    TOP_POSITIONS_MM,
    TOP_RIGHT_IDS,
    expected_file,
    nearest_end_name,
    nearest_top_ids,
)

TIME_BINS_NS = np.arange(0.0, 50.0 + 0.1, 0.1)
SUMMARY_COLUMNS = [
    "x_beam_mm", "grupo", "npe_mean", "npe_var", "npe_rms",
    "var_over_mean", "chi2ndf_poisson", "t_mean_ns", "t_rms_ns",
    "fpt_mean_ns", "fpt_rms_ns", "sigma_grupo_ps", "sigma_tavg_ps",
    "sigma_est_ps", "lambda_eff_cm",
]
KEY_POISSON_POSITIONS = (-690, -400, 0, 400, 690)
KEY_TIME_POSITIONS = (-690, -400, 0)


@dataclass
class PointResult:
    x_mm: int
    rows: list[dict]
    top_profile_mean: np.ndarray
    top_profile_sem: np.ndarray
    nearest_top_id: int
    nearest_top_3: tuple[int, ...]
    max_top_id: int
    representative_counts: dict[str, np.ndarray]
    representative_times: dict[str, np.ndarray]
    representative_fpt: dict[str, np.ndarray]
    sigma_delta_left_ps: float
    sigma_delta_right_ps: float
    sigma_tavg_ps: float
    delta_left_ns: np.ndarray
    delta_right_ns: np.ndarray
    tavg_ns: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=pathlib.Path)
    parser.add_argument(
        "--output-dir", type=pathlib.Path,
        default=pathlib.Path("analysis/exec07"),
    )
    parser.add_argument("--skip-channel-pdfs", action="store_true")
    return parser.parse_args()


def validate_inputs(data_dir: pathlib.Path) -> list[pathlib.Path]:
    errors: list[str] = []
    files: list[pathlib.Path] = []
    for x_mm in EXPECTED_POSITIONS_MM:
        path = expected_file(data_dir, x_mm)
        files.append(path)
        if not path.is_file():
            errors.append(f"missing {path}")
            continue
        events: set[int] = set()
        guns: set[float] = set()
        channels: set[int] = set()
        try:
            with uproot.open(path) as root_file:
                tree = root_file["sipm_hits"]
                for arrays in tree.iterate(
                    ["event_id", "global_id", "gun_x_mm"],
                    library="np", step_size="100 MB",
                ):
                    events.update(np.unique(arrays["event_id"]).astype(int).tolist())
                    channels.update(np.unique(arrays["global_id"]).astype(int).tolist())
                    guns.update(np.unique(arrays["gun_x_mm"]).astype(float).tolist())
        except Exception as error:  # noqa: BLE001
            errors.append(f"{path}: {type(error).__name__}: {error}")
            continue
        if events != set(range(N_EVENTS)):
            errors.append(f"{path}: expected event_id 0..1999, got {len(events)} unique")
        if guns != {float(x_mm)}:
            errors.append(f"{path}: gun_x_mm={sorted(guns)}, expected {x_mm}")
        if channels != set(range(N_CHANNELS)):
            errors.append(
                f"{path}: missing channels {sorted(set(range(N_CHANNELS)) - channels)}"
            )
    if errors:
        raise RuntimeError("blocking input validation failed:\n" + "\n".join(errors))
    return files


def spr_peak_time() -> float:
    return (
        SPR_RISE_NS * SPR_FALL_NS / (SPR_FALL_NS - SPR_RISE_NS)
        * math.log(SPR_FALL_NS / SPR_RISE_NS)
    )


def spr_norm() -> float:
    peak = spr_peak_time()
    return 1.0 / (math.exp(-peak / SPR_FALL_NS) - math.exp(-peak / SPR_RISE_NS))


def pulse(slow_state: float, fast_state: float, delta_ns: float) -> float:
    return spr_norm() * (
        slow_state * math.exp(-delta_ns / SPR_FALL_NS)
        - fast_state * math.exp(-delta_ns / SPR_RISE_NS)
    )


def leading_edge_time(arrivals: np.ndarray) -> float:
    """Port of congruent_sum4_timing.C:160-213; no walk correction (:305-315)."""
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
            if next_index < arrivals.size else math.inf
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


def poisson_chi2_ndf(counts: np.ndarray) -> float:
    mean = float(np.mean(counts))
    if mean <= 0.0:
        return math.nan
    values, observed = np.unique(counts.astype(int), return_counts=True)
    low, high = int(values.min()), int(values.max())
    support = np.arange(low, high + 1)
    obs = np.bincount(counts.astype(int) - low, minlength=len(support)).astype(float)
    expected = poisson.pmf(support, mean) * len(counts)
    expected[0] += poisson.cdf(low - 1, mean) * len(counts)
    expected[-1] += poisson.sf(high, mean) * len(counts)
    valid = expected >= 5.0
    if np.count_nonzero(valid) < 3:
        return math.nan
    return float(np.sum((obs[valid] - expected[valid]) ** 2 / expected[valid])
                 / (np.count_nonzero(valid) - 1))


def gaussian_chi2_ndf(counts: np.ndarray) -> float:
    mean = float(np.mean(counts))
    sigma = float(np.std(counts))
    if sigma <= 0.0:
        return math.nan
    low = int(np.min(counts))
    high = int(np.max(counts))
    support = np.arange(low, high + 1)
    obs = np.bincount(counts.astype(int) - low, minlength=len(support)).astype(float)
    expected = (
        norm.cdf(support + 0.5, mean, sigma)
        - norm.cdf(support - 0.5, mean, sigma)
    ) * len(counts)
    valid = expected >= 5.0
    if np.count_nonzero(valid) < 4:
        return math.nan
    return float(np.sum((obs[valid] - expected[valid]) ** 2 / expected[valid])
                 / (np.count_nonzero(valid) - 2))


def group_statistics(
    x_mm: int,
    name: str,
    ids: tuple[int, ...],
    counts_matrix: np.ndarray,
    hit_count: np.ndarray,
    time_sum: np.ndarray,
    time_sum_sq: np.ndarray,
    fpt_matrix: np.ndarray,
) -> tuple[dict, np.ndarray, np.ndarray]:
    counts = counts_matrix[:, ids].sum(axis=1)
    npe_mean = float(np.mean(counts))
    npe_var = float(np.var(counts))
    n_hits = int(np.sum(hit_count[list(ids)]))
    total_time = float(np.sum(time_sum[list(ids)]))
    total_time_sq = float(np.sum(time_sum_sq[list(ids)]))
    t_mean = total_time / n_hits if n_hits else math.nan
    t_rms = (
        math.sqrt(max(total_time_sq / n_hits - t_mean * t_mean, 0.0))
        if n_hits else math.nan
    )
    fpt = np.min(fpt_matrix[:, ids], axis=1)
    fpt[~np.isfinite(fpt)] = np.nan
    valid_fpt = fpt[np.isfinite(fpt)]
    is_top = bool(set(ids) & set(TOP_IDS)) and not bool(set(ids) & set(range(16)))
    row = {
        "x_beam_mm": x_mm,
        "grupo": name,
        "npe_mean": npe_mean,
        "npe_var": npe_var,
        "npe_rms": math.sqrt(npe_var),
        "var_over_mean": npe_var / npe_mean if npe_mean > 0 else math.nan,
        "chi2ndf_poisson": poisson_chi2_ndf(counts),
        "t_mean_ns": t_mean,
        "t_rms_ns": t_rms,
        "fpt_mean_ns": float(np.mean(valid_fpt)) if valid_fpt.size else math.nan,
        "fpt_rms_ns": float(np.std(valid_fpt)) if valid_fpt.size else math.nan,
        "sigma_grupo_ps": math.nan,
        "sigma_tavg_ps": math.nan,
        "sigma_est_ps": (
            1000.0 * SIGMA_NUMERATOR_NS / math.sqrt(npe_mean)
            if is_top and npe_mean > 0 else math.nan
        ),
        "lambda_eff_cm": math.nan,
    }
    return row, counts, fpt


def cluster_leading_edges(
    event_id: np.ndarray, global_id: np.ndarray, time_ns: np.ndarray,
    ids: tuple[int, ...],
) -> np.ndarray:
    output = np.full(N_EVENTS, np.nan)
    mask = np.isin(global_id, ids)
    events = event_id[mask]
    times = time_ns[mask]
    order = np.lexsort((times, events))
    events = events[order]
    times = times[order]
    starts = np.searchsorted(events, np.arange(N_EVENTS), side="left")
    stops = np.searchsorted(events, np.arange(N_EVENTS), side="right")
    for event in range(N_EVENTS):
        output[event] = leading_edge_time(times[starts[event]:stops[event]])
    return output


def save_channel_pages(
    npe_pdf: PdfPages,
    time_pdf: PdfPages,
    x_mm: int,
    counts_matrix: np.ndarray,
    time_hist: np.ndarray,
    fpt_matrix: np.ndarray,
) -> None:
    fig, axes = plt.subplots(9, 10, figsize=(18, 15))
    for channel, ax in enumerate(axes.flat[:N_CHANNELS]):
        counts = counts_matrix[:, channel]
        ax.hist(counts, bins=35, density=True, histtype="step", linewidth=0.7)
        mean = float(np.mean(counts))
        support = np.arange(max(0, int(np.min(counts))), int(np.max(counts)) + 1)
        ax.plot(support, poisson.pmf(support, mean), color="tab:red", linewidth=0.5)
        ax.set_title(f"ch {channel}: mean={mean:.0f}, F={np.var(counts)/mean:.1f}", fontsize=5)
        ax.tick_params(labelsize=4)
    for ax in axes.flat[N_CHANNELS:]:
        ax.axis("off")
    fig.suptitle(f"N_pe per channel, x={x_mm} mm")
    fig.tight_layout()
    npe_pdf.savefig(fig)
    plt.close(fig)

    centers = 0.5 * (TIME_BINS_NS[:-1] + TIME_BINS_NS[1:])
    fig, axes = plt.subplots(9, 10, figsize=(18, 15))
    for channel, ax in enumerate(axes.flat[:N_CHANNELS]):
        ax.step(centers, time_hist[channel], where="mid", linewidth=0.6)
        fpt = fpt_matrix[:, channel]
        fpt = fpt[np.isfinite(fpt)]
        ax.hist(fpt, bins=TIME_BINS_NS, histtype="step", color="tab:red", linewidth=0.5)
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.8)
        ax.set_title(f"ch {channel}", fontsize=5)
        ax.tick_params(labelsize=4)
    for ax in axes.flat[N_CHANNELS:]:
        ax.axis("off")
    fig.suptitle(f"All-photon time (blue) and FPT (red), x={x_mm} mm")
    fig.tight_layout()
    time_pdf.savefig(fig)
    plt.close(fig)


def analyze_file(path: pathlib.Path, x_mm: int, channel_pdfs: tuple[PdfPages, PdfPages] | None) -> PointResult:
    with uproot.open(path) as root_file:
        arrays = root_file["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns", "gun_x_mm"], library="np",
        )
    event_id = arrays["event_id"].astype(np.int32)
    global_id = arrays["global_id"].astype(np.int16)
    time_ns = arrays["time_ns"].astype(float)
    combined = event_id.astype(np.int64) * N_CHANNELS + global_id
    counts_matrix = np.bincount(
        combined, minlength=N_EVENTS * N_CHANNELS,
    ).reshape(N_EVENTS, N_CHANNELS)
    fpt_flat = np.full(N_EVENTS * N_CHANNELS, np.inf)
    np.minimum.at(fpt_flat, combined, time_ns)
    fpt_matrix = fpt_flat.reshape(N_EVENTS, N_CHANNELS)
    hit_count = np.bincount(global_id, minlength=N_CHANNELS)
    time_sum = np.bincount(global_id, weights=time_ns, minlength=N_CHANNELS)
    time_sum_sq = np.bincount(global_id, weights=time_ns * time_ns, minlength=N_CHANNELS)
    rows: list[dict] = []
    counts_by_name: dict[str, np.ndarray] = {}
    fpt_by_name: dict[str, np.ndarray] = {}
    top_profile = np.mean(counts_matrix[:, TOP_IDS], axis=0)
    top_sem = np.std(counts_matrix[:, TOP_IDS], axis=0) / math.sqrt(N_EVENTS)
    distances = np.abs(np.asarray(TOP_POSITIONS_MM) - x_mm)
    tied_nearest = np.flatnonzero(np.isclose(distances, np.min(distances)))
    nearest_local = int(tied_nearest[np.argmax(top_profile[tied_nearest])])
    ranked_locals = sorted(
        range(len(TOP_IDS)),
        key=lambda index: (
            0 if index == nearest_local else 1,
            distances[index],
            -top_profile[index],
            index,
        ),
    )
    nearest_ids = tuple(16 + index for index in ranked_locals[:3])

    groups: dict[str, tuple[int, ...]] = {
        **{f"ch_{channel:02d}": (channel,) for channel in range(N_CHANNELS)},
        **END_CLUSTERS,
        "end_left_all": END_LEFT_IDS,
        "end_right_all": END_RIGHT_IDS,
        "ends_all": END_LEFT_IDS + END_RIGHT_IDS,
        "top_left_all": TOP_LEFT_IDS,
        "top_right_all": TOP_RIGHT_IDS,
        "top_all": TOP_IDS,
    }
    groups["nearest_top"] = (nearest_ids[0],)
    groups["nearest_top_3"] = nearest_ids
    for name, ids in groups.items():
        row, counts, fpt = group_statistics(
            x_mm, name, ids, counts_matrix, hit_count, time_sum, time_sum_sq, fpt_matrix,
        )
        rows.append(row)
        counts_by_name[name] = counts
        fpt_by_name[name] = fpt

    cluster_times = {
        name: cluster_leading_edges(event_id, global_id, time_ns, ids)
        for name, ids in END_CLUSTERS.items()
    }
    sigma_delta: dict[str, float] = {}
    delta_by_side: dict[str, np.ndarray] = {}
    for side in ("left", "right"):
        a = cluster_times[f"end_{side}_A_SUM4"]
        b = cluster_times[f"end_{side}_B_SUM4"]
        valid = np.isfinite(a) & np.isfinite(b)
        delta = a[valid] - b[valid]
        delta_by_side[side] = delta
        sigma_delta[side] = float(np.std(delta) * 1000.0)
        for cluster in ("A", "B"):
            name = f"end_{side}_{cluster}_SUM4"
            next(row for row in rows if row["grupo"] == name)["sigma_grupo_ps"] = (
                sigma_delta[side] / math.sqrt(2.0)
            )
    left_a = cluster_times["end_left_A_SUM4"]
    right_a = cluster_times["end_right_A_SUM4"]
    valid = np.isfinite(left_a) & np.isfinite(right_a)
    tavg = 0.5 * (left_a[valid] + right_a[valid])
    sigma_tavg = float(np.std(tavg) * 1000.0)
    for row in rows:
        if row["grupo"] in {"end_left_all", "end_right_all", "ends_all"}:
            row["sigma_tavg_ps"] = sigma_tavg

    max_top_id = 16 + int(np.argmax(top_profile))

    representative_counts = {
        "nearest_top": counts_by_name["nearest_top"],
        "nearest_top_3": counts_by_name["nearest_top_3"],
        "nearest_end_all": counts_by_name[nearest_end_name(x_mm)],
    }
    nearest_cluster_names = (
        ("end_left_A_SUM4", "end_left_B_SUM4")
        if x_mm <= 0 else ("end_right_A_SUM4", "end_right_B_SUM4")
    )
    nearest_cluster = max(
        nearest_cluster_names, key=lambda name: np.mean(counts_by_name[name]),
    )
    representative_counts["nearest_end_SUM4"] = counts_by_name[nearest_cluster]

    representative_times: dict[str, np.ndarray] = {}
    representative_fpt: dict[str, np.ndarray] = {}
    if x_mm in KEY_TIME_POSITIONS:
        for name, ids in {
            "nearest_top": (nearest_ids[0],),
            "nearest_end_all": (
                END_LEFT_IDS if x_mm <= 0 else END_RIGHT_IDS
            ),
        }.items():
            representative_times[name] = time_ns[np.isin(global_id, ids)]
            representative_fpt[name] = np.min(fpt_matrix[:, ids], axis=1)

    if channel_pdfs is not None:
        time_hist, _, _ = np.histogram2d(
            global_id, time_ns,
            bins=(np.arange(-0.5, N_CHANNELS + 0.5, 1.0), TIME_BINS_NS),
        )
        save_channel_pages(channel_pdfs[0], channel_pdfs[1], x_mm, counts_matrix, time_hist, fpt_matrix)

    return PointResult(
        x_mm=x_mm,
        rows=rows,
        top_profile_mean=top_profile,
        top_profile_sem=top_sem,
        nearest_top_id=nearest_ids[0],
        nearest_top_3=nearest_ids,
        max_top_id=max_top_id,
        representative_counts=representative_counts,
        representative_times=representative_times,
        representative_fpt=representative_fpt,
        sigma_delta_left_ps=sigma_delta["left"],
        sigma_delta_right_ps=sigma_delta["right"],
        sigma_tavg_ps=sigma_tavg,
        delta_left_ns=delta_by_side["left"],
        delta_right_ns=delta_by_side["right"],
        tavg_ns=tavg,
    )


def validate_top_localization(points: list[PointResult], output_dir: pathlib.Path) -> None:
    """Apply the EXEC_08b physical-window exception to the localization gate."""
    rows = []
    failures = []
    for point in points:
        nearest = point.nearest_top_id
        maximum = point.max_top_id
        second = point.nearest_top_3[1]
        nearest_local = nearest - 16
        maximum_local = maximum - 16
        difference = point.top_profile_mean[maximum_local] - point.top_profile_mean[nearest_local]
        difference_error = math.hypot(
            point.top_profile_sem[maximum_local], point.top_profile_sem[nearest_local],
        )
        relative_difference = difference / point.top_profile_mean[nearest_local]
        window_track = point.x_mm % 20 == 10
        statistically_tied = difference <= difference_error
        allowed_window_dip = (
            window_track
            and maximum in {nearest, second}
            and relative_difference <= 0.15
        )
        passed = maximum == nearest or allowed_window_dip or (not window_track and statistically_tied)
        reason = (
            "strict maximum at nearest"
            if maximum == nearest else
            "known window-track dip <=15%"
            if allowed_window_dip else
            "strict statistical tie <=1 sigma"
            if not window_track and statistically_tied else
            "failed"
        )
        rows.append({
            "x_beam_mm": point.x_mm,
            "nearest_top_id": nearest,
            "second_nearest_top_id": second,
            "max_top_id": maximum,
            "max_minus_nearest_npe": difference,
            "difference_error_npe": difference_error,
            "relative_difference": relative_difference,
            "window_track_exception": window_track,
            "passed": passed,
            "reason": reason,
        })
        if not passed:
            failures.append((point.x_mm, nearest, second, maximum, relative_difference))
    pd.DataFrame(rows).to_csv(output_dir / "top_localization_gate.csv", index=False)
    if failures:
        raise RuntimeError(
            "blocking Top localization cross-check failed "
            "(x, nearest, second, max, relative difference): " + repr(failures)
        )


def fit_line(x: np.ndarray, y: np.ndarray, sigma: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    def line(value: np.ndarray, intercept: float, slope: float) -> np.ndarray:
        return intercept + slope * value
    return curve_fit(line, x, y, sigma=sigma, absolute_sigma=sigma is not None)


def add_fits(summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    fits: list[dict] = []
    for side, group in (("left", "end_left_all"), ("right", "end_right_all")):
        subset = summary[summary["grupo"] == group].sort_values("x_beam_mm").copy()
        distance_mm = (
            subset["x_beam_mm"].to_numpy() + BAR_HALF_LENGTH_MM
            if side == "left" else BAR_HALF_LENGTH_MM - subset["x_beam_mm"].to_numpy()
        )
        mean = subset["npe_mean"].to_numpy()
        sigma_log = subset["npe_rms"].to_numpy() / math.sqrt(N_EVENTS) / mean
        params, covariance = fit_line(distance_mm, np.log(mean), sigma_log)
        slope = params[1]
        slope_error = math.sqrt(covariance[1, 1])
        lambda_cm = -1.0 / slope / 10.0
        lambda_error_cm = abs(slope_error / slope * lambda_cm)
        summary.loc[summary["grupo"] == group, "lambda_eff_cm"] = lambda_cm
        time_params, time_cov = fit_line(distance_mm, subset["t_mean_ns"].to_numpy())
        time_slope = time_params[1]
        velocity_cm_ns = 1.0 / time_slope / 10.0
        velocity_error = abs(math.sqrt(time_cov[1, 1]) / time_slope * velocity_cm_ns)
        fits.append({
            "side": side,
            "lambda_eff_cm": lambda_cm,
            "lambda_eff_error_cm": lambda_error_cm,
            "v_eff_cm_ns": velocity_cm_ns,
            "v_eff_error_cm_ns": velocity_error,
        })
    return summary, pd.DataFrame(fits)


def save_geometry_figure(point: PointResult, path: pathlib.Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 3.2))
    ax.add_patch(plt.Rectangle((-700, -30), 1400, 60, facecolor="#dbeafe", edgecolor="#163b65", lw=2))
    top_x = np.asarray(TOP_POSITIONS_MM)
    ax.scatter(top_x, np.full_like(top_x, 35), marker="s", s=12, color="#777777", label="Top SiPM")
    for channel in point.nearest_top_3:
        x = TOP_POSITIONS_MM[channel - 16]
        color = "#d62728" if channel == point.nearest_top_id else "#ff9f1c"
        ax.scatter([x], [35], marker="s", s=60, color=color, zorder=4)
        ax.text(x, 43, str(channel), ha="center", va="bottom", fontsize=7, color=color)
    end_y = np.linspace(-24.5, 24.5, 8)
    ax.scatter(np.full(8, -706), end_y, marker="s", s=14, color="#1f77b4")
    ax.scatter(np.full(8, 706), end_y, marker="s", s=14, color="#1f77b4")
    ax.annotate("", xy=(point.x_mm, -28), xytext=(point.x_mm, 75),
                arrowprops={"arrowstyle": "->", "lw": 2.5, "color": "#d62728"})
    nearest_x = TOP_POSITIONS_MM[point.nearest_top_id - 16]
    subtitle = (
        f"x = {point.x_mm} mm - top nearest: ID {point.nearest_top_id} "
        f"(dx = {abs(nearest_x - point.x_mm):.0f} mm) - "
        f"distance to ends: L {point.x_mm + 700:.0f} mm / R {700 - point.x_mm:.0f} mm"
    )
    ax.set_title(subtitle, fontsize=10)
    ax.text(0, -7, "EJ-204 OPSC-101, Mylar wrapped except End faces and 70 Top windows",
            ha="center", va="center", fontsize=8)
    ax.set(xlim=(-750, 750), ylim=(-50, 85), xlabel="x [mm]")
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def save_position_profile(point: PointResult, summary: pd.DataFrame, path: pathlib.Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    channels = np.arange(16, 86)
    ax.errorbar(channels, point.top_profile_mean, yerr=point.top_profile_sem, fmt=".-", lw=0.8, ms=3)
    ax.axvline(point.nearest_top_id, color="tab:red", ls="--", label=f"nearest ID {point.nearest_top_id}")
    ax.set(xlabel="Top global channel ID", ylabel="Mean N_pe / event", title=f"Top profile, x={point.x_mm} mm")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_integrated(
    points: list[PointResult], summary: pd.DataFrame, fits: pd.DataFrame, figs: pathlib.Path,
) -> None:
    x_values = np.array([point.x_mm for point in points])
    by_group = {name: summary[summary["grupo"] == name].sort_values("x_beam_mm")
                for name in summary["grupo"].unique()}

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for name, label in [
        ("end_left_all", "End L"), ("end_right_all", "End R"),
        ("ends_all", "End L+R"), ("nearest_top", "Top nearest"), ("top_all", "Top sum"),
    ]:
        frame = by_group[name]
        sem = frame["npe_rms"].to_numpy() / math.sqrt(N_EVENTS)
        ax.errorbar(
            frame["x_beam_mm"].to_numpy(), frame["npe_mean"].to_numpy(),
            yerr=sem, marker="o", ms=3, lw=1, label=label,
        )
    ax.set(xlabel="Beam x [mm]", ylabel="Mean N_pe / event")
    ax.grid(alpha=0.3)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(figs / "P1_npe_vs_x.png", dpi=300)
    plt.close(fig)

    heatmap = np.vstack([point.top_profile_mean for point in points]).T
    fig, ax = plt.subplots(figsize=(10, 7))
    image = ax.imshow(heatmap, aspect="auto", origin="lower",
                      extent=(x_values[0], x_values[-1], 15.5, 85.5), cmap="viridis")
    ax.set(xlabel="Beam x [mm]", ylabel="Top global channel ID", title="Mean Top N_pe")
    fig.colorbar(image, ax=ax, label="Mean N_pe / event")
    fig.tight_layout()
    fig.savefig(figs / "P2_npe_heatmap_top.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for name in ("nearest_top", *END_CLUSTERS.keys(), "end_left_all", "end_right_all"):
        frame = by_group[name]
        ax.plot(
            frame["x_beam_mm"].to_numpy(), frame["var_over_mean"].to_numpy(),
            marker="o", ms=2.5, lw=0.9, label=name,
        )
    ax.axhline(1.0, color="black", ls="--", lw=1, label="Poisson")
    ax.set(xlabel="Beam x [mm]", ylabel="var(N_pe) / mean(N_pe)", yscale="log")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(figs / "P3_fano_vs_x.png", dpi=300)
    plt.close(fig)

    for point in points:
        if point.x_mm not in KEY_POISSON_POSITIONS:
            continue
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        for ax, name in zip(axes, ("nearest_top", "nearest_end_SUM4", "nearest_end_all")):
            counts = point.representative_counts[name]
            mean = float(np.mean(counts))
            sigma = float(np.std(counts))
            low, high = int(np.min(counts)), int(np.max(counts))
            bins = np.arange(low - 0.5, high + 1.5)
            ax.hist(counts, bins=bins, density=True, histtype="step", lw=1, label="data")
            support = np.arange(low, high + 1)
            ax.plot(support, poisson.pmf(support, mean), lw=1, label="Poisson")
            ax.plot(support, norm.pdf(support, mean, sigma), lw=1, label="Gaussian")
            ax.set_yscale("log")
            ax.set_ylim(bottom=1e-6)
            ax.set_title(
                f"{name}\nF={np.var(counts)/mean:.2f}, "
                f"chi2/ndf P={poisson_chi2_ndf(counts):.1f}, G={gaussian_chi2_ndf(counts):.1f}",
                fontsize=8,
            )
            ax.legend(fontsize=7)
        fig.suptitle(f"Poisson tail check, x={point.x_mm} mm")
        fig.tight_layout()
        fig.savefig(figs / f"P4_poisson_check_x{point.x_mm}.png", dpi=300)
        plt.close(fig)

    fig, axes = plt.subplots(len(KEY_TIME_POSITIONS), 2, figsize=(12, 11), sharex=True)
    for row_index, x_mm in enumerate(KEY_TIME_POSITIONS):
        point = next(point for point in points if point.x_mm == x_mm)
        for col_index, name in enumerate(("nearest_top", "nearest_end_all")):
            ax = axes[row_index, col_index]
            ax.hist(point.representative_times[name], bins=TIME_BINS_NS, histtype="step", density=True, label="all photons")
            fpt = point.representative_fpt[name]
            ax.hist(fpt[np.isfinite(fpt)], bins=TIME_BINS_NS, histtype="step", density=True, label="FPT")
            ax.set_yscale("log")
            ax.set_ylim(bottom=1e-5)
            ax.set_title(f"x={x_mm} mm, {name}")
            ax.legend(fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("Arrival time [ns]")
    fig.tight_layout()
    fig.savefig(figs / "P5_tdist_examples.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for name in ("end_left_all", "end_right_all", "nearest_top"):
        frame = by_group[name]
        axes[0].plot(
            frame["x_beam_mm"].to_numpy(), frame["t_mean_ns"].to_numpy(),
            marker="o", ms=3, label=name,
        )
        axes[1].plot(
            frame["x_beam_mm"].to_numpy(), frame["fpt_mean_ns"].to_numpy(),
            marker="o", ms=3, label=name,
        )
    axes[0].set(ylabel="Mean all-photon arrival time [ns]", xlabel="Beam x [mm]")
    axes[1].set(ylabel="Mean FPT [ns]", xlabel="Beam x [mm]")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figs / "P6_tmean_vs_x.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].plot(x_values, [p.sigma_delta_left_ps for p in points], marker="o", label="End L sigma(DeltaT)")
    axes[0].plot(x_values, [p.sigma_delta_right_ps for p in points], marker="o", label="End R sigma(DeltaT)")
    axes[0].axhline(125.0, color="tab:red", ls="--", label="hardware sigma(DeltaT)=125 ps")
    axes[0].set(xlabel="Beam x [mm]", ylabel="sigma(DeltaT) [ps]")
    axes[1].plot(x_values, [p.sigma_tavg_ps for p in points], marker="o", label="sigma((T1_LA+T1_RA)/2)")
    axes[1].set(xlabel="Beam x [mm]", ylabel="sigma(t_avg) [ps]")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Intrinsic End timing, no SPTR/electronics jitter")
    fig.tight_layout()
    fig.savefig(figs / "P7_deltaT_end.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for side, group, color in (("left", "end_left_all", "tab:blue"), ("right", "end_right_all", "tab:orange")):
        frame = by_group[group]
        distance = frame["x_beam_mm"].to_numpy() + 700 if side == "left" else 700 - frame["x_beam_mm"].to_numpy()
        fit = fits[fits["side"] == side].iloc[0]
        npe_mean = frame["npe_mean"].to_numpy()
        time_mean = frame["t_mean_ns"].to_numpy()
        axes[0].scatter(distance, np.log(npe_mean), s=18, label=f"{side}: lambda={fit.lambda_eff_cm:.1f}+/-{fit.lambda_eff_error_cm:.1f} cm", color=color)
        coefficients = np.polyfit(distance, np.log(npe_mean), 1)
        axes[0].plot(distance, np.polyval(coefficients, distance), color=color)
        axes[1].scatter(distance, time_mean, s=18, label=f"{side}: v_eff={fit.v_eff_cm_ns:.1f}+/-{fit.v_eff_error_cm_ns:.1f} cm/ns", color=color)
        coefficients = np.polyfit(distance, time_mean, 1)
        axes[1].plot(distance, np.polyval(coefficients, distance), color=color)
    axes[0].set(xlabel="Track-to-end distance [mm]", ylabel="ln(mean End N_pe)")
    axes[1].set(xlabel="Track-to-end distance [mm]", ylabel="Mean all-photon arrival time [ns]")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figs / "fits_attenuation_velocity.png", dpi=300)
    plt.close(fig)


def write_reports(
    points: list[PointResult], summary: pd.DataFrame, fits: pd.DataFrame,
    output_dir: pathlib.Path, figs: pathlib.Path,
) -> None:
    diagnostics = []
    for point in points:
        for name, counts in point.representative_counts.items():
            diagnostics.append({
                "x_beam_mm": point.x_mm,
                "group": name,
                "var_over_mean": np.var(counts) / np.mean(counts),
                "chi2ndf_poisson": poisson_chi2_ndf(counts),
                "chi2ndf_gaussian": gaussian_chi2_ndf(counts),
            })
    pd.DataFrame(diagnostics).to_csv(output_dir / "poisson_diagnostics.csv", index=False)

    with PdfPages(output_dir / "end_deltaT_distributions.pdf") as pdf:
        for point in points:
            fig, axes = plt.subplots(1, 3, figsize=(14, 4))
            for ax, values, title in (
                (axes[0], point.delta_left_ns, "End L: T1_A - T2_B"),
                (axes[1], point.delta_right_ns, "End R: T1_A - T2_B"),
                (axes[2], point.tavg_ns, "(T1_L,A + T1_R,A) / 2"),
            ):
                ax.hist(values, bins=60, histtype="step")
                ax.set_title(f"{title}\nsigma={np.std(values)*1000.0:.2f} ps", fontsize=9)
                ax.set_xlabel("Time [ns]")
                ax.grid(alpha=0.2)
            fig.suptitle(f"Intrinsic End timing distributions, x={point.x_mm} mm")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    per_position_rows = []
    for point in points:
        values = summary[summary["x_beam_mm"] == point.x_mm].set_index("grupo")
        per_position_rows.append({
            "x_beam_mm": point.x_mm,
            "nearest_top_id": point.nearest_top_id,
            "nearest_top_3": " ".join(map(str, point.nearest_top_3)),
            "max_top_id": point.max_top_id,
            "end_left_npe": values.loc["end_left_all", "npe_mean"],
            "end_right_npe": values.loc["end_right_all", "npe_mean"],
            "nearest_top_npe": values.loc["nearest_top", "npe_mean"],
            "nearest_top_t_mean_ns": values.loc["nearest_top", "t_mean_ns"],
            "nearest_top_fpt_mean_ns": values.loc["nearest_top", "fpt_mean_ns"],
        })
    pd.DataFrame(per_position_rows).to_csv(output_dir / "per_position_exec07.csv", index=False)

    sigma_group = np.array([
        point.sigma_delta_left_ps / math.sqrt(2.0) for point in points
    ] + [
        point.sigma_delta_right_ps / math.sqrt(2.0) for point in points
    ])
    sigma_tavg = np.array([point.sigma_tavg_ps for point in points])
    nearest_top = summary[summary["grupo"] == "nearest_top"]
    nearest_top_3 = summary[summary["grupo"] == "nearest_top_3"]
    dead_channels = [
        channel for channel in range(N_CHANNELS)
        if summary[summary["grupo"] == f"ch_{channel:02d}"]["npe_mean"].max() <= 0.0
    ]
    left = summary[summary["grupo"] == "end_left_all"].set_index("x_beam_mm")
    right = summary[summary["grupo"] == "end_right_all"].set_index("x_beam_mm")
    symmetry = []
    for x_mm in EXPECTED_POSITIONS_MM:
        left_value = left.loc[x_mm, "npe_mean"]
        right_value = right.loc[-x_mm, "npe_mean"]
        symmetry.append(abs(left_value - right_value) / (0.5 * (left_value + right_value)))
    symmetry = np.asarray(symmetry)
    conclusions = [
        "# EXEC_08 photon-budget conclusions",
        "",
        "Input validation: 31/31 files passed; each has event_id 0..1999, matching gun_x_mm, and IDs 0..85.",
        "",
        "## Numerical conclusions",
        "",
        f"- Effective attenuation: left {fits.iloc[0].lambda_eff_cm:.2f} +/- {fits.iloc[0].lambda_eff_error_cm:.2f} cm; right {fits.iloc[1].lambda_eff_cm:.2f} +/- {fits.iloc[1].lambda_eff_error_cm:.2f} cm. These include light extraction through 70 Top windows and are not expected to equal the 160 cm bulk ABSLENGTH.",
        f"- Effective propagation velocity from mean all-photon time: left {fits.iloc[0].v_eff_cm_ns:.2f} +/- {fits.iloc[0].v_eff_error_cm_ns:.2f} cm/ns; right {fits.iloc[1].v_eff_cm_ns:.2f} +/- {fits.iloc[1].v_eff_error_cm_ns:.2f} cm/ns.",
        f"- Typical nearest-Top var/mean: median {nearest_top.var_over_mean.median():.2f}, range {nearest_top.var_over_mean.min():.2f}-{nearest_top.var_over_mean.max():.2f}. N_pe is not forced to Poisson.",
        f"- Intrinsic End SUM4 sigma_group= sigma(DeltaT)/sqrt(2): mean {np.mean(sigma_group):.2f} ps, range {np.min(sigma_group):.2f}-{np.max(sigma_group):.2f} ps. No SPTR/electronics jitter.",
        f"- sigma(t_avg) mean {np.mean(sigma_tavg):.2f} ps, range {np.min(sigma_tavg):.2f}-{np.max(sigma_tavg):.2f} ps.",
        f"- Top analytic estimate, nearest channel: {nearest_top.sigma_est_ps.min():.2f}-{nearest_top.sigma_est_ps.max():.2f} ps; nearest-three group: {nearest_top_3.sigma_est_ps.min():.2f}-{nearest_top_3.sigma_est_ps.max():.2f} ps. This is an analytic estimate, not simulated sigma_t.",
        f"- Dead channels: {dead_channels if dead_channels else 'none'}. Mirrored End L/R asymmetry above 5%: {int(np.count_nonzero(symmetry > 0.05))}/31 positions; maximum {100.0*np.max(symmetry):.2f}%.",
        "",
        "## Method notes",
        "",
        f"- EJ-204 estimate uses sqrt(tau_r*tau_d)={SIGMA_NUMERATOR_NS:.3f} ns from tau_r={TAU_R_NS} ns and tau_d={TAU_D_NS} ns.",
        f"- SUM4 leading edge is ported from `analysis/congruent_sum4_timing.C`: normalized double exponential {SPR_RISE_NS}/{SPR_FALL_NS} ns, absolute threshold {LEADING_EDGE_THRESHOLD_PE} PE, no time-walk correction.",
        "- Gaussian overlays are diagnostic only. Poisson-tail plots use log-y and report both Poisson and Gaussian chi2/ndf.",
    ]
    (output_dir / "conclusions_exec07.md").write_text("\n".join(conclusions) + "\n", encoding="utf-8")

    with PdfPages(output_dir / "exec07_photon_budget_report.pdf") as pdf:
        fig = plt.figure(figsize=(11.69, 8.27))
        fig.text(0.06, 0.94, "EXEC_08 - EXEC_07 EndTop photon budget", fontsize=20, weight="bold")
        fig.text(0.06, 0.88, "31 positions x 2000 events, 86 channels, EJ-204 OPSC-101, Broadcom AFBR-S4N66P024M", fontsize=11)
        fig.text(0.06, 0.82, "Intrinsic simulation: /sipm/jitterSigma 0", fontsize=11, color="tab:red")
        fig.text(0.06, 0.73, "\n".join(conclusions[5:11]), fontsize=10, va="top")
        fig.text(0.06, 0.18, "Top uses 70 simulated SiPM elements and does not replicate the 32-SiPM hardware setup.", fontsize=10)
        pdf.savefig(fig)
        plt.close(fig)
        for filename in [
            "P1_npe_vs_x.png", "P2_npe_heatmap_top.png", "P3_fano_vs_x.png",
            "P4_poisson_check_x0.png", "P5_tdist_examples.png", "P6_tmean_vs_x.png",
            "P7_deltaT_end.png", "fits_attenuation_velocity.png",
        ]:
            image = plt.imread(figs / filename)
            fig, ax = plt.subplots(figsize=(11.69, 8.27))
            ax.imshow(image)
            ax.axis("off")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figs = args.output_dir / "figs"
    figs.mkdir(exist_ok=True)
    files = validate_inputs(args.data_dir)
    print("blocking validation passed: 31 files x 2000 events, gun_x and IDs 0..85")
    (args.output_dir / "validation_exec07.txt").write_text(
        "PASS: 31 files; each has event_id 0..1999, matching gun_x_mm, "
        "and aggregate global_id 0..85.\n",
        encoding="utf-8",
    )

    points: list[PointResult] = []
    if args.skip_channel_pdfs:
        for x_mm, path in zip(EXPECTED_POSITIONS_MM, files):
            print(f"analyzing x={x_mm} mm: {path.name}", flush=True)
            points.append(analyze_file(path, x_mm, None))
    else:
        with PdfPages(args.output_dir / "npe_channels.pdf") as npe_pdf, \
             PdfPages(args.output_dir / "arrival_time_channels.pdf") as time_pdf:
            for x_mm, path in zip(EXPECTED_POSITIONS_MM, files):
                print(f"analyzing x={x_mm} mm: {path.name}", flush=True)
                points.append(analyze_file(path, x_mm, (npe_pdf, time_pdf)))

    validate_top_localization(points, args.output_dir)
    summary = pd.DataFrame([row for point in points for row in point.rows])
    summary, fits = add_fits(summary)
    summary[SUMMARY_COLUMNS].to_csv(args.output_dir / "summary_exec07.csv", index=False)
    fits.to_csv(args.output_dir / "fit_results_exec07.csv", index=False)

    for point in points:
        save_geometry_figure(point, figs / f"muon_{point.x_mm}mm_geometry.png")
        save_position_profile(point, summary, figs / f"muon_{point.x_mm}mm_top_profile.png")
    plot_integrated(points, summary, fits, figs)
    write_reports(points, summary, fits, args.output_dir, figs)

    velocities = fits["v_eff_cm_ns"].to_numpy()
    if np.any(np.abs(velocities - 19.0) / 19.0 > 0.15):
        print(
            "WARNING: fitted mean-arrival v_eff is outside 19 cm/ns +/-15%; "
            "reported without forcing the expected value.",
            file=sys.stderr,
        )
    print(f"wrote EXEC_08 products to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
