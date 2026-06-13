#!/usr/bin/env python3
"""EXEC_13 (EJ-230) — time-to-threshold t_N analysis.

Mirror of exec12_tN_analysis.py, adapted for EJ-230 data (OPSC-106).
Observable, groups, thresholds, and consistency criteria are unchanged.
Only the data path, output labels, and EJ-230 scintillator constants differ.

Output (all under RESULTS_DIR):
  figures/exec13_tN_{x}mm.png   -- 6-panel figure per key position
  figures/exec13_tN_summary.png  -- sigma_fit vs x
  csv/exec13_tN_summary.csv      -- one row per (x, group, N)
  logs/exec13_tN_check.md        -- consistency check report
"""

from __future__ import annotations

import argparse
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from analysis.exec13.common13 import (  # noqa: E402
    END_CLUSTERS,
    N_EVENTS,
    TOP_IDS,
    TOP_POSITIONS_MM,
    KEY_POSITIONS,
    N_THRESHOLDS,
    results_dir,
    expected_file,
)

# HOOK_ADAPTIVE_T4_FIT: robust histogram and core-fit controls.
FIT_IQR_BIN_SCALE = 2.0
FIT_MIN_BIN_WIDTH_NS = 0.025
FIT_MIN_BINS = 6
FIT_MAX_BINS = 160
FIT_INITIAL_CORE_SIGMA = 2.0
FIT_FINAL_CORE_SIGMA = 2.5
_N_TOP = len(TOP_IDS)


def _gauss(x: np.ndarray, amp: float, mu: float, sig: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def compute_tN(
    event_id: np.ndarray,
    global_id: np.ndarray,
    times: np.ndarray,
    channel_ids: tuple[int, ...],
    N: int,
) -> tuple[np.ndarray, float]:
    mask = np.isin(global_id, channel_ids)
    ev = event_id[mask].astype(int)
    tt = times[mask]
    order = np.lexsort((tt, ev))
    ev_s = ev[order]
    tt_s = tt[order]
    starts = np.searchsorted(ev_s, np.arange(N_EVENTS), side="left")
    stops = np.searchsorted(ev_s, np.arange(N_EVENTS), side="right")
    tN = np.full(N_EVENTS, np.nan)
    have = stops - starts >= N
    tN[have] = tt_s[starts[have] + N - 1]
    frac_excl = float(np.sum(~have)) / N_EVENTS
    return tN, frac_excl


def robust_sigma(values: np.ndarray) -> float:
    q25, q75 = np.quantile(values, [0.25, 0.75])
    sigma = float((q75 - q25) / 1.3489795003921634)
    if not math.isfinite(sigma) or sigma <= 0.0:
        sigma = float(np.std(values, ddof=1))
    return sigma


def adaptive_edges(
    values: np.ndarray, minimum_width_ns: float = 0.0, use_full_range: bool = False,
) -> np.ndarray:
    """Freedman-Diaconis binning constrained by documented bin-count bounds."""
    center = float(np.median(values))
    spread = robust_sigma(values)
    if use_full_range:
        lo, hi = float(np.min(values)), float(np.max(values))
    else:
        lo = max(float(np.min(values)), center - 4.0 * spread)
        hi = min(float(np.max(values)), center + 4.0 * spread)
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return np.linspace(center - 0.5, center + 0.5, FIT_MIN_BINS + 1)
    iqr = spread * 1.3489795003921634
    width = max(
        minimum_width_ns,
        FIT_IQR_BIN_SCALE * iqr / np.cbrt(max(len(values), 1)),
    )
    requested = int(math.ceil((hi - lo) / width)) if width > 0.0 else FIT_MIN_BINS
    n_bins = min(FIT_MAX_BINS, max(1, requested))
    return np.linspace(lo, hi, n_bins + 1)


def fit_histogram(
    values: np.ndarray,
    mu0: float,
    sigma0: float,
    minimum_width_ns: float,
    range_bounds: tuple[float, float] | None = None,
) -> tuple[float, float, float, float]:
    if range_bounds is not None:
        lo, hi = range_bounds
        requested = int(math.ceil((hi - lo) / minimum_width_ns))
        edges = np.linspace(lo, hi, max(requested, 1) + 1)
    else:
        edges = adaptive_edges(values, minimum_width_ns, use_full_range=minimum_width_ns > 0.0)
    counts, _ = np.histogram(values, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    valid = np.ones_like(counts, dtype=bool) if minimum_width_ns > 0.0 else counts > 0
    if np.count_nonzero(valid) < 3:
        raise RuntimeError("fewer than three occupied adaptive bins")
    bin_width = float(edges[1] - edges[0])
    popt, pcov = curve_fit(
        _gauss,
        centers[valid],
        counts[valid].astype(float),
        p0=[float(counts.max()), mu0, sigma0],
        bounds=(
            [0.0, float(edges[0]), max(bin_width / 10.0, np.finfo(float).eps)],
            [math.inf, float(edges[-1]), float(edges[-1] - edges[0])],
        ),
        maxfev=20000,
    )
    sigma = abs(float(popt[2]))
    variance = float(pcov[2, 2])
    error = math.sqrt(variance) if variance >= 0.0 else math.nan
    counting_error = sigma / math.sqrt(2.0 * max(len(values) - 1, 1))
    if not math.isfinite(error) or error > sigma:
        error = counting_error
    return float(popt[1]), sigma, error, bin_width


def gaussian_core_fit(
    tN: np.ndarray,
) -> tuple[float, float, float, float, float, bool, str]:
    """Adaptive Gaussian core fit, stable from few-ps to broad t_N cores."""
    valid = tN[np.isfinite(tN)]
    n = len(valid)
    if n < 20:
        m = float(np.mean(valid)) if n > 0 else math.nan
        r = float(np.std(valid, ddof=1)) if n > 1 else math.nan
        return m, r, math.nan, math.nan, math.nan, False, "insufficient"
    mean0 = float(np.mean(valid))
    rms0 = float(np.std(valid, ddof=1))
    median0 = float(np.median(valid))
    robust0 = robust_sigma(valid)
    lo0 = mean0 - FIT_INITIAL_CORE_SIGMA * rms0
    hi0 = mean0 + FIT_INITIAL_CORE_SIGMA * rms0
    core0 = valid[(valid >= lo0) & (valid <= hi0)]
    try:
        mu1, sig1, err1, bin1 = fit_histogram(
            core0, mean0, rms0, FIT_MIN_BIN_WIDTH_NS, (lo0, hi0)
        )
    except Exception:
        try:
            mu1, sig1, err1, bin1 = fit_histogram(core0, median0, robust0, 0.0)
        except Exception:
            return mean0, rms0, math.nan, math.nan, math.nan, False, "failed-pass-one"
        return mean0, rms0, sig1, err1, bin1, True, "resolution-limited-adaptive"
    lo1 = mu1 - FIT_FINAL_CORE_SIGMA * sig1
    hi1 = mu1 + FIT_FINAL_CORE_SIGMA * sig1
    core1 = valid[(valid >= lo1) & (valid <= hi1)]
    try:
        _, sigma_fit, sigma_fit_err, bin_width = fit_histogram(
            core1, mu1, sig1, FIT_MIN_BIN_WIDTH_NS, (lo1, hi1)
        )
    except Exception:
        limited = bool(sig1 <= FIT_MIN_BIN_WIDTH_NS)
        return mean0, rms0, sig1, err1, bin1, limited, "adaptive-histogram-pass-one"
    limited = bool(sigma_fit <= bin_width)
    return mean0, rms0, sigma_fit, sigma_fit_err, bin_width, limited, "adaptive-histogram"


def _top_profile(event_id: np.ndarray, global_id: np.ndarray) -> np.ndarray:
    top_mask = (global_id >= 16) & (global_id <= 85)
    top_local = global_id[top_mask].astype(np.int64) - 16
    ev_top = event_id[top_mask].astype(np.int64)
    flat = np.bincount(ev_top * _N_TOP + top_local, minlength=N_EVENTS * _N_TOP)
    return flat.reshape(N_EVENTS, _N_TOP).mean(axis=0)


def _select_nearest_top(x_mm: int, profile: np.ndarray) -> int:
    distances = np.abs(np.asarray(TOP_POSITIONS_MM) - x_mm)
    tied = np.flatnonzero(np.isclose(distances, distances.min()))
    local = int(tied[int(np.argmax(profile[tied]))])
    return 16 + local


def _dominant_cluster(
    global_id: np.ndarray,
    cluster_a: tuple[int, ...],
    cluster_b: tuple[int, ...],
) -> tuple[tuple[int, ...], str]:
    npe_a = float(np.sum(np.isin(global_id, cluster_a))) / N_EVENTS
    npe_b = float(np.sum(np.isin(global_id, cluster_b))) / N_EVENTS
    if npe_a >= npe_b:
        return cluster_a, "A"
    return cluster_b, "B"


def analyze_position(x_mm: int, data_dir: pathlib.Path) -> tuple:
    path = expected_file(data_dir, x_mm)
    with uproot.open(path) as rf:
        arrays = rf["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns"], library="np")
    event_id = arrays["event_id"].astype(int)
    global_id = arrays["global_id"].astype(int)
    times = arrays["time_ns"].astype(float)

    prof = _top_profile(event_id, global_id)
    nearest_top_id = _select_nearest_top(x_mm, prof)
    left_ids, _ = _dominant_cluster(global_id,
                                    END_CLUSTERS["end_left_A_SUM4"],
                                    END_CLUSTERS["end_left_B_SUM4"])
    right_ids, _ = _dominant_cluster(global_id,
                                     END_CLUSTERS["end_right_A_SUM4"],
                                     END_CLUSTERS["end_right_B_SUM4"])
    dist_left = float(x_mm + 700)
    dist_right = float(700 - x_mm)
    if dist_left <= dist_right:
        near_ids, near_label = left_ids, "End-L"
        far_ids, far_label = right_ids, "End-R"
    else:
        near_ids, near_label = right_ids, "End-R"
        far_ids, far_label = left_ids, "End-L"
    groups = [
        ("top_nearest", (nearest_top_id,),
         f"Top nearest ID {nearest_top_id}"),
        ("end_near", near_ids,
         f"{near_label} SUM4 {{{','.join(map(str, near_ids))}}}"),
        ("end_far", far_ids,
         f"{far_label} SUM4 {{{','.join(map(str, far_ids))}}}"),
    ]

    rows: list[dict] = []
    panels: dict = {}

    for N in N_THRESHOLDS:
        for grp_key, ch_ids, label in groups:
            tN, frac_excl = compute_tN(event_id, global_id, times, ch_ids, N)
            valid = tN[np.isfinite(tN)]
            mean_v = float(np.mean(valid)) if len(valid) > 0 else math.nan
            rms_v = float(np.std(valid, ddof=1)) if len(valid) > 1 else math.nan
            if frac_excl > 0.5:
                sigma_fit = sigma_fit_err = bin_width = math.nan
                resolution_limited = False
                fit_method = "reach-below-50-percent"
            else:
                (
                    _, _, sigma_fit, sigma_fit_err, bin_width,
                    resolution_limited, fit_method,
                ) = gaussian_core_fit(tN)
            rows.append({
                "x_mm": x_mm, "group": grp_key, "N": N,
                "n_events": N_EVENTS, "frac_excluded": round(frac_excl, 4),
                "mean_ns": round(mean_v, 5) if math.isfinite(mean_v) else float("nan"),
                "rms_ns": round(rms_v, 5) if math.isfinite(rms_v) else float("nan"),
                "sigma_fit_ns": round(sigma_fit, 6) if math.isfinite(sigma_fit) else float("nan"),
                "sigma_fit_err_ns": round(sigma_fit_err, 6) if math.isfinite(sigma_fit_err) else float("nan"),
                "fit_bin_width_ns": round(bin_width, 8) if math.isfinite(bin_width) else float("nan"),
                "resolution_limited": resolution_limited, "fit_method": fit_method,
                "channels": str(ch_ids), "near_label": near_label, "far_label": far_label,
            })
            panels[(N, grp_key)] = (
                tN, label, frac_excl, sigma_fit, sigma_fit_err, bin_width,
                resolution_limited,
            )

    return rows, panels, nearest_top_id, near_label, far_label


def plot_position(x_mm: int, panels: dict, near_label: str, far_label: str,
                  figs_dir: pathlib.Path) -> None:
    col_titles = ["Top nearest", f"{near_label} SUM4", f"{far_label} SUM4"]
    row_labels = {4: r"$t_4$  [N=4 PE threshold]", 20: r"$t_{20}$  [N=20 PE]"}
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))
    for row_i, N in enumerate(N_THRESHOLDS):
        for col_i, grp_key in enumerate(("top_nearest", "end_near", "end_far")):
            ax = axes[row_i, col_i]
            (
                tN, label, frac_excl, sigma_fit, sigma_fit_err, bin_width,
                resolution_limited,
            ) = panels[(N, grp_key)]
            valid = tN[np.isfinite(tN)]
            ax.set_xlabel("Arrival time [ns]", fontsize=8)
            if col_i == 0:
                ax.set_ylabel("Events / adaptive bin", fontsize=8)
            if row_i == 0:
                ax.set_title(col_titles[col_i], fontsize=9)
            if frac_excl > 0.5:
                ax.text(0.5, 0.5,
                        f"excluded {100*frac_excl:.0f}% of events\n"
                        f"(fewer than {N} PE in this group)",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=9, color="0.45")
                ax.text(0.03, 0.97, row_labels[N], transform=ax.transAxes,
                        ha="left", va="top", fontsize=7)
                continue
            if len(valid) < 5:
                ax.text(0.5, 0.5, "insufficient data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9, color="0.45")
                continue
            mean_v = float(np.mean(valid))
            rms_v = float(np.std(valid, ddof=1))
            lo = max(0.0, mean_v - 4.0 * rms_v)
            hi = mean_v + 4.0 * rms_v
            edges = adaptive_edges(valid)
            counts, edges = np.histogram(valid, bins=edges)
            centers = 0.5 * (edges[:-1] + edges[1:])
            ax.step(centers, counts, where="mid", lw=1.3, color="tab:blue")
            if math.isfinite(sigma_fit) and math.isfinite(sigma_fit_err):
                x_fit = np.linspace(lo, hi, 400)
                amp_est = float(counts.max())
                ax.plot(x_fit, _gauss(x_fit, amp_est, mean_v, sigma_fit),
                        "r--", lw=1.1,
                        label=(rf"$\sigma_{{fit}}={sigma_fit*1000:.0f}"
                               rf"\pm{sigma_fit_err*1000:.0f}$ ps"))
                ax.legend(fontsize=7, loc="upper right")
            if resolution_limited:
                ax.text(
                    0.03, 0.82, "resolution-limited core",
                    transform=ax.transAxes, ha="left", va="top", fontsize=7,
                    color="darkorange",
                )
            excl_str = f"{100*frac_excl:.1f}%" if frac_excl > 0 else "0%"
            ax.text(0.97, 0.97,
                    f"mean = {mean_v*1000:.0f} ps\n"
                    f"RMS  = {rms_v*1000:.0f} ps\n"
                    f"bin = {bin_width*1000:.1f} ps\n"
                    f"excl. = {excl_str}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=7,
                    bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.75"})
            ax.text(0.03, 0.97, row_labels[N], transform=ax.transAxes,
                    ha="left", va="top", fontsize=7)
            ax.grid(alpha=0.22)
    fig.suptitle(rf"EJ-230 $t_N$ $|$ $x = {x_mm:+d}$ mm", fontsize=11)
    fig.tight_layout()
    out = figs_dir / f"exec13_tN_{x_mm}mm.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


def plot_synthesis(summary: pd.DataFrame, figs_dir: pathlib.Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    colors = {"top_nearest": "tab:green", "end_near": "tab:blue", "end_far": "tab:red"}
    markers = {"top_nearest": "^", "end_near": "o", "end_far": "s"}
    labels = {"top_nearest": "Top nearest", "end_near": "End near", "end_far": "End far"}
    for col_i, N in enumerate(N_THRESHOLDS):
        ax = axes[col_i]
        sub_N = summary[summary["N"] == N]
        for grp_key in ("top_nearest", "end_near", "end_far"):
            sub = sub_N[sub_N["group"] == grp_key].copy()
            x_v = sub["x_mm"].to_numpy()
            y_v = sub["sigma_fit_ns"].to_numpy() * 1000.0
            y_e = sub["sigma_fit_err_ns"].to_numpy() * 1000.0
            good = np.isfinite(y_v) & np.isfinite(y_e)
            if not np.any(good):
                good = np.isfinite(y_v)
            if np.any(good):
                ax.errorbar(x_v[good], y_v[good],
                            yerr=y_e[good] if np.any(np.isfinite(y_e[good])) else None,
                            fmt=markers[grp_key], color=colors[grp_key], ms=6,
                            capsize=3, label=labels[grp_key])
        ax.set_xlabel("Beam position x [mm]", fontsize=9)
        ax.set_ylabel(r"$\sigma(t_N)$ [ps]", fontsize=9)
        ax.set_title(rf"EJ-230 — $N = {N}$ PE threshold", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.suptitle(r"EJ-230 $\sigma(t_N)$ vs beam position (EXEC_13)", fontsize=10)
    fig.tight_layout()
    out = figs_dir / "exec13_tN_summary.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="EXEC_13: EJ-230 t_N analysis.")
    rdir = results_dir()
    parser.add_argument("--data-dir", type=pathlib.Path, default=rdir / "data")
    parser.add_argument("--figs-dir", type=pathlib.Path, default=rdir / "figures")
    parser.add_argument("--csv-dir",  type=pathlib.Path, default=rdir / "csv")
    args = parser.parse_args()

    args.figs_dir.mkdir(parents=True, exist_ok=True)
    args.csv_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []

    for x_mm in KEY_POSITIONS:
        path = expected_file(args.data_dir, x_mm)
        if not path.is_file():
            print(f"SKIP: {path} not found", flush=True)
            continue
        print(f"t_N analysis: x = {x_mm:+d} mm", flush=True)
        rows, panels, nearest_top_id, near_label, far_label = analyze_position(
            x_mm, args.data_dir)
        all_rows.extend(rows)
        plot_position(x_mm, panels, near_label, far_label, args.figs_dir)

    if not all_rows:
        print("No data found — run the scan first.", flush=True)
        return 1

    summary = pd.DataFrame(all_rows)
    csv_path = args.csv_dir / "exec13_tN_summary.csv"
    summary.to_csv(csv_path, index=False)
    print(f"  wrote {csv_path}", flush=True)
    plot_synthesis(summary, args.figs_dir)
    print("\nEXEC_13 t_N analysis complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
