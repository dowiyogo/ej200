#!/usr/bin/env python3
"""EXEC_13/14D (EJ-230) time-to-threshold t_N analysis.

The nominal thresholds remain t4 for nearest-Top/near-End and t20 for far-End.
HOOK_ADAPTIVE_TN lowers only the far-End threshold when the nominal threshold
does not reach a configurable fraction of events.

Output (all under RESULTS_DIR):
  figs/exec13_tN_{x}mm.png       -- three-panel figure per key position
  figs/exec13_tN_summary.png     -- fitted sigma vs x
  csv/exec13_tN_summary.csv      -- one row per (x, semantic group)
  csv/exec14d_adaptive_tN.csv    -- explicit adaptive-threshold audit
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
    EXPECTED_POSITIONS_MM,
    N_EVENTS,
    TOP_IDS,
    TOP_POSITIONS_MM,
    KEY_POSITIONS,
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

# HOOK_ADAPTIVE_TN: defaults are exposed through CLI arguments in main().
DEFAULT_REACH_MIN = 0.95
DEFAULT_MINIMUM_THRESHOLD = 4
NOMINAL_THRESHOLDS = {
    "top_nearest": 4,
    "end_near": 4,
    "end_far": 20,
}


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


def group_counts(
    event_id: np.ndarray,
    global_id: np.ndarray,
    channel_ids: tuple[int, ...],
) -> np.ndarray:
    """Return PE count per event for one semantic channel group."""
    mask = np.isin(global_id, channel_ids)
    return np.bincount(event_id[mask].astype(int), minlength=N_EVENTS)[:N_EVENTS]


def adaptive_threshold(
    counts: np.ndarray,
    nominal: int,
    reach_min: float,
    minimum_threshold: int,
) -> tuple[int | None, float, float, bool]:
    """Choose the greatest threshold at or below nominal with enough reach."""
    reach_nominal = float(np.mean(counts >= nominal))
    reach_floor = float(np.mean(counts >= minimum_threshold))
    if reach_floor < reach_min:
        return None, reach_nominal, reach_floor, True
    effective = max(
        threshold
        for threshold in range(minimum_threshold, nominal + 1)
        if float(np.mean(counts >= threshold)) >= reach_min
    )
    reach_effective = float(np.mean(counts >= effective))
    return effective, reach_nominal, reach_effective, False


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


def analyze_position(
    x_mm: int,
    data_dir: pathlib.Path,
    reach_min: float,
    minimum_threshold: int,
) -> tuple:
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
         f"Top nearest ID {nearest_top_id}", NOMINAL_THRESHOLDS["top_nearest"]),
        ("end_near", near_ids,
         f"{near_label} SUM4 {{{','.join(map(str, near_ids))}}}",
         NOMINAL_THRESHOLDS["end_near"]),
        ("end_far", far_ids,
         f"{far_label} SUM4 {{{','.join(map(str, far_ids))}}}",
         NOMINAL_THRESHOLDS["end_far"]),
    ]

    rows: list[dict] = []
    panels: dict = {}

    for grp_key, ch_ids, label, nominal in groups:
        counts = group_counts(event_id, global_id, ch_ids)
        effective, reach_nominal, reach_effective, starvation = adaptive_threshold(
            counts, nominal, reach_min, minimum_threshold
        )
        npe_mean = float(np.mean(counts))
        if effective is None:
            tN = np.full(N_EVENTS, np.nan)
            frac_excl = 1.0 - reach_effective
            mean_v = rms_v = sigma_fit = sigma_fit_err = bin_width = math.nan
            resolution_limited = False
            fit_method = "starvation"
        else:
            tN, frac_excl = compute_tN(event_id, global_id, times, ch_ids, effective)
            valid = tN[np.isfinite(tN)]
            mean_v = float(np.mean(valid)) if len(valid) > 0 else math.nan
            rms_v = float(np.std(valid, ddof=1)) if len(valid) > 1 else math.nan
            (
                _, _, sigma_fit, sigma_fit_err, bin_width,
                resolution_limited, fit_method,
            ) = gaussian_core_fit(tN)
        reduced = effective is not None and effective < nominal
        rows.append({
            "x_mm": x_mm,
            "group": grp_key,
            "N": effective,
            "N_nominal": nominal,
            "N_eff": effective,
            "reach_min": reach_min,
            "minimum_threshold": minimum_threshold,
            "reach_nominal": round(reach_nominal, 6),
            "reach_effective": round(reach_effective, 6),
            "starvation": starvation,
            "reduced_threshold": reduced,
            "npe_mean": round(npe_mean, 5),
            "n_events": N_EVENTS,
            "frac_excluded": round(frac_excl, 6),
            "mean_ns": round(mean_v, 5) if math.isfinite(mean_v) else float("nan"),
            "rms_ns": round(rms_v, 5) if math.isfinite(rms_v) else float("nan"),
            "sigma_fit_ns": round(sigma_fit, 6) if math.isfinite(sigma_fit) else float("nan"),
            "sigma_fit_err_ns": round(sigma_fit_err, 6) if math.isfinite(sigma_fit_err) else float("nan"),
            "fit_bin_width_ns": round(bin_width, 8) if math.isfinite(bin_width) else float("nan"),
            "resolution_limited": resolution_limited,
            "fit_method": fit_method,
            "channels": str(ch_ids),
            "near_label": near_label,
            "far_label": far_label,
        })
        panels[grp_key] = {
            "tN": tN,
            "label": label,
            "nominal": nominal,
            "effective": effective,
            "reach_nominal": reach_nominal,
            "reach_effective": reach_effective,
            "reach_min": reach_min,
            "minimum_threshold": minimum_threshold,
            "starvation": starvation,
            "reduced": reduced,
            "npe_mean": npe_mean,
            "frac_excl": frac_excl,
            "sigma_fit": sigma_fit,
            "sigma_fit_err": sigma_fit_err,
            "bin_width": bin_width,
            "resolution_limited": resolution_limited,
        }

    return rows, panels, nearest_top_id, near_label, far_label


def plot_position(x_mm: int, panels: dict, near_label: str, far_label: str,
                  figs_dir: pathlib.Path) -> None:
    col_titles = ["Top nearest", f"{near_label} SUM4", f"{far_label} SUM4"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    for col_i, grp_key in enumerate(("top_nearest", "end_near", "end_far")):
        ax = axes[col_i]
        panel = panels[grp_key]
        tN = panel["tN"]
        effective = panel["effective"]
        valid = tN[np.isfinite(tN)]
        ax.set_xlabel("Arrival time [ns]", fontsize=8)
        if col_i == 0:
            ax.set_ylabel("Events / adaptive bin", fontsize=8)
        ax.set_title(col_titles[col_i], fontsize=9)
        threshold_label = (
            rf"$t_{{{effective}}}$ [effective threshold]"
            if effective is not None else "no estimate"
        )
        ax.text(0.03, 0.97, threshold_label, transform=ax.transAxes,
                ha="left", va="top", fontsize=8)
        if panel["starvation"]:
            ax.text(
                0.5, 0.5,
                f"STARVATION\nN={panel['minimum_threshold']} reaches only "
                f"{100*panel['reach_effective']:.1f}%\n"
                f"(required {100*panel['reach_min']:.1f}%)",
                ha="center", va="center", transform=ax.transAxes, fontsize=9,
                color="darkred",
                bbox={"facecolor": "mistyrose", "edgecolor": "darkred"},
            )
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
        sigma_fit = panel["sigma_fit"]
        sigma_fit_err = panel["sigma_fit_err"]
        if math.isfinite(sigma_fit) and math.isfinite(sigma_fit_err):
            x_fit = np.linspace(lo, hi, 400)
            amp_est = float(counts.max())
            ax.plot(x_fit, _gauss(x_fit, amp_est, mean_v, sigma_fit),
                    "r--", lw=1.1,
                    label=(rf"$\sigma_{{fit}}={sigma_fit*1000:.1f}"
                           rf"\pm{sigma_fit_err*1000:.1f}$ ps"))
            ax.legend(fontsize=7, loc="upper right")
        if panel["resolution_limited"]:
            ax.text(
                0.03, 0.83,
                f"resolution-limited core\n"
                rf"$\langle N_{{pe}}\rangle={panel['npe_mean']:.1f}$",
                transform=ax.transAxes, ha="left", va="top", fontsize=7,
                color="darkorange",
            )
        if panel["reduced"]:
            ax.text(
                0.03, 0.35,
                f"WARNING: N={panel['nominal']} unreachable\n"
                rf"$\langle N_{{pe}}\rangle={panel['npe_mean']:.1f}$; "
                f"R({panel['nominal']})={100*panel['reach_nominal']:.1f}%\n"
                f"report t_N at N_eff={effective}; "
                f"R={100*panel['reach_effective']:.1f}%",
                transform=ax.transAxes, ha="left", va="top", fontsize=7,
                color="darkred",
                bbox={"facecolor": "mistyrose", "alpha": 0.92,
                      "edgecolor": "darkred", "linewidth": 1.0},
            )
        excl_str = f"{100*panel['frac_excl']:.1f}%"
        ax.text(
            0.97, 0.97,
            f"mean = {mean_v*1000:.0f} ps\n"
            f"RMS  = {rms_v*1000:.0f} ps\n"
            f"bin = {panel['bin_width']*1000:.1f} ps\n"
            f"excluded = {excl_str}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.75"},
        )
        ax.grid(alpha=0.22)
    fig.suptitle(rf"EJ-230 $t_N$ $|$ $x = {x_mm:+d}$ mm", fontsize=11)
    fig.tight_layout()
    out = figs_dir / f"exec13_tN_{x_mm}mm.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


def plot_synthesis(summary: pd.DataFrame, figs_dir: pathlib.Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    left = axes[0]
    for grp_key, color, marker, label in (
        ("top_nearest", "tab:green", "^", "Top nearest fitted t4"),
        ("end_near", "tab:blue", "o", "near-End fitted t4"),
    ):
        sub = summary[summary["group"] == grp_key]
        x_values = sub["x_mm"].to_numpy(dtype=float)
        sigma_values = sub["sigma_fit_ns"].to_numpy(dtype=float) * 1000.0
        error_values = sub["sigma_fit_err_ns"].to_numpy(dtype=float) * 1000.0
        left.errorbar(
            x_values, sigma_values, yerr=error_values,
            fmt=marker, color=color, ms=5, capsize=2, label=label,
        )
    left.set_title("Photon-counting limit: fitted sigma(t4)", fontsize=10)
    left.set_xlabel("Beam position x [mm]", fontsize=9)
    left.set_ylabel(r"$\sigma_{fit}(t_4)$ [ps]", fontsize=9)
    left.legend(fontsize=8)
    left.grid(alpha=0.25)

    right = axes[1]
    far = summary[summary["group"] == "end_far"].copy()
    regular = far[(~far["reduced_threshold"]) & (~far["starvation"])]
    reduced = far[far["reduced_threshold"]]
    right.errorbar(
        regular["x_mm"].to_numpy(dtype=float),
        regular["sigma_fit_ns"].to_numpy(dtype=float) * 1000.0,
        yerr=regular["sigma_fit_err_ns"].to_numpy(dtype=float) * 1000.0,
        fmt="s", color="tab:red", ms=5, capsize=2, label="far-End, N_eff=20",
    )
    right.errorbar(
        reduced["x_mm"].to_numpy(dtype=float),
        reduced["sigma_fit_ns"].to_numpy(dtype=float) * 1000.0,
        yerr=reduced["sigma_fit_err_ns"].to_numpy(dtype=float) * 1000.0,
        fmt="X", color="darkorange", ms=7, capsize=2,
        label="far-End, reduced N_eff",
    )
    for row in reduced.itertuples():
        right.annotate(
            f"N={int(row.N_eff)}", (row.x_mm, row.sigma_fit_ns * 1000.0),
            xytext=(0, 7), textcoords="offset points", ha="center", fontsize=6,
        )
    right.set_title("Far-End adaptive time-to-threshold", fontsize=10)
    right.set_xlabel("Beam position x [mm]", fontsize=9)
    right.set_ylabel(r"$\sigma_{fit}(t_{N_{eff}})$ [ps]", fontsize=9)
    right.legend(fontsize=8)
    right.grid(alpha=0.25)
    fig.suptitle(r"EJ-230 fitted timing metrics with HOOK_ADAPTIVE_TN", fontsize=10)
    fig.tight_layout()
    out = figs_dir / "exec13_tN_summary.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


def write_adaptive_doc(
    summary: pd.DataFrame,
    document: pathlib.Path,
    reach_min: float,
    minimum_threshold: int,
) -> None:
    far = summary[summary["group"] == "end_far"].copy()
    lines = [
        "# EXEC_14D adaptive t_N",
        "",
        "## HOOK_ADAPTIVE_TN",
        "",
        f"- `reach_min = {reach_min:.3f}` (CLI: `--reach-min`)",
        f"- `minimum_threshold = {minimum_threshold}` PE "
        "(CLI: `--minimum-threshold`)",
        "- Nominal thresholds: nearest-Top = 4 PE, near-End = 4 PE, "
        "far-End = 20 PE.",
        "- `N_eff` is the largest integer threshold at or below `N_nominal` "
        "whose event reach is at least `reach_min`.",
        "- If the minimum threshold does not reach `reach_min`, the result is "
        "marked as genuine starvation and no timing value is fabricated.",
        "",
        "## Far-End threshold audit",
        "",
        "| x [mm] | N nominal | N eff | mean Npe | R(nominal) | R(effective) | status |",
        "|---:|---:|---:|---:|---:|---:|:---|",
    ]
    for row in far.itertuples():
        effective = "-" if pd.isna(row.N_eff) else str(int(row.N_eff))
        if row.starvation:
            status = "starvation"
        elif row.reduced_threshold:
            status = "reduced"
        else:
            status = "nominal"
        lines.append(
            f"| {int(row.x_mm)} | {int(row.N_nominal)} | {effective} | "
            f"{row.npe_mean:.2f} | {100*row.reach_nominal:.1f}% | "
            f"{100*row.reach_effective:.1f}% | {status} |"
        )
    lines.extend([
        "",
        "## Forward note",
        "",
        "A future EJ-204/EJ-230 comparison must rerun this same reusable hook "
        "on EJ-204 data. Comparing a fixed t20 EJ-204 series with this adaptive "
        "EJ-230 series would not be apples-to-apples.",
        "",
    ])
    document.parent.mkdir(parents=True, exist_ok=True)
    document.write_text("\n".join(lines), encoding="utf-8")
    print(f"  wrote {document}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="EXEC_13: EJ-230 t_N analysis.")
    rdir = results_dir()
    parser.add_argument("--data-dir", type=pathlib.Path, default=rdir / "data")
    parser.add_argument("--figs-dir", type=pathlib.Path, default=rdir / "figures")
    parser.add_argument("--csv-dir",  type=pathlib.Path, default=rdir / "csv")
    parser.add_argument("--reach-min", type=float, default=DEFAULT_REACH_MIN)
    parser.add_argument(
        "--minimum-threshold", type=int, default=DEFAULT_MINIMUM_THRESHOLD,
    )
    parser.add_argument(
        "--document", type=pathlib.Path,
        default=pathlib.Path("docs/exec14d_adaptive_tN.md"),
    )
    args = parser.parse_args()
    if not 0.0 < args.reach_min <= 1.0:
        parser.error("--reach-min must be in (0, 1]")
    if args.minimum_threshold < 1:
        parser.error("--minimum-threshold must be positive")

    args.figs_dir.mkdir(parents=True, exist_ok=True)
    args.csv_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []

    for x_mm in EXPECTED_POSITIONS_MM:
        path = expected_file(args.data_dir, x_mm)
        if not path.is_file():
            print(f"SKIP: {path} not found", flush=True)
            continue
        print(f"t_N analysis: x = {x_mm:+d} mm", flush=True)
        rows, panels, nearest_top_id, near_label, far_label = analyze_position(
            x_mm, args.data_dir, args.reach_min, args.minimum_threshold)
        all_rows.extend(rows)
        if x_mm in KEY_POSITIONS:
            plot_position(x_mm, panels, near_label, far_label, args.figs_dir)

    if not all_rows:
        print("No data found — run the scan first.", flush=True)
        return 1

    summary = pd.DataFrame(all_rows).sort_values(["x_mm", "group"]).reset_index(drop=True)
    csv_path = args.csv_dir / "exec13_tN_summary.csv"
    summary.to_csv(csv_path, index=False)
    print(f"  wrote {csv_path}", flush=True)
    adaptive_path = args.csv_dir / "exec14d_adaptive_tN.csv"
    adaptive_columns = [
        "x_mm", "group", "N_nominal", "N_eff", "reach_min",
        "minimum_threshold", "reach_nominal", "reach_effective",
        "npe_mean", "reduced_threshold", "starvation",
        "sigma_fit_ns", "sigma_fit_err_ns", "resolution_limited",
    ]
    summary[adaptive_columns].to_csv(adaptive_path, index=False)
    print(f"  wrote {adaptive_path}", flush=True)
    plot_synthesis(summary, args.figs_dir)
    write_adaptive_doc(
        summary, args.document, args.reach_min, args.minimum_threshold,
    )
    print("\nEXEC_13 t_N analysis complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
