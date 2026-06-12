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

DT_FINE_NS = 0.025
DT_COARSE_NS = 0.100
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


def gaussian_core_fit(
    tN: np.ndarray, dt_ns: float
) -> tuple[float, float, float, float]:
    valid = tN[np.isfinite(tN)]
    n = len(valid)
    if n < 20:
        m = float(np.mean(valid)) if n > 0 else math.nan
        r = float(np.std(valid, ddof=1)) if n > 1 else math.nan
        return m, r, math.nan, math.nan
    mean0 = float(np.mean(valid))
    rms0 = float(np.std(valid, ddof=1))
    lo0, hi0 = mean0 - 2.0 * rms0, mean0 + 2.0 * rms0
    core0 = valid[(valid >= lo0) & (valid <= hi0)]
    if len(core0) < 10:
        return mean0, rms0, math.nan, math.nan
    edges0 = np.arange(lo0, hi0 + dt_ns * 0.5, dt_ns)
    if len(edges0) < 4:
        return mean0, rms0, math.nan, math.nan
    counts0, _ = np.histogram(core0, bins=edges0)
    centers0 = 0.5 * (edges0[:-1] + edges0[1:])
    try:
        popt0, _ = curve_fit(_gauss, centers0, counts0.astype(float),
                             p0=[float(counts0.max()), mean0, rms0], maxfev=5000)
        mu1, sig1 = float(popt0[1]), abs(float(popt0[2]))
    except Exception:
        return mean0, rms0, math.nan, math.nan
    lo1, hi1 = mu1 - 2.0 * sig1, mu1 + 2.0 * sig1
    core1 = valid[(valid >= lo1) & (valid <= hi1)]
    if len(core1) < 10:
        return mean0, rms0, sig1, math.nan
    edges1 = np.arange(lo1, hi1 + dt_ns * 0.5, dt_ns)
    if len(edges1) < 4:
        return mean0, rms0, sig1, math.nan
    counts1, _ = np.histogram(core1, bins=edges1)
    centers1 = 0.5 * (edges1[:-1] + edges1[1:])
    try:
        popt1, pcov1 = curve_fit(_gauss, centers1, counts1.astype(float),
                                 p0=[float(counts1.max()), mu1, sig1], maxfev=5000)
        sigma_fit = abs(float(popt1[2]))
        var_sig = float(pcov1[2, 2])
        sigma_fit_err = math.sqrt(var_sig) if var_sig >= 0.0 else math.nan
    except Exception:
        return mean0, rms0, sig1, math.nan
    return mean0, rms0, sigma_fit, sigma_fit_err


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
    far_dt = DT_FINE_NS if abs(dist_left - dist_right) < 50.0 else DT_COARSE_NS

    groups = [
        ("top_nearest", (nearest_top_id,), DT_FINE_NS,
         f"Top nearest ID {nearest_top_id}"),
        ("end_near", near_ids, DT_FINE_NS,
         f"{near_label} SUM4 {{{','.join(map(str, near_ids))}}}"),
        ("end_far", far_ids, far_dt,
         f"{far_label} SUM4 {{{','.join(map(str, far_ids))}}}"),
    ]

    rows: list[dict] = []
    panels: dict = {}

    for N in N_THRESHOLDS:
        for grp_key, ch_ids, dt_ns, label in groups:
            tN, frac_excl = compute_tN(event_id, global_id, times, ch_ids, N)
            valid = tN[np.isfinite(tN)]
            mean_v = float(np.mean(valid)) if len(valid) > 0 else math.nan
            rms_v = float(np.std(valid, ddof=1)) if len(valid) > 1 else math.nan
            if frac_excl > 0.5:
                sigma_fit = sigma_fit_err = math.nan
            else:
                _, _, sigma_fit, sigma_fit_err = gaussian_core_fit(tN, dt_ns)
            rows.append({
                "x_mm": x_mm, "group": grp_key, "N": N,
                "n_events": N_EVENTS, "frac_excluded": round(frac_excl, 4),
                "mean_ns": round(mean_v, 5) if math.isfinite(mean_v) else float("nan"),
                "rms_ns": round(rms_v, 5) if math.isfinite(rms_v) else float("nan"),
                "sigma_fit_ns": round(sigma_fit, 6) if math.isfinite(sigma_fit) else float("nan"),
                "sigma_fit_err_ns": round(sigma_fit_err, 6) if math.isfinite(sigma_fit_err) else float("nan"),
                "channels": str(ch_ids), "near_label": near_label, "far_label": far_label,
            })
            panels[(N, grp_key)] = (tN, dt_ns, label, frac_excl, sigma_fit, sigma_fit_err)

    return rows, panels, nearest_top_id, near_label, far_label


def plot_position(x_mm: int, panels: dict, near_label: str, far_label: str,
                  figs_dir: pathlib.Path) -> None:
    col_titles = ["Top nearest", f"{near_label} SUM4", f"{far_label} SUM4"]
    row_labels = {4: r"$t_4$  [N=4 PE threshold]", 20: r"$t_{20}$  [N=20 PE]"}
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))
    for row_i, N in enumerate(N_THRESHOLDS):
        for col_i, grp_key in enumerate(("top_nearest", "end_near", "end_far")):
            ax = axes[row_i, col_i]
            tN, dt_ns, label, frac_excl, sigma_fit, sigma_fit_err = panels[(N, grp_key)]
            valid = tN[np.isfinite(tN)]
            ax.set_xlabel("Arrival time [ns]", fontsize=8)
            if col_i == 0:
                ax.set_ylabel(f"Events / {dt_ns*1000:.0f} ps", fontsize=8)
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
            edges = np.arange(lo, hi + dt_ns * 0.5, dt_ns)
            if len(edges) < 3:
                edges = np.linspace(lo, hi, 50)
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
            excl_str = f"{100*frac_excl:.1f}%" if frac_excl > 0 else "0%"
            ax.text(0.97, 0.97,
                    f"mean = {mean_v*1000:.0f} ps\n"
                    f"RMS  = {rms_v*1000:.0f} ps\n"
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
