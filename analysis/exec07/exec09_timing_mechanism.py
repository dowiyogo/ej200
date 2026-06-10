#!/usr/bin/env python3
"""Test whether Top windows suppress the late End-photon tail.

The comparison uses the same End channel groups in both campaigns. Statistical
errors are delete-one-block jackknife errors over event IDs, preserving the
within-event correlation among detected photons.
"""

from __future__ import annotations

import argparse
import math
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

POSITIONS = (-400, 0, 400)
ENDTOP_FILES = {x: f"photon_hits_x{x}mm.root" for x in POSITIONS}
END_ONLY_FILES = {
    -400: "photon_hits_run007.root",
    0: "photon_hits_run015.root",
    400: "photon_hits_run023.root",
}
GROUPS = {
    "all_end": tuple(range(16)),
    "end_left_A_SUM4": (0, 1, 2, 3),
    "end_left_B_SUM4": (4, 5, 6, 7),
    "end_right_A_SUM4": (8, 9, 10, 11),
    "end_right_B_SUM4": (12, 13, 14, 15),
}
METRICS = ("mean_ns", "rms_ns", "t50_ns", "t90_ns", "t99_ns", "frac_gt10", "frac_gt20")
N_BLOCKS = 20


def metrics(times: np.ndarray) -> dict[str, float]:
    quantiles = np.quantile(times, [0.5, 0.9, 0.99])
    return {
        "mean_ns": float(np.mean(times)),
        "rms_ns": float(np.std(times)),
        "t50_ns": float(quantiles[0]),
        "t90_ns": float(quantiles[1]),
        "t99_ns": float(quantiles[2]),
        "frac_gt10": float(np.mean(times > 10.0)),
        "frac_gt20": float(np.mean(times > 20.0)),
    }


def jackknife_metrics(event_id: np.ndarray, times: np.ndarray) -> tuple[dict[str, float], dict[str, float]]:
    central = metrics(times)
    replicas = {
        name: [] for name in METRICS
    }
    block = event_id % N_BLOCKS
    for omitted in range(N_BLOCKS):
        replica = metrics(times[block != omitted])
        for name in METRICS:
            replicas[name].append(replica[name])
    errors = {}
    for name in METRICS:
        values = np.asarray(replicas[name])
        errors[name] = float(math.sqrt((N_BLOCKS - 1) / N_BLOCKS * np.sum((values - values.mean()) ** 2)))
    return central, errors


def load_end_hits(path: pathlib.Path, expected_x: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with uproot.open(path) as root_file:
        arrays = root_file["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns", "gun_x_mm"],
            cut="global_id < 16",
            library="np",
        )
    gun_x = np.unique(arrays["gun_x_mm"])
    if gun_x.tolist() != [float(expected_x)]:
        raise RuntimeError(f"{path}: gun_x_mm={gun_x}, expected {expected_x}")
    return (
        arrays["event_id"].astype(np.int32),
        arrays["global_id"].astype(np.int16),
        arrays["time_ns"].astype(float),
    )


def analyze_file(path: pathlib.Path, campaign: str, x_mm: int) -> tuple[list[dict], np.ndarray, int]:
    event_id, global_id, time_ns = load_end_hits(path, x_mm)
    n_events = int(np.max(event_id)) + 1
    rows = []
    for group, ids in GROUPS.items():
        selected = np.isin(global_id, ids)
        central, errors = jackknife_metrics(event_id[selected], time_ns[selected])
        row = {
            "campaign": campaign,
            "x_beam_mm": x_mm,
            "group": group,
            "n_photons": int(np.count_nonzero(selected)),
        }
        for name in METRICS:
            row[name] = central[name]
            row[f"{name}_error"] = errors[name]
        rows.append(row)
    return rows, time_ns, n_events


def comparisons(raw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for x_mm in POSITIONS:
        for group in GROUPS:
            subset = raw[(raw.x_beam_mm == x_mm) & (raw.group == group)].set_index("campaign")
            top = subset.loc["EndTop"]
            pure = subset.loc["End-only"]
            row = {"x_beam_mm": x_mm, "group": group}
            for name in METRICS:
                difference = top[name] - pure[name]
                error = math.hypot(top[f"{name}_error"], pure[f"{name}_error"])
                row[f"{name}_difference"] = difference
                row[f"{name}_difference_error"] = error
                row[f"{name}_difference_sigma"] = difference / error
            rows.append(row)
    return pd.DataFrame(rows)


def plot_overlays(
    all_end: dict[tuple[str, int], tuple[np.ndarray, int]], raw: pd.DataFrame, output: pathlib.Path,
) -> None:
    bins = np.linspace(0.0, 50.0, 251)
    bin_width = bins[1] - bins[0]
    fig, axes = plt.subplots(3, 2, figsize=(13.5, 12), sharex=True)
    for row_index, x_mm in enumerate(POSITIONS):
        for campaign, color in (("End-only", "tab:blue"), ("EndTop", "tab:red")):
            times, n_events = all_end[(campaign, x_mm)]
            axes[row_index, 0].hist(
                times, bins=bins,
                weights=np.full(len(times), 1.0 / (n_events * bin_width)),
                histtype="step", linewidth=1.5, color=color, label=campaign,
            )
            axes[row_index, 1].hist(
                times, bins=bins, density=True, histtype="step",
                linewidth=1.5, color=color, label=campaign,
            )
        for axis in axes[row_index]:
            axis.set_yscale("log")
            axis.set_xlabel("End photon arrival time [ns]")
            axis.grid(alpha=0.25)
        axes[row_index, 0].set_title(f"x beam = {x_mm:+d} mm: absolute rate")
        axes[row_index, 1].set_title(f"x beam = {x_mm:+d} mm: shape only")
        metrics_at_x = raw[(raw.x_beam_mm == x_mm) & (raw.group == "all_end")].set_index("campaign")
        annotation = []
        for campaign in ("End-only", "EndTop"):
            row = metrics_at_x.loc[campaign]
            annotation.append(
                f"{campaign}: t99={row.t99_ns:.3f}+/-{row.t99_ns_error:.3f} ns\n"
                f"  f(t>10)={row.frac_gt10:.5f}+/-{row.frac_gt10_error:.5f}\n"
                f"  f(t>20)={row.frac_gt20:.6f}+/-{row.frac_gt20_error:.6f}"
            )
        axes[row_index, 1].text(
            0.98, 0.97, "\n".join(annotation), transform=axes[row_index, 1].transAxes,
            ha="right", va="top", fontsize=7,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "0.75"},
        )
    for axis in axes[:, 0]:
        axis.set_ylabel("Mean End photons / event / ns")
    for axis in axes[:, 1]:
        axis.set_ylabel("Area-normalized density [1/ns]")
    axes[0, 0].legend()
    axes[0, 1].legend()
    fig.suptitle("EXEC_10: End arrival times - absolute photon rate and normalized shape")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endtop-dir", required=True, type=pathlib.Path)
    parser.add_argument("--end-only-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("analysis/exec07"))
    args = parser.parse_args()

    rows = []
    all_end = {}
    for x_mm in POSITIONS:
        for campaign, directory, files in (
            ("EndTop", args.endtop_dir, ENDTOP_FILES),
            ("End-only", args.end_only_dir, END_ONLY_FILES),
        ):
            path = directory / files[x_mm]
            print(f"{campaign} x={x_mm}: {path}", flush=True)
            file_rows, all_times, n_events = analyze_file(path, campaign, x_mm)
            rows.extend(file_rows)
            all_end[(campaign, x_mm)] = (all_times, n_events)

    raw = pd.DataFrame(rows)
    compared = comparisons(raw)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw.to_csv(args.output_dir / "exec09_tail_metrics.csv", index=False)
    compared.to_csv(args.output_dir / "exec09_tail_comparison.csv", index=False)
    plot_overlays(all_end, raw, args.output_dir / "figs" / "exec09_tail_comparison.png")

    aggregate = compared[compared.group == "all_end"].copy()
    late_suppressed = (
        (aggregate["t99_ns_difference_sigma"] < -3.0)
        & (aggregate["frac_gt10_difference_sigma"] < -3.0)
        & (aggregate["frac_gt20_difference_sigma"] < -3.0)
    )
    early_similar = aggregate["t50_ns_difference"].abs() < 0.5
    confirmed = bool(late_suppressed.all() and early_similar.all())
    verdict = [
        f"mechanism_confirmed={str(confirmed).lower()}",
        "criterion=all three positions: t99, frac_gt10, frac_gt20 lower by >3 sigma; |delta t50| < 0.5 ns",
        aggregate.to_string(index=False),
    ]
    (args.output_dir / "exec09_timing_verdict.txt").write_text("\n".join(verdict) + "\n")
    print("\n".join(verdict))
    return 0 if confirmed else 2


if __name__ == "__main__":
    raise SystemExit(main())
