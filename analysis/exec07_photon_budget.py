#!/usr/bin/env python3
"""EXEC_07 photon-budget and arrival-time analysis for 86-channel EndTop data."""

from __future__ import annotations

import argparse
import math
import pathlib
import re
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import poisson

N_CHANNELS = 86
TOP_IDS = tuple(range(16, 86))
TOP_LEFT_IDS = tuple(range(16, 51))
TOP_RIGHT_IDS = tuple(range(51, 86))
END_CLUSTERS = {
    "end_left_A_SUM4": (0, 1, 2, 3),
    "end_left_B_SUM4": (4, 5, 6, 7),
    "end_right_A_SUM4": (8, 9, 10, 11),
    "end_right_B_SUM4": (12, 13, 14, 15),
}
TOP_POSITIONS_MM = np.array(
    [-692.0 + 20.0 * i for i in range(35)] +
    [12.0 + 20.0 * i for i in range(35)]
)
TAU_R_NS = 0.7
TAU_D_NS = 1.8


@dataclass
class GroupStats:
    name: str
    ids: tuple[int, ...]
    counts: np.ndarray
    arrival_times: np.ndarray
    fpt: np.ndarray


def discover_files(input_dir: pathlib.Path) -> list[pathlib.Path]:
    files = sorted(input_dir.rglob("photon_hits_x*mm.root"))
    return files or sorted(input_dir.rglob("*.root"))


def x_from_file(path: pathlib.Path, gun_x: np.ndarray) -> float:
    match = re.search(r"_x(-?\d+)mm", path.name)
    if match:
        return float(match.group(1))
    unique = np.unique(gun_x)
    if len(unique) != 1:
        raise ValueError(f"{path}: expected one gun_x_mm value, found {unique}")
    return float(unique[0])


def cluster_time(event_id: np.ndarray, time_ns: np.ndarray, mask: np.ndarray,
                 n_events: int, threshold_pe: int) -> np.ndarray:
    frame = pd.DataFrame({"event": event_id[mask], "time": time_ns[mask]})
    if frame.empty:
        return np.full(n_events, np.nan)
    ranked = frame.sort_values(["event", "time"])
    ranked = ranked[ranked.groupby("event").cumcount() == threshold_pe - 1]
    result = np.full(n_events, np.nan)
    result[ranked["event"].to_numpy(dtype=int)] = ranked["time"].to_numpy()
    return result


def make_group(name: str, ids: tuple[int, ...], event_id: np.ndarray,
               global_id: np.ndarray, time_ns: np.ndarray, n_events: int) -> GroupStats:
    mask = np.isin(global_id, ids)
    counts = np.bincount(event_id[mask], minlength=n_events)
    fpt = np.full(n_events, np.inf)
    np.minimum.at(fpt, event_id[mask], time_ns[mask])
    fpt[~np.isfinite(fpt)] = np.nan
    return GroupStats(name, ids, counts, time_ns[mask], fpt)


def summary_row(x_beam: float, group: GroupStats, sigma_dt_ps: float = math.nan) -> dict:
    mean = float(np.mean(group.counts))
    variance = float(np.var(group.counts))
    is_top = bool(set(group.ids) & set(TOP_IDS)) and not bool(set(group.ids) & set(range(16)))
    sigma_est = (
        1000.0 * math.sqrt(TAU_R_NS * TAU_D_NS) / math.sqrt(mean)
        if is_top and mean > 0 else math.nan
    )
    return {
        "x_beam": x_beam,
        "channel/group": group.name,
        "N_pe_mean": mean,
        "N_pe_var": variance,
        "var_over_mean": variance / mean if mean > 0 else math.nan,
        "t_mean_ns": float(np.mean(group.arrival_times)) if group.arrival_times.size else math.nan,
        "FPT_mean_ns": float(np.nanmean(group.fpt)) if np.any(np.isfinite(group.fpt)) else math.nan,
        "sigma_est_ps": sigma_est,
        "sigma_DT_ps": sigma_dt_ps,
    }


def overlay_poisson(ax: plt.Axes, counts: np.ndarray, title: str) -> None:
    mean = float(np.mean(counts))
    maximum = max(int(np.max(counts)), 1)
    bins = np.arange(-0.5, maximum + 1.5, 1.0)
    ax.hist(counts, bins=bins, density=True, histtype="step", linewidth=0.8)
    x = np.arange(maximum + 1)
    ax.plot(x, poisson.pmf(x, mean), color="tab:red", linewidth=0.7)
    fano = np.var(counts) / mean if mean > 0 else math.nan
    ax.set_title(f"{title}\nmean={mean:.1f}, var/mean={fano:.2f}", fontsize=6)
    ax.tick_params(labelsize=5)


def channel_pages(pdf: PdfPages, x_beam: float, groups: list[GroupStats],
                  mode: str) -> None:
    fig, axes = plt.subplots(9, 10, figsize=(18, 15))
    for ax, group in zip(axes.flat, groups):
        if mode == "npe":
            overlay_poisson(ax, group.counts, group.name)
        else:
            ax.hist(group.arrival_times, bins=50, histtype="step", density=True,
                    label="all photons", linewidth=0.8)
            valid_fpt = group.fpt[np.isfinite(group.fpt)]
            ax.hist(valid_fpt, bins=40, histtype="step", density=True,
                    label="FPT", linewidth=0.8)
            ax.set_title(
                f"{group.name}\n<t>={np.mean(group.arrival_times):.2f}, "
                f"<FPT>={np.mean(valid_fpt):.2f} ns",
                fontsize=6,
            )
            ax.tick_params(labelsize=5)
    for ax in axes.flat[len(groups):]:
        ax.axis("off")
    fig.suptitle(f"EXEC_07 {mode}: x_beam={x_beam:g} mm")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def group_page(pdf: PdfPages, x_beam: float, groups: list[GroupStats], mode: str) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    for ax, group in zip(axes.flat, groups):
        if mode == "npe":
            overlay_poisson(ax, group.counts, group.name)
        else:
            ax.hist(group.arrival_times, bins=60, histtype="step", density=True,
                    label="all photons")
            valid_fpt = group.fpt[np.isfinite(group.fpt)]
            ax.hist(valid_fpt, bins=50, histtype="step", density=True, label="FPT")
            ax.set_title(
                f"{group.name}\n<t>={np.mean(group.arrival_times):.2f}, "
                f"<FPT>={np.mean(valid_fpt):.2f} ns",
                fontsize=8,
            )
            ax.legend(fontsize=6)
    for ax in axes.flat[len(groups):]:
        ax.axis("off")
    fig.suptitle(f"EXEC_07 grouped {mode}: x_beam={x_beam:g} mm")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_trends(summary: pd.DataFrame, output_dir: pathlib.Path) -> None:
    selected = [
        "end_left_A_SUM4", "end_left_B_SUM4",
        "end_right_A_SUM4", "end_right_B_SUM4",
        "top_left_all", "top_right_all", "nearest_top",
    ]
    for column, ylabel, filename in [
        ("N_pe_mean", "Mean detected photoelectrons / event", "npe_vs_x.png"),
        ("t_mean_ns", "Mean arrival time [ns]", "mean_t_vs_x.png"),
        ("FPT_mean_ns", "Mean first-photon time [ns]", "fpt_mean_vs_x.png"),
    ]:
        fig, ax = plt.subplots(figsize=(9, 5))
        for name in selected:
            subset = summary[summary["channel/group"] == name].sort_values("x_beam")
            if not subset.empty:
                ax.plot(subset["x_beam"].to_numpy(), subset[column].to_numpy(),
                        marker="o", label=name)
        ax.set(xlabel="Beam x [mm]", ylabel=ylabel)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(output_dir / filename, dpi=180)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for name in ("end_left_delta_AB", "end_right_delta_AB"):
        subset = summary[summary["channel/group"] == name].sort_values("x_beam")
        ax.plot(subset["x_beam"].to_numpy(), subset["sigma_DT_ps"].to_numpy(),
                marker="o", label=name)
    ax.set(xlabel="Beam x [mm]", ylabel="sigma(DeltaT) [ps]")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "end_delta_t.png", dpi=180)
    plt.close(fig)


def analyze_file(path: pathlib.Path, threshold_pe: int, channel_npe_pdf: PdfPages,
                 channel_time_pdf: PdfPages, group_npe_pdf: PdfPages,
                 group_time_pdf: PdfPages, dt_pdf: PdfPages) -> list[dict]:
    with uproot.open(path) as root_file:
        arrays = root_file["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns", "gun_x_mm"], library="np"
        )
    event_id = arrays["event_id"].astype(int)
    global_id = arrays["global_id"].astype(int)
    time_ns = arrays["time_ns"].astype(float)
    x_beam = x_from_file(path, arrays["gun_x_mm"])
    n_events = int(np.max(event_id)) + 1

    channels = [
        make_group(f"ch_{channel:02d}", (channel,), event_id, global_id, time_ns, n_events)
        for channel in range(N_CHANNELS)
    ]
    channel_pages(channel_npe_pdf, x_beam, channels, "N_pe with equal-mean Poisson")
    channel_pages(channel_time_pdf, x_beam, channels, "arrival time and FPT")

    groups = [
        make_group(name, ids, event_id, global_id, time_ns, n_events)
        for name, ids in END_CLUSTERS.items()
    ]
    groups.extend([
        make_group("top_left_all", TOP_LEFT_IDS, event_id, global_id, time_ns, n_events),
        make_group("top_right_all", TOP_RIGHT_IDS, event_id, global_id, time_ns, n_events),
    ])
    nearest_local = int(np.argmin(np.abs(TOP_POSITIONS_MM - x_beam)))
    nearest = make_group(
        "nearest_top", (16 + nearest_local,), event_id, global_id, time_ns, n_events
    )
    groups.append(nearest)
    group_page(group_npe_pdf, x_beam, groups, "N_pe with equal-mean Poisson")
    group_page(group_time_pdf, x_beam, groups, "arrival time and FPT")

    rows = [summary_row(x_beam, group) for group in channels + groups]
    cluster_times = {
        name: cluster_time(event_id, time_ns, np.isin(global_id, ids), n_events, threshold_pe)
        for name, ids in END_CLUSTERS.items()
    }
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, side in zip(axes, ("left", "right")):
        a = cluster_times[f"end_{side}_A_SUM4"]
        b = cluster_times[f"end_{side}_B_SUM4"]
        valid = np.isfinite(a) & np.isfinite(b)
        delta = a[valid] - b[valid]
        sigma_ps = float(np.std(delta) * 1000.0)
        ax.hist(delta, bins=60, histtype="step")
        ax.set_title(
            f"x={x_beam:g} mm, {side}: sigma(DeltaT)={sigma_ps:.1f} ps, "
            f"sigma_group={sigma_ps / math.sqrt(2.0):.1f} ps"
        )
        ax.set_xlabel("DeltaT [ns]")
        rows.append({
            "x_beam": x_beam,
            "channel/group": f"end_{side}_delta_AB",
            "N_pe_mean": math.nan,
            "N_pe_var": math.nan,
            "var_over_mean": math.nan,
            "t_mean_ns": float(np.mean(delta)) if delta.size else math.nan,
            "FPT_mean_ns": math.nan,
            "sigma_est_ps": math.nan,
            "sigma_DT_ps": sigma_ps,
        })
    fig.tight_layout()
    dt_pdf.savefig(fig)
    plt.close(fig)

    left_t1 = cluster_times["end_left_A_SUM4"]
    right_t1 = cluster_times["end_right_A_SUM4"]
    valid = np.isfinite(left_t1) & np.isfinite(right_t1)
    average = 0.5 * (left_t1[valid] + right_t1[valid])
    rows.append({
        "x_beam": x_beam,
        "channel/group": "end_T1_left_right_average",
        "N_pe_mean": math.nan,
        "N_pe_var": math.nan,
        "var_over_mean": math.nan,
        "t_mean_ns": float(np.mean(average)) if average.size else math.nan,
        "FPT_mean_ns": math.nan,
        "sigma_est_ps": math.nan,
        "sigma_DT_ps": float(np.std(average) * 1000.0) if average.size else math.nan,
    })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path("results/exec07_analysis"))
    parser.add_argument("--leading-edge-pe", type=int, default=4)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    files = discover_files(args.input_dir)
    if not files:
        raise SystemExit(f"no ROOT files found under {args.input_dir}")

    rows: list[dict] = []
    with PdfPages(args.output_dir / "npe_poisson_channels.pdf") as npe_pdf, \
         PdfPages(args.output_dir / "arrival_time_channels.pdf") as time_pdf, \
         PdfPages(args.output_dir / "npe_poisson_groups.pdf") as group_npe_pdf, \
         PdfPages(args.output_dir / "arrival_time_groups.pdf") as group_time_pdf, \
         PdfPages(args.output_dir / "end_delta_t_distributions.pdf") as dt_pdf:
        for path in files:
            rows.extend(analyze_file(
                path, args.leading_edge_pe, npe_pdf, time_pdf,
                group_npe_pdf, group_time_pdf, dt_pdf
            ))

    summary = pd.DataFrame(rows).sort_values(["x_beam", "channel/group"])
    summary.to_csv(args.output_dir / "summary_exec07.csv", index=False)
    plot_trends(summary, args.output_dir)
    print(f"wrote {len(summary)} summary rows and plots to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
