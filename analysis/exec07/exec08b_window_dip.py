#!/usr/bin/env python3
"""Quantify the EXEC_08b local Top-window dip from existing and short-run ROOT files."""

from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from exec07.common import N_EVENTS, TOP_POSITIONS_MM  # noqa: E402

LOCAL_IDS = tuple(range(16, 27))


def profile(path: pathlib.Path, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    event_counts = np.zeros((N_EVENTS, 86), dtype=np.int32)
    events: set[int] = set()
    gun_x: set[float] = set()
    impact_rows: list[pd.DataFrame] = []
    with uproot.open(path) as root_file:
        for arrays in root_file["sipm_hits"].iterate(
            ["event_id", "global_id", "gun_x_mm", "x_mm", "z_mm"],
            library="np", step_size="200 MB",
        ):
            event_id = arrays["event_id"].astype(int)
            global_id = arrays["global_id"].astype(int)
            combined = event_id * 86 + global_id
            event_counts += np.bincount(
                combined, minlength=N_EVENTS * 86,
            ).reshape(N_EVENTS, 86)
            events.update(np.unique(event_id).tolist())
            gun_x.update(np.unique(arrays["gun_x_mm"]).astype(float).tolist())
            selected = global_id == 18
            if np.any(selected):
                impact_rows.append(pd.DataFrame({
                    "run": label,
                    "x_local_mm": arrays["x_mm"][selected] - TOP_POSITIONS_MM[2],
                    "z_local_mm": arrays["z_mm"][selected],
                }))
    if events != set(range(N_EVENTS)) or len(gun_x) != 1:
        raise RuntimeError(f"invalid run {path}: events={len(events)}, gun_x={gun_x}")
    n_events = len(events)
    rows = []
    for channel in LOCAL_IDS:
        counts = event_counts[:, channel].astype(float)
        mean = float(np.mean(counts))
        rms = float(np.std(counts))
        rows.append({
            "run": label,
            "gun_x_mm": next(iter(gun_x)),
            "channel": channel,
            "channel_x_mm": TOP_POSITIONS_MM[channel - 16],
            "npe_mean": mean,
            "npe_rms": rms,
            "npe_sem": rms / np.sqrt(n_events),
        })
    impacts = pd.concat(impact_rows, ignore_index=True) if impact_rows else pd.DataFrame()
    return pd.DataFrame(rows), impacts


def difference(frame: pd.DataFrame, first: int, second: int) -> tuple[float, float, float]:
    a = frame.set_index("channel").loc[first]
    b = frame.set_index("channel").loc[second]
    value = float(a.npe_mean - b.npe_mean)
    error = float(np.hypot(a.npe_sem, b.npe_sem))
    return value, error, value / error if error > 0 else np.nan


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--existing", required=True, type=pathlib.Path)
    parser.add_argument("--run-a", required=True, type=pathlib.Path)
    parser.add_argument("--run-b", required=True, type=pathlib.Path)
    parser.add_argument("--run-c1", required=True, type=pathlib.Path)
    parser.add_argument("--run-c2", type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("analysis/exec07"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    analyzed = [
        profile(args.existing, "existing_x-650"),
        profile(args.run_a, "run_A_x-652"),
        profile(args.run_b, "run_B_x-642"),
        profile(args.run_c1, "run_C1_x-648"),
    ]
    if args.run_c2:
        analyzed.append(profile(args.run_c2, "run_C2_x-654_exact_mirror"))
    frames = [item[0] for item in analyzed]
    impacts = [item[1] for item in analyzed]
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(args.output_dir / "exec08b_window_dip_profiles.csv", index=False)

    n_columns = 3 if len(frames) > 4 else 2
    n_rows = int(np.ceil(len(frames) / n_columns))
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(6.5 * n_columns, 4.5 * n_rows), sharey=True)
    for ax, frame in zip(axes.flat, frames):
        ax.errorbar(
            frame["channel_x_mm"].to_numpy(),
            frame["npe_mean"].to_numpy(),
            yerr=frame["npe_sem"].to_numpy(),
            fmt="o-", capsize=2,
        )
        gun_x = frame["gun_x_mm"].iloc[0]
        ax.axvline(gun_x, color="tab:red", ls="--", label=f"track x={gun_x:g} mm")
        for row in frame.itertuples():
            ax.annotate(str(row.channel), (row.channel_x_mm, row.npe_mean),
                        xytext=(0, 5), textcoords="offset points", ha="center", fontsize=7)
        ax.set(xlabel="Top SiPM center x [mm]", title=frame["run"].iloc[0])
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    for ax in axes.flat[len(frames):]:
        ax.axis("off")
    axes.flat[0].set_ylabel("Mean detected N_pe / event")
    fig.tight_layout()
    fig.savefig(args.output_dir / "figs" / "exec08b_window_dip_profiles.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(n_rows, n_columns, figsize=(6 * n_columns, 4.5 * n_rows))
    for ax, frame, impact in zip(axes.flat, frames, impacts):
        image = ax.hist2d(
            impact["x_local_mm"], impact["z_local_mm"],
            bins=(60, 60), range=((-3, 3), (-3, 3)), cmap="viridis",
        )
        ax.set(
            xlabel="ID 18 hit x - center [mm]", ylabel="ID 18 hit z [mm]",
            title=frame["run"].iloc[0],
        )
        fig.colorbar(image[3], ax=ax, label="Detected photons")
    for ax in axes.flat[len(frames):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(args.output_dir / "figs" / "exec08b_id18_impact_maps.png", dpi=300)
    plt.close(fig)

    for frame in frames:
        gun_x = frame["gun_x_mm"].iloc[0]
        nearest = frame.iloc[np.argmin(np.abs(frame["channel_x_mm"] - gun_x))]
        maximum = frame.iloc[np.argmax(frame["npe_mean"])]
        ratio = maximum.npe_mean / nearest.npe_mean
        difference_value, difference_error, significance = difference(
            frame, int(maximum.channel), int(nearest.channel),
        )
        print(
            f"{frame['run'].iloc[0]}: nearest ID {nearest.channel} "
            f"Npe={nearest.npe_mean:.3f}+/-{nearest.npe_sem:.3f}; "
            f"max ID {maximum.channel} Npe={maximum.npe_mean:.3f}+/-{maximum.npe_sem:.3f}; "
            f"max/nearest={ratio:.5f}; delta={difference_value:.3f}+/-"
            f"{difference_error:.3f} ({significance:.2f} sigma)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
