#!/usr/bin/env python3
"""Generate absolute and normalized arrival-time examples for EXEC_10."""

from __future__ import annotations

import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

POSITIONS = (-690, -400, 0)
BINS = np.linspace(0.0, 30.0, 301)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("analysis/exec07"))
    args = parser.parse_args()
    positions = pd.read_csv(args.output_dir / "per_position_exec07.csv").set_index("x_beam_mm")
    fig, axes = plt.subplots(3, 2, figsize=(13.5, 11.5), sharex=True)
    bin_width = BINS[1] - BINS[0]
    for row_index, x_mm in enumerate(POSITIONS):
        nearest_id = int(positions.loc[x_mm, "nearest_top_id"])
        with uproot.open(args.data_dir / f"photon_hits_x{x_mm}mm.root") as root_file:
            arrays = root_file["sipm_hits"].arrays(["event_id", "global_id", "time_ns"], library="np")
        selected = arrays["global_id"] == nearest_id
        events = arrays["event_id"][selected].astype(np.int32)
        times = arrays["time_ns"][selected].astype(float)
        n_events = int(np.max(arrays["event_id"])) + 1
        fpt = np.full(n_events, np.inf)
        np.minimum.at(fpt, events, times)
        fpt = fpt[np.isfinite(fpt)]
        axes[row_index, 0].hist(
            times, bins=BINS, weights=np.full(len(times), 1.0 / (n_events * bin_width)),
            histtype="step", lw=1.4, label="all detected photons",
        )
        axes[row_index, 0].hist(
            fpt, bins=BINS, weights=np.full(len(fpt), 1.0 / (n_events * bin_width)),
            histtype="step", lw=1.4, label="FPT per event",
        )
        axes[row_index, 1].hist(times, bins=BINS, density=True, histtype="step", lw=1.4, label="all detected photons")
        axes[row_index, 1].hist(fpt, bins=BINS, density=True, histtype="step", lw=1.4, label="FPT per event")
        for axis in axes[row_index]:
            axis.set_yscale("log")
            axis.grid(alpha=0.25)
            axis.set_xlabel("Arrival time [ns]")
        axes[row_index, 0].set_title(f"x={x_mm:+d} mm, nearest Top ID {nearest_id}: absolute rate")
        axes[row_index, 1].set_title(f"x={x_mm:+d} mm, nearest Top ID {nearest_id}: shape only")
    for axis in axes[:, 0]:
        axis.set_ylabel("Entries / event / ns")
    for axis in axes[:, 1]:
        axis.set_ylabel("Area-normalized density [1/ns]")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].legend(fontsize=8)
    fig.suptitle("Arrival-time examples; FPT = first detected photon time in each event")
    fig.tight_layout()
    path = args.output_dir / "figs" / "P5_tdist_examples.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
