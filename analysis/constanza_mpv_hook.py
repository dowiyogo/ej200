#!/usr/bin/env python3
"""Build the provisional Constanza ToT-MPV versus simulated-PE hook table."""

import argparse
import csv
from pathlib import Path

import numpy as np
import uproot


# Presentation TB April 2026, slides 22/23. ASIC0/channel5 is explicitly N/A.
DATA_SUM4 = [
    ("ASIC0_CH1", 241.54, 5.89),
    ("ASIC1_CH1", 207.39, 5.04),
    ("ASIC1_CH5", 198.58, 4.84),
]


def stats(values):
    histogram = np.bincount(values)
    return float(values.mean()), float(np.median(values)), int(histogram.argmax())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("center_root", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    with uproot.open(args.center_root) as root_file:
        tree = root_file["sipm_hits"]
        arrays = tree.arrays(["event_id", "global_id", "gun_x_mm"], library="np")

    positions = np.unique(arrays["gun_x_mm"])
    if len(positions) != 1 or positions[0] != 0:
        raise RuntimeError(f"Expected the center run, got positions {positions}")

    n_events = int(arrays["event_id"].max()) + 1
    counts = np.zeros((n_events, 16), dtype=np.int32)
    mask = arrays["global_id"] < 16
    np.add.at(counts, (arrays["event_id"][mask], arrays["global_id"][mask]), 1)

    clusters = [
        ("sim_L_0_3", counts[:, 0:4].sum(axis=1)),
        ("sim_L_4_7", counts[:, 4:8].sum(axis=1)),
        ("sim_R_8_11", counts[:, 8:12].sum(axis=1)),
        ("sim_R_12_15", counts[:, 12:16].sum(axis=1)),
    ]
    cluster_stats = [(name, *stats(values)) for name, values in clusters]
    representative_mean = float(np.mean([row[1] for row in cluster_stats]))
    representative_mode = int(round(np.median([row[3] for row in cluster_stats])))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="") as handle:
        fields = [
            "kind", "label", "sim_mean_pe", "sim_median_pe", "sim_mode_pe",
            "data_mpv_lsb", "data_mpv_ns", "ps_per_lsb",
            "implied_lsb_per_sim_mean_pe", "implied_lsb_per_sim_mode_pe",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for label, mean, median, mode in cluster_stats:
            writer.writerow({
                "kind": "simulation_center_sum4",
                "label": label,
                "sim_mean_pe": f"{mean:.6f}",
                "sim_median_pe": f"{median:.6f}",
                "sim_mode_pe": mode,
            })
        for label, lsb, width_ns in DATA_SUM4:
            writer.writerow({
                "kind": "constanza_sum4_mpv",
                "label": label,
                "sim_mean_pe": f"{representative_mean:.6f}",
                "sim_mode_pe": representative_mode,
                "data_mpv_lsb": f"{lsb:.6f}",
                "data_mpv_ns": f"{width_ns:.6f}",
                "ps_per_lsb": f"{1000.0 * width_ns / lsb:.6f}",
                "implied_lsb_per_sim_mean_pe": f"{lsb / representative_mean:.6f}",
                "implied_lsb_per_sim_mode_pe": f"{lsb / representative_mode:.6f}",
            })

    print(f"Representative simulated center SUM4: mean={representative_mean:.3f} PE, "
          f"mode={representative_mode} PE")
    print("The implied LSB/PE values are a provisional ToT hook, not a charge calibration.")


if __name__ == "__main__":
    main()
