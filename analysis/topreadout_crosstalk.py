#!/usr/bin/env python3
"""
topreadout_crosstalk.py — Photon cross-talk between top SiPMs.

For a muon at the midpoint between top SiPMs i and j, quantify photons reaching:
  - the two flanking SiPMs (i, j)
  - the next neighbours (i-1, j+1)
  - farther neighbours (i-2, j+2), (i-3, j+3), ...

Output:
  topreadout_crosstalk.pdf — occupancy per top-SiPM index relative to midpoint,
                             normalized by event.
"""

import argparse
import glob
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

TOP_FACE_TYPE = 2
TOP_GLOBAL_OFFSET = 16


def default_root_files():
    files = sorted(glob.glob("photon_hits_run*.root"))
    if files:
        return files
    for name in ("photon_hits_merged.root", "photon_hits.root"):
        if pathlib.Path(name).exists():
            return [name]
    return []


def load_hits(root_files, step_size="100 MB"):
    chunks = []
    event_offset = 0
    for path in root_files:
        with uproot.open(path) as f:
            tree = f["sipm_hits"]
            cols = ["event_id", "face_type", "global_id"]
            file_chunks = []
            max_event = -1
            for chunk in tree.iterate(cols, library="pd", step_size=step_size):
                if len(chunk) == 0:
                    continue
                max_event = max(max_event, int(chunk["event_id"].max()))
                chunk = chunk[chunk["face_type"] == TOP_FACE_TYPE].copy()
                if len(chunk) == 0:
                    continue
                chunk["event_id"] = chunk["event_id"].astype(np.int64) + event_offset
                file_chunks.append(chunk[["event_id", "global_id"]])
            if file_chunks:
                chunks.extend(file_chunks)
            event_offset += max_event + 1 if max_event >= 0 else 0
    if not chunks:
        raise RuntimeError("No top-SiPM hits found in input files.")
    return pd.concat(chunks, ignore_index=True)


def distance_from_pair(local_id, left_idx, right_idx):
    return min(abs(local_id - left_idx), abs(local_id - right_idx))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_files", nargs="*", help="ROOT files (default: photon_hits_run*.root)")
    parser.add_argument("--left", type=int, default=9, help="Left top-SiPM local index")
    parser.add_argument("--right", type=int, default=10, help="Right top-SiPM local index")
    parser.add_argument("--out", default="topreadout_crosstalk.pdf")
    parser.add_argument("--step-size", default="100 MB")
    args = parser.parse_args()

    root_files = args.root_files or default_root_files()
    if not root_files:
        print("ERROR: no ROOT files found.")
        sys.exit(1)

    df = load_hits(root_files, step_size=args.step_size)
    df["local_id"] = df["global_id"] - TOP_GLOBAL_OFFSET
    midpoint = 0.5 * (args.left + args.right)
    df["relative_index"] = df["local_id"] - midpoint
    df["distance"] = df["local_id"].apply(
        lambda v: distance_from_pair(int(v), args.left, args.right))

    n_events = df["event_id"].nunique()
    occ = (df.groupby(["local_id", "relative_index"])
             .size()
             .reset_index(name="n_hits")
             .sort_values("local_id"))
    occ["hits_per_event"] = occ["n_hits"] / float(n_events)

    near = df[df["distance"].isin([0, 1])]
    far = df[df["distance"] >= 2]
    total = len(df)
    frac_far = len(far) / total if total else np.nan
    frac_near = len(near) / total if total else np.nan

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.bar(occ["relative_index"].to_numpy(), occ["hits_per_event"].to_numpy(), width=0.8,
           color="#4dac26", edgecolor="black", linewidth=0.5)
    ax.axvline(args.left - midpoint, color="#2166ac", ls="--", lw=1.2,
               label=f"flanking {args.left}")
    ax.axvline(args.right - midpoint, color="#2166ac", ls="--", lw=1.2,
               label=f"flanking {args.right}")
    ax.set_xlabel("Top SiPM index relative to midpoint")
    ax.set_ylabel("Detected photons per event")
    ax.set_title("Top-readout optical cross-talk")
    ax.grid(axis="y", alpha=0.3)
    ax.text(0.99, 0.95,
            f"events: {n_events}\n"
            f"near fraction (distance 0 or 1): {frac_near:.3f}\n"
            f"far fraction (distance >= 2): {frac_far:.3f}",
            transform=ax.transAxes, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8"))
    ax.legend(framealpha=0.9)
    fig.tight_layout()
    fig.savefig(args.out)

    print(f"Events with top hits: {n_events}")
    print(f"Photons in flanking+next neighbours (distance 0 or 1): {len(near)} ({frac_near:.3%})")
    print(f"Photons in farther SiPMs (distance >= 2): {len(far)} ({frac_far:.3%})")
    print(f"  -> {args.out}")


if __name__ == "__main__":
    main()
