#!/usr/bin/env python3
"""
compare_edge_wraps.py — overlay edge timing and light yield for wrap modes.

Usage:
    python analysis/compare_edge_wraps.py mylar air black

Each directory is expected to contain edge_summary.csv from edge_resolution.py.
Output:
    edge_wrap_comparison.pdf
"""

import argparse
import pathlib

import matplotlib.pyplot as plt
import pandas as pd


def load_summary(path):
    p = pathlib.Path(path)
    csv = p / "edge_summary.csv" if p.is_dir() else p
    df = pd.read_csv(csv)
    label = p.name if p.is_dir() else p.stem.replace("edge_summary_", "")
    return label, df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summaries", nargs=3,
                        help="Three folders or edge_summary.csv files: mylar air black")
    parser.add_argument("--out", default="edge_wrap_comparison.pdf")
    args = parser.parse_args()

    loaded = [load_summary(p) for p in args.summaries]
    colors = {
        "mylar": "#542788",
        "air": "#2166ac",
        "black": "#b2182b",
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True)
    ax_sig, ax_npe = axes

    for label, df in loaded:
        key = label.lower()
        color = colors.get(key, None)
        ax_sig.plot(df["x_mm"], df["sigma_end_ps"], marker="o", lw=1.8,
                    color=color, label=f"{label} end")
        if "sigma_top_ps" in df:
            ax_sig.plot(df["x_mm"], df["sigma_top_ps"], marker="s", lw=1.2,
                        ls="--", color=color, alpha=0.75, label=f"{label} top")
        ax_npe.plot(df["x_mm"], df["npe_left"] + df["npe_right"],
                    marker="o", lw=1.8, color=color, label=f"{label} end L+R")
        ax_npe.plot(df["x_mm"], df["npe_top"],
                    marker="s", lw=1.2, ls="--", color=color, alpha=0.75,
                    label=f"{label} top")

    ax_sig.set_xlabel("Muon x [mm]")
    ax_sig.set_ylabel(r"$\sigma_t$ [ps]")
    ax_sig.set_title("Timing resolution")
    ax_sig.grid(alpha=0.3)
    ax_sig.legend(fontsize=8, framealpha=0.9)

    ax_npe.set_xlabel("Muon x [mm]")
    ax_npe.set_ylabel(r"$\langle N_{pe}\rangle$ per event")
    ax_npe.set_title("Light yield")
    ax_npe.grid(alpha=0.3)
    ax_npe.legend(fontsize=8, framealpha=0.9)

    fig.suptitle("Edge-wrap comparison")
    fig.tight_layout()
    fig.savefig(args.out)
    print(f"  -> {args.out}")


if __name__ == "__main__":
    main()
