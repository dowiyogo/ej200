"""
Shared figure style for EXEC_12T-B Beamer.
Import and call setup_style() before any figure.
"""

import matplotlib
import matplotlib.pyplot as plt
import subprocess
import os

FOOTER = (
    "EJ-204 | Pair (28,29) | Top readout | intrinsic"
    " | data f431c01"
)
ANA_SHA = "EXEC_12T-B"

COLOR_4PE = "#1f77b4"   # blue
COLOR_20PE = "#d62728"  # red

LABEL_4PE = "4th-hit (4PE)"
LABEL_20PE = "20th-hit (20PE)"


def setup_style():
    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "lines.linewidth": 1.6,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "figure.figsize": (7.0, 4.2),
    })


def add_footer(fig, extra=""):
    text = FOOTER
    if extra:
        text += " | " + extra
    fig.text(0.5, -0.02, text,
             ha="center", va="bottom",
             fontsize=6.5, color="#555555",
             transform=fig.transFigure)


def save(fig, path, extra=""):
    add_footer(fig, extra)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
