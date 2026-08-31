#!/usr/bin/env python3
"""
fig_gen.py  —  Generate the 4 ROOT-data and analytic figures for napkin_first_principles deck.

Output files (in figs/):
  fig_npe_x.pdf        — Npe/end vs x for EJ-200/204/230 (Vikuiti, END-only, phase_ab_optimal.csv)
  fig_sigma_t_x.pdf    — σ(T0) vs x for same three materials
  fig_survival_RN.pdf  — R^N vs N  (analytic, R=0.98)
  fig_bulk_survival.pdf— bar chart exp(−70cm/λ_att) for the three materials
"""

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────
DECK     = Path(__file__).parent
FIGS     = DECK / "figs"
FIGS.mkdir(exist_ok=True)
DATA_CSV = Path("/home/rrios/ej200/analysis/optim/phase_ab_optimal.csv")

# ── Style ──────────────────────────────────────────────────────────────────
NAVY   = "#0d1b2e"
TEAL   = "#1ea8c8"
ORANGE = "#d4a850"
BLUE   = "#3b82f6"
WHITE  = "#f0f4f8"
GRAY   = "#8aa0b8"

COLORS = {"EJ-200": BLUE, "EJ-204": TEAL, "EJ-230": ORANGE}
LABELS = {"EJ-200": "EJ-200  (λ = 380 cm)", "EJ-204": "EJ-204  (λ = 160 cm)", "EJ-230": "EJ-230  (λ = 120 cm)"}

def dark_fig(w=7, h=4.5):
    fig, ax = plt.subplots(figsize=(w, h), facecolor=NAVY)
    ax.set_facecolor(NAVY)
    ax.tick_params(colors=WHITE, which="both", length=4)
    ax.spines["bottom"].set_color(GRAY)
    ax.spines["left"].set_color(GRAY)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.label.set_color(WHITE)
    ax.yaxis.label.set_color(WHITE)
    ax.title.set_color(WHITE)
    return fig, ax

def save(fig, name, ax=None):
    if ax:
        ax.legend(framealpha=0.15, facecolor=NAVY, edgecolor=GRAY,
                  labelcolor=WHITE, fontsize=10)
    fig.tight_layout(pad=1.2)
    path = FIGS / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=NAVY)
    print(f"  saved {path}")
    plt.close(fig)

# ── Load phase_ab_optimal.csv ─────────────────────────────────────────────
data = {}
with open(DATA_CSV) as f:
    for row in csv.DictReader(f):
        mat = row["scint"]  # ej200, ej204, ej230
        x   = int(row["x"])
        npe = 0.5 * (float(row["npe_L"]) + float(row["npe_R"]))   # per-end average
        sig = float(row["sig_m_opt"])
        mat_key = mat.upper().replace("EJ", "EJ-")
        if mat_key not in data:
            data[mat_key] = {"x": [], "npe": [], "sigma": []}
        data[mat_key]["x"].append(x)
        data[mat_key]["npe"].append(npe)
        data[mat_key]["sigma"].append(sig)

# Sort by x
for mat in data:
    idx = np.argsort(data[mat]["x"])
    for key in ("x", "npe", "sigma"):
        data[mat][key] = np.array(data[mat][key])[idx]

# ════════════════════════════════════════════════════════════════════════════
# Figure 1: Npe vs x
# ════════════════════════════════════════════════════════════════════════════
fig, ax = dark_fig()
for mat in ["EJ-200", "EJ-204", "EJ-230"]:
    d = data[mat]
    ax.plot(d["x"], d["npe"], "o-", color=COLORS[mat], linewidth=2,
            markersize=5, label=LABELS[mat])
ax.set_xlabel("Muon position  x  [mm]", fontsize=11)
ax.set_ylabel("$N_{\\mathrm{pe}}$ per END face", fontsize=11)
ax.set_xlim(-650, 650)
ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
ax.yaxis.set_major_locator(ticker.MultipleLocator(200))
ax.grid(color=GRAY, alpha=0.25, linewidth=0.7, linestyle="--")
ax.axvline(0, color=GRAY, alpha=0.4, linewidth=0.8, linestyle=":")
# Ordering annotation
ax.text(640, data["EJ-200"]["npe"][data["EJ-200"]["x"].tolist().index(0)] + 30,
        "EJ-200", color=COLORS["EJ-200"], ha="right", va="bottom", fontsize=9, fontweight="bold")
ax.text(640, data["EJ-204"]["npe"][data["EJ-204"]["x"].tolist().index(0)] + 30,
        "EJ-204", color=COLORS["EJ-204"], ha="right", va="bottom", fontsize=9, fontweight="bold")
ax.text(640, data["EJ-230"]["npe"][data["EJ-230"]["x"].tolist().index(0)] + 30,
        "EJ-230", color=COLORS["EJ-230"], ha="right", va="bottom", fontsize=9, fontweight="bold")
save(fig, "fig_npe_x.pdf", ax)

# ════════════════════════════════════════════════════════════════════════════
# Figure 2: σ(T0) vs x
# ════════════════════════════════════════════════════════════════════════════
fig, ax = dark_fig()
for mat in ["EJ-200", "EJ-204", "EJ-230"]:
    d = data[mat]
    ax.plot(d["x"], d["sigma"], "o-", color=COLORS[mat], linewidth=2,
            markersize=5, label=LABELS[mat])
ax.set_xlabel("Muon position  x  [mm]", fontsize=11)
ax.set_ylabel(r"$\sigma(T_0)$  [ps]", fontsize=11)
ax.set_xlim(-650, 650)
ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax.grid(color=GRAY, alpha=0.25, linewidth=0.7, linestyle="--")
ax.axvline(0, color=GRAY, alpha=0.4, linewidth=0.8, linestyle=":")
# Ordering annotation (inverted vs Npe)
idx0 = data["EJ-230"]["x"].tolist().index(0)
ax.text(640, data["EJ-230"]["sigma"][idx0] - 1.5,
        "EJ-230", color=COLORS["EJ-230"], ha="right", va="top", fontsize=9, fontweight="bold")
ax.text(640, data["EJ-204"]["sigma"][data["EJ-204"]["x"].tolist().index(0)] - 1.5,
        "EJ-204", color=COLORS["EJ-204"], ha="right", va="top", fontsize=9, fontweight="bold")
ax.text(640, data["EJ-200"]["sigma"][data["EJ-200"]["x"].tolist().index(0)] - 1.5,
        "EJ-200", color=COLORS["EJ-200"], ha="right", va="top", fontsize=9, fontweight="bold")
save(fig, "fig_sigma_t_x.pdf", ax)

# ════════════════════════════════════════════════════════════════════════════
# Figure 3: R^N vs N  (analytic, R = 0.98, napkin value)
# ════════════════════════════════════════════════════════════════════════════
fig, ax = dark_fig(w=6, h=4)
R    = 0.98
Nmax = 100
N_pts = np.arange(0, Nmax + 1, 10)
surv  = R ** N_pts
N_fine = np.linspace(0, Nmax, 300)
ax.plot(N_fine, R ** N_fine, color=TEAL, linewidth=2, zorder=2)
ax.plot(N_pts, surv, "o", color=TEAL, markersize=6, zorder=3)

# Key value labels
for Nk, label in [(10, "0.82"), (50, "0.36"), (100, "0.13")]:
    val = R**Nk
    ax.annotate(f"$0.98^{{{Nk}}}={label}$",
                xy=(Nk, val), xytext=(Nk + 3, val + 0.06),
                color=WHITE, fontsize=8.5,
                arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.8))

ax.set_xlabel("Number of reflections  $N$", fontsize=11)
ax.set_ylabel("Survival  $R^{\\,N}$", fontsize=11)
ax.set_xlim(0, Nmax)
ax.set_ylim(0, 1.05)
ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.grid(color=GRAY, alpha=0.25, linewidth=0.7, linestyle="--")
ax.text(50, 0.95, "$R = 0.98$  (Vikuiti 3M ESR)", color=WHITE, fontsize=10,
        ha="center", va="center",
        bbox=dict(facecolor=NAVY, edgecolor=GRAY, alpha=0.7, boxstyle="round,pad=0.3"))
fig.tight_layout(pad=1.2)
fig.savefig(FIGS / "fig_survival_RN.pdf", dpi=150, bbox_inches="tight", facecolor=NAVY)
print(f"  saved {FIGS}/fig_survival_RN.pdf")
plt.close(fig)

# ════════════════════════════════════════════════════════════════════════════
# Figure 4: Bulk survival bar chart (axial best case, exp(−70 cm / λ_att))
# ════════════════════════════════════════════════════════════════════════════
fig, ax = dark_fig(w=5.5, h=4)
mats  = ["EJ-200", "EJ-204", "EJ-230"]
latts = [380, 160, 120]   # cm
d_cm  = 70.0              # cm (bar half-length)
survs = [math.exp(-d_cm / la) for la in latts]
xlabs = [f"EJ-200\n(380 cm)", f"EJ-204\n(160 cm)", f"EJ-230\n(120 cm)"]
clrs  = [COLORS[m] for m in mats]

bars = ax.bar(xlabs, survs, color=clrs, width=0.55, edgecolor=NAVY, linewidth=0.8, zorder=3)
ax.set_ylim(0, 1.05)
ax.set_ylabel("Bulk survival  $\\exp(-70\\,\\mathrm{cm}/\\lambda_{\\mathrm{att}})$", fontsize=10)
ax.grid(axis="y", color=GRAY, alpha=0.25, linewidth=0.7, linestyle="--", zorder=1)
ax.tick_params(axis="x", colors=WHITE, labelsize=10)
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
for bar, val in zip(bars, survs):
    ax.text(bar.get_x() + bar.get_width() / 2.0, val + 0.015,
            f"{val:.2f}", ha="center", va="bottom", color=WHITE,
            fontsize=12, fontweight="bold")
ax.set_title("Axial best case:  centre → one end", color=WHITE, fontsize=10, pad=6)
fig.tight_layout(pad=1.2)
fig.savefig(FIGS / "fig_bulk_survival.pdf", dpi=150, bbox_inches="tight", facecolor=NAVY)
print(f"  saved {FIGS}/fig_bulk_survival.pdf")
plt.close(fig)

print("\nAll figures generated.")
