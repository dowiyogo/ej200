"""
EXEC_12T-B — Generate all Beamer figures from cached CSV/NPZ data.
Usage:
    python analysis/exec12tb_figures.py <outdir>

<outdir>/figures/ must exist.  All figures written as PDF.
"""

import sys
import os
import json
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from scipy.stats import gaussian_kde

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO = Path("/home/reriosto/SHiP/ej200")
R12T = REPO / "results" / "exec12t_20260612_195426"
R11  = REPO / "results" / "exec11_20260612_182454"
R12  = REPO / "results" / "exec12_20260612_191000"

ANALYSIS12T = R12T / "analysis"
DERIVED12T  = R12T / "derived" / "order_statistics"
TABLES12T   = R12T / "tables"

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

sys.path.insert(0, str(REPO / "analysis"))
from exec12tb_figstyle import (
    setup_style, add_footer, save,
    COLOR_4PE, COLOR_20PE, LABEL_4PE, LABEL_20PE,
)

setup_style()

FOOTER_EXTRA = ""   # can add e.g. analysis commit SHA later

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

tpos = pd.read_csv(ANALYSIS12T / "temporal_position_summary.csv")
tpos4 = tpos[tpos.threshold == 4].sort_values("x_true_mm").reset_index(drop=True)
tpos20 = tpos[tpos.threshold == 20].sort_values("x_true_mm").reset_index(drop=True)

sweep = pd.read_csv(ANALYSIS12T / "threshold_sweep_summary.csv")

cal4  = pd.read_csv(ANALYSIS12T / "calibration_4pe.csv").iloc[0]
cal20 = pd.read_csv(ANALYSIS12T / "calibration_20pe.csv").iloc[0]

pos_rec_4  = pd.read_csv(ANALYSIS12T / "position_reconstruction_4pe.csv")
pos_rec_20 = pd.read_csv(ANALYSIS12T / "position_reconstruction_20pe.csv")

win_dip = pd.read_csv(ANALYSIS12T / "window_dip_summary.csv")

glob = pd.read_csv(TABLES12T / "global_context.csv")

ref_pos = pd.read_csv(R11 / "analysis" / "reference_positions.csv")
REF1 = float(ref_pos[ref_pos.reference == "POS_REF_1"].x_true_mm.values[0])
REF2 = float(ref_pos[ref_pos.reference == "POS_REF_2"].x_true_mm.values[0])

# LOO predictions
loo4_npz  = np.load(ANALYSIS12T / "loo_predictions_4pe.npz")
loo20_npz = np.load(ANALYSIS12T / "loo_predictions_20pe.npz")

x_all_4   = loo4_npz["x_true"]
xrec_all_4  = loo4_npz["x_rec"]
x_all_20  = loo20_npz["x_true"]
xrec_all_20 = loo20_npz["x_rec"]

# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def rms68(arr):
    """Half-width of the central 68% of arr."""
    q = np.nanpercentile(np.abs(arr - np.nanmedian(arr)), 68.0)
    return float(q)


def load_npz(x_mm):
    fname = DERIVED12T / f"pair_x{x_mm:+.1f}mm.npz"
    return np.load(fname)


def dt_ps(npz, k):
    """Δt in ps for the k-th hit (k=4 or 20)."""
    col = k - 1
    tA = npz["t_A"][:, col]
    tB = npz["t_B"][:, col]
    mask = np.isfinite(tA) & np.isfinite(tB)
    return (tA[mask] - tB[mask]) * 1000.0


def xrec_at(x_mm, threshold):
    arr = loo4_npz if threshold == 4 else loo20_npz
    xall = arr["x_true"]
    xrec = arr["x_rec"]
    mask = np.abs(xall - x_mm) < 0.1
    return xrec[mask]


# ---------------------------------------------------------------------------
# Output directory (set by caller)
# ---------------------------------------------------------------------------

OUTDIR = None


def fig_path(name):
    return OUTDIR / "figures" / name


# ===========================================================================
# FIGURE 1 — Fine scan geometry
# ===========================================================================

def fig_fine_scan_geometry():
    fig, ax = plt.subplots(figsize=(8.5, 3.2))

    # Bar (scaled: 1400 mm long, 60 mm tall)
    bar_len = 1400.0
    bar_h   = 60.0
    bar_x0  = -700.0

    bar_rect = mpatches.FancyBboxPatch(
        (bar_x0, -bar_h / 2), bar_len, bar_h,
        boxstyle="round,pad=2", linewidth=1.2,
        edgecolor="#333", facecolor="#d4e6f1", alpha=0.8,
    )
    ax.add_patch(bar_rect)
    ax.text(0, 0, "EJ-204 bar  1400 × 60 × 10 mm",
            ha="center", va="center", fontsize=9, color="#1a1a6e")

    # Fine scan zone
    x_scan_min, x_scan_max = -462.0, -422.0
    scan_rect = mpatches.Rectangle(
        (x_scan_min, -bar_h / 2 - 8), x_scan_max - x_scan_min, bar_h + 16,
        linewidth=0, facecolor="#f9e79f", alpha=0.55, zorder=2,
    )
    ax.add_patch(scan_rect)

    # 41 scan positions
    for xi in range(-462, -421):
        ax.axvline(xi, ymin=0.08, ymax=0.92, color="#e67e22", lw=0.4, alpha=0.6, zorder=3)

    # SiPM A and B
    sipm_y = -bar_h / 2 - 22
    ax.annotate("SiPM A\n(−452 mm)", xy=(-452, -bar_h / 2),
                xytext=(-452, sipm_y), ha="center", fontsize=8,
                arrowprops=dict(arrowstyle="-|>", color="#1f77b4", lw=1.2),
                color="#1f77b4")
    ax.annotate("SiPM B\n(−432 mm)", xy=(-432, -bar_h / 2),
                xytext=(-432, sipm_y), ha="center", fontsize=8,
                arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=1.2),
                color="#d62728")

    # Pair centre
    ax.axvline(-442, color="#666", lw=0.8, ls="--", alpha=0.6)
    ax.text(-442, bar_h / 2 + 6, "pair centre\n−442 mm",
            ha="center", va="bottom", fontsize=7, color="#666")

    # Bar ends
    ax.text(bar_x0 + 5, 0, "End L\n(−700)", ha="left", va="center",
            fontsize=7.5, color="#555")
    ax.text(bar_x0 + bar_len - 5, 0, "End R\n(+700)", ha="right",
            va="center", fontsize=7.5, color="#555")

    ax.set_xlim(-720, 720)
    ax.set_ylim(-75, 55)
    ax.set_xlabel("x (mm)", fontsize=10)
    ax.set_yticks([])
    ax.set_title("EJ-204 bar layout — fine scan of pair (28, 29)")
    ax.grid(False)

    save(fig, fig_path("fine_scan_geometry.pdf"))


# ===========================================================================
# FIGURE 2 — Order-statistic schematic (from real event)
# ===========================================================================

def fig_order_statistic_schematic():
    npz = load_npz(REF1)
    # pick a real event with many hits
    npe_A = npz["npe_A"]
    evt_idx = np.argmax(npe_A)   # event with most photons in A
    tA = npz["t_A"][evt_idx]     # sorted order stats, shape (30,)

    fig, ax = plt.subplots(figsize=(7.5, 3.0))

    # Draw hit ticks
    valid = np.isfinite(tA)
    tA_valid = tA[valid] * 1000.0   # ns → ps (relative to t1)
    t0 = tA_valid[0]
    tA_rel = tA_valid - t0

    for i, t in enumerate(tA_rel[:25]):
        color = "#aaaaaa"
        lw = 0.8
        if i == 3:
            color = COLOR_4PE
            lw = 2.2
        elif i == 19:
            color = COLOR_20PE
            lw = 2.2
        ax.axvline(t, ymin=0.1, ymax=0.9, color=color, lw=lw, alpha=0.85)

    # Labels
    t4  = tA_rel[3]
    t20 = tA_rel[19]
    t1  = tA_rel[0]

    ax.annotate("$t_1$ (1st hit)", xy=(t1, 0.85), xytext=(t1 + 30, 0.95),
                fontsize=8, ha="left", color="#555",
                arrowprops=dict(arrowstyle="->", color="#555", lw=0.8),
                xycoords=("data", "axes fraction"),
                textcoords=("data", "axes fraction"))

    ax.annotate("$t_4$ (4th hit)", xy=(t4, 0.85), xytext=(t4 + 40, 0.78),
                fontsize=9, ha="left", color=COLOR_4PE, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=COLOR_4PE, lw=1.2),
                xycoords=("data", "axes fraction"),
                textcoords=("data", "axes fraction"))

    ax.annotate("$t_{20}$ (20th hit)", xy=(t20, 0.85),
                xytext=(t20 - 80, 0.62),
                fontsize=9, ha="right", color=COLOR_20PE, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=COLOR_20PE, lw=1.2),
                xycoords=("data", "axes fraction"),
                textcoords=("data", "axes fraction"))

    ax.set_xlabel(r"$t - t_1$ (ps)", fontsize=10)
    ax.set_yticks([])
    ax.set_title("Detected-hit sequence for one event — order statistics")
    ax.set_xlim(-50, min(tA_rel[24] + 100, 1600))

    save(fig, fig_path("order_statistic_schematic.pdf"))


# ===========================================================================
# FIGURE 3 — CFD vs order-statistic distinction
# ===========================================================================

def fig_cfd_vs_orderstat():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.5))

    # --- Left: CFD waveform (synthetic) ---
    t = np.linspace(0, 10, 500)
    # rising edge sigmoid, then slow fall
    wf = np.where(t < 5,
                  1 / (1 + np.exp(-3 * (t - 3))),
                  1 / (1 + np.exp(-3 * (t - 3))) * np.exp(-0.4 * (t - 5)))
    wf = wf / wf.max()

    ax1.plot(t, wf, color="#333", lw=2)
    ax1.axhline(0.20, color="#e67e22", lw=1.4, ls="--", label="20% of peak")

    # find crossing
    from scipy.interpolate import interp1d
    rising = t[t < 3.5]
    wf_r   = wf[t < 3.5]
    mask   = wf_r > 0.01
    f_interp = interp1d(wf_r[mask], rising[mask], bounds_error=False)
    t_cfd = float(f_interp(0.20))
    ax1.axvline(t_cfd, color="#e67e22", lw=1.6)
    ax1.annotate("$t_\\mathrm{CFD}$\n(20% crossing)", xy=(t_cfd, 0.20),
                 xytext=(t_cfd + 1.2, 0.38), fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color="#e67e22"),
                 color="#e67e22")

    ax1.set_xlabel("Time (ns)", fontsize=10)
    ax1.set_ylabel("Normalised amplitude", fontsize=10)
    ax1.set_title("Lv et al. method\n20% constant-fraction discriminator",
                  fontsize=9.5)
    ax1.legend(fontsize=8)
    ax1.set_ylim(-0.05, 1.15)

    # --- Right: order-statistic sequence ---
    np.random.seed(42)
    hits = np.sort(np.random.exponential(0.35, 30))
    hits += 1.0  # offset

    for i, h in enumerate(hits[:25]):
        col = "#aaaaaa"
        lw  = 0.8
        if i == 19:
            col = COLOR_20PE
            lw  = 2.4
        ax2.axvline(h, ymin=0.05, ymax=0.92, color=col, lw=lw)

    ax2.annotate("$t_{20}$: 20th detected hit\n(order statistic)",
                 xy=(hits[19], 0.85), xytext=(hits[19] + 0.3, 0.72),
                 fontsize=8.5, ha="left", color=COLOR_20PE, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=COLOR_20PE, lw=1.2),
                 xycoords=("data", "axes fraction"),
                 textcoords=("data", "axes fraction"))

    ax2.set_xlabel("Arrival time (ns)", fontsize=10)
    ax2.set_yticks([])
    ax2.set_title("This analysis\n20th detected-hit — order statistic",
                  fontsize=9.5)
    ax2.set_xlim(0.6, 3.0)

    fig.suptitle(r"$\mathbf{t_{20} \neq 20\%\ CFD}$  — conceptually different estimators",
                 fontsize=11, y=1.01)
    fig.tight_layout()
    save(fig, fig_path("cfd_vs_orderstat_schematic.pdf"))


# ===========================================================================
# FIGURE 4 — μ_Δt(x) with calibration fits and residuals
# ===========================================================================

def fig_mean_dt_vs_x():
    x4  = tpos4.x_true_mm.values
    y4  = tpos4.mean_dt.values        # ps
    x20 = tpos20.x_true_mm.values
    y20 = tpos20.mean_dt.values

    m4, b4   = cal4.slope_ps_per_mm,  cal4.intercept_ps
    m20, b20 = cal20.slope_ps_per_mm, cal20.intercept_ps

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(7.5, 5.5),
                                          gridspec_kw={"height_ratios": [3, 1.2]},
                                          sharex=True)

    ax_top.plot(x4,  y4,  "o", ms=4, color=COLOR_4PE,
                label=f"{LABEL_4PE}  m={m4:.2f} ps/mm")
    ax_top.plot(x20, y20, "s", ms=4, color=COLOR_20PE,
                label=f"{LABEL_20PE}  m={m20:.2f} ps/mm")

    xx = np.linspace(x4.min(), x4.max(), 200)
    ax_top.plot(xx, m4  * xx + b4,  "-", color=COLOR_4PE,  lw=1.2, alpha=0.7)
    ax_top.plot(xx, m20 * xx + b20, "-", color=COLOR_20PE, lw=1.2, alpha=0.7)

    ax_top.set_ylabel(r"$\mu_{\Delta t}$ (ps)")
    ax_top.set_title(r"Mean $\Delta t$ vs. beam position — linear calibrations")
    ax_top.legend(fontsize=8.5)

    # Residuals
    res4  = y4  - (m4  * x4  + b4)
    res20 = y20 - (m20 * x20 + b20)
    ax_bot.plot(x4,  res4,  "o", ms=3.5, color=COLOR_4PE)
    ax_bot.plot(x20, res20, "s", ms=3.5, color=COLOR_20PE)
    ax_bot.axhline(0, color="#888", lw=0.8)
    ax_bot.set_ylabel("Residual (ps)")
    ax_bot.set_xlabel("x (mm)")

    fig.tight_layout()
    save(fig, fig_path("mean_dt_vs_x_with_residuals.pdf"))


# ===========================================================================
# FIGURE 5 — RMS68 Δt vs x
# ===========================================================================

def fig_rms68_dt_vs_x():
    fig, ax = plt.subplots()

    ax.plot(tpos4.x_true_mm.values,  tpos4.rms68_dt.values,  "o-", ms=3.5,
            color=COLOR_4PE,  label=f"{LABEL_4PE}")
    ax.plot(tpos20.x_true_mm.values, tpos20.rms68_dt.values, "s-", ms=3.5,
            color=COLOR_20PE, label=f"{LABEL_20PE}")

    ax.set_xlabel("x (mm)")
    ax.set_ylabel(r"RMS$_{68}(\Delta t)$ (ps)")
    ax.set_title(r"Timing width $\Delta t$ along the fine scan")
    ax.legend()

    save(fig, fig_path("rms68_dt_vs_x_4_20.pdf"))


# ===========================================================================
# FIGURE 6 — ρ(A, B) vs x
# ===========================================================================

def fig_rho_ab_vs_x():
    fig, ax = plt.subplots()

    ax.plot(tpos4.x_true_mm.values,  tpos4.rho_ab.values,  "o-", ms=3.5,
            color=COLOR_4PE,  label=f"{LABEL_4PE}")
    ax.plot(tpos20.x_true_mm.values, tpos20.rho_ab.values, "s-", ms=3.5,
            color=COLOR_20PE, label=f"{LABEL_20PE}")

    ax.axhline(0, color="#888", lw=0.8)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel(r"$\rho(t_A^{(k)},\, t_B^{(k)})$")
    ax.set_title("A–B correlation of timing estimators vs. position")
    ax.legend()

    save(fig, fig_path("rho_ab_vs_x_4_20.pdf"))


# ===========================================================================
# FIGURE 7/8 — Δt overlays at REF1 and REF2
# ===========================================================================

def _dt_overlay(x_mm, fname):
    npz = load_npz(x_mm)
    dt4  = dt_ps(npz, 4)
    dt20 = dt_ps(npz, 20)

    fig, ax = plt.subplots()

    bins4  = np.linspace(dt4.min(),  dt4.max(),  80)
    bins20 = np.linspace(dt20.min(), dt20.max(), 80)

    ax.hist(dt4,  bins=bins4,  histtype="step", lw=1.8, color=COLOR_4PE,
            label=f"{LABEL_4PE}  μ={dt4.mean():.0f} ps  RMS₆₈={rms68(dt4):.0f} ps",
            density=True)
    ax.hist(dt20, bins=bins20, histtype="step", lw=1.8, color=COLOR_20PE,
            label=f"{LABEL_20PE}  μ={dt20.mean():.0f} ps  RMS₆₈={rms68(dt20):.0f} ps",
            density=True)

    ax.set_xlabel(r"$\Delta t = t_A^{(k)} - t_B^{(k)}$ (ps)")
    ax.set_ylabel("Density (arb. units)")
    ax.set_title(f"$\\Delta t$ distribution — $x = {x_mm:.0f}$ mm")
    ax.legend(fontsize=8)

    save(fig, fig_path(fname))


def fig_ref1_dt_overlay():
    _dt_overlay(REF1, "ref1_dt_overlay.pdf")


def fig_ref2_dt_overlay():
    _dt_overlay(REF2, "ref2_dt_overlay.pdf")


# ===========================================================================
# FIGURE 9/10 — x_rec overlays at REF1 and REF2
# ===========================================================================

def _xrec_overlay(x_mm, fname):
    xrec4  = xrec_at(x_mm, 4)
    xrec20 = xrec_at(x_mm, 20)

    fig, ax = plt.subplots()

    # remove NaN/inf
    xrec4  = xrec4[np.isfinite(xrec4)]
    xrec20 = xrec20[np.isfinite(xrec20)]

    all_vals = np.concatenate([xrec4, xrec20])
    lo, hi = np.percentile(all_vals, 0.5), np.percentile(all_vals, 99.5)
    bins = np.linspace(lo, hi, 80)

    ax.hist(xrec4,  bins=bins, histtype="step", lw=1.8, color=COLOR_4PE,
            label=(f"{LABEL_4PE}\n"
                   f"bias={xrec4.mean()-x_mm:.1f} mm  "
                   f"RMS₆₈={rms68(xrec4-x_mm):.1f} mm"),
            density=True)
    ax.hist(xrec20, bins=bins, histtype="step", lw=1.8, color=COLOR_20PE,
            label=(f"{LABEL_20PE}\n"
                   f"bias={xrec20.mean()-x_mm:.1f} mm  "
                   f"RMS₆₈={rms68(xrec20-x_mm):.1f} mm"),
            density=True)

    ax.axvline(x_mm, color="#000", lw=1.2, ls="--", label=f"$x_{{true}}={x_mm:.0f}$ mm")
    ax.set_xlabel("$x_\\mathrm{rec}$ (mm)")
    ax.set_ylabel("Density (arb. units)")
    ax.set_title(f"Reconstructed position — LOO-CV — $x = {x_mm:.0f}$ mm")
    ax.legend(fontsize=8)

    save(fig, fig_path(fname))


def fig_ref1_xrec_overlay():
    _xrec_overlay(REF1, "ref1_xrec_overlay.pdf")


def fig_ref2_xrec_overlay():
    _xrec_overlay(REF2, "ref2_xrec_overlay.pdf")


# ===========================================================================
# FIGURE 11 — σ(t+) vs x
# ===========================================================================

def fig_tplus_rms68_vs_x():
    fig, ax = plt.subplots()

    ax.plot(tpos4.x_true_mm.values,  tpos4.sigma_tplus_core.values,  "o-", ms=3.5,
            color=COLOR_4PE,  label=f"{LABEL_4PE}")
    ax.plot(tpos20.x_true_mm.values, tpos20.sigma_tplus_core.values, "s-", ms=3.5,
            color=COLOR_20PE, label=f"{LABEL_20PE}")

    ax.set_xlabel("x (mm)")
    ax.set_ylabel(r"$\sigma_\mathrm{core}(t_+)$ (ps)")
    ax.set_title(r"Common-time estimator $t_+ = (t_A+t_B)/2$ width vs. position")
    ax.legend()

    save(fig, fig_path("tplus_rms68_vs_x_4_20.pdf"))


# ===========================================================================
# FIGURES 12a/b — Threshold sweep panels
# ===========================================================================

def _sweep_highlight(ax, k_vals=(4, 20)):
    for k, c, lbl in zip(k_vals, [COLOR_4PE, COLOR_20PE],
                          ["k=4", "k=20"]):
        ax.axvline(k, color=c, lw=1.2, ls="--", alpha=0.7, label=lbl)


def fig_sweep_rms68x():
    fig, ax = plt.subplots()
    ax.plot(sweep.threshold.values, sweep.mean_rms68_x_cv.values, "o-", ms=3.5, color="#2c3e50")
    _sweep_highlight(ax)
    # annotate k=4 and k=20
    for k, c in [(4, COLOR_4PE), (20, COLOR_20PE)]:
        row = sweep[sweep.threshold == k].iloc[0]
        ax.annotate(f"k={k}\n{row.mean_rms68_x_cv:.1f} mm",
                    xy=(k, row.mean_rms68_x_cv),
                    xytext=(k + 2, row.mean_rms68_x_cv + 0.5),
                    fontsize=8, color=c,
                    arrowprops=dict(arrowstyle="->", color=c, lw=0.9))
    ax.set_xlabel("Threshold k (detected hits)")
    ax.set_ylabel(r"Mean RMS$_{68}(x_\mathrm{rec})$ (mm)")
    ax.set_title("Position resolution vs. threshold (LOO-CV)")
    ax.legend(fontsize=8)
    save(fig, fig_path("sweep_rms68x_vs_k.pdf"))


def fig_sweep_slope():
    fig, ax = plt.subplots()
    ax.plot(sweep.threshold.values, sweep.slope_ps_per_mm.values, "o-", ms=3.5, color="#2c3e50")
    _sweep_highlight(ax)
    for k, c in [(4, COLOR_4PE), (20, COLOR_20PE)]:
        row = sweep[sweep.threshold == k].iloc[0]
        ax.annotate(f"k={k}\n{row.slope_ps_per_mm:.2f} ps/mm",
                    xy=(k, row.slope_ps_per_mm),
                    xytext=(k + 2, row.slope_ps_per_mm + 0.3),
                    fontsize=8, color=c,
                    arrowprops=dict(arrowstyle="->", color=c, lw=0.9))
    ax.set_xlabel("Threshold k")
    ax.set_ylabel("Calibration slope (ps/mm)")
    ax.set_title("Position sensitivity (slope) vs. threshold")
    ax.legend(fontsize=8)
    save(fig, fig_path("sweep_slope_vs_k.pdf"))


def fig_sweep_chi2():
    fig, ax = plt.subplots()
    ax.plot(sweep.threshold.values, sweep.calibration_chi2_ndf.values, "o-", ms=3.5, color="#2c3e50")
    _sweep_highlight(ax)
    ax.axhline(1.0, color="#888", lw=0.9, ls=":", label="χ²/ndf = 1")
    ax.set_xlabel("Threshold k")
    ax.set_ylabel(r"$\chi^2/\mathrm{ndf}$ (calibration)")
    ax.set_title("Calibration goodness-of-fit vs. threshold")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    save(fig, fig_path("sweep_chi2_vs_k.pdf"))


def fig_sweep_bias():
    fig, ax = plt.subplots()
    ax.plot(sweep.threshold.values, sweep.max_abs_bias.values, "o-", ms=3.5, color="#2c3e50")
    _sweep_highlight(ax)
    for k, c in [(4, COLOR_4PE), (20, COLOR_20PE)]:
        row = sweep[sweep.threshold == k].iloc[0]
        ax.annotate(f"k={k}: {row.max_abs_bias:.2f} mm",
                    xy=(k, row.max_abs_bias),
                    xytext=(k + 2.5, row.max_abs_bias + 0.3),
                    fontsize=8, color=c,
                    arrowprops=dict(arrowstyle="->", color=c, lw=0.9))
    ax.set_xlabel("Threshold k")
    ax.set_ylabel("|max bias| (mm)")
    ax.set_title("Maximum absolute reconstruction bias vs. threshold")
    ax.legend(fontsize=8)
    save(fig, fig_path("sweep_bias_vs_k.pdf"))


# ===========================================================================
# FIGURE 16 — Trade-off: RMS68_x vs |bias| (replaces degenerate Pareto)
# ===========================================================================

def fig_tradeoff_resolution_bias():
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    x = sweep.max_abs_bias.values
    y = sweep.mean_rms68_x_cv.values
    k = sweep.threshold.values

    sc = ax.scatter(x, y, c=k, cmap="viridis_r", s=48, zorder=3, alpha=0.85)
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("k (threshold)", fontsize=9)

    # Highlight k=4 and k=20
    for ki, color, lbl in [(4, COLOR_4PE, "k=4"), (20, COLOR_20PE, "k=20")]:
        row = sweep[sweep.threshold == ki].iloc[0]
        ax.scatter(row.max_abs_bias, row.mean_rms68_x_cv,
                   s=160, color=color, zorder=5,
                   marker="*", label=lbl, edgecolors="#000", lw=0.6)
        ax.annotate(lbl,
                    xy=(row.max_abs_bias, row.mean_rms68_x_cv),
                    xytext=(row.max_abs_bias + 0.15, row.mean_rms68_x_cv + 0.3),
                    fontsize=9, color=color, fontweight="bold")

    ax.set_xlabel("|max bias| across scan (mm)")
    ax.set_ylabel(r"Mean RMS$_{68}(x_\mathrm{rec})$ (mm)")
    ax.set_title("Resolution–bias trade-off  [efficiency = 100% for all k]")
    ax.legend(fontsize=9)

    save(fig, fig_path("tradeoff_resolution_vs_bias.pdf"))


# ===========================================================================
# FIGURES 17/18 — Window-dip
# ===========================================================================

def fig_window_dip_counts():
    run_labels = win_dip.run.apply(
        lambda s: s.replace("photon_hits_run_", "Run ").replace("_x-", "\nx=−")
                   .replace("mm", " mm")
    ).unique()

    runs_4  = win_dip[win_dip.threshold == 4].sort_values("x_true_mm")
    runs_20 = win_dip[win_dip.threshold == 20].sort_values("x_true_mm")

    x_vals = runs_4.x_true_mm.values
    x_str  = [f"x={xi:.0f}" for xi in x_vals]

    fig, ax = plt.subplots(figsize=(7, 3.8))
    w = 0.35
    xs = np.arange(len(x_vals))
    ax.bar(xs - w/2, runs_4.mean_npe.values, width=w, color=COLOR_4PE,
           alpha=0.8, label="4PE (raw NPE)", edgecolor="#333", lw=0.6)
    ax.bar(xs + w/2, runs_20.mean_npe.values, width=w, color=COLOR_20PE,
           alpha=0.8, label="20PE (same events)", edgecolor="#333", lw=0.6)

    ax.set_xticks(xs)
    ax.set_xticklabels(x_str, fontsize=9)
    ax.set_ylabel("Mean NPE")
    ax.set_title("Mean photon yield — window-dip runs")
    ax.legend()

    save(fig, fig_path("window_dip_counts.pdf"))


def fig_window_dip_timing():
    runs_4  = win_dip[win_dip.threshold == 4].sort_values("x_true_mm")
    runs_20 = win_dip[win_dip.threshold == 20].sort_values("x_true_mm")

    x_vals = runs_4.x_true_mm.values
    x_str  = [f"x={xi:.0f}" for xi in x_vals]

    t4_ns  = runs_4.mean_time_ns.values * 1000   # → ps
    t20_ns = runs_20.mean_time_ns.values * 1000

    fig, ax = plt.subplots(figsize=(7, 3.8))
    xs = np.arange(len(x_vals))
    w  = 0.35
    ax.bar(xs - w/2, t4_ns,  width=w, color=COLOR_4PE,  alpha=0.8,
           label=LABEL_4PE, edgecolor="#333", lw=0.6)
    ax.bar(xs + w/2, t20_ns, width=w, color=COLOR_20PE, alpha=0.8,
           label=LABEL_20PE, edgecolor="#333", lw=0.6)

    ax.set_xticks(xs)
    ax.set_xticklabels(x_str, fontsize=9)
    ax.set_ylabel(r"$\langle t_k \rangle$ (ps)")
    ax.set_title("Mean arrival time — window-dip runs")
    ax.legend()

    save(fig, fig_path("window_dip_t4_t20.pdf"))


# ===========================================================================
# FIGURE 19 — Global context bar chart
# ===========================================================================

def fig_global_context():
    method_labels = {
        "A_end":            "Amplitude\n(End)",
        "R_end":            "Ratio\n(End)",
        "dt_end_pool":      "Δt pooled\n(End)",
        "dt_end_weighted":  "Δt weighted\n(End)",
        "local_R":          "Ratio local\n(Pair 28/29)",
        "x_top_centroid_raw": "Centroid\n(Top)",
    }

    rows = glob.copy()
    rows["label"] = rows.method.map(method_labels).fillna(rows.method)

    fig, ax = plt.subplots(figsize=(8.5, 4.0))
    xs = np.arange(len(rows))
    w  = 0.35
    bars_s = ax.bar(xs - w/2, rows.sigma_core_mean, width=w,
                    color="#2980b9", alpha=0.85, label=r"$\sigma_\mathrm{core}$ (mm)",
                    edgecolor="#1a5276", lw=0.7)
    bars_b = ax.bar(xs + w/2, rows.max_abs_bias,    width=w,
                    color="#e74c3c", alpha=0.85, label="|max bias| (mm)",
                    edgecolor="#922b21", lw=0.7)

    ax.set_xticks(xs)
    ax.set_xticklabels(rows.label, fontsize=8.5)
    ax.set_ylabel("Position (mm)")
    ax.set_title("Global X reconstruction summary — EXEC_12 context")
    ax.legend()

    for bar in bars_s:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 0.5,
                f"{h:.1f}", ha="center", va="bottom", fontsize=7)
    for bar in bars_b:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 0.5,
                f"{h:.1f}", ha="center", va="bottom", fontsize=7)

    save(fig, fig_path("global_context_sigma_bias.pdf"))


# ===========================================================================
# FIGURE 20 — LOO position resolution vs x (both thresholds)
# ===========================================================================

def fig_position_resolution_vs_x():
    pr4  = pos_rec_4.sort_values("x_true_mm")
    pr20 = pos_rec_20.sort_values("x_true_mm")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 6.0), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 1]})

    ax1.plot(pr4.x_true_mm.values,  pr4.rms68.values,  "o-", ms=3.5, color=COLOR_4PE,  label=LABEL_4PE)
    ax1.plot(pr20.x_true_mm.values, pr20.rms68.values, "s-", ms=3.5, color=COLOR_20PE, label=LABEL_20PE)
    ax1.set_ylabel(r"RMS$_{68}(x_\mathrm{rec} - x_\mathrm{true})$ (mm)")
    ax1.set_title("LOO-CV position resolution vs. beam position")
    ax1.legend()

    ax2.plot(pr4.x_true_mm.values,  pr4["mean"].values,  "o-", ms=3.5, color=COLOR_4PE)
    ax2.plot(pr20.x_true_mm.values, pr20["mean"].values, "s-", ms=3.5, color=COLOR_20PE)
    ax2.axhline(0, color="#888", lw=0.8)
    ax2.set_ylabel("Bias: $\\langle x_\\mathrm{rec}\\rangle - x_\\mathrm{true}$ (mm)")
    ax2.set_xlabel("x (mm)")

    fig.tight_layout()
    save(fig, fig_path("position_resolution_4_20_vs_x.pdf"))


# ===========================================================================
# Run all
# ===========================================================================

FIGURES = [
    ("fine_scan_geometry.pdf",            fig_fine_scan_geometry),
    ("order_statistic_schematic.pdf",     fig_order_statistic_schematic),
    ("cfd_vs_orderstat_schematic.pdf",    fig_cfd_vs_orderstat),
    ("mean_dt_vs_x_with_residuals.pdf",   fig_mean_dt_vs_x),
    ("rms68_dt_vs_x_4_20.pdf",           fig_rms68_dt_vs_x),
    ("rho_ab_vs_x_4_20.pdf",             fig_rho_ab_vs_x),
    ("ref1_dt_overlay.pdf",              fig_ref1_dt_overlay),
    ("ref2_dt_overlay.pdf",              fig_ref2_dt_overlay),
    ("ref1_xrec_overlay.pdf",            fig_ref1_xrec_overlay),
    ("ref2_xrec_overlay.pdf",            fig_ref2_xrec_overlay),
    ("tplus_rms68_vs_x_4_20.pdf",        fig_tplus_rms68_vs_x),
    ("sweep_rms68x_vs_k.pdf",            fig_sweep_rms68x),
    ("sweep_slope_vs_k.pdf",             fig_sweep_slope),
    ("sweep_chi2_vs_k.pdf",              fig_sweep_chi2),
    ("sweep_bias_vs_k.pdf",              fig_sweep_bias),
    ("tradeoff_resolution_vs_bias.pdf",  fig_tradeoff_resolution_bias),
    ("window_dip_counts.pdf",            fig_window_dip_counts),
    ("window_dip_t4_t20.pdf",            fig_window_dip_timing),
    ("global_context_sigma_bias.pdf",    fig_global_context),
    ("position_resolution_4_20_vs_x.pdf", fig_position_resolution_vs_x),
]


def main(outdir):
    global OUTDIR
    OUTDIR = Path(outdir)
    (OUTDIR / "figures").mkdir(parents=True, exist_ok=True)

    ok = []
    fail = []
    for name, fn in FIGURES:
        try:
            fn()
            size = (OUTDIR / "figures" / name).stat().st_size
            print(f"  OK  {name}  ({size//1024} kB)")
            ok.append(name)
        except Exception as e:
            print(f"  FAIL {name}: {e}")
            fail.append((name, str(e)))

    print(f"\n{len(ok)}/{len(FIGURES)} figures generated.")
    if fail:
        print("FAILURES:")
        for n, e in fail:
            print(f"  {n}: {e}")
        sys.exit(1)

    # Write SHA manifest
    sha_map = {}
    for name, _ in FIGURES:
        p = OUTDIR / "figures" / name
        if p.exists():
            sha_map[name] = hashlib.sha256(p.read_bytes()).hexdigest()[:12]
    (OUTDIR / "figures" / "sha_manifest.json").write_text(
        json.dumps(sha_map, indent=2)
    )
    print("SHA manifest written.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python exec12tb_figures.py <outdir>")
        sys.exit(1)
    main(sys.argv[1])
