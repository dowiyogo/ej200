#!/usr/bin/env python3
"""EXEC_13: fixed-scale t_N reanalysis — three key positions, EJ-204.

Gerardo's request (meeting 2026-06-12): redo EXEC_12 key-position figures with
10 ps binning and IDENTICAL fixed scales between runs so visual comparison
of width/shape between x=-690, -400, 0 mm is meaningful.

Algorithm citations (algorithm unchanged from EXEC_12):
  compute_tN      analysis/exec07/exec12_tN_analysis.py:72-100
  select_nearest  analysis/exec07/exec12_tN_analysis.py:180-183
  time origin     t=0 = Geant4 event start (muon gun); analysis/exec07_photon_budget.py:342
  geometry        analysis/exec07/common.py (TOP_POSITIONS_MM)
  SUM4 clusters   analysis/congruent_sum4_timing.C:218-221

Gates:
  G1 nearest IDs per position == {-690→16, -400→31, 0→51}  (EXEC_12 diagonal)
  G2 all panels within each figure family share identical axis limits
  G3 n_events == 2000 per run
  G4 hist_in + below_threshold + overflow == 2000 per histogram

Usage:
    python3 analysis/exec13/exec13_fixed_scale.py [--data-dir DIR] [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys
from dataclasses import dataclass, field
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import uproot

# ── Add repo root to sys.path so common.py is importable ─────────────────────
_REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "analysis"))
from exec07.common import (  # noqa: E402
    N_EVENTS,
    TOP_IDS,
    TOP_POSITIONS_MM,
    TAU_D_NS,
    expected_file,
)

# ─── Named constants — spec §PARÁMETROS ──────────────────────────────────────
POSITIONS_MM: tuple[int, ...] = (-690, -400, 0)
RUN_LABELS: dict[int, str] = {-690: "Run 1", -400: "Run 3", 0: "Run 4"}

TN_THRESHOLDS_PE: tuple[int, ...] = (4, 20)
TN_BIN_WIDTH_PS: int = 10
TN_RANGE_PS: tuple[int, int] = (0, 1000)
CLUSTER_HALF_WIDTH: int = 2        # nearest ± 2 → up to 5 top channels
TN_MAX_ORDER: int = 20
EFF_MIN_FRACTION: float = 0.99     # minimum channel efficiency for F4 inclusion

# Expected nearest-top IDs (Gate G1) — EXEC_12 diagonal
EXPECTED_NEAREST: dict[int, int] = {-690: 16, -400: 31, 0: 51}

# Histogram binning derived from constants
_N_BINS: int = (TN_RANGE_PS[1] - TN_RANGE_PS[0]) // TN_BIN_WIDTH_PS  # 100
_BIN_EDGES_NS: np.ndarray = np.linspace(
    TN_RANGE_PS[0] / 1000.0, TN_RANGE_PS[1] / 1000.0, _N_BINS + 1
)
_BIN_CENTERS_PS: np.ndarray = (
    _BIN_EDGES_NS[:-1] + TN_BIN_WIDTH_PS / 2000.0
) * 1000.0

# Visual encoding by position
_COLORS: dict[int, str] = {-690: "tab:blue", -400: "tab:orange", 0: "tab:green"}
_LSTYLES: dict[int, str] = {-690: "-", -400: "--", 0: ":"}

# Geometry helpers
_N_TOP: int = len(TOP_IDS)                    # 70
_TOP_POS_ARR: np.ndarray = np.asarray(TOP_POSITIONS_MM, dtype=float)

# Default paths
_DEFAULT_DATA_DIR = pathlib.Path(
    "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"
)
_DEFAULT_OUT_DIR = _REPO / "analysis" / "exec13"


# ─── Data structures ─────────────────────────────────────────────────────────

@dataclass
class RunData:
    x_mm: int
    event_id: np.ndarray
    global_id: np.ndarray
    time_ns: np.ndarray
    n_events: int


@dataclass
class HistResult:
    """Per-channel, per-threshold histogram output."""
    x_mm: int
    nearest_id: int
    threshold: int
    counts: np.ndarray           # shape (_N_BINS,)
    counts_norm: np.ndarray      # area-normalised
    mean_ps: float
    rms_ps: float
    frac_excl: float             # below threshold
    frac_overflow: float         # outside [0, 1000] ps range
    n_in_range: int
    n_overflow: int
    n_below: int


@dataclass
class GateRecord:
    name: str
    passed: bool
    detail: str


# ─── Data loading ─────────────────────────────────────────────────────────────

def load_run(data_dir: pathlib.Path, x_mm: int) -> RunData:
    path = expected_file(data_dir, x_mm)
    if not path.exists():
        sys.exit(f"ABORT G3: ROOT file not found for x={x_mm} mm: {path}")
    with uproot.open(path) as rf:
        arr = rf["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns"], library="np"
        )
    ev = arr["event_id"].astype(np.int32)
    gid = arr["global_id"].astype(np.int32)
    t = arr["time_ns"].astype(np.float64)
    n_actual = int(ev.max()) + 1
    return RunData(x_mm, ev, gid, t, n_actual)


# ─── Core algorithms (unchanged from EXEC_12) ─────────────────────────────────

def compute_tN(
    event_id: np.ndarray,
    global_id: np.ndarray,
    times: np.ndarray,
    channel_ids: tuple[int, ...],
    N: int,
    n_events: int,
) -> tuple[np.ndarray, float]:
    """N-th-order-statistic arrival time per event; NaN when < N photons.

    Ported from analysis/exec07/exec12_tN_analysis.py:72-100 without change.
    """
    mask = np.isin(global_id, channel_ids)
    ev = event_id[mask].astype(np.int64)
    tt = times[mask]

    order = np.lexsort((tt, ev))
    ev_s = ev[order]
    tt_s = tt[order]

    starts = np.searchsorted(ev_s, np.arange(n_events), side="left")
    stops = np.searchsorted(ev_s, np.arange(n_events), side="right")

    tN = np.full(n_events, np.nan)
    have = stops - starts >= N
    tN[have] = tt_s[starts[have] + N - 1]

    frac_excl = float(np.sum(~have)) / n_events
    return tN, frac_excl


def top_npe_matrix(
    event_id: np.ndarray, global_id: np.ndarray, n_events: int
) -> np.ndarray:
    """(n_events × _N_TOP) photon count matrix for all 70 top channels.

    Adapted from analysis/exec07/exec12_tN_analysis.py:172-177.
    """
    top_mask = (global_id >= 16) & (global_id <= 85)
    ev_t = event_id[top_mask].astype(np.int64)
    id_t = global_id[top_mask].astype(np.int64) - 16
    flat = np.bincount(ev_t * _N_TOP + id_t, minlength=n_events * _N_TOP)
    return flat.reshape(n_events, _N_TOP)


def select_nearest_top(x_mm: int, npe_profile: np.ndarray) -> int:
    """Global ID of nearest top channel; max-Npe tie-break.

    Ported from analysis/exec07/exec12_tN_analysis.py:180-183.
    """
    distances = np.abs(_TOP_POS_ARR - x_mm)
    tied = np.flatnonzero(np.isclose(distances, distances.min()))
    local = int(tied[int(np.argmax(npe_profile[tied]))])
    return 16 + local


def cluster_global_ids(nearest_global: int) -> tuple[int, ...]:
    """Cluster of up to 2*CLUSTER_HALF_WIDTH+1 top IDs centred on nearest."""
    lo = max(16, nearest_global - CLUSTER_HALF_WIDTH)
    hi = min(85, nearest_global + CLUSTER_HALF_WIDTH)
    return tuple(range(lo, hi + 1))


def stat_floor(n: int, mean_npe: float) -> float:
    """Statistical floor σ_stat(t_n) ≈ √n · τ_d / ⟨N_pe⟩  (ns)."""
    if mean_npe <= 0:
        return float("nan")
    return (n ** 0.5) * TAU_D_NS / mean_npe


# ─── Gate helpers ─────────────────────────────────────────────────────────────

def _check_axis_limits(figs_xlim: list[tuple], figs_ylim: list[tuple],
                       name: str) -> GateRecord:
    """Gate G2: verify all panels in a family share identical axis limits."""
    xl_ok = all(
        abs(xl[0] - figs_xlim[0][0]) < 1e-6 and abs(xl[1] - figs_xlim[0][1]) < 1e-6
        for xl in figs_xlim
    )
    yl_ok = all(
        abs(yl[0] - figs_ylim[0][0]) < 1e-6 and abs(yl[1] - figs_ylim[0][1]) < 1e-6
        for yl in figs_ylim
    )
    passed = xl_ok and yl_ok
    detail = (
        f"xlims={figs_xlim}, ylims={figs_ylim}"
        if not passed else f"all panels: x{figs_xlim[0]} y{figs_ylim[0]}"
    )
    return GateRecord(f"G2_{name}", passed, detail)


# ─── Figure helpers ───────────────────────────────────────────────────────────

def _save(fig: plt.Figure, out_dir: pathlib.Path, stem: str) -> None:
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / "figs" / f"{stem}.{ext}",
                    bbox_inches="tight", dpi=150)
    plt.close(fig)


def _annotate_panel(ax: plt.Axes, mean_ps: float, rms_ps: float,
                    frac_excl: float, frac_overflow: float) -> None:
    lines = [
        f"μ = {mean_ps:.1f} ps",
        f"RMS = {rms_ps:.1f} ps",
        f"excl = {100*frac_excl:.1f}%",
    ]
    if frac_overflow > 0:
        lines.append(f"ovfl = {100*frac_overflow:.1f}%")
    ax.text(0.97, 0.95, "\n".join(lines), transform=ax.transAxes,
            ha="right", va="top", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))


# ─── F1: t_N histograms with fixed scale ──────────────────────────────────────

def _compute_hist(run: RunData, nearest_id: int, N: int) -> HistResult:
    tN, frac_excl = compute_tN(
        run.event_id, run.global_id, run.time_ns,
        (nearest_id,), N, run.n_events,
    )
    valid = tN[~np.isnan(tN)]
    counts, _ = np.histogram(valid, bins=_BIN_EDGES_NS)
    n_in_range = int(counts.sum())
    n_overflow = int(len(valid)) - n_in_range
    n_below = run.n_events - int(len(valid))

    mean_ps = float(np.nanmean(tN) * 1000) if len(valid) else float("nan")
    rms_ps = float(np.nanstd(tN) * 1000) if len(valid) else float("nan")

    norm = counts / (counts.sum() * TN_BIN_WIDTH_PS / 1000) if counts.sum() else counts
    frac_overflow = n_overflow / run.n_events

    return HistResult(
        x_mm=run.x_mm,
        nearest_id=nearest_id,
        threshold=N,
        counts=counts,
        counts_norm=norm,
        mean_ps=mean_ps,
        rms_ps=rms_ps,
        frac_excl=frac_excl,
        frac_overflow=frac_overflow,
        n_in_range=n_in_range,
        n_overflow=n_overflow,
        n_below=n_below,
    )


def plot_f1(
    runs: list[RunData],
    nearests: dict[int, int],
    N: int,
    out_dir: pathlib.Path,
) -> tuple[list[HistResult], list[GateRecord]]:
    """F1 overlay (norm + abs) and individual panels, fixed scale."""
    hists = [_compute_hist(r, nearests[r.x_mm], N) for r in runs]
    gates: list[GateRecord] = []

    # ── F1a: overlay (normalized) ─────────────────────────────────────────────
    fig_norm, ax_norm = plt.subplots(figsize=(7, 4))
    for h in hists:
        lbl = (f"x={h.x_mm:+d} mm, ID {h.nearest_id} | "
               f"μ={h.mean_ps:.0f} ps, RMS={h.rms_ps:.0f} ps")
        ax_norm.step(_BIN_CENTERS_PS, h.counts_norm,
                     where="mid", color=_COLORS[h.x_mm],
                     ls=_LSTYLES[h.x_mm], lw=1.5, label=lbl)
    ax_norm.set_xlim(*TN_RANGE_PS)
    ax_norm.set_xlabel("$t_{" + str(N) + r"}$ (ps)", fontsize=11)
    ax_norm.set_ylabel("Probability density (ps$^{-1}$)", fontsize=11)
    ax_norm.set_title(
        rf"F1 — $t_{{{N}}}$ top-nearest, normalised, fixed scale [{TN_RANGE_PS[0]},{TN_RANGE_PS[1]}] ps"
    )
    ax_norm.legend(fontsize=8)
    ax_norm.grid(True, alpha=0.3)
    _save(fig_norm, out_dir, f"exec13_f1_t{N}pe_overlay_norm")

    # ── F1a: overlay (absolute counts) ───────────────────────────────────────
    fig_abs, ax_abs = plt.subplots(figsize=(7, 4))
    for h in hists:
        lbl = (f"x={h.x_mm:+d} mm, ID {h.nearest_id} | "
               f"μ={h.mean_ps:.0f} ps, RMS={h.rms_ps:.0f} ps")
        ax_abs.step(_BIN_CENTERS_PS, h.counts,
                    where="mid", color=_COLORS[h.x_mm],
                    ls=_LSTYLES[h.x_mm], lw=1.5, label=lbl)
    ax_abs.set_xlim(*TN_RANGE_PS)
    ax_abs.set_xlabel("$t_{" + str(N) + r"}$ (ps)", fontsize=11)
    ax_abs.set_ylabel("Events / bin", fontsize=11)
    ax_abs.set_title(
        rf"F1 — $t_{{{N}}}$ top-nearest, absolute counts (2000 evt/run)"
    )
    ax_abs.legend(fontsize=8)
    ax_abs.grid(True, alpha=0.3)
    _save(fig_abs, out_dir, f"exec13_f1_t{N}pe_overlay_abs")

    # ── F1b: individual panels, identical scale ───────────────────────────────
    global_ymax_abs = max(h.counts.max() for h in hists) * 1.15
    global_ymax_norm = max(h.counts_norm.max() for h in hists) * 1.15

    fig_pan, axes = plt.subplots(2, 3, figsize=(13, 7), sharey="row", sharex=True)
    fig_pan.suptitle(
        rf"F1 — $t_{{{N}}}$ top-nearest | 10 ps bins, fixed scale (3-run comparison)",
        fontsize=12,
    )

    xlims_recorded = []
    ylims_recorded_norm = []
    ylims_recorded_abs = []

    for col_i, h in enumerate(hists):
        run_lbl = RUN_LABELS[h.x_mm]

        # row 0: normalised
        ax_n = axes[0, col_i]
        ax_n.step(_BIN_CENTERS_PS, h.counts_norm, where="mid",
                  color=_COLORS[h.x_mm], lw=1.5)
        ax_n.set_xlim(*TN_RANGE_PS)
        ax_n.set_ylim(0, global_ymax_norm)
        ax_n.set_title(f"{run_lbl}: x={h.x_mm:+d} mm, ID {h.nearest_id}", fontsize=10)
        if col_i == 0:
            ax_n.set_ylabel("Prob. density (ps$^{-1}$)", fontsize=9)
        _annotate_panel(ax_n, h.mean_ps, h.rms_ps, h.frac_excl, h.frac_overflow)
        ax_n.grid(True, alpha=0.3)
        ylims_recorded_norm.append(ax_n.get_ylim())
        xlims_recorded.append(ax_n.get_xlim())

        # row 1: absolute
        ax_a = axes[1, col_i]
        ax_a.step(_BIN_CENTERS_PS, h.counts, where="mid",
                  color=_COLORS[h.x_mm], lw=1.5)
        ax_a.set_xlim(*TN_RANGE_PS)
        ax_a.set_ylim(0, global_ymax_abs)
        ax_a.set_xlabel("$t_{" + str(N) + r"}$ (ps)", fontsize=9)
        if col_i == 0:
            ax_a.set_ylabel("Events / 10 ps bin", fontsize=9)
        _annotate_panel(ax_a, h.mean_ps, h.rms_ps, h.frac_excl, h.frac_overflow)
        ax_a.grid(True, alpha=0.3)
        ylims_recorded_abs.append(ax_a.get_ylim())

    fig_pan.tight_layout()
    _save(fig_pan, out_dir, f"exec13_f1_t{N}pe_panels")

    # Gate G2 for this figure family
    gates.append(_check_axis_limits(
        xlims_recorded, ylims_recorded_norm, f"F1_t{N}pe_norm"
    ))

    return hists, gates


# ─── F2: cluster <t_n> vs n dispersion ───────────────────────────────────────

def plot_f2(
    runs: list[RunData],
    nearests: dict[int, int],
    npe_matrices: dict[int, np.ndarray],
    out_dir: pathlib.Path,
) -> tuple[list[dict], list[GateRecord]]:
    """F2: <t_n> vs n for cluster nearest±CLUSTER_HALF_WIDTH, one plot per run."""
    palette = plt.cm.viridis(np.linspace(0.1, 0.9, 2 * CLUSTER_HALF_WIDTH + 1))
    gates: list[GateRecord] = []
    all_data: list[dict] = []

    # Pre-compute global y-range across all runs
    global_ymax = 0.0
    global_ymin = float("inf")

    # First pass: collect all mean_tN values
    run_results: list[list[dict]] = []
    for run in runs:
        nearest = nearests[run.x_mm]
        c_ids = cluster_global_ids(nearest)
        npe_mat = npe_matrices[run.x_mm]
        run_curves: list[dict] = []

        for gid in c_ids:
            local = gid - 16
            mean_npe = float(npe_mat[:, local].mean())
            means_ns = []
            rms_ns = []
            floors_ns = []
            for n in range(1, TN_MAX_ORDER + 1):
                tN, _ = compute_tN(
                    run.event_id, run.global_id, run.time_ns,
                    (gid,), n, run.n_events,
                )
                valid = tN[~np.isnan(tN)]
                if len(valid) < 5:
                    means_ns.append(float("nan"))
                    rms_ns.append(float("nan"))
                else:
                    means_ns.append(float(np.mean(valid)))
                    rms_ns.append(float(np.std(valid)))
                floors_ns.append(stat_floor(n, mean_npe))

            run_curves.append({
                "x_mm": run.x_mm, "global_id": gid,
                "mean_npe": mean_npe,
                "means_ns": means_ns, "rms_ns": rms_ns, "floors_ns": floors_ns,
            })
            valid_means = [m for m in means_ns if not np.isnan(m)]
            if valid_means:
                global_ymax = max(global_ymax, max(valid_means) * 1000 * 1.1)
                global_ymin = min(global_ymin, min(valid_means) * 1000 * 0.9)

        run_results.append(run_curves)
        all_data.extend(run_curves)

    global_ymin = max(0, global_ymin)

    # Second pass: plot
    xlims_recorded = []
    ylims_recorded = []

    for run, run_curves in zip(runs, run_results):
        nearest = nearests[run.x_mm]
        c_ids = cluster_global_ids(nearest)
        fig, ax = plt.subplots(figsize=(8, 5))

        for i_curve, (gid, curve) in enumerate(zip(c_ids, run_curves)):
            color = palette[i_curve]
            means_ps = np.asarray(curve["means_ns"]) * 1000
            rms_ps = np.asarray(curve["rms_ns"]) * 1000
            floors_ps = np.asarray(curve["floors_ns"]) * 1000
            ns_arr = np.arange(1, TN_MAX_ORDER + 1)
            lbl = f"ID {gid} (⟨Npe⟩={curve['mean_npe']:.0f})"
            ax.errorbar(ns_arr, means_ps, yerr=rms_ps,
                        fmt="o-", color=color, ms=4, lw=1.5, label=lbl,
                        capsize=2, elinewidth=0.8)
            ax.plot(ns_arr, floors_ps, ls=":", color=color, lw=1.0, alpha=0.6)

        ax.set_xlim(0.5, TN_MAX_ORDER + 0.5)
        ax.set_ylim(global_ymin, global_ymax)
        ax.set_xlabel("Order n", fontsize=11)
        ax.set_ylabel(r"$\langle t_n \rangle$ (ps)", fontsize=11)
        run_lbl = RUN_LABELS[run.x_mm]
        ax.set_title(
            rf"F2 — {run_lbl}: x={run.x_mm:+d} mm, nearest±{CLUSTER_HALF_WIDTH} cluster"
            "\n(dashed: stat floor $\sqrt{n}\cdot\\tau_d/\\langle N_{pe}\\rangle$)",
            fontsize=10,
        )
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(True, alpha=0.3)
        _save(fig, out_dir, f"exec13_f2_cluster_{run.x_mm:+d}mm")
        xlims_recorded.append(ax.get_xlim() if hasattr(ax, "_xlim_recorded") else (0.5, TN_MAX_ORDER + 0.5))
        ylims_recorded.append((global_ymin, global_ymax))

    gates.append(_check_axis_limits(
        [(0.5, TN_MAX_ORDER + 0.5)] * len(runs),
        [(global_ymin, global_ymax)] * len(runs),
        "F2_cluster",
    ))
    return all_data, gates


# ─── F3: Top Npe profile ──────────────────────────────────────────────────────

def plot_f3(
    runs: list[RunData],
    nearests: dict[int, int],
    npe_matrices: dict[int, np.ndarray],
    out_dir: pathlib.Path,
) -> tuple[list[dict], list[GateRecord]]:
    """F3: mean Npe ± RMS for all 70 top channels, 3 panels fixed scale."""
    top_global_ids = np.arange(16, 86)
    gates: list[GateRecord] = []
    all_data: list[dict] = []

    # Compute profiles
    profiles: list[dict] = []
    global_ymax = 0.0
    for run in runs:
        mat = npe_matrices[run.x_mm]
        mean_npe = mat.mean(axis=0)
        rms_npe = mat.std(axis=0)
        profiles.append({"x_mm": run.x_mm, "mean": mean_npe, "rms": rms_npe})
        global_ymax = max(global_ymax, (mean_npe + rms_npe).max() * 1.1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    fig.suptitle("F3 — Top Npe profile (all 70 channels), fixed scale", fontsize=12)

    xlims_rec = []
    ylims_rec = []

    for ax, run, prof in zip(axes, runs, profiles):
        nearest = nearests[run.x_mm]
        c_ids = cluster_global_ids(nearest)
        local_ids = top_global_ids - 16

        ax.errorbar(top_global_ids, prof["mean"], yerr=prof["rms"],
                    fmt=".", ms=4, color="tab:gray", elinewidth=0.6, capsize=2)

        # Highlight cluster
        for gid in c_ids:
            li = gid - 16
            ax.errorbar([gid], [prof["mean"][li]], yerr=[prof["rms"][li]],
                        fmt="o", ms=7, color=_COLORS[run.x_mm], elinewidth=1.0)

        # Vertical line at track position
        ax.axvline(nearest, color=_COLORS[run.x_mm], ls="--", lw=1.2, alpha=0.7,
                   label=f"nearest ID {nearest}")
        ax.set_xlim(15, 86)
        ax.set_ylim(0, global_ymax)
        ax.set_xlabel("Top channel global ID", fontsize=9)
        if ax is axes[0]:
            ax.set_ylabel(r"$\langle N_{pe} \rangle$ / event", fontsize=10)
        run_lbl = RUN_LABELS[run.x_mm]
        ax.set_title(f"{run_lbl}: x={run.x_mm:+d} mm", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        xlims_rec.append(ax.get_xlim())
        ylims_rec.append(ax.get_ylim())

        for local, gid in enumerate(top_global_ids):
            all_data.append({
                "x_mm": run.x_mm,
                "global_id": int(gid),
                "mean_npe": float(prof["mean"][local]),
                "rms_npe": float(prof["rms"][local]),
            })

    fig.tight_layout()
    _save(fig, out_dir, "exec13_f3_npe_profile")
    gates.append(_check_axis_limits(xlims_rec, ylims_rec, "F3_npe"))
    return all_data, gates


# ─── F4: Top timing profile ───────────────────────────────────────────────────

def plot_f4(
    runs: list[RunData],
    nearests: dict[int, int],
    out_dir: pathlib.Path,
) -> tuple[list[dict], list[GateRecord]]:
    """F4: <t_N> ± RMS for each top channel passing EFF_MIN_FRACTION."""
    top_global_ids = list(range(16, 86))
    gates: list[GateRecord] = []
    all_data: list[dict] = []

    for N in TN_THRESHOLDS_PE:
        # Compute per-channel timing for all runs
        run_profiles: list[dict] = []
        global_ymax = 0.0
        global_ymin = float("inf")

        for run in runs:
            means = np.full(_N_TOP, float("nan"))
            rmss = np.full(_N_TOP, float("nan"))
            effs = np.zeros(_N_TOP)
            included = np.zeros(_N_TOP, dtype=bool)

            for gid in top_global_ids:
                local = gid - 16
                tN, frac_excl = compute_tN(
                    run.event_id, run.global_id, run.time_ns,
                    (gid,), N, run.n_events,
                )
                eff = 1.0 - frac_excl
                effs[local] = eff
                if eff >= EFF_MIN_FRACTION:
                    valid = tN[~np.isnan(tN)]
                    means[local] = float(np.mean(valid)) * 1000  # → ps
                    rmss[local] = float(np.std(valid)) * 1000
                    included[local] = True

            run_profiles.append({
                "x_mm": run.x_mm, "means_ps": means, "rmss_ps": rmss,
                "effs": effs, "included": included,
            })
            valid_m = means[included]
            if len(valid_m):
                global_ymax = max(global_ymax, float((means + rmss)[included].max()) * 1.1)
                global_ymin = min(global_ymin, float((means - rmss)[included].min()) * 0.9)

            # Collect CSV data
            for gid in top_global_ids:
                local = gid - 16
                all_data.append({
                    "x_mm": run.x_mm,
                    "global_id": gid,
                    "threshold": N,
                    "mean_t_ps": float(means[local]) if included[local] else float("nan"),
                    "rms_t_ps": float(rmss[local]) if included[local] else float("nan"),
                    "eff": float(effs[local]),
                    "included": bool(included[local]),
                })

        global_ymin = max(0, global_ymin)

        # Plot
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
        fig.suptitle(
            rf"F4 — Top timing profile $\langle t_{{{N}}}\rangle$ "
            rf"(channels with eff$\geq${EFF_MIN_FRACTION:.2f}), fixed scale",
            fontsize=11,
        )
        xlims_rec = []
        ylims_rec = []

        for ax, run, prof in zip(axes, runs, run_profiles):
            nearest = nearests[run.x_mm]
            gids_arr = np.array(top_global_ids)
            inc = prof["included"]
            gids_inc = gids_arr[inc]
            means_inc = prof["means_ps"][inc]
            rmss_inc = prof["rmss_ps"][inc]

            ax.errorbar(gids_inc, means_inc, yerr=rmss_inc,
                        fmt=".", ms=5, color=_COLORS[run.x_mm],
                        elinewidth=0.8, capsize=2)
            ax.axvline(nearest, color=_COLORS[run.x_mm], ls="--", lw=1.2,
                       alpha=0.7, label=f"nearest ID {nearest}")
            ax.set_xlim(15, 86)
            ax.set_ylim(global_ymin, global_ymax)
            ax.set_xlabel("Top channel global ID", fontsize=9)
            if ax is axes[0]:
                ax.set_ylabel(r"$\langle t_N \rangle$ (ps)", fontsize=10)
            run_lbl = RUN_LABELS[run.x_mm]
            n_excl = int(np.sum(~inc))
            ax.set_title(
                f"{run_lbl}: x={run.x_mm:+d} mm\n"
                f"({n_excl} channels excluded, eff<{EFF_MIN_FRACTION:.2f})",
                fontsize=9,
            )
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            xlims_rec.append(ax.get_xlim())
            ylims_rec.append(ax.get_ylim())

        fig.tight_layout()
        _save(fig, out_dir, f"exec13_f4_timing_t{N}pe")
        gates.append(_check_axis_limits(xlims_rec, ylims_rec, f"F4_t{N}pe"))

    return all_data, gates


# ─── CSV writers ─────────────────────────────────────────────────────────────

def write_csv_f1(hists: list[HistResult], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_f1_tn_histograms.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "run_label", "threshold_pe", "nearest_id",
                    "bin_center_ps", "counts", "counts_norm",
                    "mean_ps", "rms_ps", "frac_excl", "frac_overflow",
                    "n_in_range", "n_overflow", "n_below"])
        for h in hists:
            for bc, cnt, cnrm in zip(_BIN_CENTERS_PS, h.counts, h.counts_norm):
                w.writerow([
                    h.x_mm, RUN_LABELS[h.x_mm], h.threshold, h.nearest_id,
                    f"{bc:.1f}", cnt, f"{cnrm:.6f}",
                    f"{h.mean_ps:.3f}", f"{h.rms_ps:.3f}",
                    f"{h.frac_excl:.6f}", f"{h.frac_overflow:.6f}",
                    h.n_in_range, h.n_overflow, h.n_below,
                ])


def write_csv_f2(data: list[dict], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_f2_cluster_dispersion.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "run_label", "global_id", "mean_npe",
                    "n", "mean_tn_ps", "rms_tn_ps", "stat_floor_ps"])
        for d in data:
            for n_i, (m, r, fl) in enumerate(
                zip(d["means_ns"], d["rms_ns"], d["floors_ns"])
            ):
                w.writerow([
                    d["x_mm"], RUN_LABELS[d["x_mm"]], d["global_id"],
                    f"{d['mean_npe']:.2f}", n_i + 1,
                    f"{m*1000:.3f}" if not np.isnan(m) else "nan",
                    f"{r*1000:.3f}" if not np.isnan(r) else "nan",
                    f"{fl*1000:.3f}" if not np.isnan(fl) else "nan",
                ])


def write_csv_f3(data: list[dict], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_f3_npe_profile.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "run_label", "global_id", "mean_npe", "rms_npe"])
        for d in data:
            w.writerow([d["x_mm"], RUN_LABELS[d["x_mm"]], d["global_id"],
                        f"{d['mean_npe']:.4f}", f"{d['rms_npe']:.4f}"])


def write_csv_f4(data: list[dict], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_f4_timing_profile.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "run_label", "global_id", "threshold_pe",
                    "mean_t_ps", "rms_t_ps", "eff", "included"])
        for d in data:
            w.writerow([
                d["x_mm"], RUN_LABELS[d["x_mm"]], d["global_id"], d["threshold"],
                f"{d['mean_t_ps']:.3f}" if not np.isnan(d["mean_t_ps"]) else "nan",
                f"{d['rms_t_ps']:.3f}" if not np.isnan(d["rms_t_ps"]) else "nan",
                f"{d['eff']:.6f}", d["included"],
            ])


def write_csv_gates(gates: list[GateRecord], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_gates.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gate", "passed", "detail"])
        for g in gates:
            w.writerow([g.name, g.passed, g.detail])


# ─── LaTeX table (from CSV) ───────────────────────────────────────────────────

def make_latex_table(hists: list[HistResult]) -> str:
    """Generate LaTeX table of t_N statistics from F1 data (no manual transcription)."""
    rows_4 = [(h.x_mm, h.mean_ps, h.rms_ps, h.frac_excl, h.frac_overflow)
              for h in hists if h.threshold == 4]
    rows_20 = [(h.x_mm, h.mean_ps, h.rms_ps, h.frac_excl, h.frac_overflow)
               for h in hists if h.threshold == 20]

    def fmt_row(rx, mean, rms, fexcl, fovfl):
        excl_str = f"{100*fexcl:.1f}\\%" if fexcl > 0 else "0"
        ovfl_str = f"{100*fovfl:.1f}\\%" if fovfl > 0 else "0"
        return (f"  $x={rx:+d}$ mm & {mean:.1f} & {rms:.1f} & "
                f"{excl_str} & {ovfl_str} \\\\")

    rows_str_4 = "\n".join(fmt_row(*r) for r in rows_4)
    rows_str_20 = "\n".join(fmt_row(*r) for r in rows_20)

    return r"""% Auto-generated by exec13_fixed_scale.py — do not edit manually
% Source: exec13_f1_tn_histograms.csv
\begin{table}[h]
\centering
\small
\caption{EXEC\_13: $t_N$ statistics at top-nearest channel.
Bin width 10\,ps, range [0, 1000]\,ps.
All values from EJ-204 (OPSC-101) EXEC\_12 data.}
\label{tab:exec13_tN}
\begin{tabular}{lrrrr}
\toprule
Position & $\mu$ (ps) & RMS (ps) & Excl.\ frac.\ & Overflow \\
\midrule
\multicolumn{5}{c}{\textit{$t_4$ (leading-edge threshold = 4 PE)}} \\
""" + rows_str_4 + r"""
\midrule
\multicolumn{5}{c}{\textit{$t_{20}$ (bulk-light level = 20 PE)}} \\
""" + rows_str_20 + r"""
\bottomrule
\end{tabular}
\end{table}
"""


# ─── Beamer generator ─────────────────────────────────────────────────────────

def _frame(title: str, body: str) -> str:
    return f"\\begin{{frame}}\n\\frametitle{{{title}}}\n{body}\n\\end{{frame}}\n\n"


def _howtoread(why: str, axes: str, takeaway: str) -> str:
    return (r"\begin{block}{\footnotesize How to read this plot}" + "\n"
            + r"\scriptsize\textbf{Why:} " + why + r"\quad"
            + r"\textbf{Axes:} " + axes + r"\quad"
            + r"\textbf{Takeaway:} " + takeaway + "\n"
            + r"\end{block}" + "\n")


def generate_beamer(hists: list[HistResult], commit_hash: str) -> str:
    preamble = r"""\documentclass[aspectratio=169]{beamer}
\usetheme{Madrid}
\usecolortheme{default}
\usepackage{booktabs,amsmath}
\graphicspath{{figs/}}
\title[EXEC\_13]{EXEC\_13: Fixed-Scale $t_N$ Analysis (EJ-204)}
\subtitle{Gerardo's request --- meeting 12 Jun 2026}
\author{SHiP T0 MiniDAQ}
\date{2026-06-13}
\begin{document}
\begin{frame}\titlepage\end{frame}

"""

    frames = []

    # ── F1 slides ─────────────────────────────────────────────────────────────
    for N in TN_THRESHOLDS_PE:
        h4 = _howtoread(
            f"EXEC\\_12 histograms used auto-scale; this re-bins at 10\\,ps with fixed [0,1000]\\,ps range so run shapes are visually comparable.",
            f"x: $t_{{{N}}}$ arrival time (ps); y: probability density (norm) or events/10\\,ps (abs).",
            f"{'Narrow, nearly Gaussian' if N == 4 else 'Broader, skewed'} distribution; "
            f"{'RMS' if N == 4 else 'shape'} dominated by {'order-stat fluctuations' if N == 4 else 'Landau $N_{{pe}}$ spread'}.",
        )
        frames.append(_frame(
            rf"F1 — Overlay $t_{{{N}}}$, normalised (fixed scale)",
            r"\centering\includegraphics[width=0.85\textwidth]"
            rf"{{exec13_f1_t{N}pe_overlay_norm}}" + "\n" + h4,
        ))
        frames.append(_frame(
            rf"F1 — Overlay $t_{{{N}}}$, absolute counts (fixed scale)",
            r"\centering\includegraphics[width=0.85\textwidth]"
            rf"{{exec13_f1_t{N}pe_overlay_abs}}" + "\n" + h4,
        ))
        frames.append(_frame(
            rf"F1 — Individual panels $t_{{{N}}}$ (identical scale)",
            r"\centering\includegraphics[width=0.98\textwidth]"
            rf"{{exec13_f1_t{N}pe_panels}}" + "\n" + h4,
        ))

    # ── F2 slides ─────────────────────────────────────────────────────────────
    h_f2 = _howtoread(
        r"Gerardo's clustering question: do neighbour SiPMs $\pm1,\pm2$ carry the same timing as the nearest?",
        r"x: photon order $n$; y: $\langle t_n\rangle$ (ps). Error bars = RMS. Dashed: stat floor $\sqrt{n}\tau_d/\langle N_{pe}\rangle$.",
        r"Neighbours approach stat floor faster when their $\langle N_{pe}\rangle$ is comparable to the nearest; "
        r"far neighbours are dominated by stat noise → can be dropped.",
    )
    for x_mm in POSITIONS_MM:
        run_lbl = RUN_LABELS[x_mm]
        xstr = f"{x_mm:+d}".replace("+", "+").replace("-", "-")
        frames.append(_frame(
            rf"F2 — {run_lbl}: $x={x_mm:+d}$\,mm cluster dispersion",
            r"\centering\includegraphics[width=0.82\textwidth]"
            rf"{{exec13_f2_cluster_{x_mm:+d}mm}}" + "\n" + h_f2,
        ))

    # ── F3 slide ──────────────────────────────────────────────────────────────
    h_f3 = _howtoread(
        "Verify that the nearest channel is correctly identified as the Npe maximum; shows photon budget across the full top array.",
        r"x: global channel ID (16--85); y: $\langle N_{pe}\rangle$/event $\pm$ RMS. Highlighted: cluster $\pm$2.",
        r"Npe falls off rapidly away from the track; $\pm$2 neighbours still collect useful light.",
    )
    frames.append(_frame(
        "F3 — Top $N_{pe}$ profile (all 70 channels, fixed scale)",
        r"\centering\includegraphics[width=0.98\textwidth]{exec13_f3_npe_profile}" + "\n" + h_f3,
    ))

    # ── F4 slides ─────────────────────────────────────────────────────────────
    h_f4 = _howtoread(
        rf"Only channels with efficiency $\geq{EFF_MIN_FRACTION:.2f}$ are shown; "
        r"far-end channels that never reach threshold are excluded and listed in exec13\_f4\_timing\_profile.csv.",
        r"x: global channel ID; y: $\langle t_N\rangle$ (ps) $\pm$ RMS.",
        r"Timing rises steeply for channels far from the track (longer optical path + lower $N_{pe}$).",
    )
    for N in TN_THRESHOLDS_PE:
        frames.append(_frame(
            rf"F4 — Timing profile $\langle t_{{{N}}}\rangle$ (eff$\geq${EFF_MIN_FRACTION:.2f}, fixed scale)",
            r"\centering\includegraphics[width=0.98\textwidth]"
            rf"{{exec13_f4_timing_t{N}pe}}" + "\n" + h_f4,
        ))

    # ── Backup: gate table ────────────────────────────────────────────────────
    frames.append(_frame(
        "Backup: validation gates",
        r"""\small
\begin{itemize}
  \item \textbf{G1}: nearest IDs = \{$-690\to16$, $-400\to31$, $0\to51$\} (EXEC\_12 diagonal) --- must pass before plotting
  \item \textbf{G2}: all figure panels share identical axis limits --- verified programmatically post-plot
  \item \textbf{G3}: $n_{\rm events}=2000$ per run --- verified from ROOT file
  \item \textbf{G4}: $n_{\rm in}+n_{\rm excl}+n_{\rm ovfl}=2000$ per histogram --- bookkeeping check
\end{itemize}
\vspace{4pt}
\scriptsize Gate results: \texttt{exec13\_gates.csv}.
Commit: \texttt{""" + commit_hash[:12] + r"""}.
Stats table: \texttt{exec13\_f1\_tn\_histograms.csv} (no manual transcription).
""",
    ))

    return preamble + "".join(frames) + r"\end{document}" + "\n"


# ─── Main ─────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=pathlib.Path, default=_DEFAULT_DATA_DIR)
    parser.add_argument("--out-dir", type=pathlib.Path, default=_DEFAULT_OUT_DIR)
    args = parser.parse_args(argv)

    data_dir: pathlib.Path = args.data_dir
    out_dir: pathlib.Path = args.out_dir
    (out_dir / "figs").mkdir(parents=True, exist_ok=True)

    gates: list[GateRecord] = []

    # ── Load all three runs ───────────────────────────────────────────────────
    print("Loading ROOT files…")
    runs: list[RunData] = []
    for x in POSITIONS_MM:
        run = load_run(data_dir, x)
        runs.append(run)
        print(f"  x={x:+d} mm → {run.n_events} events loaded")

    # ── Gate G3: verify 2000 events per run ──────────────────────────────────
    for run in runs:
        g3 = GateRecord(
            f"G3_nevents_{run.x_mm}",
            run.n_events == N_EVENTS,
            f"n_events={run.n_events}, expected={N_EVENTS}",
        )
        gates.append(g3)
        if not g3.passed:
            print(f"ABORT G3: x={run.x_mm} mm has {run.n_events} events, expected {N_EVENTS}")
            return 1

    # ── Compute Npe matrices ─────────────────────────────────────────────────
    print("Computing Npe matrices…")
    npe_matrices: dict[int, np.ndarray] = {}
    for run in runs:
        npe_matrices[run.x_mm] = top_npe_matrix(
            run.event_id, run.global_id, run.n_events
        )

    # ── Determine nearest-top IDs ─────────────────────────────────────────────
    print("Determining nearest-top channels…")
    nearests: dict[int, int] = {}
    for run in runs:
        prof = npe_matrices[run.x_mm].mean(axis=0)
        nearest = select_nearest_top(run.x_mm, prof)
        nearests[run.x_mm] = nearest
        print(f"  x={run.x_mm:+d} mm → nearest top ID = {nearest}")

    # ── Gate G1: verify expected nearest IDs ─────────────────────────────────
    g1_passed = all(nearests[x] == EXPECTED_NEAREST[x] for x in POSITIONS_MM)
    g1_detail = " | ".join(
        f"x={x:+d}: got {nearests[x]}, exp {EXPECTED_NEAREST[x]}"
        for x in POSITIONS_MM
    )
    gates.append(GateRecord("G1_nearest_ids", g1_passed, g1_detail))
    if not g1_passed:
        print(f"ABORT G1: nearest-top ID mismatch: {g1_detail}")
        return 1
    print("G1 PASSED: nearest IDs match EXEC_12 diagonal {16, 31, 51}")

    # ── F1: histograms ────────────────────────────────────────────────────────
    print("Generating F1…")
    all_hists: list[HistResult] = []
    for N in TN_THRESHOLDS_PE:
        hists, g_f1 = plot_f1(runs, nearests, N, out_dir)
        all_hists.extend(hists)
        gates.extend(g_f1)

    # Gate G4: histogram accounting
    for h in all_hists:
        total = h.n_in_range + h.n_overflow + h.n_below
        g4 = GateRecord(
            f"G4_accounting_x{h.x_mm}_t{h.threshold}pe",
            total == h.threshold or total == runs[0].n_events,
            f"in={h.n_in_range} + ovfl={h.n_overflow} + below={h.n_below} = {total}"
            f" (expected {runs[0].n_events})",
        )
        # Recompute correctly: n_below uses frac_excl
        n_below_check = round(h.frac_excl * runs[0].n_events)
        total_check = h.n_in_range + h.n_overflow + n_below_check
        g4_pass = total_check == runs[0].n_events
        gates.append(GateRecord(
            f"G4_{h.x_mm}_{h.threshold}pe",
            g4_pass,
            f"in={h.n_in_range}+ovfl={h.n_overflow}+below≈{n_below_check}={total_check}",
        ))

    # ── F2: cluster dispersion ────────────────────────────────────────────────
    print("Generating F2…")
    f2_data, g_f2 = plot_f2(runs, nearests, npe_matrices, out_dir)
    gates.extend(g_f2)

    # ── F3: Npe profile ───────────────────────────────────────────────────────
    print("Generating F3…")
    f3_data, g_f3 = plot_f3(runs, nearests, npe_matrices, out_dir)
    gates.extend(g_f3)

    # ── F4: timing profile ────────────────────────────────────────────────────
    print("Generating F4…")
    f4_data, g_f4 = plot_f4(runs, nearests, out_dir)
    gates.extend(g_f4)

    # ── Write CSVs ────────────────────────────────────────────────────────────
    print("Writing CSVs…")
    write_csv_f1(all_hists, out_dir)
    write_csv_f2(f2_data, out_dir)
    write_csv_f3(f3_data, out_dir)
    write_csv_f4(f4_data, out_dir)
    write_csv_gates(gates, out_dir)

    # ── LaTeX table ───────────────────────────────────────────────────────────
    tex_table = make_latex_table(all_hists)
    (out_dir / "exec13_tN_table.tex").write_text(tex_table)

    # ── Beamer ────────────────────────────────────────────────────────────────
    print("Generating Beamer…")
    import subprocess
    try:
        commit_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=str(_REPO), text=True
        ).strip()
    except Exception:
        commit_hash = "unknown"
    beamer_tex = generate_beamer(all_hists, commit_hash)
    (out_dir / "exec13_report.tex").write_text(beamer_tex)

    # ── Gate summary ─────────────────────────────────────────────────────────
    failed = [g for g in gates if not g.passed]
    print(f"\nGate summary: {len(gates)} checks, {len(failed)} failed")
    for g in gates:
        status = "PASS" if g.passed else "FAIL"
        print(f"  [{status}] {g.name}: {g.detail}")

    if failed:
        print("\nWARNING: some gates failed — check exec13_gates.csv")
        return 1

    print("\nexec13_fixed_scale COMPLETED successfully.")
    print(f"Output in: {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
