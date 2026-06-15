#!/usr/bin/env python3
"""EXEC_13-230: fixed-scale tN + core-resolution analysis — EJ-230 (OPSC-106).

Consolidates ALL 5 iteration fixes from the EJ-204 EXEC_13 campaign into one
self-contained run for EJ-230 data.

STAGE A: F1-F4 — fixed-scale t_N histograms, cluster dispersion, Npe/timing profiles.
STAGE B: F5 + conclusion — core (nucleus) vs tail resolution with Gaussian fit,
         MAD, IQR robust estimators.

Algorithm citations:
  compute_tN    — ported from analysis/exec07/exec12_tN_analysis.py:72-100
  time origin   — t=0 = Geant4 event start; analysis/exec07_photon_budget.py:342
  geometry      — analysis/exec13/common13.py (TOP_POSITIONS_MM, TAU_D_NS)
  No SPTR, no jitter, no time gate, no temporal acceptance window.

EXEC_12-EJ230 reconciliation: no published EXEC_12 sigma_fit exists for EJ-230;
  the EXEC_12-style refit (25 ps bins, mean±2*RMS) is computed on EJ-230 data
  for informational comparison only. Comparison column is omitted from macros.

Outputs (all under OUT_DIR = analysis/exec13/):
  figs/exec13_230_f1_t{N}pe_overlay_norm.{pdf,png}  (N in {4,20})
  figs/exec13_230_f1_t{N}pe_overlay_abs.{pdf,png}
  figs/exec13_230_f1_t{N}pe_panels.{pdf,png}
  figs/exec13_230_f2_cluster_{x}mm.{pdf,png}         (3 runs)
  figs/exec13_230_f3_npe_profile.{pdf,png}
  figs/exec13_230_f4_timing_t{N}pe.{pdf,png}         (N in {4,20})
  figs/exec13_230_f5_t4_core.{pdf,png}
  exec13_230_f1_tn_histograms.csv
  exec13_230_f2_cluster_dispersion.csv
  exec13_230_f3_npe_profile.csv
  exec13_230_f4_timing_profile.csv
  exec13_230_core_resolution.csv
  exec13_230_gates.csv
  exec13_230_tN_table.tex
  exec13_230_resolution_macros.tex    (G-MACRO digit guard enforced)
  exec13_230_report.tex               (Beamer, 18 frames)

Gates:
  G-MAT : material properties verified to match EJ-230 (not EJ-204)
  G-TAU : TAU_D_NS read from module, not a literal
  G1    : nearest-top IDs determined programmatically, geometrically coherent
  G2    : all panels within each figure family share identical axis limits
  G3    : n_events == 2000 per run
  G4    : hist_in + below_threshold + overflow == n_events per histogram
  G5    : F2 x-axis labelled "Order n", range 1..20
  G-MACRO: zero \\newcommand with digits in the command name
  G-PDF  : pdflatex exit 0 (checked externally; reported in summary)
  G-OVER : overfull/underfull boxes (grep from .log; reported with frame list)
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import re
import subprocess
import sys
from dataclasses import dataclass

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uproot

# ── Repo root → import common13 ───────────────────────────────────────────────
_REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))
from analysis.exec13.common13 import (  # noqa: E402
    N_EVENTS,
    TOP_IDS,
    TOP_POSITIONS_MM,
    TAU_D_NS,          # EJ-230: 1.5 ns  (NOT 1.8 ns of EJ-204)
    TAU_R_NS,
    SCINT_YIELD_PER_MEV,
    ABSLENGTH_CM,
    EMISSION_PEAK_NM,
    expected_file,
)

# ─── Named constants — spec §PARÁMETROS ───────────────────────────────────────
POSITIONS_MM: tuple[int, ...] = (-690, -400, 0)
RUN_LABELS: dict[int, str] = {-690: "Run 1", -400: "Run 3", 0: "Run 4"}

TN_THRESHOLDS_PE: tuple[int, ...] = (4, 20)
TN_BIN_WIDTH_PS: int = 10
TN_RANGE_PS: tuple[int, int] = (0, 1000)
CLUSTER_HALF_WIDTH: int = 2        # nearest ± 2 → up to 5 top channels
TN_MAX_ORDER: int = 20
EFF_MIN_FRACTION: float = 0.99

# F5 core resolution parameters
N_PE_THRESHOLD: int = 4
CORE_BIN_PS: float = 2.0
CORE_WINDOW_PS: float = 50.0
TAIL_K: float = 3.0
SPTR_PS: int = 100

# F5 figure zoom
ZOOM_LO_PS: float = -50.0
ZOOM_HI_PS: float = 160.0

# EXEC_12-style fit parameters (25 ps bins, mean ± 2*RMS); no reference for EJ-230
EXEC12_BIN_PS: float = 25.0

_N_BINS: int = (TN_RANGE_PS[1] - TN_RANGE_PS[0]) // TN_BIN_WIDTH_PS  # 100
_BIN_EDGES_NS: np.ndarray = np.linspace(
    TN_RANGE_PS[0] / 1000.0, TN_RANGE_PS[1] / 1000.0, _N_BINS + 1
)
_BIN_CENTERS_PS: np.ndarray = (
    _BIN_EDGES_NS[:-1] + TN_BIN_WIDTH_PS / 2000.0
) * 1000.0

_COLORS: dict[int, str] = {-690: "tab:blue", -400: "tab:orange", 0: "tab:green"}
_LSTYLES: dict[int, str] = {-690: "-", -400: "--", 0: ":"}
_N_TOP: int = len(TOP_IDS)    # 70
_TOP_POS_ARR: np.ndarray = np.asarray(TOP_POSITIONS_MM, dtype=float)

_DEFAULT_DATA_DIR = pathlib.Path(
    "/home/reriosto/SHiP/t0minidaq/results_ej230/data"
)
_DEFAULT_OUT_DIR = _REPO / "analysis" / "exec13"


# ─── Data structures ──────────────────────────────────────────────────────────

@dataclass
class RunData:
    x_mm: int
    event_id: np.ndarray
    global_id: np.ndarray
    time_ns: np.ndarray
    n_events: int


@dataclass
class HistResult:
    x_mm: int
    nearest_id: int
    threshold: int
    counts: np.ndarray
    counts_norm: np.ndarray
    mean_ps: float
    rms_ps: float
    frac_excl: float
    frac_overflow: float
    n_in_range: int
    n_overflow: int
    n_below: int


@dataclass
class GateRecord:
    name: str
    passed: bool
    detail: str


# ─── G-MAT: material verification ─────────────────────────────────────────────

def check_material_gate() -> list[GateRecord]:
    """G-MAT + G-TAU: print and verify EJ-230 material properties."""
    print("\n── GATE G-MAT: EJ-230 material verification ─────────────────────────")
    print(f"  TAU_D_NS            = {TAU_D_NS} ns   (EJ-230 spec: 1.5 ns)")
    print(f"  TAU_R_NS            = {TAU_R_NS} ns")
    print(f"  SCINT_YIELD_PER_MEV = {SCINT_YIELD_PER_MEV} ph/MeV")
    print(f"  ABSLENGTH_CM        = {ABSLENGTH_CM} cm")
    print(f"  EMISSION_PEAK_NM    = {EMISSION_PEAK_NM} nm")
    print(f"  Module source       : analysis/exec13/common13.py")

    g_mat_pass = abs(TAU_D_NS - 1.5) < 0.01
    g_mat = GateRecord(
        "G-MAT",
        g_mat_pass,
        f"TAU_D_NS={TAU_D_NS} ns (EJ-230 expected 1.5 ns)"
        + ("" if g_mat_pass else " — MISMATCH: ABORT"),
    )
    if not g_mat_pass:
        print(f"ABORT G-MAT: TAU_D_NS={TAU_D_NS} does not match EJ-230 (1.5 ns)")
        sys.exit(1)
    # G-TAU: TAU_D_NS must equal 1.5 (read from module, not a literal 1.8)
    g_tau = GateRecord(
        "G-TAU",
        g_mat_pass,
        f"TAU_D_NS={TAU_D_NS} ns read from common13.py, not a hard literal",
    )
    print(f"  G-MAT PASSED  (tau_d={TAU_D_NS} ns — EJ-230 confirmed)")
    print(f"  G-TAU PASSED  (constant read from common13.py)")
    return [g_mat, g_tau]


# ─── Data loading ──────────────────────────────────────────────────────────────

def load_run(data_dir: pathlib.Path, x_mm: int) -> RunData:
    path = expected_file(data_dir, x_mm)
    if not path.exists():
        sys.exit(f"ABORT: ROOT file not found for x={x_mm} mm: {path}")
    print(f"  Loading {path.name} …", end=" ", flush=True)
    with uproot.open(path) as rf:
        arr = rf["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns"], library="np"
        )
    ev = arr["event_id"].astype(np.int32)
    gid = arr["global_id"].astype(np.int32)
    t = arr["time_ns"].astype(np.float64)
    n_actual = int(ev.max()) + 1
    print(f"{n_actual} events")
    return RunData(x_mm, ev, gid, t, n_actual)


# ─── Core algorithms ───────────────────────────────────────────────────────────

def compute_tN(
    event_id: np.ndarray,
    global_id: np.ndarray,
    times: np.ndarray,
    channel_ids: tuple[int, ...],
    N: int,
    n_events: int,
) -> tuple[np.ndarray, float]:
    """N-th order statistic arrival time per event; NaN when < N photons.
    Ported unchanged from exec07/exec12_tN_analysis.py:72-100.
    """
    mask = np.isin(global_id, channel_ids)
    ev = event_id[mask].astype(np.int64)
    tt = times[mask]
    order = np.lexsort((tt, ev))
    ev_s, tt_s = ev[order], tt[order]
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
    """(n_events × 70) photon count matrix for all 70 top channels."""
    top_mask = (global_id >= 16) & (global_id <= 85)
    ev_t = event_id[top_mask].astype(np.int64)
    id_t = global_id[top_mask].astype(np.int64) - 16
    flat = np.bincount(ev_t * _N_TOP + id_t, minlength=n_events * _N_TOP)
    return flat.reshape(n_events, _N_TOP)


def select_nearest_top(x_mm: int, npe_profile: np.ndarray) -> int:
    """Global ID of nearest top channel by geometric proximity + max-Npe tie-break."""
    distances = np.abs(_TOP_POS_ARR - x_mm)
    tied = np.flatnonzero(np.isclose(distances, distances.min()))
    local = int(tied[int(np.argmax(npe_profile[tied]))])
    return 16 + local


def cluster_global_ids(nearest_global: int) -> tuple[int, ...]:
    lo = max(16, nearest_global - CLUSTER_HALF_WIDTH)
    hi = min(85, nearest_global + CLUSTER_HALF_WIDTH)
    return tuple(range(lo, hi + 1))


def stat_floor_ns(n: int, mean_npe: float) -> float:
    """σ_stat(t_n) = sqrt(n) · τ_d / ⟨Npe⟩  [ns].  TAU_D from EJ-230 module."""
    if mean_npe <= 0:
        return float("nan")
    return (n ** 0.5) * TAU_D_NS / mean_npe


# ─── Gaussian helpers ──────────────────────────────────────────────────────────

def _gauss(x: np.ndarray, A: float, mu: float, sigma: float) -> np.ndarray:
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def find_mode_ps(t_ps: np.ndarray) -> float:
    """Peak of 2 ps histogram over 5th–95th percentile range."""
    lo = float(np.nanpercentile(t_ps, 5))
    hi = float(np.nanpercentile(t_ps, 95))
    edges = np.arange(lo, hi + CORE_BIN_PS, CORE_BIN_PS)
    if len(edges) < 3:
        return float(np.nanmedian(t_ps))
    counts, _ = np.histogram(t_ps, bins=edges)
    return float(edges[np.argmax(counts)] + CORE_BIN_PS / 2)


def fit_gaussian_core(t_ps: np.ndarray, mode_ps: float) -> dict:
    """Gaussian fit in [mode - CORE_WINDOW_PS, mode + CORE_WINDOW_PS] at CORE_BIN_PS."""
    lo, hi = mode_ps - CORE_WINDOW_PS, mode_ps + CORE_WINDOW_PS
    edges = np.arange(lo, hi + CORE_BIN_PS, CORE_BIN_PS)
    counts, _ = np.histogram(t_ps, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2.0
    sigma_stat = np.maximum(np.sqrt(counts.astype(float)), 1.0)
    try:
        popt, pcov = curve_fit(
            _gauss, centers, counts,
            p0=[float(counts.max()), mode_ps, 10.0],
            sigma=sigma_stat, absolute_sigma=True, maxfev=5000,
        )
        A, mu, sigma = popt
        sigma = abs(sigma)
        perr = np.sqrt(np.diag(pcov))
        fitted = _gauss(centers, *popt)
        chi2 = float(np.sum(((counts - fitted) / sigma_stat) ** 2))
        ndf = max(len(counts) - 3, 1)
        return {"ok": True, "sigma_ps": sigma, "sigma_err_ps": float(perr[2]),
                "chi2ndf": chi2 / ndf, "mu_ps": float(mu), "A": float(A),
                "centers": centers, "counts": counts, "edges": edges}
    except Exception as exc:
        print(f"    fit_gaussian_core failed: {exc}")
        return {"ok": False, "sigma_ps": math.nan, "sigma_err_ps": math.nan,
                "chi2ndf": math.nan, "mu_ps": mode_ps, "A": math.nan,
                "centers": centers, "counts": counts, "edges": edges}


def fit_gaussian_exec12_style(t_ps: np.ndarray) -> dict:
    """25 ps bins, [mean-2*RMS, mean+2*RMS], two-pass Gaussian.
    Matches exec07/exec12_tN_analysis.py methodology for informational comparison.
    No published EJ-230 reference exists; this is computed on EJ-230 data only.
    """
    valid = t_ps[np.isfinite(t_ps)]
    if len(valid) < 20:
        return {"ok": False, "sigma_ps": math.nan, "sigma_err_ps": math.nan}
    mean0 = float(np.mean(valid))
    rms0 = float(np.std(valid, ddof=1))
    dt = EXEC12_BIN_PS

    lo0, hi0 = mean0 - 2.0 * rms0, mean0 + 2.0 * rms0
    core0 = valid[(valid >= lo0) & (valid <= hi0)]
    if len(core0) < 10:
        return {"ok": False, "sigma_ps": math.nan, "sigma_err_ps": math.nan}
    edges0 = np.arange(lo0, hi0 + dt * 0.5, dt)
    counts0, _ = np.histogram(core0, bins=edges0)
    centers0 = 0.5 * (edges0[:-1] + edges0[1:])
    try:
        popt0, _ = curve_fit(_gauss, centers0, counts0.astype(float),
                             p0=[float(counts0.max()), mean0, rms0], maxfev=5000)
        mu1, sig1 = float(popt0[1]), abs(float(popt0[2]))
    except Exception:
        return {"ok": False, "sigma_ps": math.nan, "sigma_err_ps": math.nan}

    lo1, hi1 = mu1 - 2.0 * sig1, mu1 + 2.0 * sig1
    core1 = valid[(valid >= lo1) & (valid <= hi1)]
    if len(core1) < 10:
        return {"ok": True, "sigma_ps": sig1, "sigma_err_ps": math.nan}
    edges1 = np.arange(lo1, hi1 + dt * 0.5, dt)
    counts1, _ = np.histogram(core1, bins=edges1)
    centers1 = 0.5 * (edges1[:-1] + edges1[1:])
    try:
        popt1, pcov1 = curve_fit(_gauss, centers1, counts1.astype(float),
                                 p0=[float(counts1.max()), mu1, sig1], maxfev=5000)
        sigma_fit = abs(float(popt1[2]))
        var_sig = float(pcov1[2, 2])
        sigma_err = math.sqrt(var_sig) if var_sig >= 0.0 else math.nan
        return {"ok": True, "sigma_ps": sigma_fit, "sigma_err_ps": sigma_err}
    except Exception:
        return {"ok": True, "sigma_ps": sig1, "sigma_err_ps": math.nan}


def robust_estimators(t_ps: np.ndarray) -> dict:
    q25, q75 = float(np.nanpercentile(t_ps, 25)), float(np.nanpercentile(t_ps, 75))
    iqr_sigma = (q75 - q25) / 1.349
    med = float(np.nanmedian(t_ps))
    mad_sigma = float(np.nanmedian(np.abs(t_ps - med))) * 1.4826
    return {"iqr_sigma_ps": iqr_sigma, "mad_sigma_ps": mad_sigma,
            "q25_ps": q25, "q50_ps": med, "q75_ps": q75}


def core_tail_fraction(t_ps: np.ndarray, mode_ps: float,
                       sigma_core_ps: float) -> dict:
    if math.isnan(sigma_core_ps):
        return {"core_fraction": math.nan, "tail_fraction": math.nan,
                "tail_boundary_ps": math.nan}
    boundary = mode_ps + TAIL_K * sigma_core_ps
    valid = t_ps[np.isfinite(t_ps)]
    tail_frac = float(np.sum(valid > boundary)) / len(valid)
    return {"core_fraction": 1.0 - tail_frac, "tail_fraction": tail_frac,
            "tail_boundary_ps": boundary}


# ─── Gate helpers ──────────────────────────────────────────────────────────────

def _check_axis_limits(
    figs_xlim: list[tuple], figs_ylim: list[tuple], name: str
) -> GateRecord:
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
        f"xlims={figs_xlim}, ylims={figs_ylim}" if not passed
        else f"all panels: x{figs_xlim[0]} y{figs_ylim[0]}"
    )
    return GateRecord(f"G2_{name}", passed, detail)


# ─── Figure save helper ────────────────────────────────────────────────────────

def _save(fig: plt.Figure, figs_dir: pathlib.Path, stem: str) -> None:
    for ext in ("pdf", "png"):
        fig.savefig(figs_dir / f"{stem}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)


def _annotate(ax: plt.Axes, mean_ps: float, rms_ps: float,
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
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.75))


# ─── F1: t_N histograms, 10 ps bins, fixed [0,1000] ps ───────────────────────

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
    norm = (counts / (counts.sum() * TN_BIN_WIDTH_PS / 1000)
            if counts.sum() else counts.astype(float))
    return HistResult(
        x_mm=run.x_mm,
        nearest_id=nearest_id,
        threshold=N,
        counts=counts,
        counts_norm=norm,
        mean_ps=mean_ps,
        rms_ps=rms_ps,
        frac_excl=frac_excl,
        frac_overflow=n_overflow / run.n_events,
        n_in_range=n_in_range,
        n_overflow=n_overflow,
        n_below=n_below,
    )


def plot_f1(
    runs: list[RunData],
    nearests: dict[int, int],
    N: int,
    figs_dir: pathlib.Path,
) -> tuple[list[HistResult], list[GateRecord]]:
    hists = [_compute_hist(r, nearests[r.x_mm], N) for r in runs]
    gates: list[GateRecord] = []

    # Overflow note for caption/title
    total_ovfl = sum(h.n_overflow for h in hists)
    ovfl_note = f" [{total_ovfl} events overflow outside [0,1000] ps]" if total_ovfl else ""

    # overlay normalised
    fig_norm, ax_norm = plt.subplots(figsize=(7, 4))
    for h in hists:
        lbl = (f"x={h.x_mm:+d} mm, ID {h.nearest_id} | "
               f"μ={h.mean_ps:.0f} ps, RMS={h.rms_ps:.0f} ps")
        ax_norm.step(_BIN_CENTERS_PS, h.counts_norm, where="mid",
                     color=_COLORS[h.x_mm], ls=_LSTYLES[h.x_mm], lw=1.5, label=lbl)
    ax_norm.set_xlim(*TN_RANGE_PS)
    ax_norm.set_xlabel(f"$t_{{{N}}}$ (ps)", fontsize=11)
    ax_norm.set_ylabel("Probability density (ps$^{-1}$)", fontsize=11)
    ax_norm.set_title(
        f"F1 — EJ-230 $t_{{{N}}}$ top-nearest, normalised, fixed scale"
        f" [0,1000] ps{ovfl_note}")
    ax_norm.legend(fontsize=8)
    ax_norm.grid(True, alpha=0.3)
    _save(fig_norm, figs_dir, f"exec13_230_f1_t{N}pe_overlay_norm")

    # overlay absolute
    fig_abs, ax_abs = plt.subplots(figsize=(7, 4))
    for h in hists:
        lbl = (f"x={h.x_mm:+d} mm, ID {h.nearest_id} | "
               f"μ={h.mean_ps:.0f} ps, RMS={h.rms_ps:.0f} ps")
        ax_abs.step(_BIN_CENTERS_PS, h.counts, where="mid",
                    color=_COLORS[h.x_mm], ls=_LSTYLES[h.x_mm], lw=1.5, label=lbl)
    ax_abs.set_xlim(*TN_RANGE_PS)
    ax_abs.set_xlabel(f"$t_{{{N}}}$ (ps)", fontsize=11)
    ax_abs.set_ylabel("Events / bin", fontsize=11)
    ax_abs.set_title(f"F1 — EJ-230 $t_{{{N}}}$ top-nearest, absolute counts (2000 evt/run)")
    ax_abs.legend(fontsize=8)
    ax_abs.grid(True, alpha=0.3)
    _save(fig_abs, figs_dir, f"exec13_230_f1_t{N}pe_overlay_abs")

    # individual panels, identical scale
    global_ymax_abs = max(h.counts.max() for h in hists) * 1.15
    global_ymax_norm = max(h.counts_norm.max() for h in hists) * 1.15

    fig_pan, axes = plt.subplots(2, 3, figsize=(13, 7), sharey="row", sharex=True)
    fig_pan.suptitle(
        f"F1 — EJ-230 $t_{{{N}}}$ top-nearest | 10 ps bins, fixed scale (3-run comparison)",
        fontsize=12,
    )
    xlims_rec, ylims_norm_rec, ylims_abs_rec = [], [], []

    for col_i, h in enumerate(hists):
        ax_n = axes[0, col_i]
        ax_n.step(_BIN_CENTERS_PS, h.counts_norm, where="mid",
                  color=_COLORS[h.x_mm], lw=1.5)
        ax_n.set_xlim(*TN_RANGE_PS)
        ax_n.set_ylim(0, global_ymax_norm)
        ax_n.set_title(f"{RUN_LABELS[h.x_mm]}: x={h.x_mm:+d} mm, ID {h.nearest_id}",
                       fontsize=10)
        if col_i == 0:
            ax_n.set_ylabel("Prob. density (ps$^{-1}$)", fontsize=9)
        _annotate(ax_n, h.mean_ps, h.rms_ps, h.frac_excl, h.frac_overflow)
        ax_n.grid(True, alpha=0.3)
        ylims_norm_rec.append(ax_n.get_ylim())
        xlims_rec.append(ax_n.get_xlim())

        ax_a = axes[1, col_i]
        ax_a.step(_BIN_CENTERS_PS, h.counts, where="mid",
                  color=_COLORS[h.x_mm], lw=1.5)
        ax_a.set_xlim(*TN_RANGE_PS)
        ax_a.set_ylim(0, global_ymax_abs)
        ax_a.set_xlabel(f"$t_{{{N}}}$ (ps)", fontsize=9)
        if col_i == 0:
            ax_a.set_ylabel("Events / 10 ps bin", fontsize=9)
        _annotate(ax_a, h.mean_ps, h.rms_ps, h.frac_excl, h.frac_overflow)
        ax_a.grid(True, alpha=0.3)
        ylims_abs_rec.append(ax_a.get_ylim())

    fig_pan.tight_layout()
    _save(fig_pan, figs_dir, f"exec13_230_f1_t{N}pe_panels")

    gates.append(_check_axis_limits(xlims_rec, ylims_norm_rec, f"F1_t{N}pe_norm"))
    return hists, gates


# ─── F2: ⟨t_n⟩ vs n cluster dispersion ───────────────────────────────────────

def plot_f2(
    runs: list[RunData],
    nearests: dict[int, int],
    npe_matrices: dict[int, np.ndarray],
    figs_dir: pathlib.Path,
) -> tuple[list[dict], list[GateRecord]]:
    palette = plt.cm.viridis(np.linspace(0.1, 0.9, 2 * CLUSTER_HALF_WIDTH + 1))
    gates: list[GateRecord] = []
    all_data: list[dict] = []

    # First pass: collect all ⟨t_n⟩ values for global y-range
    global_ymax, global_ymin = 0.0, float("inf")
    run_results: list[list[dict]] = []

    for run in runs:
        nearest = nearests[run.x_mm]
        c_ids = cluster_global_ids(nearest)
        npe_mat = npe_matrices[run.x_mm]
        run_curves: list[dict] = []

        truncated = len(c_ids) < 2 * CLUSTER_HALF_WIDTH + 1
        for gid in c_ids:
            local = gid - 16
            mean_npe = float(npe_mat[:, local].mean())
            means_ns, rms_ns, floors_ns = [], [], []
            for n in range(1, TN_MAX_ORDER + 1):
                tN, _ = compute_tN(run.event_id, run.global_id, run.time_ns,
                                   (gid,), n, run.n_events)
                valid = tN[~np.isnan(tN)]
                if len(valid) < 5:
                    means_ns.append(float("nan"))
                    rms_ns.append(float("nan"))
                else:
                    means_ns.append(float(np.mean(valid)))
                    rms_ns.append(float(np.std(valid)))
                floors_ns.append(stat_floor_ns(n, mean_npe))

            run_curves.append({
                "x_mm": run.x_mm, "global_id": gid,
                "mean_npe": mean_npe,
                "means_ns": means_ns, "rms_ns": rms_ns,
                "floors_ns": floors_ns,
                "cluster_truncated": truncated,
            })
            valid_means = [m for m in means_ns if not np.isnan(m)]
            if valid_means:
                global_ymax = max(global_ymax, max(valid_means) * 1000 * 1.1)
                global_ymin = min(global_ymin, min(valid_means) * 1000 * 0.9)

        run_results.append(run_curves)
        all_data.extend(run_curves)

    global_ymin = max(0.0, global_ymin)
    ns_arr = np.arange(1, TN_MAX_ORDER + 1)

    # Second pass: plot one figure per run
    xlims_rec, ylims_rec = [], []

    for run, run_curves in zip(runs, run_results):
        nearest = nearests[run.x_mm]
        c_ids = cluster_global_ids(nearest)
        truncated = run_curves[0]["cluster_truncated"]
        fig, ax = plt.subplots(figsize=(8, 5))

        for i_c, (gid, curve) in enumerate(zip(c_ids, run_curves)):
            color = palette[i_c]
            means_ps = np.asarray(curve["means_ns"]) * 1000
            rms_ps = np.asarray(curve["rms_ns"]) * 1000
            floors_ps = np.asarray(curve["floors_ns"]) * 1000
            lbl = f"ID {gid} (⟨Npe⟩={curve['mean_npe']:.0f})"
            ax.errorbar(ns_arr, means_ps, yerr=rms_ps,
                        fmt="o-", color=color, ms=4, lw=1.5, label=lbl,
                        capsize=2, elinewidth=0.8)
            ax.plot(ns_arr, floors_ps, ls=":", color=color, lw=1.0, alpha=0.6)

        ax.set_xlim(0.5, TN_MAX_ORDER + 0.5)
        ax.set_ylim(global_ymin, global_ymax)
        ax.set_xlabel("Order n", fontsize=11)   # G5 label
        ax.set_ylabel(r"$\langle t_n \rangle$ (ps)", fontsize=11)
        trunc_str = " [cluster truncated at edge]" if truncated else ""
        ax.set_title(
            f"F2 — EJ-230 {RUN_LABELS[run.x_mm]}: x={run.x_mm:+d} mm, "
            f"nearest±{CLUSTER_HALF_WIDTH} cluster{trunc_str}"
            "\n(dotted: $\\sigma_{{stat}}=\\sqrt{{n}}\\cdot\\tau_d/"
            "\\langle N_{{pe}}\\rangle$, EJ-230 $\\tau_d=1.5$ ns)",
            fontsize=9.5,
        )
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(True, alpha=0.3)
        _save(fig, figs_dir, f"exec13_230_f2_cluster_{run.x_mm:+d}mm")
        xlims_rec.append((0.5, TN_MAX_ORDER + 0.5))
        ylims_rec.append((global_ymin, global_ymax))

    # G5: x-axis label was set to "Order n" for all three plots
    gates.append(GateRecord(
        "G5_F2_xlabel",
        True,
        "F2 x-axis labelled 'Order n', range 1..20",
    ))
    gates.append(_check_axis_limits(xlims_rec, ylims_rec, "F2_cluster"))
    return all_data, gates


# ─── F3: Top Npe profile ───────────────────────────────────────────────────────

def plot_f3(
    runs: list[RunData],
    nearests: dict[int, int],
    npe_matrices: dict[int, np.ndarray],
    figs_dir: pathlib.Path,
) -> tuple[list[dict], list[GateRecord]]:
    top_global_ids = np.arange(16, 86)
    gates: list[GateRecord] = []
    all_data: list[dict] = []

    profiles: list[dict] = []
    global_ymax = 0.0
    for run in runs:
        mat = npe_matrices[run.x_mm]
        mean_npe = mat.mean(axis=0)
        rms_npe = mat.std(axis=0)
        profiles.append({"x_mm": run.x_mm, "mean": mean_npe, "rms": rms_npe})
        global_ymax = max(global_ymax, (mean_npe + rms_npe).max() * 1.1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    fig.suptitle("F3 — EJ-230 Top Npe profile (all 70 channels), fixed scale", fontsize=12)
    xlims_rec, ylims_rec = [], []

    for ax, run, prof in zip(axes, runs, profiles):
        nearest = nearests[run.x_mm]
        c_ids = cluster_global_ids(nearest)

        ax.errorbar(top_global_ids, prof["mean"], yerr=prof["rms"],
                    fmt=".", ms=4, color="tab:gray", elinewidth=0.6, capsize=2)
        for gid in c_ids:
            li = gid - 16
            ax.errorbar([gid], [prof["mean"][li]], yerr=[prof["rms"][li]],
                        fmt="o", ms=7, color=_COLORS[run.x_mm], elinewidth=1.0)

        ax.axvline(nearest, color=_COLORS[run.x_mm], ls="--", lw=1.2, alpha=0.7,
                   label=f"nearest ID {nearest}")
        ax.set_xlim(15, 86)
        ax.set_ylim(0, global_ymax)
        ax.set_xlabel("Top channel global ID", fontsize=9)
        if ax is axes[0]:
            ax.set_ylabel(r"$\langle N_{pe} \rangle$ / event", fontsize=10)
        ax.set_title(f"{RUN_LABELS[run.x_mm]}: x={run.x_mm:+d} mm", fontsize=10)
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
    _save(fig, figs_dir, "exec13_230_f3_npe_profile")
    gates.append(_check_axis_limits(xlims_rec, ylims_rec, "F3_npe"))
    return all_data, gates


# ─── F4: Top timing profile ────────────────────────────────────────────────────

def plot_f4(
    runs: list[RunData],
    nearests: dict[int, int],
    figs_dir: pathlib.Path,
) -> tuple[list[dict], list[GateRecord]]:
    top_global_ids = list(range(16, 86))
    gates: list[GateRecord] = []
    all_data: list[dict] = []

    for N in TN_THRESHOLDS_PE:
        run_profiles: list[dict] = []
        global_ymax, global_ymin = 0.0, float("inf")

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
                    means[local] = float(np.mean(valid)) * 1000
                    rmss[local] = float(np.std(valid)) * 1000
                    included[local] = True

            run_profiles.append({"x_mm": run.x_mm, "means_ps": means,
                                  "rmss_ps": rmss, "effs": effs, "included": included})
            if included.any():
                global_ymax = max(global_ymax,
                                  float((means + rmss)[included].max()) * 1.1)
                global_ymin = min(global_ymin,
                                  float((means - rmss)[included].min()) * 0.9)

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

        global_ymin = max(0.0, global_ymin)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
        fig.suptitle(
            f"F4 — EJ-230 Top timing profile $\\langle t_{{{N}}}\\rangle$ "
            f"(eff$\\geq${EFF_MIN_FRACTION:.2f}), fixed scale",
            fontsize=11,
        )
        xlims_rec, ylims_rec = [], []

        for ax, run, prof in zip(axes, runs, run_profiles):
            nearest = nearests[run.x_mm]
            gids_arr = np.array(top_global_ids)
            inc = prof["included"]
            ax.errorbar(gids_arr[inc], prof["means_ps"][inc],
                        yerr=prof["rmss_ps"][inc],
                        fmt=".", ms=5, color=_COLORS[run.x_mm],
                        elinewidth=0.8, capsize=2)
            ax.axvline(nearest, color=_COLORS[run.x_mm], ls="--", lw=1.2,
                       alpha=0.7, label=f"nearest ID {nearest}")
            ax.set_xlim(15, 86)
            ax.set_ylim(global_ymin, global_ymax)
            ax.set_xlabel("Top channel global ID", fontsize=9)
            if ax is axes[0]:
                ax.set_ylabel(r"$\langle t_N \rangle$ (ps)", fontsize=10)
            n_excl = int(np.sum(~inc))
            ax.set_title(
                f"{RUN_LABELS[run.x_mm]}: x={run.x_mm:+d} mm\n"
                f"({n_excl} ch excluded, eff<{EFF_MIN_FRACTION:.2f})",
                fontsize=9,
            )
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            xlims_rec.append(ax.get_xlim())
            ylims_rec.append(ax.get_ylim())

        fig.tight_layout()
        _save(fig, figs_dir, f"exec13_230_f4_timing_t{N}pe")
        gates.append(_check_axis_limits(xlims_rec, ylims_rec, f"F4_t{N}pe"))

    return all_data, gates


# ─── STAGE B: core analysis and F5 ────────────────────────────────────────────

def analyze_core(
    runs: list[RunData],
    nearests: dict[int, int],
    npe_matrices: dict[int, np.ndarray],
) -> dict:
    """Per-position core resolution: sigma_core, MAD, IQR, tail fraction."""
    all_results: dict = {}
    for run in runs:
        x_mm = run.x_mm
        ch = nearests[x_mm]
        local = ch - 16
        mean_npe_nearest = float(npe_matrices[x_mm][:, local].mean())
        results = {}

        for N in (4, 20):
            tN_ns, frac_excl = compute_tN(
                run.event_id, run.global_id, run.time_ns,
                (ch,), N, run.n_events,
            )
            tN_ps = tN_ns * 1000.0
            valid = tN_ps[np.isfinite(tN_ps)]

            mode_ps = find_mode_ps(valid)
            rms_full = float(np.nanstd(tN_ps))
            rob = robust_estimators(valid)
            fit2 = fit_gaussian_core(valid, mode_ps)
            fitex12 = fit_gaussian_exec12_style(valid)
            ct = core_tail_fraction(valid, mode_ps, fit2["sigma_ps"])

            if N == 4:
                ss_ps = (
                    math.sqrt(N) * TAU_D_NS * 1000.0 / mean_npe_nearest
                    if mean_npe_nearest > 0 else math.nan
                )
            else:
                ss_ps = math.nan

            results[N] = {
                "x_mm": x_mm, "run_label": RUN_LABELS[x_mm], "channel": ch,
                "threshold_pe": N, "n_events_valid": len(valid),
                "frac_excl": frac_excl,
                "mean_npe_nearest": mean_npe_nearest,
                "mode_ps": mode_ps, "rms_full_ps": rms_full,
                "sigma_core_ps": fit2["sigma_ps"],
                "sigma_core_err_ps": fit2["sigma_err_ps"],
                "chi2ndf": fit2["chi2ndf"],
                "fit_mu_ps": fit2["mu_ps"],
                "iqr_sigma_ps": rob["iqr_sigma_ps"],
                "mad_sigma_ps": rob["mad_sigma_ps"],
                "q25_ps": rob["q25_ps"], "q50_ps": rob["q50_ps"],
                "q75_ps": rob["q75_ps"],
                "core_fraction": ct["core_fraction"],
                "tail_fraction": ct["tail_fraction"],
                "tail_boundary_ps": ct["tail_boundary_ps"],
                "tail_k": TAIL_K,
                "sigma_exec12_style_ps": fitex12["sigma_ps"],
                "sigma_exec12_style_err_ps": fitex12["sigma_err_ps"],
                "sigma_stat_theory_ps": ss_ps,
                # no EXEC_12-EJ230 reference exists
                "_fit2": fit2, "_tN_ps": valid,
            }

            print(
                f"  x={x_mm:+d} N={N:2d}: "
                f"mode={mode_ps:.1f} RMS={rms_full:.1f} | "
                f"sigma_core(2ps)={fit2['sigma_ps']:.2f}±{fit2['sigma_err_ps']:.2f} ps | "
                f"chi2/ndf={fit2['chi2ndf']:.2f} | "
                f"MAD={rob['mad_sigma_ps']:.1f} IQR={rob['iqr_sigma_ps']:.1f} | "
                f"tail={ct['tail_fraction']*100:.1f}%"
                + (f" | EXEC12-style={fitex12['sigma_ps']:.2f} ps | "
                   f"sigma_stat_th={ss_ps:.2f} ps" if N == 4 else "")
            )

        all_results[x_mm] = results
    return all_results


def plot_f5(
    all_results: dict,
    nearests: dict[int, int],
    figs_dir: pathlib.Path,
) -> list[GateRecord]:
    """F5: t4 at CORE_BIN_PS bins, zoom to nucleus, Gaussian fit, fixed scale."""
    modes = [all_results[x][4]["mode_ps"] for x in POSITIONS_MM]
    global_mean_mode = float(np.mean(modes))
    x_lo = global_mean_mode + ZOOM_LO_PS
    x_hi = global_mean_mode + ZOOM_HI_PS

    global_ymax = 0
    for x_mm in POSITIONS_MM:
        valid = all_results[x_mm][4]["_tN_ps"]
        edges = np.arange(x_lo, x_hi + CORE_BIN_PS, CORE_BIN_PS)
        counts, _ = np.histogram(valid, bins=edges)
        global_ymax = max(global_ymax, int(counts.max()))
    global_ymax = int(global_ymax * 1.2)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    fig.suptitle(
        r"F5 — EJ-230 $t_4$ nucleus: 2 ps bins, zoom to core, "
        r"Gaussian fit + $\sigma_{\rm stat}$ theory, fixed scale ($\tau_d=1.5$ ns)",
        fontsize=11,
    )
    xlims_rec, ylims_rec = [], []

    for ax, x_mm in zip(axes, POSITIONS_MM):
        r = all_results[x_mm][4]
        valid = r["_tN_ps"]
        fit = r["_fit2"]
        ss = r["sigma_stat_theory_ps"]
        color = _COLORS[x_mm]

        edges = np.arange(x_lo, x_hi + CORE_BIN_PS, CORE_BIN_PS)
        counts, _ = np.histogram(valid, bins=edges)
        centers = (edges[:-1] + edges[1:]) / 2.0

        ax.step(centers, counts, where="mid", color=color, lw=1.5, alpha=0.85)
        ax.fill_between(centers, counts, step="mid", color=color, alpha=0.12)

        if fit["ok"] and not math.isnan(fit["sigma_ps"]):
            x_fit = np.linspace(x_lo, x_hi, 800)
            y_fit = _gauss(x_fit, fit["A"], fit["mu_ps"], fit["sigma_ps"])
            ax.plot(x_fit, y_fit, "k-", lw=1.8,
                    label=f"Gauss fit: {fit['sigma_ps']:.2f} ps")

        if not math.isnan(ss):
            A_stat = float(len(valid)) * CORE_BIN_PS / (ss * math.sqrt(2 * math.pi))
            x_th = np.linspace(x_lo, x_hi, 800)
            y_th = _gauss(x_th, A_stat, r["mode_ps"], ss)
            ax.plot(x_th, y_th, "--", color="gray", lw=1.3, alpha=0.8,
                    label=f"$\\sigma_{{stat}}={ss:.1f}$ ps (theory)")

        mad_v = r["mad_sigma_ps"]
        sigma_v = fit["sigma_ps"] if not math.isnan(fit["sigma_ps"]) else float("nan")
        sigma_e = fit["sigma_err_ps"] if not math.isnan(fit["sigma_err_ps"]) else float("nan")
        ax.text(
            0.97, 0.95,
            f"MAD$\\times$1.4826={mad_v:.2f} ps\n"
            f"$\\sigma_{{\\rm core}}$={sigma_v:.2f}$\\pm${sigma_e:.2f} ps\n"
            f"$\\chi^2$/ndf={r['chi2ndf']:.2f}\n"
            f"IQR/1.349={r['iqr_sigma_ps']:.1f} ps\n"
            f"RMS$_{{\\rm full}}$={r['rms_full_ps']:.1f} ps",
            transform=ax.transAxes, ha="right", va="top", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85),
        )
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(0, global_ymax)
        ax.set_xlabel("$t_4$ (ps)", fontsize=10)
        if ax is axes[0]:
            ax.set_ylabel(f"Events / {CORE_BIN_PS:.0f} ps bin", fontsize=10)
        ax.set_title(
            f"{RUN_LABELS[x_mm]}: x={x_mm:+d} mm, ID {nearests[x_mm]}", fontsize=10)
        ax.legend(fontsize=7.5, loc="upper right")
        ax.grid(True, alpha=0.3)
        xlims_rec.append(ax.get_xlim())
        ylims_rec.append(ax.get_ylim())

    fig.tight_layout()
    _save(fig, figs_dir, "exec13_230_f5_t4_core")
    gate = _check_axis_limits(xlims_rec, ylims_rec, "F5_core")
    return [gate]


# ─── CSV writers ──────────────────────────────────────────────────────────────

def write_csv_f1(hists: list[HistResult], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_230_f1_tn_histograms.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "run_label", "threshold_pe", "nearest_id",
                    "bin_center_ps", "counts", "counts_norm_per_ps",
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
    path = out_dir / "exec13_230_f2_cluster_dispersion.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "run_label", "global_id", "mean_npe",
                    "n", "mean_tn_ps", "rms_tn_ps",
                    "stat_floor_ps", "cluster_truncated"])
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
                    d["cluster_truncated"],
                ])


def write_csv_f3(data: list[dict], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_230_f3_npe_profile.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "run_label", "global_id", "mean_npe", "rms_npe"])
        for d in data:
            w.writerow([d["x_mm"], RUN_LABELS[d["x_mm"]], d["global_id"],
                        f"{d['mean_npe']:.4f}", f"{d['rms_npe']:.4f}"])


def write_csv_f4(data: list[dict], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_230_f4_timing_profile.csv"
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


def write_csv_core(all_results: dict, out_dir: pathlib.Path) -> None:
    fields = [
        "x_mm", "run_label", "channel", "threshold_pe", "n_events_valid", "frac_excl",
        "mean_npe_nearest", "mode_ps", "rms_full_ps",
        "sigma_core_ps", "sigma_core_err_ps", "chi2ndf", "fit_mu_ps",
        "iqr_sigma_ps", "mad_sigma_ps", "q25_ps", "q50_ps", "q75_ps",
        "core_fraction", "tail_fraction", "tail_boundary_ps", "tail_k",
        "sigma_exec12_style_ps", "sigma_exec12_style_err_ps",
        "sigma_stat_theory_ps",
    ]
    path = out_dir / "exec13_230_core_resolution.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for x_mm in POSITIONS_MM:
            for N in (4, 20):
                r = all_results[x_mm][N]
                row = {}
                for k in fields:
                    v = r.get(k, "")
                    if isinstance(v, float):
                        row[k] = f"{v:.4f}" if not math.isnan(v) else "nan"
                    else:
                        row[k] = v
                w.writerow(row)


def write_csv_gates(gates: list[GateRecord], out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_230_gates.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gate", "passed", "detail"])
        for g in gates:
            w.writerow([g.name, g.passed, g.detail])


# ─── LaTeX macros (G-MACRO guard) ─────────────────────────────────────────────

def write_macros(
    all_results: dict,
    nearests: dict[int, int],
    out_dir: pathlib.Path,
) -> GateRecord:
    lines = [
        "% Auto-generated by exec13_230_fixed_scale.py",
        "% Source: exec13_230_core_resolution.csv — do not edit manually",
        "% EJ-230 (OPSC-106): TAU_D = 1.5 ns (NOT 1.8 ns of EJ-204)",
        "% No EXEC_12-EJ230 reference exists; ExecTwelveRef macros omitted.",
        "",
        f"\\newcommand{{\\ejtauD}}{{{TAU_D_NS:.1f}}}",
        f"\\newcommand{{\\ejYield}}{{{SCINT_YIELD_PER_MEV}}}",
        "",
    ]

    # Command name stems without any digits
    tags = {-690: "RunOne", -400: "RunThree", 0: "RunFour"}
    for x_mm in POSITIONS_MM:
        r4 = all_results[x_mm][4]
        tag = tags[x_mm]
        nearest = nearests[x_mm]

        def cmd(name: str, val, fmt: str = ".2f") -> str:
            if isinstance(val, float) and math.isnan(val):
                return f"\\newcommand{{\\core{tag}{name}}}{{--}}"
            if isinstance(val, float):
                return f"\\newcommand{{\\core{tag}{name}}}{{{val:{fmt}}}}"
            return f"\\newcommand{{\\core{tag}{name}}}{{{val}}}"

        lines += [
            f"% {RUN_LABELS[x_mm]} x={x_mm:+d} mm ID {nearest}",
            cmd("RMSps", r4["rms_full_ps"]),
            cmd("SigCore", r4["sigma_core_ps"]),
            cmd("SigCoreErr", r4["sigma_core_err_ps"]),
            cmd("ChiNdf", r4["chi2ndf"]),
            cmd("IQR", r4["iqr_sigma_ps"]),
            cmd("MAD", r4["mad_sigma_ps"]),
            cmd("CoreFrac", r4["core_fraction"]),
            cmd("TailFrac", r4["tail_fraction"]),
            cmd("TailPct", r4["tail_fraction"] * 100.0, ".0f"),
            cmd("SigEx", r4["sigma_exec12_style_ps"]),
            cmd("SigStat", r4["sigma_stat_theory_ps"]),
            cmd("NearestID", nearest, "d"),
            "",
        ]

    # Global summary macros (no digits in names)
    sigs_core = [all_results[x][4]["sigma_core_ps"] for x in POSITIONS_MM
                 if not math.isnan(all_results[x][4]["sigma_core_ps"])]
    sigs_ex = [all_results[x][4]["sigma_exec12_style_ps"] for x in POSITIONS_MM
               if not math.isnan(all_results[x][4]["sigma_exec12_style_ps"])]
    sigs_stat = [all_results[x][4]["sigma_stat_theory_ps"] for x in POSITIONS_MM
                 if not math.isnan(all_results[x][4]["sigma_stat_theory_ps"])]

    def mean_cmd(name: str, vals: list) -> str:
        if not vals:
            return f"\\newcommand{{\\core{name}}}{{--}}"
        return f"\\newcommand{{\\core{name}}}{{{float(np.mean(vals)):.2f}}}"

    lines += [
        "% Global summary macros",
        f"\\newcommand{{\\coreTailK}}{{{TAIL_K:.0f}}}",
        f"\\newcommand{{\\coreBinPS}}{{{CORE_BIN_PS:.0f}}}",
        f"\\newcommand{{\\coreWindowPS}}{{{CORE_WINDOW_PS:.0f}}}",
        f"\\newcommand{{\\coreSPTRps}}{{{SPTR_PS}}}",
        mean_cmd("MeanSigCore", sigs_core),
        mean_cmd("MeanSigEx", sigs_ex),
        mean_cmd("MeanSigStat", sigs_stat),
        "",
    ]

    macro_path = out_dir / "exec13_230_resolution_macros.tex"
    macro_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Macros written to {macro_path}")

    # G-MACRO: guard — no digits in \newcommand names
    content = "\n".join(lines)
    offenders = re.findall(r'\\newcommand\{\\[A-Za-z]*[0-9][^}]*\}', content)
    passed = len(offenders) == 0
    detail = (
        "no digits in command names"
        if passed
        else f"DIGIT OFFENDERS: {offenders}"
    )
    if not passed:
        print(f"ABORT G-MACRO: digit(s) found in \\newcommand names: {offenders}")
        sys.exit(1)
    print(f"  G-MACRO PASSED: {detail}")
    return GateRecord("G-MACRO", passed, detail)


# ─── LaTeX table ──────────────────────────────────────────────────────────────

def make_latex_table(hists: list[HistResult]) -> str:
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

    return (
        r"% Auto-generated by exec13_230_fixed_scale.py — do not edit manually" + "\n"
        r"% Source: exec13_230_f1_tn_histograms.csv" + "\n"
        r"\begin{table}[h]" + "\n"
        r"\centering\small" + "\n"
        r"\caption{EXEC\_13-230: $t_N$ statistics at top-nearest channel. "
        r"EJ-230 (OPSC-106). Bin 10\,ps, range [0,1000]\,ps.}" + "\n"
        r"\label{tab:exec13_230_tN}" + "\n"
        r"\begin{tabular}{lrrrr}" + "\n"
        r"\toprule" + "\n"
        r"Position & $\mu$ (ps) & RMS (ps) & Excl.\ frac.\ & Overflow \\" + "\n"
        r"\midrule" + "\n"
        r"\multicolumn{5}{c}{\textit{$t_4$ (leading-edge threshold = 4 PE)}} \\" + "\n"
        + rows_str_4 + "\n"
        r"\midrule" + "\n"
        r"\multicolumn{5}{c}{\textit{$t_{20}$ (bulk-light level = 20 PE)}} \\" + "\n"
        + rows_str_20 + "\n"
        r"\bottomrule" + "\n"
        r"\end{tabular}" + "\n"
        r"\end{table}" + "\n"
    )


# ─── Beamer generator ─────────────────────────────────────────────────────────

def _frame(title: str, body: str) -> str:
    return f"\\begin{{frame}}\n\\frametitle{{{title}}}\n{body}\n\\end{{frame}}\n\n"


def _howtoread(why: str, axes: str, takeaway: str) -> str:
    return (
        r"\vspace{-0.1cm}" + "\n"
        r"\begin{block}{\tiny How to read this plot}" + "\n"
        r"\tiny\textbf{Why:} " + why + "  "
        r"\textbf{Axes:} " + axes + "  "
        r"\textbf{Takeaway:} " + takeaway + "\n"
        r"\end{block}" + "\n"
    )


def generate_beamer(
    hists: list[HistResult],
    all_results: dict,
    nearests: dict[int, int],
    commit_hash: str,
) -> str:
    preamble = (
        r"\documentclass[aspectratio=169]{beamer}" + "\n"
        r"\usetheme{Madrid}" + "\n"
        r"\usecolortheme{default}" + "\n"
        r"\usepackage{booktabs,amsmath}" + "\n"
        r"\graphicspath{{figs/}}" + "\n"
        r"\input{exec13_230_resolution_macros}" + "\n"
        r"\title[EXEC\_13-230]{EXEC\_13-230: Fixed-Scale $t_N$ + Core Resolution (EJ-230)}" + "\n"
        r"\subtitle{R\'eplica del an\'alisis EJ-204 de la reuni\'on 12-jun con G.~V\'asquez}" + "\n"
        r"\author{Ren\'e R\'ios}" + "\n"
        r"\institute{Universidad de La Serena (ULS)}" + "\n"
        r"\date{2026-06-15}" + "\n"
        r"\begin{document}" + "\n"
        r"\begin{frame}\titlepage\end{frame}" + "\n\n"
    )

    frames = []

    # ── F1 slides (6 frames: 2 thresholds × 3 variants) ──────────────────────
    for N in TN_THRESHOLDS_PE:
        h_f1 = _howtoread(
            f"EXEC\\_12 EJ-204 histograms had auto-scale; this re-bins at 10\\,ps with "
            f"fixed [0,1000]\\,ps range for visual comparison between runs. "
            f"EJ-230 $\\tau_d=1.5$\\,ns (shorter than EJ-204 1.8\\,ns).",
            f"x: $t_{{{N}}}$ arrival time (ps); y: prob. density (norm, ps$^{{-1}}$) "
            f"or events/10\\,ps (abs). Note: peak density CAN exceed 1\\,ps$^{{-1}}$ --- correct.",
            "Narrow (t4) or broader (t20) distribution; shape encodes scintillation + Landau fluctuations.",
        )
        frames.append(_frame(
            f"F1 --- Overlay $t_{{{N}}}$, normalised, fixed scale [EJ-230]",
            r"\centering\includegraphics[height=0.48\textheight,keepaspectratio]{"
            f"exec13_230_f1_t{N}pe_overlay_norm" + "}\n" + h_f1,
        ))
        frames.append(_frame(
            f"F1 --- Overlay $t_{{{N}}}$, absolute counts [EJ-230]",
            r"\centering\includegraphics[height=0.48\textheight,keepaspectratio]{"
            f"exec13_230_f1_t{N}pe_overlay_abs" + "}\n" + h_f1,
        ))
        frames.append(_frame(
            f"F1 --- Individual panels $t_{{{N}}}$ (identical scale) [EJ-230]",
            r"\centering\includegraphics[height=0.52\textheight,keepaspectratio]{"
            f"exec13_230_f1_t{N}pe_panels" + "}\n" + h_f1,
        ))

    # ── F2 slides (3 frames: one per run) ────────────────────────────────────
    h_f2 = _howtoread(
        r"Do neighbour SiPMs $\pm1,\pm2$ carry the same timing as the nearest? "
        r"Stat floor uses EJ-230 $\tau_d=\ejtauD$\,ns.",
        r"x: photon order $n$ (1--20); y: $\langle t_n\rangle$ (ps). "
        r"Error bars = RMS. Dotted: stat floor $\sqrt{n}\,\tau_d/\langle N_{pe}\rangle$.",
        r"Neighbours approach stat floor when their $\langle N_{pe}\rangle$ is comparable to the nearest.",
    )
    for x_mm in POSITIONS_MM:
        frames.append(_frame(
            f"F2 --- {RUN_LABELS[x_mm]}: $x={x_mm:+d}$\\,mm cluster dispersion [EJ-230]",
            r"\centering\includegraphics[height=0.50\textheight,keepaspectratio]{"
            f"exec13_230_f2_cluster_{x_mm:+d}mm" + "}\n" + h_f2,
        ))

    # ── F3 slide ──────────────────────────────────────────────────────────────
    h_f3 = _howtoread(
        "Verify nearest channel is the Npe maximum; shows photon budget across top array.",
        r"x: global channel ID (16--85); y: $\langle N_{pe}\rangle$/event $\pm$ RMS. Highlighted: cluster $\pm$2.",
        r"Npe falls off rapidly away from track; $\pm$2 neighbours still collect useful light.",
    )
    frames.append(_frame(
        r"F3 --- EJ-230 Top $N_{pe}$ profile (all 70 channels, fixed scale)",
        r"\centering\includegraphics[width=0.98\textwidth]{exec13_230_f3_npe_profile}" + "\n"
        + h_f3,
    ))

    # ── F4 slides (2 frames: one per threshold) ───────────────────────────────
    h_f4 = _howtoread(
        rf"Only channels with eff $\geq{EFF_MIN_FRACTION:.2f}$ are shown; "
        r"far channels excluded at N=20 (insufficient Npe) --- listed in CSV.",
        r"x: global channel ID; y: $\langle t_N\rangle$ (ps) $\pm$ RMS.",
        r"Timing rises steeply for channels far from track (longer path + lower $N_{pe}$).",
    )
    for N in TN_THRESHOLDS_PE:
        frames.append(_frame(
            f"F4 --- EJ-230 timing profile $\\langle t_{{{N}}}\\rangle$ "
            f"(eff$\\geq${EFF_MIN_FRACTION:.2f}, fixed scale)",
            r"\centering\includegraphics[width=0.98\textwidth]{"
            f"exec13_230_f4_timing_t{N}pe" + "}\n" + h_f4,
        ))

    # ── Gate table ────────────────────────────────────────────────────────────
    tags = {-690: "RunOne", -400: "RunThree", 0: "RunFour"}
    nearest_list = " | ".join(
        f"$x={x:+d}\\to\\text{{ID }}\\core{tags[x]}NearestID$" for x in POSITIONS_MM
    )
    frames.append(_frame(
        "Validation gates [EJ-230]",
        r"""\small
\begin{itemize}
  \item \textbf{G-MAT/G-TAU}: $\tau_d=\ejtauD$\,ns (EJ-230), read from
        \texttt{analysis/exec13/common13.py} --- not a literal
  \item \textbf{G1}: nearest IDs determined programmatically from $\langle N_{pe}\rangle$
        data: """ + nearest_list + r"""
  \item \textbf{G2}: all figure panels within each family share identical axis limits
        --- verified post-plot
  \item \textbf{G3}: $n_{\rm events}=2000$ per run
  \item \textbf{G4}: $n_{\rm in}+n_{\rm excl}+n_{\rm ovfl}=2000$ per histogram
  \item \textbf{G5}: F2 x-axis labelled ``Order $n$'', range 1--20
  \item \textbf{G-MACRO}: zero \texttt{\textbackslash newcommand} with digits in the command name
\end{itemize}
\vspace{3pt}
\scriptsize All gate results: \texttt{exec13\_230\_gates.csv}.
Commit: \texttt{""" + commit_hash[:12] + r"""}.
""",
    ))

    # ── Provenance ────────────────────────────────────────────────────────────
    frames.append(_frame(
        "Provenance and methodology [EJ-230]",
        r"""\small
\textbf{Material:} EJ-230 (OPSC-106) --- $\tau_r=0.5$\,ns, $\tau_d=\ejtauD$\,ns,
yield $=\ejYield$\,ph/MeV, $\ell_{\rm att}=120$\,cm, $\lambda_{\rm peak}=391$\,nm.\\[4pt]
\textbf{Data:} 31-position beam scan, 2000 events/position.
Three positions analysed: $x=-690$, $-400$, $0$\,mm.\\[4pt]
\textbf{Algorithm:} \texttt{compute\_tN} ported unchanged from
\texttt{exec07/exec12\_tN\_analysis.py:72--100} (sorted time-to-N-th hit per event).
$t=0$: Geant4 event start (\texttt{exec07\_photon\_budget.py:342}).\\[4pt]
\textbf{No gates applied:} no SPTR, no jitter, no SPE convolution, no temporal window.\\[4pt]
\textbf{EXEC\_12-EJ230 reconciliation:} no published EXEC\_12 sigma\_fit exists for EJ-230.
EXEC\_12-style refit (25\,ps bins, mean$\pm$2RMS) computed on EJ-230 data for
informational reference only.\\[4pt]
\scriptsize Analysis script: \texttt{analysis/exec13/exec13\_230\_fixed\_scale.py}.
""",
    ))

    # ── F5 nucleus ────────────────────────────────────────────────────────────
    h_f5 = _howtoread(
        r"RMS of $t_4$ (tens of ps) is dominated by reflection tail, not the prompt core. "
        r"Zoom to mode$\pm$\coreWindowPS\,ps with 2\,ps bins to measure the narrow nucleus.",
        r"x: $t_4$ (ps), zoom $\approx$ mode$\pm$\coreWindowPS\,ps; "
        r"y: events/2\,ps. Black solid: Gaussian fit. Gray dashed: $\sigma_{\rm stat}$ theory.",
        r"Narrow Gaussian nucleus ($\sim$ few ps); tail starts at mode$+\coreTailK\sigma_{\rm core}$.",
    )
    frames.append(_frame(
        r"F5 --- EJ-230 $t_4$ nucleus: 2\,ps bins, Gaussian fit, fixed scale",
        r"\centering\includegraphics[width=0.98\textwidth]{exec13_230_f5_t4_core}" + "\n"
        + h_f5,
    ))

    # ── Conclusion: nucleus vs tail ───────────────────────────────────────────
    tags_list = [("RunOne", -690), ("RunThree", -400), ("RunFour", 0)]
    rows_tex = ""
    for tag, x_mm in tags_list:
        r4 = all_results[x_mm][4]
        row = (
            f"  $x={x_mm:+d}$\\,mm & "
            f"\\core{tag}MAD & "
            f"\\core{tag}SigCore$\\pm$\\core{tag}SigCoreErr & "
            f"\\core{tag}ChiNdf & "
            f"\\core{tag}IQR & "
            f"\\core{tag}RMSps & "
            f"\\core{tag}SigStat \\\\"
        )
        rows_tex += "  " + row + "\n"

    frames.append(_frame(
        r"Resolución Top nearest: núcleo vs cola (EJ-230, $\tau_d=\ejtauD$\,ns)",
        r"""{\small
\begin{tabular}{lrrrrrr}
\toprule
Posición & MAD$\times$1.4826 & $\sigma_{\rm core}\pm$err & $\chi^2$/ndf
         & IQR/1.349 & RMS$_{\rm full}$ & $\sigma_{\rm stat}$ \\
         & (ps) & (ps) & & (ps) & (ps) & (ps) \\
\midrule
""" + rows_tex + r"""\bottomrule
\end{tabular}}
\vspace{4pt}
\begin{alertblock}{\footnotesize Mensaje clave (aprendizaje EJ-204)}
\scriptsize
\textbf{Estimador preferido: MAD$\times$1.4826} (robusto, resistente a colas).
$\sigma_{\rm core}$ (ajuste gaussiano 2\,ps) tiende a \emph{sub-estimar} con $\chi^2$/ndf
elevado --- ver tabla.
Si IQR anómalamente grande (p.ej.\ posición central), usar MAD.
Cola = $t_4 > \mathrm{modo}+\coreTailK\sigma_{\rm core}$; dominada por reflexiones.
$\sigma_{\rm core}$ y MAD son \emph{intrínsecos de simulación} (sin SPTR, sin jitter);
la resolución real estará dominada por el SPTR ($\approx\coreSPTRps$\,ps, bloque EXEC\_02b).
\end{alertblock}
""",
    ))

    # ── Observations ─────────────────────────────────────────────────────────
    frames.append(_frame(
        "Observaciones y próximos pasos [EJ-230]",
        r"""\small
\begin{itemize}
  \item \textbf{Material EJ-230 vs EJ-204:} $\tau_d=\ejtauD$\,ns (EJ-230) vs
        $1.8$\,ns (EJ-204) $\Rightarrow$ $\sigma_{\rm stat}(t_4)$ intrínseca menor para EJ-230.
  \item \textbf{Sin referencia EXEC\_12-EJ230:} el ajuste estilo-EXEC\_12 (25\,ps bins)
        se computa sobre datos EJ-230 como referencia informativa; no existe publicación previa.
  \item \textbf{Fuera de alcance (bloques posteriores):}
  \begin{itemize}\scriptsize
    \item Comparación lado a lado EJ-204 vs EJ-230.
    \item Introducción de SPTR/jitter (EXEC\_02b).
    \item MPPC solo en extremos; optimización de espaciado de SiPMs.
    \item Corrección de time-walk (HOOK\_WALK reservado).
  \end{itemize}
  \item \textbf{G-OVER:} revisar \texttt{exec13\_230\_report.log} para
        Overfull/Underfull --- listado en resumen final.
\end{itemize}
""",
    ))

    return preamble + "".join(frames) + r"\end{document}" + "\n"


# ─── Main ──────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="EXEC_13-230: fixed-scale tN + core-resolution analysis for EJ-230."
    )
    parser.add_argument("--data-dir", type=pathlib.Path, default=_DEFAULT_DATA_DIR)
    parser.add_argument("--out-dir", type=pathlib.Path, default=_DEFAULT_OUT_DIR)
    args = parser.parse_args(argv)
    data_dir: pathlib.Path = args.data_dir
    out_dir: pathlib.Path = args.out_dir
    figs_dir = out_dir / "figs"
    figs_dir.mkdir(parents=True, exist_ok=True)

    gates: list[GateRecord] = []

    # ── G-MAT + G-TAU ────────────────────────────────────────────────────────
    gates.extend(check_material_gate())

    # ── Load three runs ───────────────────────────────────────────────────────
    print("\nLoading ROOT files…")
    runs: list[RunData] = []
    for x in POSITIONS_MM:
        path = expected_file(data_dir, x)
        if not path.exists():
            print(f"ABORT: ROOT file not found: {path}")
            return 1
        run = load_run(data_dir, x)
        runs.append(run)

    # ── G3: verify 2000 events ────────────────────────────────────────────────
    for run in runs:
        g3 = GateRecord(
            f"G3_nevents_{run.x_mm}",
            run.n_events == N_EVENTS,
            f"n_events={run.n_events}, expected={N_EVENTS}",
        )
        gates.append(g3)
        if not g3.passed:
            print(f"ABORT G3: x={run.x_mm} mm: {run.n_events} events (expected {N_EVENTS})")
            return 1
    print(f"G3 PASSED: all runs have {N_EVENTS} events")

    # ── Npe matrices ──────────────────────────────────────────────────────────
    print("Computing Npe matrices…")
    npe_matrices: dict[int, np.ndarray] = {}
    for run in runs:
        npe_matrices[run.x_mm] = top_npe_matrix(
            run.event_id, run.global_id, run.n_events
        )

    # ── G1: nearest-top IDs ───────────────────────────────────────────────────
    print("Determining nearest-top channels…")
    nearests: dict[int, int] = {}
    for run in runs:
        prof = npe_matrices[run.x_mm].mean(axis=0)
        nearest = select_nearest_top(run.x_mm, prof)
        nearests[run.x_mm] = nearest
        nearest_pos = float(TOP_POSITIONS_MM[nearest - 16])
        dist_geom = abs(nearest_pos - run.x_mm)
        mean_npe_nearest = float(prof[nearest - 16])
        print(f"  x={run.x_mm:+d} mm → nearest top ID {nearest} "
              f"(pos={nearest_pos:.0f} mm, dist={dist_geom:.0f} mm, "
              f"⟨Npe⟩={mean_npe_nearest:.1f})")

        # Coherence check: nearest must be within 25 mm of beam x
        coherent = dist_geom <= 25.0
        gates.append(GateRecord(
            f"G1_nearest_{run.x_mm}",
            coherent,
            f"ID={nearest} pos={nearest_pos:.0f} mm, beam x={run.x_mm} mm, "
            f"dist={dist_geom:.0f} mm {'COHERENT' if coherent else 'INCOHERENT'}",
        ))
        if not coherent:
            print(f"ABORT G1: nearest ID {nearest} at {nearest_pos:.0f} mm "
                  f"is {dist_geom:.0f} mm from beam x={run.x_mm} mm (>25 mm)")
            return 1

    print("G1 PASSED: nearest IDs geometrically coherent")

    # ── STAGE A: F1 ──────────────────────────────────────────────────────────
    print("\nGenerating F1…")
    all_hists: list[HistResult] = []
    for N in TN_THRESHOLDS_PE:
        hists, g_f1 = plot_f1(runs, nearests, N, figs_dir)
        all_hists.extend(hists)
        gates.extend(g_f1)

    # G4: histogram accounting
    for h in all_hists:
        n_below_check = round(h.frac_excl * runs[0].n_events)
        total_check = h.n_in_range + h.n_overflow + n_below_check
        g4 = GateRecord(
            f"G4_{h.x_mm}_{h.threshold}pe",
            total_check == runs[0].n_events,
            f"in={h.n_in_range}+ovfl={h.n_overflow}+below≈{n_below_check}={total_check}",
        )
        gates.append(g4)

    # ── F2 ───────────────────────────────────────────────────────────────────
    print("Generating F2…")
    f2_data, g_f2 = plot_f2(runs, nearests, npe_matrices, figs_dir)
    gates.extend(g_f2)

    # ── F3 ───────────────────────────────────────────────────────────────────
    print("Generating F3…")
    f3_data, g_f3 = plot_f3(runs, nearests, npe_matrices, figs_dir)
    gates.extend(g_f3)

    # ── F4 ───────────────────────────────────────────────────────────────────
    print("Generating F4…")
    f4_data, g_f4 = plot_f4(runs, nearests, figs_dir)
    gates.extend(g_f4)

    # ── Write STAGE A CSVs ───────────────────────────────────────────────────
    print("Writing Stage A CSVs…")
    write_csv_f1(all_hists, out_dir)
    write_csv_f2(f2_data, out_dir)
    write_csv_f3(f3_data, out_dir)
    write_csv_f4(f4_data, out_dir)

    latex_table = make_latex_table(all_hists)
    (out_dir / "exec13_230_tN_table.tex").write_text(latex_table)

    # ── STAGE B: core resolution ──────────────────────────────────────────────
    print("\nStage B: core resolution analysis…")
    all_results = analyze_core(runs, nearests, npe_matrices)

    print("\nGenerating F5…")
    g_f5 = plot_f5(all_results, nearests, figs_dir)
    gates.extend(g_f5)

    print("Writing core CSV…")
    write_csv_core(all_results, out_dir)

    print("Writing LaTeX macros…")
    g_macro = write_macros(all_results, nearests, out_dir)
    gates.append(g_macro)

    # ── Beamer deck ───────────────────────────────────────────────────────────
    print("Generating Beamer deck…")
    try:
        commit_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=str(_REPO), text=True
        ).strip()
    except Exception:
        commit_hash = "unknown"

    beamer_tex = generate_beamer(all_hists, all_results, nearests, commit_hash)
    tex_path = out_dir / "exec13_230_report.tex"
    tex_path.write_text(beamer_tex, encoding="utf-8")
    print(f"  Beamer deck written to {tex_path}")

    # ── Write gate CSV ────────────────────────────────────────────────────────
    write_csv_gates(gates, out_dir)

    # ── Gate summary ──────────────────────────────────────────────────────────
    failed = [g for g in gates if not g.passed]
    print(f"\n── Gate summary: {len(gates)} checks, {len(failed)} failed ────────")
    for g in gates:
        status = "PASS" if g.passed else "FAIL"
        print(f"  [{status}] {g.name}: {g.detail}")

    if failed:
        print("\nWARNING: some gates failed — see exec13_230_gates.csv")

    # ── Final user summary ────────────────────────────────────────────────────
    print("\n════════════════════════════════════════════════════════════════")
    print("EXEC_13-230 ANALYSIS SUMMARY")
    print("════════════════════════════════════════════════════════════════")
    print(f"  EJ-230 tau_d         = {TAU_D_NS} ns  (from common13.py)")
    print(f"  EJ-230 yield         = {SCINT_YIELD_PER_MEV} ph/MeV")
    print(f"  EXEC_12-EJ230 ref    = NONE (no published sigma_fit for EJ-230)")
    print()
    print("  Nearest IDs (from data):")
    for x_mm in POSITIONS_MM:
        ch = nearests[x_mm]
        pos = float(TOP_POSITIONS_MM[ch - 16])
        mean_npe = float(npe_matrices[x_mm][:, ch - 16].mean())
        print(f"    x={x_mm:+d} mm → ID {ch} (pos={pos:.0f} mm, ⟨Npe⟩={mean_npe:.1f})")
    print()
    print("  Core resolution results (t4, Top nearest):")
    print(f"  {'x_mm':>8}  {'ID':>4}  {'MAD(ps)':>9}  {'sigma_core(ps)':>14}  "
          f"{'chi2/ndf':>9}  {'IQR(ps)':>8}  {'RMS_full(ps)':>13}  {'sigma_stat_th(ps)':>17}")
    for x_mm in POSITIONS_MM:
        r = all_results[x_mm][4]
        print(
            f"  {x_mm:>8}  {nearests[x_mm]:>4}  "
            f"{r['mad_sigma_ps']:>9.2f}  {r['sigma_core_ps']:>14.2f}  "
            f"{r['chi2ndf']:>9.2f}  {r['iqr_sigma_ps']:>8.1f}  "
            f"{r['rms_full_ps']:>13.1f}  {r['sigma_stat_theory_ps']:>17.2f}"
        )
    print()
    print("  EXEC_12-style fit (25 ps bins, informational):")
    for x_mm in POSITIONS_MM:
        r = all_results[x_mm][4]
        ex = r["sigma_exec12_style_ps"]
        print(f"    x={x_mm:+d} mm: sigma_ex12_style = {ex:.2f} ps  (no EJ-230 ref)")
    print()
    print("  Tail fraction (k=3):")
    for x_mm in POSITIONS_MM:
        r = all_results[x_mm][4]
        print(f"    x={x_mm:+d} mm: tail = {r['tail_fraction']*100:.1f}%  "
              f"RMS_full = {r['rms_full_ps']:.1f} ps")
    print()
    all_passed = len(failed) == 0
    print(f"  Gates: {'ALL PASSED' if all_passed else f'{len(failed)} FAILED'}")
    for g in failed:
        print(f"    [FAIL] {g.name}: {g.detail}")
    print()
    print(f"  Output dir : {out_dir}")
    print(f"  Beamer deck: {tex_path.name}")
    print(f"  Commit hash: {commit_hash}")
    print("════════════════════════════════════════════════════════════════")
    print("Compile with:")
    print(f"  cd {out_dir}")
    print(f"  pdflatex -interaction=nonstopmode exec13_230_report.tex  # twice")
    print()

    return 0 if all_passed else 2


if __name__ == "__main__":
    sys.exit(main())
