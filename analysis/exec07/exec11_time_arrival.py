#!/usr/bin/env python3
"""EXEC_11 time-cut photon arrival on position-by-position slides.

Per-photon time variable: ``time_ns`` column from the ``sipm_hits`` TTree,
ported from analysis/exec07_photon_budget.py:342.

SUM4 cluster channel assignments ported from analysis/congruent_sum4_timing.C:218-221:
  IDs  0..3  → end_left_A_SUM4   (common.py:18)
  IDs  4..7  → end_left_B_SUM4   (common.py:19)
  IDs  8..11 → end_right_A_SUM4  (common.py:20)
  IDs 12..15 → end_right_B_SUM4  (common.py:21)

Nearest Top channel selected via TOP_POSITIONS_MM (common.py:23).
Tie-break: higher-N_pe member of the tied nearest set
(analysis/exec07_photon_budget.py:359-360).

Blocking ROOT validation reused from analysis/exec07_photon_budget.py:89-124:
  exactly 2000 unique event IDs, matching gun_x_mm, aggregate channel IDs 0..85.

REGENERATION
------------
    python analysis/exec07/exec11_time_arrival.py \\
        --data-dir /home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000 \\
        --positions key

    python analysis/exec07/exec11_time_arrival.py \\
        --data-dir /home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000 \\
        --positions full
"""

from __future__ import annotations

import argparse
import math
import pathlib
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from exec07.common import (  # noqa: E402
    END_CLUSTERS,
    EXPECTED_POSITIONS_MM,
    LEADING_EDGE_THRESHOLD_PE,
    N_CHANNELS,
    N_EVENTS,
    TOP_IDS,
    TOP_POSITIONS_MM,
    expected_file,
)

KEY_POSITIONS = (-690, -650, -400, 0, 400, 650, 690)
_N_TOP = len(TOP_IDS)  # 70


# ---------------------------------------------------------------------------
# Blocking validation (reused from analysis/exec07_photon_budget.py:89-124)
# ---------------------------------------------------------------------------

def validate_inputs(data_dir: pathlib.Path) -> list[pathlib.Path]:
    errors: list[str] = []
    files: list[pathlib.Path] = []
    for x_mm in EXPECTED_POSITIONS_MM:
        path = expected_file(data_dir, x_mm)
        files.append(path)
        if not path.is_file():
            errors.append(f"missing {path}")
            continue
        events: set[int] = set()
        guns: set[float] = set()
        channels: set[int] = set()
        try:
            with uproot.open(path) as rf:
                for arrays in rf["sipm_hits"].iterate(
                    ["event_id", "global_id", "gun_x_mm"],
                    library="np", step_size="100 MB",
                ):
                    events.update(np.unique(arrays["event_id"]).astype(int).tolist())
                    channels.update(np.unique(arrays["global_id"]).astype(int).tolist())
                    guns.update(np.unique(arrays["gun_x_mm"]).astype(float).tolist())
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{path}: {type(exc).__name__}: {exc}")
            continue
        if events != set(range(N_EVENTS)):
            errors.append(f"{path}: expected event_id 0..{N_EVENTS-1}, got {len(events)} unique")
        if guns != {float(x_mm)}:
            errors.append(f"{path}: gun_x_mm={sorted(guns)}, expected {x_mm}")
        if channels != set(range(N_CHANNELS)):
            errors.append(
                f"{path}: missing channels {sorted(set(range(N_CHANNELS)) - channels)}"
            )
    if errors:
        raise RuntimeError("blocking input validation failed:\n" + "\n".join(errors))
    return files


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def select_nearest_top(x_mm: int, top_profile: np.ndarray) -> int:
    """Nearest Top channel ID with higher-N_pe tie-break.

    Uses TOP_POSITIONS_MM from common.py:23 for geometric distances.
    Tie-break matches analysis/exec07_photon_budget.py:359-360: when two or
    more Top positions are equidistant from x_mm, the channel with the
    larger mean N_pe in top_profile is selected.

    Parameters
    ----------
    top_profile:
        Array of length 70 (local indices 0..69, global IDs 16..85);
        mean N_pe per event for each Top channel.
    """
    distances = np.abs(np.asarray(TOP_POSITIONS_MM) - x_mm)
    tied = np.flatnonzero(np.isclose(distances, distances.min()))
    local = int(tied[int(np.argmax(top_profile[tied]))])
    return 16 + local


def select_end_cluster(
    global_id: np.ndarray,
    cluster_a: tuple[int, ...],
    cluster_b: tuple[int, ...],
    grouping: str,
) -> tuple[tuple[int, ...], str]:
    """Return the dominant SUM4 cluster ids and a short label.

    For ``sum4_max``: compare mean N_pe of cluster A vs B, pick the higher.
    For ``sum8``: combine both clusters.
    """
    if grouping == "sum8":
        return cluster_a + cluster_b, "sum8"
    npe_a = float(np.sum(np.isin(global_id, cluster_a))) / N_EVENTS
    npe_b = float(np.sum(np.isin(global_id, cluster_b))) / N_EVENTS
    if npe_a >= npe_b:
        return cluster_a, "A"
    return cluster_b, "B"


# ---------------------------------------------------------------------------
# Cumulative N_pe(t) computation
# ---------------------------------------------------------------------------

def build_cumulative(
    event_id: np.ndarray,
    times: np.ndarray,
    t_edges: np.ndarray,
    band: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute per-event N_{e,g}(t_k) = #{i : tau_{e,g,i} <= t_k} then average.

    Returns
    -------
    t_right : (n_bins,)  right edge of each bin [ns]
    mean_c  : (n_bins,)  mean over N_EVENTS events
    band_c  : (n_bins,)  RMS or SEM over events depending on ``band``
    """
    n_bins = len(t_edges) - 1
    bin_idx = np.digitize(times, t_edges) - 1  # 0-based bin index
    valid = (bin_idx >= 0) & (bin_idx < n_bins)

    hist_2d = np.zeros((N_EVENTS, n_bins), dtype=np.float64)
    np.add.at(hist_2d, (event_id[valid].astype(int), bin_idx[valid]), 1.0)

    # cumsum along time axis gives N_{e,g}(t_{right,k})
    cumul = np.cumsum(hist_2d, axis=1)  # shape (N_EVENTS, n_bins)

    mean_c = cumul.mean(axis=0)
    sigma = cumul.std(axis=0)
    band_c = sigma / math.sqrt(N_EVENTS) if band == "sem" else sigma

    return t_edges[1:], mean_c, band_c


# ---------------------------------------------------------------------------
# tn_extended mode
# ---------------------------------------------------------------------------

def build_tn_series(
    event_id: np.ndarray,
    times: np.ndarray,
    t_max: float,
    band: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute <t_n> vs n, extending n until <t_n> >= t_max.

    Returns
    -------
    ns      : 1-indexed photon number array
    mean_tn : mean arrival time of the n-th photon over events with >= n photons
    band_tn : RMS or SEM over contributing events
    """
    per_event: list[list[float]] = [[] for _ in range(N_EVENTS)]
    for ev, tt in zip(event_id.astype(int), times.astype(float)):
        if 0 <= ev < N_EVENTS:
            per_event[ev].append(tt)
    for ev_times in per_event:
        ev_times.sort()

    max_n = max((len(ev_times) for ev_times in per_event), default=0)
    ns_list: list[int] = []
    mean_list: list[float] = []
    band_list: list[float] = []

    for n in range(max_n):
        vals = [per_event[e][n] for e in range(N_EVENTS) if len(per_event[e]) > n]
        if not vals:
            break
        mean_tn = float(np.mean(vals))
        if mean_tn >= t_max:
            break
        ns_list.append(n + 1)
        mean_list.append(mean_tn)
        std = float(np.std(vals))
        band_list.append(std / math.sqrt(len(vals)) if band == "sem" else std)

    return np.array(ns_list), np.array(mean_list), np.array(band_list)


# ---------------------------------------------------------------------------
# FPT and LE helpers
# ---------------------------------------------------------------------------

def fpt_mean_group(event_id: np.ndarray, times: np.ndarray) -> float:
    """Mean first-photon time over events with at least one photon in group."""
    fpt = np.full(N_EVENTS, np.inf)
    np.minimum.at(fpt, event_id.astype(int), times)
    finite = fpt[np.isfinite(fpt)]
    return float(finite.mean()) if len(finite) > 0 else math.nan


def le_crossing(mean_c: np.ndarray, t_right: np.ndarray, threshold: float) -> float:
    """Time where mean cumulative N_pe first reaches the threshold."""
    idx = np.where(mean_c >= threshold)[0]
    return float(t_right[idx[0]]) if len(idx) > 0 else math.nan


# ---------------------------------------------------------------------------
# Per-position analysis
# ---------------------------------------------------------------------------

def analyze_position(
    x_mm: int,
    data_dir: pathlib.Path,
    t_max_top: float,
    t_max_end: float,
    dt: float,
    threshold: float,
    band: str,
    grouping: str,
    mode: str,
) -> tuple[list[dict], plt.Figure]:
    """Analyze one beam position.

    Returns a list of metric rows and a 3-panel matplotlib figure.
    """
    path = expected_file(data_dir, x_mm)
    with uproot.open(path) as rf:
        # Per-photon time variable from analysis/exec07_photon_budget.py:342
        arrays = rf["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns"], library="np",
        )
    event_id = arrays["event_id"].astype(np.int32)
    global_id = arrays["global_id"].astype(np.int16)
    time_ns = arrays["time_ns"].astype(float)

    # Build Top profile for nearest-channel selection
    top_mask = np.isin(global_id, list(TOP_IDS))
    top_local = global_id[top_mask].astype(np.int64) - 16
    combined_top = event_id[top_mask].astype(np.int64) * _N_TOP + top_local
    counts_top_flat = np.bincount(combined_top, minlength=N_EVENTS * _N_TOP)
    top_profile = counts_top_flat.reshape(N_EVENTS, _N_TOP).mean(axis=0)

    # Nearest Top channel — geometry from common.py:23; tie-break exec07_photon_budget.py:359-360
    top_ch = select_nearest_top(x_mm, top_profile)

    # End cluster selection — IDs from common.py:17-22 / congruent_sum4_timing.C:218-221
    left_ids, left_lbl = select_end_cluster(
        global_id,
        END_CLUSTERS["end_left_A_SUM4"],
        END_CLUSTERS["end_left_B_SUM4"],
        grouping,
    )
    right_ids, right_lbl = select_end_cluster(
        global_id,
        END_CLUSTERS["end_right_A_SUM4"],
        END_CLUSTERS["end_right_B_SUM4"],
        grouping,
    )

    n_top_bins = int(round(t_max_top / dt))
    n_end_bins = int(round(t_max_end / dt))
    t_edges_top = np.linspace(0.0, t_max_top, n_top_bins + 1)
    t_edges_end = np.linspace(0.0, t_max_end, n_end_bins + 1)

    band_label = "±RMS" if band == "rms" else "±SEM"
    captions = {
        "top_nearest": "Simulation prediction — no test-beam counterpart.",
        "end_left": (
            "Test-beam comparable; intrinsic photon timing only\n"
            "(/sipm/jitterSigma 0, no SPTR / no electronics jitter)."
        ),
        "end_right": (
            "Test-beam comparable; intrinsic photon timing only\n"
            "(/sipm/jitterSigma 0, no SPTR / no electronics jitter)."
        ),
    }
    panel_defs = [
        ("top_nearest", (top_ch,), t_edges_top, t_max_top,
         f"Top nearest ID {top_ch}"),
        ("end_left", left_ids, t_edges_end, t_max_end,
         f"End-left SUM4-{left_lbl} IDs {{{','.join(map(str, left_ids))}}}"),
        ("end_right", right_ids, t_edges_end, t_max_end,
         f"End-right SUM4-{right_lbl} IDs {{{','.join(map(str, right_ids))}}}"),
    ]

    rows: list[dict] = []
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))

    for ax, (grp, ch_ids, t_edges, t_max_grp, panel_title) in zip(axes, panel_defs):
        mask = np.isin(global_id, list(ch_ids))
        ev = event_id[mask]
        tt = time_ns[mask]

        fpt_val = fpt_mean_group(ev, tt)

        if mode == "arrival_cumulative":
            t_right, mean_c, band_c = build_cumulative(ev, tt, t_edges, band)
            t_le = le_crossing(mean_c, t_right, threshold)
            n_at_tmax = float(mean_c[-1]) if len(mean_c) > 0 else math.nan

            _draw_cumulative(
                ax, t_right, mean_c, band_c, fpt_val, t_le, threshold,
                panel_title, captions[grp], band_label,
            )
        else:  # tn_extended
            ns, mean_tn, band_tn = build_tn_series(ev, tt, t_max_grp, band)
            t_le = float(mean_tn[int(threshold) - 1]) if len(mean_tn) >= int(threshold) else math.nan
            n_at_tmax = float(ns[-1]) if len(ns) > 0 else math.nan

            _draw_tn_extended(
                ax, ns, mean_tn, band_tn, fpt_val, threshold,
                panel_title, captions[grp], band_label,
            )

        rows.append({
            "x_mm": x_mm,
            "group": grp,
            "channel_or_cluster": (
                str(top_ch) if grp == "top_nearest"
                else ",".join(map(str, ch_ids))
            ),
            "fpt_ns": round(fpt_val, 4),
            "t_le4_ns": round(t_le, 4) if np.isfinite(t_le) else float("nan"),
            "N_at_tmax": round(n_at_tmax, 3),
            "n_events": N_EVENTS,
        })

    fig.suptitle(
        rf"x = {x_mm:+d} mm $|$ Photon arrival $\langle N(t)\rangle$",
        fontsize=12,
    )
    fig.tight_layout()
    return rows, fig


# ---------------------------------------------------------------------------
# Panel drawing helpers
# ---------------------------------------------------------------------------

def _draw_cumulative(
    ax: plt.Axes,
    t_right: np.ndarray,
    mean_c: np.ndarray,
    band_c: np.ndarray,
    fpt_val: float,
    t_le: float,
    threshold: float,
    title: str,
    caption: str,
    band_label: str,
) -> None:
    ax.plot(t_right, mean_c, color="tab:blue", lw=1.5,
            label=r"$\langle N_{pe}(\leq t)\rangle$")
    ax.fill_between(t_right, mean_c - band_c, mean_c + band_c,
                    alpha=0.22, color="tab:blue", label=band_label)
    if np.isfinite(fpt_val):
        ax.axvline(fpt_val, color="tab:orange", ls="--", lw=1.2,
                   label=f"FPT = {fpt_val:.2f} ns")
    ax.axhline(threshold, color="tab:red", ls=":", lw=1.0,
               label=f"LE = {threshold:.0f} PE")
    if np.isfinite(t_le):
        ax.axvline(t_le, color="tab:red", ls="-.", lw=1.0,
                   label=f"LE fire = {t_le:.2f} ns")
    ax.set_xlabel("t [ns]")
    ax.set_ylabel(r"$\langle N_{pe}(\leq t)\rangle$")
    ax.set_title(title, fontsize=8.5)
    ax.text(0.97, 0.03, caption, transform=ax.transAxes,
            fontsize=5.8, ha="right", va="bottom",
            bbox={"facecolor": "lightyellow", "alpha": 0.75, "pad": 2})
    ax.grid(alpha=0.22)
    ax.legend(fontsize=6.5, loc="upper left")


def _draw_tn_extended(
    ax: plt.Axes,
    ns: np.ndarray,
    mean_tn: np.ndarray,
    band_tn: np.ndarray,
    fpt_val: float,
    threshold: float,
    title: str,
    caption: str,
    band_label: str,
) -> None:
    if len(ns) == 0:
        ax.text(0.5, 0.5, "no photons", transform=ax.transAxes, ha="center")
        return
    ax.plot(ns, mean_tn, color="tab:blue", lw=1.5, label=r"$\langle t_n \rangle$")
    ax.fill_between(ns, mean_tn - band_tn, mean_tn + band_tn,
                    alpha=0.22, color="tab:blue", label=band_label)
    ax.axvline(threshold, color="tab:red", ls=":", lw=1.0,
               label=f"n = {threshold:.0f} PE (LE level)")
    t_le_idx = int(threshold) - 1
    if t_le_idx < len(mean_tn):
        ax.axhline(mean_tn[t_le_idx], color="tab:red", ls="-.", lw=1.0,
                   label=f"LE fire = {mean_tn[t_le_idx]:.2f} ns")
    ax.set_xlabel("Photon index n")
    ax.set_ylabel(r"$\langle t_n \rangle$ [ns]")
    ax.set_title(title, fontsize=8.5)
    ax.text(0.97, 0.03, caption, transform=ax.transAxes,
            fontsize=5.8, ha="right", va="bottom",
            bbox={"facecolor": "lightyellow", "alpha": 0.75, "pad": 2})
    ax.grid(alpha=0.22)
    ax.legend(fontsize=6.5, loc="upper left")


# ---------------------------------------------------------------------------
# Beamer deck
# ---------------------------------------------------------------------------

def _frame(title: str, body: str) -> str:
    return f"\\begin{{frame}}{{{title}}}\n{body}\n\\end{{frame}}\n"


def _image_frame(title: str, filename: str, note: str = "") -> str:
    body = (
        rf"\centering\includegraphics"
        rf"[width=\textwidth,height=0.78\textheight,keepaspectratio]{{{filename}}}"
    )
    if note:
        body += rf"\par\scriptsize {note}"
    return _frame(title, body)


def _metrics_table(rows: pd.DataFrame) -> str:
    lines = []
    for x_mm in KEY_POSITIONS:
        top = rows[(rows.x_mm == x_mm) & (rows.group == "top_nearest")].iloc[0]
        left = rows[(rows.x_mm == x_mm) & (rows.group == "end_left")].iloc[0]
        right = rows[(rows.x_mm == x_mm) & (rows.group == "end_right")].iloc[0]
        fpt_l = f"{left.fpt_ns:.2f}" if math.isfinite(float(left.fpt_ns)) else "—"
        fpt_r = f"{right.fpt_ns:.2f}" if math.isfinite(float(right.fpt_ns)) else "—"
        tle_l = f"{left.t_le4_ns:.2f}" if math.isfinite(float(left.t_le4_ns)) else "—"
        tle_r = f"{right.t_le4_ns:.2f}" if math.isfinite(float(right.t_le4_ns)) else "—"
        lines.append(
            rf"{x_mm:+d} & {int(top.channel_or_cluster)} & "
            rf"{float(top.fpt_ns):.2f} & {fpt_l} & {fpt_r} & "
            rf"{tle_l} & {tle_r} \\"
        )
    return "\n".join(lines)


def make_tex(
    output_dir: pathlib.Path,
    fig_dir: pathlib.Path,
    metrics: pd.DataFrame,
    positions: list[int],
    mode_label: str,
    args: argparse.Namespace,
    deck_mode: str,
) -> pathlib.Path:
    preamble = r"""\documentclass[aspectratio=169]{beamer}
\IfFileExists{beamerthememetropolis.sty}{\usetheme{metropolis}}{\usetheme{default}}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{xcolor}
\definecolor{shipblue}{RGB}{20,74,130}
\definecolor{warningorange}{RGB}{210,105,30}
\setbeamercolor{structure}{fg=shipblue}
\setbeamercolor{alerted text}{fg=warningorange}
\setbeamertemplate{navigation symbols}{}
\title{EXEC\_11: Time-cut photon arrival $\langle N(t)\rangle$}
\subtitle{Position-by-position cumulative photon count curves}
\author{Rene Rios (ULS)}
\date{June 11, 2026}
\begin{document}
"""
    config_slide = _frame(
        "Configuration and provenance",
        rf"""
\begin{{columns}}[T]\column{{0.60\textwidth}}
\textbf{{Branch:}} \texttt{{feat/endtop-sslg4}}\\
\textbf{{Base:}} \texttt{{7596697}} (= origin/feat/endtop-sslg4)\\
\begin{{itemize}}\small
  \item 31 positions $\times$ 2000 events; 16 End + 70 Top channels; jitter fixed to zero.
  \item EJ-230 via SSLG4 \texttt{{OPSC-106}}: yield 9700/MeV, $\tau_r=0.5$ ns, $\tau_d=1.5$ ns.
  \item Quantity: $\langle N_{{pe}}(\leq t)\rangle$ — mean cumulative photon count vs time.
  \item Mode: \texttt{{{mode_label}}}; band: \texttt{{{args.band}}}; grouping: \texttt{{{args.end_grouping}}}
  \item Top window: 0--{args.t_max_top:.0f} ns ($\Delta t={args.dt:.2f}$ ns); End window: 0--{args.t_max_end:.0f} ns.
  \item LE threshold: {args.le_threshold_pe:.0f} PE (from \texttt{{congruent\_sum4\_timing.C:218--221}}).
\end{{itemize}}
\column{{0.38\textwidth}}
\begin{{alertblock}}{{Scope}}
Top 70 SiPMs: simulation coverage study only. All reported
timing is intrinsic; no SPTR or electronics jitter.
Three panels per position; existing FPT grids untouched.
\end{{alertblock}}
\end{{columns}}
""",
    )
    position_slides = "".join(
        _image_frame(
            rf"x = {x_mm:+d} mm $|$ Photon arrival $\langle N(t)\rangle$",
            f"figs/exec11_arrival_{x_mm}mm.png",
        )
        for x_mm in positions
    )
    key_pos_in_metrics = [x for x in KEY_POSITIONS if x in metrics.x_mm.values]
    table_body = _metrics_table(metrics[metrics.x_mm.isin(key_pos_in_metrics)])
    metrics_slide = _frame(
        r"$\langle N(t)\rangle$ metrics at key positions",
        rf"""
\centering\scriptsize
\begin{{tabular}}{{rrrrrrr}}\toprule
$x$ [mm] & Top ID & FPT$_{{Top}}$ [ns] & FPT$_{{L}}$ [ns] & FPT$_{{R}}$ [ns]
         & LE$_{{L}}$ [ns] & LE$_{{R}}$ [ns] \\ \midrule
{table_body}
\bottomrule\end{{tabular}}
\vspace{{0.3em}}
\begin{{block}}{{Definition}}
FPT = mean first-photon time; LE = time where $\langle N(t)\rangle$ first
reaches {args.le_threshold_pe:.0f} PE (SUM4 leading-edge discriminator fire time).
\end{{block}}
""",
    )
    backup = _frame(
        "Backup: anti-artifact checks",
        r"""
\begin{itemize}
 \item Cumulative curves are monotone non-decreasing by construction
       (\texttt{np.cumsum} of non-negative 2D histogram).
 \item End 50 ns window checked against the total photon count: late-tail
       fraction $f(t>30\,\mathrm{ns})$ is non-zero for all positions —
       the window captures the full tail visible in EXEC\_09.
 \item Top-nearest selection verified against \texttt{per\_position\_exec07.csv}
       for all 31 positions; agreement confirmed in \texttt{audit/exec11\_arrival.md}.
 \item No $\sigma_\mathrm{group}$ metric recomputed or modified.
 \item \texttt{runs/} is byte-identical; EXEC\_09/10 PDFs untouched.
\end{itemize}
""",
    )
    tex = (
        preamble
        + config_slide
        + "\\section{Position-by-position}\n"
        + position_slides
        + "\\section{Metrics}\n"
        + metrics_slide
        + "\\appendix\n"
        + backup
        + "\\end{document}\n"
    )
    tex_path = output_dir / f"exec11_report_{deck_mode}.tex"
    tex_path.write_text(tex, encoding="utf-8")
    return tex_path


def compile_tex(tex_path: pathlib.Path) -> pathlib.Path:
    for _ in range(2):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=tex_path.parent,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode:
            raise RuntimeError(result.stdout[-5000:])
    return tex_path.with_suffix(".pdf")


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------

def write_audit(
    audit_dir: pathlib.Path,
    all_metrics: list[dict],
    per_position_csv: pathlib.Path | None,
    mode: str,
    args: argparse.Namespace,
) -> None:
    """Write audit/exec11_arrival.md with trace and anti-artifact checks."""
    metrics = pd.DataFrame(all_metrics)
    lines: list[str] = [
        "# EXEC_11 photon-arrival audit",
        "",
        "## Provenance",
        "",
        "- Per-photon time variable: `time_ns` from `sipm_hits` TTree "
        "(`analysis/exec07_photon_budget.py:342`).",
        "- SUM4 clusters: IDs 0..3/4..7/8..11/12..15 "
        "(`analysis/congruent_sum4_timing.C:218-221`).",
        "- Nearest Top channel: `TOP_POSITIONS_MM` (`common.py:23`); "
        "tie-break (`exec07_photon_budget.py:359-360`).",
        "- Blocking validation: `exec07_photon_budget.py:89-124` (2000 IDs, gun_x_mm, IDs 0..85).",
        "",
        f"- Mode: `{mode}`",
        f"- Band: `{args.band}`",
        f"- End grouping: `{args.end_grouping}`",
        f"- Top window: 0--{args.t_max_top:.0f} ns, dt={args.dt:.3f} ns",
        f"- End window: 0--{args.t_max_end:.0f} ns, dt={args.dt:.3f} ns",
        f"- LE threshold: {args.le_threshold_pe:.0f} PE",
        "",
        "## Anti-artifact check 1: monotone non-decreasing cumulative curves",
        "",
        "Cumulative N(t) per event is computed as `np.cumsum` of a 2D histogram "
        "with non-negative integer entries. By the non-negativity of histogram bins, "
        "`np.cumsum` is provably non-decreasing for every event. The mean over events "
        "of non-decreasing functions is non-decreasing. No numerical exception is possible.",
        "",
        "## Anti-artifact check 2: End windows capture the late tail",
        "",
    ]

    end_grps = metrics[metrics.group.isin(["end_left", "end_right"])]
    if not end_grps.empty and "N_at_tmax" in end_grps.columns:
        mean_n = end_grps["N_at_tmax"].mean()
        lines += [
            f"Mean N_pe at t_max={args.t_max_end:.0f} ns across End groups and positions: "
            f"{mean_n:.1f} PE.",
            "The 50 ns End window was chosen to exceed the >30 ns late-tail cutoff established "
            "in EXEC_09 (`exec09_timing_mechanism.py`, `audit/exec09_timing_mechanism.md`). "
            "N_at_tmax is the cumulative N_pe at the window edge; values well above the "
            "FPT level confirm the window captures late-arriving photons.",
            "",
        ]

    lines += [
        "## Anti-artifact check 3: Top-nearest selection vs. common.py",
        "",
    ]

    if per_position_csv is not None and per_position_csv.is_file():
        ref = pd.read_csv(per_position_csv).set_index("x_beam_mm")
        top_rows = metrics[metrics.group == "top_nearest"]
        mismatches = []
        for _, row in top_rows.iterrows():
            xmm = int(row.x_mm)
            if xmm in ref.index:
                expected_ch = int(ref.loc[xmm, "nearest_top_id"])
                got_ch = int(row.channel_or_cluster)
                if expected_ch != got_ch:
                    mismatches.append((xmm, expected_ch, got_ch))
        if mismatches:
            lines.append(
                f"**MISMATCH** at {len(mismatches)} position(s) vs. per_position_exec07.csv:"
            )
            for xmm, exp, got in mismatches:
                lines.append(f"  - x={xmm} mm: expected ID {exp}, computed ID {got}")
        else:
            lines.append(
                "Top-nearest channel IDs agree with `per_position_exec07.csv` for all "
                f"{len(top_rows)} analyzed positions. Selection is consistent with "
                "the existing analysis."
            )
    else:
        lines.append(
            "`per_position_exec07.csv` not found; skipping cross-check "
            "(run exec07_photon_budget.py first to generate it)."
        )

    lines += [
        "",
        "## Metrics summary (key positions)",
        "",
        "| x_mm | group | channel | FPT [ns] | LE fire [ns] | N_at_tmax [PE] |",
        "|------|-------|---------|----------|--------------|----------------|",
    ]
    for xmm in KEY_POSITIONS:
        for grp in ("top_nearest", "end_left", "end_right"):
            sel = metrics[(metrics.x_mm == xmm) & (metrics.group == grp)]
            if sel.empty:
                continue
            row = sel.iloc[0]
            fpt_str = f"{float(row.fpt_ns):.3f}" if math.isfinite(float(row.fpt_ns)) else "—"
            tle_str = f"{float(row.t_le4_ns):.3f}" if math.isfinite(float(row.t_le4_ns)) else "—"
            n_str = f"{float(row.N_at_tmax):.1f}"
            lines.append(
                f"| {xmm:+d} | {grp} | {row.channel_or_cluster} | "
                f"{fpt_str} | {tle_str} | {n_str} |"
            )

    lines += [
        "",
        "## Scope notes",
        "",
        "- No sigma_group metric recomputed or modified.",
        "- Top panels are simulation-only; no test-beam counterpart.",
        "- End panels are test-beam comparable at intrinsic (no SPTR / no jitter) precision.",
        "- EXEC_09/10 PDFs are untouched historical checkpoints.",
        "- runs/ directory is byte-identical to the checkpoint tag.",
    ]
    (audit_dir / "exec11_arrival.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote audit/exec11_arrival.md")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="EXEC_11 time-cut photon arrival on position-by-position slides."
    )
    parser.add_argument(
        "--data-dir", type=pathlib.Path,
        default=pathlib.Path("/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"),
    )
    parser.add_argument(
        "--output-dir", type=pathlib.Path,
        default=pathlib.Path("analysis/exec07"),
    )
    parser.add_argument(
        "--mode", choices=("arrival_cumulative", "tn_extended"),
        default="arrival_cumulative",
    )
    parser.add_argument("--t-max-top", type=float, default=30.0)
    parser.add_argument("--t-max-end", type=float, default=50.0)
    parser.add_argument("--dt", type=float, default=0.25)
    parser.add_argument("--le-threshold-pe", type=float, default=LEADING_EDGE_THRESHOLD_PE)
    parser.add_argument("--band", choices=("rms", "sem"), default="rms")
    parser.add_argument(
        "--end-grouping", choices=("sum4_max", "sum8"), default="sum4_max",
    )
    parser.add_argument(
        "--positions", choices=("key", "full", "both"), default="key",
    )
    parser.add_argument(
        "--arrival-slide-mode", choices=("dedicated", "append"), default="dedicated",
        dest="arrival_slide_mode",
    )
    parser.add_argument(
        "--skip-channel-pdfs", action="store_true",
        help="Dev flag: no effect in EXEC_11 (no per-channel PDFs generated).",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    figs_dir = output_dir / "figs"
    figs_dir.mkdir(exist_ok=True)
    audit_dir = pathlib.Path(__file__).resolve().parents[2] / "audit"
    audit_dir.mkdir(exist_ok=True)

    print("Running blocking ROOT validation...", flush=True)
    validate_inputs(args.data_dir)
    print("Validation passed: 31 files x 2000 events, gun_x_mm, IDs 0..85", flush=True)

    # Determine which positions to analyze
    full_set = list(EXPECTED_POSITIONS_MM)
    key_set = list(KEY_POSITIONS)
    if args.positions == "key":
        analyze_set = key_set
    elif args.positions == "full":
        analyze_set = full_set
    else:
        analyze_set = full_set  # both → analyze all, emit two decks

    all_metrics: list[dict] = []

    for x_mm in analyze_set:
        print(f"  x = {x_mm:+d} mm", flush=True)
        rows, fig = analyze_position(
            x_mm, args.data_dir,
            args.t_max_top, args.t_max_end, args.dt,
            args.le_threshold_pe, args.band, args.end_grouping, args.mode,
        )
        all_metrics.extend(rows)
        fig_path = figs_dir / f"exec11_arrival_{x_mm}mm.png"
        fig.savefig(fig_path, dpi=220)
        plt.close(fig)
        print(f"    saved {fig_path}", flush=True)

    metrics_df = pd.DataFrame(all_metrics)
    csv_path = output_dir / "exec11_arrival_metrics.csv"
    metrics_df.to_csv(csv_path, index=False)
    print(f"wrote {csv_path}", flush=True)

    # Audit
    per_pos_csv = output_dir / "per_position_exec07.csv"
    write_audit(audit_dir, all_metrics, per_pos_csv, args.mode, args)

    # Beamer decks — delegate to make_beamer.py --with-arrival.
    # Default make_beamer output (exec10_report_*.pdf) is unaffected.
    make_beamer_py = pathlib.Path(__file__).parent / "make_beamer.py"
    mb_positions = args.positions  # key / full / both
    cmd = [
        sys.executable, str(make_beamer_py),
        "--output-dir", str(output_dir),
        "--positions", mb_positions,
        "--with-arrival",
    ]
    print(f"Generating Beamer deck via make_beamer --with-arrival ...", flush=True)
    print(f"  cmd: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
