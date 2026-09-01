#!/usr/bin/env python3
"""EXEC_12b T4: per-channel <t_n> vs n order-statistics dispersion figures.

For each of the 7 key beam positions, produce 3 multi-panel PNG canvases:
  figs/tn_{x}mm_endL.png  -- End L channels 0-7  (2x4 grid, blue)
  figs/tn_{x}mm_endR.png  -- End R channels 8-15 (2x4 grid, red)
  figs/tn_{x}mm_top.png   -- 9 nearest Top channels (3x3 grid, green)

Each subplot shows:
  y-axis: <t_n> [ns] (mean arrival time of the n-th detected photon per event)
  x-axis: n = 1..20 (n-ésimo fotón)
  Error bars: RMS over events (event-to-event dispersion)
  Dashed band: statistical floor sigma_stat(t_n) = sqrt(n)*tau_d/<Npe>

Requirement: >=10 events per (channel, n) point, else the point is omitted.

Data source:
  TTree "sipm_hits", branches event_id/global_id/time_ns
  (analysis/exec07_photon_budget.py:337-342)
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import uproot

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from exec07.common import (  # noqa: E402
    N_EVENTS,
    TAU_D_NS,
    TOP_POSITIONS_MM,
    expected_file,
)

KEY_POSITIONS = (-690, -650, -400, 0, 400, 650, 690)
N_MAX = 20
MIN_ENTRIES = 10
_DPI = 150

# ── Geometry helpers ──────────────────────────────────────────────────────────

def _top_9_neighbors(x_mm: float) -> list[int]:
    """Return global IDs of the 9 nearest Top channels (±4 by position index)."""
    top_pos = list(TOP_POSITIONS_MM)
    nearest_idx = min(range(70), key=lambda i: (abs(top_pos[i] - x_mm), i))
    half = 4
    start, end = nearest_idx - half, nearest_idx + half
    if start < 0:
        end = min(69, end - start)
        start = 0
    elif end > 69:
        start = max(0, start - (end - 69))
        end = 69
    return [16 + i for i in range(start, end + 1)]

# ── Core computation ─────────────────────────────────────────────────────────

def _load_arrays(data_dir: pathlib.Path, x_mm: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read event_id, global_id, time_ns from the sipm_hits TTree."""
    path = expected_file(data_dir, x_mm)
    with uproot.open(path) as f:
        arr = f["sipm_hits"].arrays(["event_id", "global_id", "time_ns"], library="np")
    return (
        arr["event_id"].astype(np.int32),
        arr["global_id"].astype(np.int32),
        arr["time_ns"].astype(np.float64),
    )


def _per_channel_stats(
    event_id: np.ndarray,
    global_id: np.ndarray,
    time_ns: np.ndarray,
    ch_id: int,
) -> tuple[list[int], list[float], list[float], float]:
    """Compute <t_n> and RMS(t_n) for n=1..N_MAX for a single channel.

    Returns (ns, means_ns, rmss_ns, mean_npe).
    Points with <MIN_ENTRIES events are omitted.
    """
    mask = global_id == ch_id
    ev = event_id[mask]
    tt = time_ns[mask]

    mean_npe = float(len(ev)) / N_EVENTS

    if len(ev) == 0:
        return [], [], [], 0.0

    order = np.lexsort((tt, ev))
    ev_s = ev[order].astype(np.int64)
    tt_s = tt[order]

    all_ev = np.arange(N_EVENTS, dtype=np.int64)
    starts = np.searchsorted(ev_s, all_ev, side="left")
    stops = np.searchsorted(ev_s, all_ev, side="right")

    ns_out: list[int] = []
    means: list[float] = []
    rmss: list[float] = []

    for n in range(1, N_MAX + 1):
        have = stops - starts >= n
        tn_vals = tt_s[starts[have] + n - 1]
        if len(tn_vals) >= MIN_ENTRIES:
            ns_out.append(n)
            means.append(float(np.mean(tn_vals)))
            rmss.append(float(np.std(tn_vals, ddof=1)))

    return ns_out, means, rmss, mean_npe


def _stat_floor(ns: list[int], mean_npe: float) -> list[float]:
    """σ_stat(t_n) ≈ sqrt(n)*tau_d/<Npe> (pure order-statistics floor)."""
    if mean_npe <= 0:
        return [0.0] * len(ns)
    return [float(np.sqrt(n)) * TAU_D_NS / mean_npe for n in ns]

# ── Plotting helpers ──────────────────────────────────────────────────────────

def _draw_pad(
    ax: plt.Axes,
    ns: list[int],
    means: list[float],
    rmss: list[float],
    mean_npe: float,
    color: str,
    title: str,
) -> None:
    """Populate one subplot with <t_n> ± RMS and the statistical floor band."""
    if not ns:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="gray")
        ax.set_title(title, fontsize=6, pad=2)
        ax.set_visible(True)
        return

    ns_arr = np.array(ns)
    means_arr = np.array(means)
    rmss_arr = np.array(rmss)
    floor = np.array(_stat_floor(ns, mean_npe))

    # Main data: mean ± RMS
    ax.errorbar(
        ns_arr, means_arr, yerr=rmss_arr,
        fmt="o-", color=color, markersize=2.5, linewidth=1.0,
        elinewidth=0.8, capsize=2, label=r"$\langle t_n\rangle\pm\mathrm{RMS}$",
    )

    # Statistical floor band around mean
    ax.fill_between(
        ns_arr,
        means_arr - floor,
        means_arr + floor,
        alpha=0.20, color="gray",
        label=r"piso estadístico ($\sqrt{n}\tau_d/N_{pe}$)",
    )
    ax.plot(ns_arr, means_arr - floor, color="gray", linewidth=0.7,
            linestyle="--", alpha=0.7)
    ax.plot(ns_arr, means_arr + floor, color="gray", linewidth=0.7,
            linestyle="--", alpha=0.7)

    ax.set_xlim(0.5, N_MAX + 0.5)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(1))
    ax.tick_params(axis="both", labelsize=6)
    ax.set_title(title, fontsize=6.5, pad=2)

    npe_str = f"$\\langle N_{{pe}}\\rangle={mean_npe:.1f}$"
    ax.annotate(npe_str, xy=(0.97, 0.03), xycoords="axes fraction",
                ha="right", va="bottom", fontsize=5.5)


def _make_end_figure(
    event_id: np.ndarray,
    global_id: np.ndarray,
    time_ns: np.ndarray,
    ids: list[int],
    x_mm: int,
    color: str,
    face_label: str,
    output_path: pathlib.Path,
) -> None:
    """2×4 grid for 8 End channels."""
    fig, axes = plt.subplots(2, 4, figsize=(14, 5.5), dpi=_DPI)
    fig.subplots_adjust(hspace=0.55, wspace=0.38,
                        left=0.06, right=0.98, top=0.86, bottom=0.12)

    fig.suptitle(
        f"$\\langle t_n\\rangle$ per channel — {face_label}  |  "
        f"$\\mu$ @ $x = {x_mm:+d}$ mm\n"
        r"Error bars $=$ RMS (event-to-event).  Dashed band: piso estadístico $\sqrt{n}\,\tau_d/\langle N_{pe}\rangle$",
        fontsize=8, y=0.97,
    )

    for ax, ch in zip(axes.flat, ids):
        ns, means, rmss, npe = _per_channel_stats(event_id, global_id, time_ns, ch)
        _draw_pad(ax, ns, means, rmss, npe, color, f"ID {ch}")

    # Shared axis labels
    for ax in axes[1, :]:
        ax.set_xlabel("$n$-ésimo fotón", fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\langle t_n\rangle$ [ns]", fontsize=7)

    fig.savefig(output_path, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {output_path}")


def _make_top_figure(
    event_id: np.ndarray,
    global_id: np.ndarray,
    time_ns: np.ndarray,
    ids: list[int],
    x_mm: int,
    output_path: pathlib.Path,
) -> None:
    """3×3 grid for 9 nearest Top channels."""
    n_ch = len(ids)
    n_rows, n_cols = 3, 3
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 8), dpi=_DPI)
    fig.subplots_adjust(hspace=0.55, wspace=0.38,
                        left=0.07, right=0.98, top=0.88, bottom=0.08)

    fig.suptitle(
        f"$\\langle t_n\\rangle$ per Top channel (nearest $\\pm4$)  |  "
        f"$\\mu$ @ $x = {x_mm:+d}$ mm\n"
        r"Error bars $=$ RMS.  Dashed band: piso estadístico $\sqrt{n}\,\tau_d/\langle N_{pe}\rangle$",
        fontsize=8, y=0.97,
    )

    color = "#228B22"  # forest green
    top_pos = list(TOP_POSITIONS_MM)

    for i, ax in enumerate(axes.flat):
        if i >= n_ch:
            ax.set_visible(False)
            continue
        ch = ids[i]
        pos_mm = top_pos[ch - 16]  # global_id = 16 + index
        ns, means, rmss, npe = _per_channel_stats(event_id, global_id, time_ns, ch)
        title = f"ID {ch} | $x_{{w}}$={pos_mm:+.0f} mm"
        _draw_pad(ax, ns, means, rmss, npe, color, title)

    for ax in axes[n_rows - 1, :]:
        ax.set_xlabel("$n$-ésimo fotón", fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\langle t_n\rangle$ [ns]", fontsize=7)

    fig.savefig(output_path, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {output_path}")

# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default="/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000",
        type=pathlib.Path,
    )
    parser.add_argument(
        "--output-dir",
        default="analysis/exec07",
        type=pathlib.Path,
    )
    args = parser.parse_args(argv)
    fig_dir = args.output_dir / "figs"
    fig_dir.mkdir(parents=True, exist_ok=True)

    end_l_ids = list(range(0, 8))
    end_r_ids = list(range(8, 16))

    for x_mm in KEY_POSITIONS:
        print(f"Processing x = {x_mm:+d} mm …")
        event_id, global_id, time_ns = _load_arrays(args.data_dir, x_mm)
        top_ids = _top_9_neighbors(x_mm)

        _make_end_figure(
            event_id, global_id, time_ns,
            end_l_ids, x_mm, "steelblue", "End L (IDs 0–7)",
            fig_dir / f"tn_{x_mm}mm_endL.png",
        )
        _make_end_figure(
            event_id, global_id, time_ns,
            end_r_ids, x_mm, "firebrick", "End R (IDs 8–15)",
            fig_dir / f"tn_{x_mm}mm_endR.png",
        )
        _make_top_figure(
            event_id, global_id, time_ns,
            top_ids, x_mm,
            fig_dir / f"tn_{x_mm}mm_top.png",
        )

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
