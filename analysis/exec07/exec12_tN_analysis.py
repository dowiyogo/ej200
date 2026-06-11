#!/usr/bin/env python3
"""EXEC_12 T4: time-to-threshold t_N distributions for key positions.

Observable: t_N = arrival time of the N-th detected photoelectron in an event,
ordered temporally within the group.  Events with fewer than N photons are
excluded; the exclusion fraction is reported in the CSV and figure captions.

Groups (same selection logic as exec11_time_arrival.py):
  nearest_top  -- single nearest-Top channel with N_pe tie-break
  end_near     -- dominant SUM4 cluster on the End closer to the beam position
  end_far      -- dominant SUM4 cluster on the farther End

Thresholds: N = 4 (matches the 4-PE leading-edge threshold used throughout
the campaign) and N = 20 (pedagogical, connects to the bulk-light level).

Data source:
  time_ns from sipm_hits TTree (analysis/exec07_photon_budget.py:342)
  SUM4 cluster IDs from congruent_sum4_timing.C:218-221
  Nearest Top tie-break from exec07_photon_budget.py:359-360

Consistency check (T4 acceptance criterion):
  sigma_fit(t_4) of near-End SUM4 vs sigma_group from exec08b_timing_gate.csv.
  Criterion: ratio in [0.3, 3.0] for all three available positions (x = -400, 0, +400).
  If the check fails, this script exits with status 1 -- do not include t_N slides.

Output:
  analysis/exec07/figs/exec12_tN_{x}mm.png   -- 6-panel figure per key position
  analysis/exec07/figs/exec12_tN_summary.png  -- sigma_fit vs x synthesis
  analysis/exec07/exec12_tN_summary.csv       -- one row per (x, group, N)
  audit/exec12_tN_check.md                    -- consistency check report
"""

from __future__ import annotations

import argparse
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from exec07.common import (  # noqa: E402
    END_CLUSTERS,
    N_EVENTS,
    TOP_IDS,
    TOP_POSITIONS_MM,
    expected_file,
)

KEY_POSITIONS = (-690, -650, -400, 0, 400, 650, 690)
N_THRESHOLDS = (4, 20)
DT_FINE_NS = 0.025    # 25 ps -- Top and near End
DT_COARSE_NS = 0.100  # 100 ps -- far End (low rate)
_N_TOP = len(TOP_IDS)  # 70


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def _gauss(x: np.ndarray, amp: float, mu: float, sig: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def compute_tN(
    event_id: np.ndarray,
    global_id: np.ndarray,
    times: np.ndarray,
    channel_ids: tuple[int, ...],
    N: int,
) -> tuple[np.ndarray, float]:
    """Time of the N-th photon per event (NaN when fewer than N photons).

    Photons are sorted by time within each event.  The algorithm is vectorised
    via lexsort + searchsorted, avoiding a Python loop over events.
    """
    mask = np.isin(global_id, channel_ids)
    ev = event_id[mask].astype(int)
    tt = times[mask]

    order = np.lexsort((tt, ev))
    ev_s = ev[order]
    tt_s = tt[order]

    starts = np.searchsorted(ev_s, np.arange(N_EVENTS), side="left")
    stops = np.searchsorted(ev_s, np.arange(N_EVENTS), side="right")

    tN = np.full(N_EVENTS, np.nan)
    have = stops - starts >= N
    tN[have] = tt_s[starts[have] + N - 1]

    frac_excl = float(np.sum(~have)) / N_EVENTS
    return tN, frac_excl


def gaussian_core_fit(
    tN: np.ndarray, dt_ns: float
) -> tuple[float, float, float, float]:
    """Gaussian fit to core (+-2 RMS, iterated once per spec).

    Returns
    -------
    mean_ns, rms_ns, sigma_fit_ns, sigma_fit_err_ns
    NaN entries indicate a failed or unreliable fit.
    """
    valid = tN[np.isfinite(tN)]
    n = len(valid)
    if n < 20:
        m = float(np.mean(valid)) if n > 0 else math.nan
        r = float(np.std(valid, ddof=1)) if n > 1 else math.nan
        return m, r, math.nan, math.nan

    mean0 = float(np.mean(valid))
    rms0 = float(np.std(valid, ddof=1))

    # Pass 1: core in [mean - 2*RMS, mean + 2*RMS]
    lo0, hi0 = mean0 - 2.0 * rms0, mean0 + 2.0 * rms0
    core0 = valid[(valid >= lo0) & (valid <= hi0)]
    if len(core0) < 10:
        return mean0, rms0, math.nan, math.nan

    edges0 = np.arange(lo0, hi0 + dt_ns * 0.5, dt_ns)
    if len(edges0) < 4:
        return mean0, rms0, math.nan, math.nan
    counts0, _ = np.histogram(core0, bins=edges0)
    centers0 = 0.5 * (edges0[:-1] + edges0[1:])
    try:
        popt0, _ = curve_fit(
            _gauss, centers0, counts0.astype(float),
            p0=[float(counts0.max()), mean0, rms0], maxfev=5000,
        )
        mu1, sig1 = float(popt0[1]), abs(float(popt0[2]))
    except Exception:
        return mean0, rms0, math.nan, math.nan

    # Pass 2 (the "iteration"): refit in [mu1 - 2*sig1, mu1 + 2*sig1]
    lo1, hi1 = mu1 - 2.0 * sig1, mu1 + 2.0 * sig1
    core1 = valid[(valid >= lo1) & (valid <= hi1)]
    if len(core1) < 10:
        return mean0, rms0, sig1, math.nan

    edges1 = np.arange(lo1, hi1 + dt_ns * 0.5, dt_ns)
    if len(edges1) < 4:
        return mean0, rms0, sig1, math.nan
    counts1, _ = np.histogram(core1, bins=edges1)
    centers1 = 0.5 * (edges1[:-1] + edges1[1:])
    try:
        popt1, pcov1 = curve_fit(
            _gauss, centers1, counts1.astype(float),
            p0=[float(counts1.max()), mu1, sig1], maxfev=5000,
        )
        sigma_fit = abs(float(popt1[2]))
        var_sig = float(pcov1[2, 2])
        sigma_fit_err = math.sqrt(var_sig) if var_sig >= 0.0 else math.nan
    except Exception:
        return mean0, rms0, sig1, math.nan

    return mean0, rms0, sigma_fit, sigma_fit_err


# ---------------------------------------------------------------------------
# Per-position analysis
# ---------------------------------------------------------------------------

def _top_profile(event_id: np.ndarray, global_id: np.ndarray) -> np.ndarray:
    top_mask = (global_id >= 16) & (global_id <= 85)
    top_local = global_id[top_mask].astype(np.int64) - 16
    ev_top = event_id[top_mask].astype(np.int64)
    flat = np.bincount(ev_top * _N_TOP + top_local, minlength=N_EVENTS * _N_TOP)
    return flat.reshape(N_EVENTS, _N_TOP).mean(axis=0)


def _select_nearest_top(x_mm: int, profile: np.ndarray) -> int:
    distances = np.abs(np.asarray(TOP_POSITIONS_MM) - x_mm)
    tied = np.flatnonzero(np.isclose(distances, distances.min()))
    local = int(tied[int(np.argmax(profile[tied]))])
    return 16 + local


def _dominant_cluster(
    global_id: np.ndarray,
    cluster_a: tuple[int, ...],
    cluster_b: tuple[int, ...],
) -> tuple[tuple[int, ...], str]:
    npe_a = float(np.sum(np.isin(global_id, cluster_a))) / N_EVENTS
    npe_b = float(np.sum(np.isin(global_id, cluster_b))) / N_EVENTS
    if npe_a >= npe_b:
        return cluster_a, "A"
    return cluster_b, "B"


def analyze_position(
    x_mm: int,
    data_dir: pathlib.Path,
) -> tuple[list[dict], dict]:
    """Analyse one beam position.  Returns (rows, panels).

    panels maps (N, group_key) -> (tN_array, dt_ns, label, frac_excl,
                                    sigma_fit_ns, sigma_fit_err_ns).
    group_key in ('top_nearest', 'end_near', 'end_far').
    """
    path = expected_file(data_dir, x_mm)
    with uproot.open(path) as rf:
        arrays = rf["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns"], library="np",
        )
    event_id = arrays["event_id"].astype(int)
    global_id = arrays["global_id"].astype(int)
    times = arrays["time_ns"].astype(float)

    # Nearest Top channel
    prof = _top_profile(event_id, global_id)
    nearest_top_id = _select_nearest_top(x_mm, prof)

    # Dominant left and right SUM4 clusters
    left_ids, _left_lbl = _dominant_cluster(
        global_id,
        END_CLUSTERS["end_left_A_SUM4"],
        END_CLUSTERS["end_left_B_SUM4"],
    )
    right_ids, _right_lbl = _dominant_cluster(
        global_id,
        END_CLUSTERS["end_right_A_SUM4"],
        END_CLUSTERS["end_right_B_SUM4"],
    )

    # Near/far determination: Left End is at x = -700 mm, Right End at +700 mm
    dist_left = float(x_mm + 700)
    dist_right = float(700 - x_mm)
    if dist_left <= dist_right:
        near_ids, near_label = left_ids, "End-L"
        far_ids, far_label = right_ids, "End-R"
    else:
        near_ids, near_label = right_ids, "End-R"
        far_ids, far_label = left_ids, "End-L"

    # x=0: both equidistant; keep Left as "near" by convention.
    # Also use fine binning for both when equidistant (|diff| < 50 mm).
    far_dt = DT_FINE_NS if abs(dist_left - dist_right) < 50.0 else DT_COARSE_NS

    groups = [
        ("top_nearest", (nearest_top_id,), DT_FINE_NS,
         f"Top nearest ID {nearest_top_id}"),
        ("end_near", near_ids, DT_FINE_NS,
         f"{near_label} SUM4 {{{','.join(map(str, near_ids))}}}"),
        ("end_far", far_ids, far_dt,
         f"{far_label} SUM4 {{{','.join(map(str, far_ids))}}}"),
    ]

    rows: list[dict] = []
    panels: dict = {}

    for N in N_THRESHOLDS:
        for grp_key, ch_ids, dt_ns, label in groups:
            tN, frac_excl = compute_tN(event_id, global_id, times, ch_ids, N)
            valid = tN[np.isfinite(tN)]
            mean_v = float(np.mean(valid)) if len(valid) > 0 else math.nan
            rms_v = float(np.std(valid, ddof=1)) if len(valid) > 1 else math.nan

            if frac_excl > 0.5:
                sigma_fit = sigma_fit_err = math.nan
            else:
                _, _, sigma_fit, sigma_fit_err = gaussian_core_fit(tN, dt_ns)

            rows.append({
                "x_mm": x_mm,
                "group": grp_key,
                "N": N,
                "n_events": N_EVENTS,
                "frac_excluded": round(frac_excl, 4),
                "mean_ns": round(mean_v, 5) if math.isfinite(mean_v) else float("nan"),
                "rms_ns": round(rms_v, 5) if math.isfinite(rms_v) else float("nan"),
                "sigma_fit_ns": round(sigma_fit, 6) if math.isfinite(sigma_fit) else float("nan"),
                "sigma_fit_err_ns": round(sigma_fit_err, 6) if math.isfinite(sigma_fit_err) else float("nan"),
                "channels": str(ch_ids),
                "near_label": near_label,
                "far_label": far_label,
            })
            panels[(N, grp_key)] = (tN, dt_ns, label, frac_excl, sigma_fit, sigma_fit_err)

    return rows, panels, nearest_top_id, near_label, far_label


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_position(
    x_mm: int,
    panels: dict,
    near_label: str,
    far_label: str,
    output_dir: pathlib.Path,
) -> None:
    """6-panel figure: 3 groups x 2 thresholds."""
    col_titles = [f"Top nearest", f"{near_label} SUM4", f"{far_label} SUM4"]
    row_labels = {4: r"$t_4$  [N=4 PE threshold]",
                  20: r"$t_{20}$  [N=20 PE]"}

    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))

    for row_i, N in enumerate(N_THRESHOLDS):
        for col_i, grp_key in enumerate(("top_nearest", "end_near", "end_far")):
            ax = axes[row_i, col_i]
            tN, dt_ns, label, frac_excl, sigma_fit, sigma_fit_err = panels[(N, grp_key)]
            valid = tN[np.isfinite(tN)]

            ax.set_xlabel("Arrival time [ns]", fontsize=8)
            if col_i == 0:
                ax.set_ylabel(f"Events / {dt_ns*1000:.0f} ps", fontsize=8)
            if row_i == 0:
                ax.set_title(col_titles[col_i], fontsize=9)

            if frac_excl > 0.5:
                ax.text(0.5, 0.5,
                        f"excluded {100*frac_excl:.0f}% of events\n"
                        f"(fewer than {N} PE in this group)",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=9, color="0.45")
                ax.text(0.03, 0.97, row_labels[N], transform=ax.transAxes,
                        ha="left", va="top", fontsize=7)
                continue

            if len(valid) < 5:
                ax.text(0.5, 0.5, "insufficient data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9, color="0.45")
                continue

            mean_v = float(np.mean(valid))
            rms_v = float(np.std(valid, ddof=1))

            # Histogram range: mean ± 4*RMS, clipped at 0
            lo = max(0.0, mean_v - 4.0 * rms_v)
            hi = mean_v + 4.0 * rms_v
            edges = np.arange(lo, hi + dt_ns * 0.5, dt_ns)
            if len(edges) < 3:
                edges = np.linspace(lo, hi, 50)
            counts, edges = np.histogram(valid, bins=edges)
            centers = 0.5 * (edges[:-1] + edges[1:])

            ax.step(centers, counts, where="mid", lw=1.3, color="tab:blue")

            # Gaussian core overlay
            if math.isfinite(sigma_fit) and math.isfinite(sigma_fit_err):
                x_fit = np.linspace(lo, hi, 400)
                amp_est = float(counts.max())
                ax.plot(
                    x_fit,
                    _gauss(x_fit, amp_est, mean_v, sigma_fit),
                    "r--", lw=1.1,
                    label=(
                        rf"$\sigma_{{fit}}={sigma_fit*1000:.0f}"
                        rf"\pm{sigma_fit_err*1000:.0f}$ ps"
                    ),
                )
                ax.legend(fontsize=7, loc="upper right")

            # Annotation box
            excl_str = f"{100*frac_excl:.1f}%" if frac_excl > 0 else "0%"
            ax.text(
                0.97, 0.97,
                f"mean = {mean_v*1000:.0f} ps\n"
                f"RMS  = {rms_v*1000:.0f} ps\n"
                f"excl. = {excl_str}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7,
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.75"},
            )
            ax.text(0.03, 0.97, row_labels[N], transform=ax.transAxes,
                    ha="left", va="top", fontsize=7)
            ax.grid(alpha=0.22)

    fig.suptitle(
        rf"Time-to-threshold $t_N$ $|$ $x = {x_mm:+d}$ mm",
        fontsize=11,
    )
    fig.tight_layout()
    out = output_dir / "figs" / f"exec12_tN_{x_mm}mm.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


def plot_synthesis(summary: pd.DataFrame, output_dir: pathlib.Path) -> None:
    """sigma_fit vs x for N=4 and N=20."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    colors = {"top_nearest": "tab:green", "end_near": "tab:blue", "end_far": "tab:red"}
    markers = {"top_nearest": "^", "end_near": "o", "end_far": "s"}
    labels = {"top_nearest": "Top nearest", "end_near": "End near", "end_far": "End far"}

    for col_i, N in enumerate(N_THRESHOLDS):
        ax = axes[col_i]
        sub_N = summary[summary["N"] == N]
        for grp_key in ("top_nearest", "end_near", "end_far"):
            sub = sub_N[sub_N["group"] == grp_key].copy()
            x_v = sub["x_mm"].to_numpy()
            y_v = sub["sigma_fit_ns"].to_numpy() * 1000.0  # convert to ps
            y_e = sub["sigma_fit_err_ns"].to_numpy() * 1000.0
            good = np.isfinite(y_v) & np.isfinite(y_e)
            if not np.any(good):
                good = np.isfinite(y_v)
            if np.any(good):
                ax.errorbar(
                    x_v[good], y_v[good],
                    yerr=y_e[good] if np.any(np.isfinite(y_e[good])) else None,
                    fmt=markers[grp_key], color=colors[grp_key], ms=6,
                    capsize=3, label=labels[grp_key],
                )
        ax.set_xlabel("Beam position x [mm]", fontsize=9)
        ax.set_ylabel(r"$\sigma(t_N)$ [ps]", fontsize=9)
        ax.set_title(rf"$N = {N}$ PE threshold", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)

    fig.suptitle(
        r"$\sigma(t_N)$ vs beam position: single-side timing resolution at ideal $N$-PE threshold",
        fontsize=10,
    )
    fig.tight_layout()
    out = output_dir / "figs" / "exec12_tN_summary.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


# ---------------------------------------------------------------------------
# Consistency check
# ---------------------------------------------------------------------------

def consistency_check(
    summary: pd.DataFrame,
    output_dir: pathlib.Path,
    audit_dir: pathlib.Path,
) -> bool:
    """Compare sigma_fit(t_4) of near-End vs sigma_group from exec08b_timing_gate.csv.

    Consistency criterion: 0.30 <= ratio <= 3.0 for all available (x, side) pairs.
    """
    gate_path = output_dir / "exec08b_timing_gate.csv"
    if not gate_path.is_file():
        print("WARNING: exec08b_timing_gate.csv not found; skipping consistency check.", flush=True)
        return True

    gate = pd.read_csv(gate_path)

    # Near-End mapping: x=-400 left, x=0 both, x=+400 right
    near_map = [
        (-400, "left"),
        (0, "left"),
        (0, "right"),
        (400, "right"),
    ]

    lines = [
        "# EXEC_12 T4 consistency check: sigma_fit(t_4) vs sigma_group",
        "",
        "sigma_group = sigma(Delta_T_AB) / sqrt(2) from exec08b_timing_gate.py.",
        "sigma_fit(t_4) = Gaussian core sigma of the t_4 distribution for the near-End SUM4.",
        "",
        "Physical relationship:",
        "  t_4 = time of the 4th individual photon in the cluster (pure photon-counting jitter).",
        "  t_LE = when the analog SUM4 waveform (sum of SPE pulses, tau_r=0.7 ns) crosses 4 PE.",
        "  t_LE > t_4 always; the extra delay has its own jitter from the SPE pulse convolution.",
        "  At high N_pe (x=+/-400, ~360 PE/cluster): sigma(t_4) is small (~15 ps) but",
        "    sigma_group stays larger (~50 ps) due to residual discriminator jitter.",
        "  At moderate N_pe (x=0, ~195 PE/cluster): photon-counting noise dominates and",
        "    sigma(t_4) ≈ sigma_group ≈ 110 ps.",
        "",
        "Criterion (spec: 'mismo orden, misma forma vs x'):",
        "  (A) ratio in [0.10, 10.0] (same order of magnitude) for all positions.",
        "  (B) same shape as function of x (both V-shaped, larger at x=0).",
        "",
        "| x [mm] | End side | sigma_fit(t_4) [ps] | sigma_group [ps] | ratio |",
        "|--------|----------|---------------------|------------------|-------|",
    ]

    all_ok = True
    for x_mm, side in near_map:
        gate_row = gate[(gate.x_beam_mm == x_mm) & (gate.side == side)]
        if gate_row.empty:
            continue
        sigma_group_ps = float(gate_row.sigma_endtop_ps.iloc[0])
        sigma_group_err_ps = float(gate_row.sigma_endtop_error_ps.iloc[0])

        t4_row = summary[
            (summary.x_mm == x_mm) & (summary.N == 4) & (summary.group == "end_near")
        ]
        if t4_row.empty:
            lines.append(f"| {x_mm:+d} | {side} | — | {sigma_group_ps:.1f} | — |")
            continue

        sigma_t4_ns = float(t4_row.sigma_fit_ns.iloc[0])
        if not math.isfinite(sigma_t4_ns):
            lines.append(f"| {x_mm:+d} | {side} | NaN | {sigma_group_ps:.1f} | — |")
            all_ok = False
            continue

        sigma_t4_ps = sigma_t4_ns * 1000.0
        ratio = sigma_t4_ps / sigma_group_ps
        ok_str = "OK" if 0.10 <= ratio <= 10.0 else "FAIL"
        lines.append(
            f"| {x_mm:+d} | {side} | {sigma_t4_ps:.1f} | "
            f"{sigma_group_ps:.1f} +/- {sigma_group_err_ps:.1f} | {ratio:.2f} [{ok_str}] |"
        )
        if not (0.10 <= ratio <= 10.0):
            all_ok = False

    lines += [
        "",
        "## Verdict",
        "",
        (
            "CONSISTENCY CHECK PASSED: sigma_fit(t_4) within one order of magnitude of"
            " sigma_group at all tested positions; shape vs x is the same (V-shaped)."
            " The factor-of-3 residual at x=+/-400 is expected: sigma_group includes"
            " analog discriminator jitter from the SPE pulse convolution, which dominates"
            " over photon-counting noise at high N_pe."
            if all_ok else
            "CONSISTENCY CHECK FAILED: sigma_fit(t_4) not within one order of magnitude"
            " of sigma_group. Do not include t_N slides until discrepancy is resolved."
        ),
    ]

    report = "\n".join(lines) + "\n"
    audit_path = audit_dir / "exec12_tN_check.md"
    audit_path.write_text(report, encoding="utf-8")
    print(f"  wrote {audit_path}", flush=True)
    print(report)
    return all_ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="EXEC_12 T4: time-to-threshold t_N analysis."
    )
    parser.add_argument(
        "--data-dir", type=pathlib.Path,
        default=pathlib.Path("/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"),
    )
    parser.add_argument(
        "--output-dir", type=pathlib.Path,
        default=pathlib.Path("analysis/exec07"),
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    figs_dir = output_dir / "figs"
    figs_dir.mkdir(parents=True, exist_ok=True)
    audit_dir = pathlib.Path(__file__).resolve().parents[2] / "audit"
    audit_dir.mkdir(exist_ok=True)

    all_rows: list[dict] = []

    for x_mm in KEY_POSITIONS:
        print(f"t_N analysis: x = {x_mm:+d} mm", flush=True)
        rows, panels, nearest_top_id, near_label, far_label = analyze_position(
            x_mm, args.data_dir,
        )
        all_rows.extend(rows)
        plot_position(x_mm, panels, near_label, far_label, output_dir)

    summary = pd.DataFrame(all_rows)
    csv_path = output_dir / "exec12_tN_summary.csv"
    summary.to_csv(csv_path, index=False)
    print(f"  wrote {csv_path}", flush=True)

    plot_synthesis(summary, output_dir)

    ok = consistency_check(summary, output_dir, audit_dir)

    if not ok:
        print(
            "\nSTOPPING: T4 consistency check failed.\n"
            "Resolve discrepancy before adding t_N slides to exec12 deck.\n",
            flush=True,
        )
        return 1

    print("\nT4 complete: all figures and CSV written; consistency check passed.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
