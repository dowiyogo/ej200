#!/usr/bin/env python3
"""
analysis_endonly_mylar.py — End-only + Mylar-wrapped EJ-204 bar simulation analysis.

Reads from DATA_DIR (default: /home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614):
  analysis/attenuation_curve.csv, attenuation_fits.csv, sigma_t_sum4.csv
  run_x*.log  (boundary census)
  run_metadata.txt
  photon_hits_x*.root  (timing histograms)

Writes to OUT_DIR (default: ./):
  figures/*.png   — DPI >= 150
  tables/*.csv
  tables/*.tex
  deck_values.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import uproot
from scipy.optimize import curve_fit

# ─── constants (physical, EJ-204, from literature / code) ───────────────────
BAR_HALF_MM   = 700.0        # bar half-length  (mm); from GDML/endonly_sum4.py
LAMBDA_BULK_CM = 160.0       # EJ-204 bulk attenuation length (cm), datasheet
SCINT_YIELD   = 10000        # EJ-204 photons/MeV (OPSC-101 from SSLG4)
TAU_RISE_NS   = 0.7          # EJ-204 rise time (ns)
TAU_DECAY_NS  = 1.8          # EJ-204 decay time (ns)
PDE_PEAK      = 0.63         # SiPM PDE at peak (420 nm)  AFBR-S4N66P014M
SIPM_MODEL    = "AFBR-S4N66P014M"
SIPM_MODEL_PDE = "AFBR-S4N66P024M"  # used for PDE curve only
DPI           = 150
C_MM_NS       = 299.792      # speed of light in mm/ns

# Representative positions for timing histograms
REPR_MM = [-690, -400, 0, 400, 690]

# ─── matplotlib style ────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.size":      10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi":     DPI,
    "axes.grid":      True,
    "grid.alpha":     0.3,
})


# ─── helpers ─────────────────────────────────────────────────────────────────

def _pos_tag(x: float) -> str:
    sign = "+" if x >= 0 else ""
    return f"x={sign}{int(x):d} mm"


def _read_metadata(data_dir: pathlib.Path) -> dict:
    meta: dict = {}
    path = data_dir / "run_metadata.txt"
    for line in path.read_text().splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            meta[k.strip()] = v.strip()
    return meta


def _parse_position_log(path: pathlib.Path) -> dict:
    """Return boundary census and final run summary from a per-position log."""
    text = path.read_text(errors="replace")

    # Boundary census: appears once after all threads merge
    census: dict[str, int] = {}
    m = re.search(r"Boundary Census.*?\n(.*?)={10}", text, re.DOTALL)
    if m:
        for line in m.group(1).splitlines():
            if ":" in line:
                key, _, val_str = line.strip().partition(":")
                try:
                    census[key.strip()] = int(val_str.strip().replace(",", ""))
                except ValueError:
                    pass

    # Find the final (aggregate) run summary = one with the most events
    summaries = re.findall(
        r"=== EJ Scintillator Bar Run Summary ===(.*?)={5,}", text, re.DOTALL
    )
    final: dict[str, float] = {}
    max_events = 0
    for block in summaries:
        fields: dict[str, float] = {}
        for line in block.splitlines():
            if ":" in line:
                k, _, v = line.strip().partition(":")
                try:
                    fields[k.strip()] = float(v.strip().replace(",", ""))
                except ValueError:
                    pass
        n = int(fields.get("Events run", 0))
        if n > max_events:
            max_events = n
            final = fields
    return {"census": census, "summary": final}


def _load_csv(path: pathlib.Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# ─── ROOT timing extraction ───────────────────────────────────────────────────

SPR_RISE_NS   = 0.5
SPR_FALL_NS   = 5.0
LE_THRESHOLD  = 4.0   # PE


def _spr_pulse(slow: float, fast: float, delta: float) -> float:
    return (slow * math.exp(-delta / SPR_FALL_NS)
            - fast * math.exp(-delta / SPR_RISE_NS))


def _leading_edge_time(arrivals: np.ndarray) -> float:
    """SUM4 leading-edge time from sorted photon arrivals. Port of endonly_sum4.py."""
    if arrivals.size == 0:
        return math.nan
    arrivals = np.sort(arrivals)
    slow = fast = 0.0
    i = 0
    while i < arrivals.size:
        cur = float(arrivals[i])
        j = i
        while j < arrivals.size and arrivals[j] == cur:
            slow += 1.0
            fast += 1.0
            j += 1
        interval = float(arrivals[j] - cur) if j < arrivals.size else math.inf
        d0 = fast / SPR_RISE_NS - slow / SPR_FALL_NS
        if d0 > 0.0:
            peak_d = math.log(fast * SPR_FALL_NS / (slow * SPR_RISE_NS)) / (
                1.0 / SPR_RISE_NS - 1.0 / SPR_FALL_NS
            )
            reach = min(peak_d, interval)
            if reach >= 0.0 and _spr_pulse(slow, fast, reach) >= LE_THRESHOLD:
                lo, hi = 0.0, reach
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if _spr_pulse(slow, fast, mid) >= LE_THRESHOLD:
                        hi = mid
                    else:
                        lo = mid
                return cur + hi
        if j >= arrivals.size:
            break
        slow *= math.exp(-interval / SPR_FALL_NS)
        fast *= math.exp(-interval / SPR_RISE_NS)
        i = j
    return math.nan


def _earliest(a: float, b: float) -> float:
    if not math.isfinite(a):
        return b
    if not math.isfinite(b):
        return a
    return min(a, b)


def _per_event_timing(root_path: pathlib.Path, n_events: int) -> dict:
    """
    Load one ROOT file and compute per-event timing quantities:
      fpt_left, fpt_right       — first photon time per end
      mean_left, mean_right     — mean photon arrival time per end
      t50_left, t50_right       — median photon arrival time per end
      sum4_left, sum4_right     — SUM4 leading-edge trigger time per end

    Returns dict with 1-D arrays of length n_events (NaN when no photons).
    """
    with uproot.open(root_path) as f:
        arr = f["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns"], library="np"
        )

    ev   = arr["event_id"].astype(np.int64)
    gid  = arr["global_id"].astype(np.int64)
    tns  = arr["time_ns"].astype(np.float64)

    n = n_events
    fpt_l  = np.full(n, np.nan); fpt_r  = np.full(n, np.nan)
    mean_l = np.full(n, np.nan); mean_r = np.full(n, np.nan)
    t50_l  = np.full(n, np.nan); t50_r  = np.full(n, np.nan)
    sum4_l = np.full(n, np.nan); sum4_r = np.full(n, np.nan)

    left_mask  = gid < 8
    right_mask = gid >= 8

    # FPT: min per event per side  (fast with np.minimum.at)
    fpt_l_work = np.full(n, np.inf); fpt_r_work = np.full(n, np.inf)
    np.minimum.at(fpt_l_work, ev[left_mask],  tns[left_mask])
    np.minimum.at(fpt_r_work, ev[right_mask], tns[right_mask])
    fpt_l = np.where(fpt_l_work < np.inf, fpt_l_work, np.nan)
    fpt_r = np.where(fpt_r_work < np.inf, fpt_r_work, np.nan)

    # Mean and t50: need to group by event_id (one pass over events with hits)
    sum_l  = np.zeros(n); cnt_l  = np.zeros(n, dtype=np.int64)
    sum_r  = np.zeros(n); cnt_r  = np.zeros(n, dtype=np.int64)
    np.add.at(sum_l, ev[left_mask],  tns[left_mask])
    np.add.at(cnt_l, ev[left_mask],  1)
    np.add.at(sum_r, ev[right_mask], tns[right_mask])
    np.add.at(cnt_r, ev[right_mask], 1)
    has_l = cnt_l > 0; has_r = cnt_r > 0
    mean_l = np.where(has_l, sum_l / cnt_l, np.nan)
    mean_r = np.where(has_r, sum_r / cnt_r, np.nan)

    # t50: median per event per side — requires sorting; do efficiently
    # Sort hits by (side<<1 | left/right, event_id) then walk groups
    for side_mask, t50_arr in ((left_mask, t50_l), (right_mask, t50_r)):
        ev_s  = ev[side_mask]
        tns_s = tns[side_mask]
        if ev_s.size == 0:
            continue
        order = np.argsort(ev_s, kind="stable")
        ev_s  = ev_s[order]
        tns_s = tns_s[order]
        unique_ev, starts = np.unique(ev_s, return_index=True)
        stops = np.r_[starts[1:], ev_s.size]
        for ue, st, en in zip(unique_ev, starts, stops):
            if 0 <= ue < n:
                t50_arr[ue] = float(np.median(tns_s[st:en]))

    # SUM4 trigger: groups {0-3},{4-7} = left; {8-11},{12-15} = right
    # Replicate endonly_sum4.py logic: earliest of the two group triggers
    triggers = np.full((n, 4), np.nan)
    combined = ev * 4 + gid // 4   # channel → group index 0..3
    order = np.argsort(combined, kind="stable")
    combined_s = combined[order]
    tns_s      = tns[order]
    unique_g, g_starts = np.unique(combined_s, return_index=True)
    g_stops = np.r_[g_starts[1:], combined_s.size]
    for ug, gs, ge in zip(unique_g, g_starts, g_stops):
        ue  = int(ug // 4)
        grp = int(ug % 4)
        if 0 <= ue < n and 0 <= grp < 4:
            triggers[ue, grp] = _leading_edge_time(tns_s[gs:ge])

    for i in range(n):
        sum4_l[i] = _earliest(triggers[i, 0], triggers[i, 1])
        sum4_r[i] = _earliest(triggers[i, 2], triggers[i, 3])

    return {
        "fpt_l": fpt_l, "fpt_r": fpt_r,
        "mean_l": mean_l, "mean_r": mean_r,
        "t50_l": t50_l, "t50_r": t50_r,
        "sum4_l": sum4_l, "sum4_r": sum4_r,
    }


# ─── figure helpers ───────────────────────────────────────────────────────────

def _save(fig: plt.Figure, path: pathlib.Path, caption: str) -> None:
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    caption_path = path.with_suffix(".txt")
    caption_path.write_text(caption)
    print(f"  [fig] {path.name}")


def _exp_model(d: np.ndarray, n0: float, lam: float) -> np.ndarray:
    return n0 * np.exp(-d / lam)


# ─── main analysis ────────────────────────────────────────────────────────────

def run(data_dir: pathlib.Path, out_dir: pathlib.Path) -> None:
    fig_dir = out_dir / "figures"
    tab_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    print("Reading metadata...")
    meta = _read_metadata(data_dir)
    n_events     = int(meta["events_per_position"])
    mylar_r      = float(meta["mylar_reflectivity"])
    mylar_lobe   = float(meta["mylar_specular_lobe"])
    mylar_sigma  = float(meta["mylar_sigma_alpha_deg"])
    positions_str = meta.get("positions", "")
    expected_pos = [float(x) for x in positions_str.split()]

    print("Reading CSVs...")
    att_rows = _load_csv(data_dir / "analysis" / "attenuation_curve.csv")
    fit_rows = _load_csv(data_dir / "analysis" / "attenuation_fits.csv")
    sig_rows = _load_csv(data_dir / "analysis" / "sigma_t_sum4.csv")

    # Convert to arrays
    x_mm       = np.array([float(r["x_mm"])          for r in att_rows])
    npe_l      = np.array([float(r["npe_left_mean"])  for r in att_rows])
    npe_l_err  = np.array([float(r["npe_left_sem"])   for r in att_rows])
    npe_r      = np.array([float(r["npe_right_mean"]) for r in att_rows])
    npe_r_err  = np.array([float(r["npe_right_sem"])  for r in att_rows])
    dist_l     = np.array([float(r["distance_left_mm"])  for r in att_rows])
    dist_r     = np.array([float(r["distance_right_mm"]) for r in att_rows])
    n_end_hits = np.array([float(r["n_end_hits"]) for r in att_rows])

    npe_sum     = npe_l + npe_r
    npe_sum_err = np.sqrt(npe_l_err**2 + npe_r_err**2)

    sig_x  = np.array([float(r["x_mm"])                for r in sig_rows])
    sig_s  = np.array([float(r["sigma_single_ps"])      for r in sig_rows])
    sig_se = np.array([float(r["sigma_single_error_ps"])for r in sig_rows])
    trig_eff = np.array([float(r["trigger_efficiency"]) for r in sig_rows])
    mean_delta = np.array([float(r["mean_delta_lr_ns"]) for r in sig_rows])

    # Attenuation fit results
    fits: dict[str, dict] = {r["side"]: r for r in fit_rows}

    print("Parsing boundary census from logs...")
    census_list: list[dict] = []
    for xv in expected_pos:
        tag = f"{int(xv):+d}mm" if xv != 0 else "x0mm"
        log_path = data_dir / f"run_x{tag}.log"
        if not log_path.exists():
            # try without sign
            log_path = data_dir / f"run_x{int(xv)}mm.log"
        if log_path.exists():
            res = _parse_position_log(log_path)
            res["x_mm"] = xv
            census_list.append(res)

    # Aggregate census across all positions
    total_census: dict[str, int] = {}
    total_summary: dict[str, float] = {}
    for cl in census_list:
        for k, v in cl["census"].items():
            total_census[k] = total_census.get(k, 0) + v
        for k, v in cl.get("summary", {}).items():
            total_summary[k] = total_summary.get(k, 0.0) + v

    print("Loading ROOT files for timing analysis...")
    # For ALL positions: compute mean FPT and mean t50
    fpt_l_mean  = np.full(len(x_mm), np.nan)
    fpt_r_mean  = np.full(len(x_mm), np.nan)
    mean_l_avg  = np.full(len(x_mm), np.nan)
    mean_r_avg  = np.full(len(x_mm), np.nan)
    t50_l_mean  = np.full(len(x_mm), np.nan)
    t50_r_mean  = np.full(len(x_mm), np.nan)

    # Per-repr-position distributions for timing histograms
    repr_timing: dict[int, dict] = {}

    root_files = sorted(data_dir.glob("photon_hits_x*.root"))
    if not root_files:
        print("  WARNING: no ROOT files found — timing histograms will be skipped.")
    else:
        for idx, xv in enumerate(x_mm):
            tag = f"{int(xv):+d}mm" if xv != 0 else "x0mm"
            rpath = data_dir / f"photon_hits_x{tag}.root"
            if not rpath.exists():
                rpath = data_dir / f"photon_hits_x{int(xv)}mm.root"
            if not rpath.exists():
                print(f"  WARNING: ROOT file missing for x={xv} mm")
                continue
            timing = _per_event_timing(rpath, n_events)

            fpt_l_mean[idx]  = float(np.nanmean(timing["fpt_l"]))
            fpt_r_mean[idx]  = float(np.nanmean(timing["fpt_r"]))
            mean_l_avg[idx]  = float(np.nanmean(timing["mean_l"]))
            mean_r_avg[idx]  = float(np.nanmean(timing["mean_r"]))
            t50_l_mean[idx]  = float(np.nanmean(timing["t50_l"]))
            t50_r_mean[idx]  = float(np.nanmean(timing["t50_r"]))

            if int(xv) in REPR_MM:
                repr_timing[int(xv)] = timing
                print(f"  Loaded {rpath.name} [repr]")
            else:
                print(f"  Loaded {rpath.name}")

    # ─── FIG 1: Npe vs position ───────────────────────────────────────────────
    print("Generating figures...")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(x_mm, npe_l, yerr=npe_l_err, fmt="o-", ms=3, lw=1.2,
                color="C0", label="End-Left (IDs 0–7)")
    ax.errorbar(x_mm, npe_r, yerr=npe_r_err, fmt="s-", ms=3, lw=1.2,
                color="C1", label="End-Right (IDs 8–15)")
    ax.errorbar(x_mm, npe_sum, yerr=npe_sum_err, fmt="^-", ms=3, lw=1.2,
                color="C2", label="Sum L+R")
    ax.set_xlabel("Beam position $x$ (mm)")
    ax.set_ylabel("Mean photon count per event (PE)")
    ax.set_title("Photon count at ends vs. beam position — End-only + Mylar, EJ-204")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(-720, 720)
    _save(fig, fig_dir / "fig_npe_vs_x.png",
          f"Mean photoelectron count at End-Left (IDs 0–7), End-Right (IDs 8–15), "
          f"and their sum vs. beam position $x$ (longitudinal axis X). "
          f"{n_events} muon events per position. Log scale. "
          f"Mylar wrapping on ±Y, ±Z faces: R={mylar_r:.2f}, specular lobe={mylar_lobe}.")

    # ─── FIG 2: Attenuation curves ────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    for ax, side, dist, npe, npe_err, color, label in [
        (axes[0], "left",  dist_l, npe_l, npe_l_err, "C0", "End-Left"),
        (axes[1], "right", dist_r, npe_r, npe_r_err, "C1", "End-Right"),
    ]:
        ax.errorbar(dist/10, npe, yerr=npe_err, fmt="o", ms=3, color=color,
                    label=label, zorder=3)
        fd = fits.get(side, {})
        if fd.get("status") == "ok":
            n0 = float(fd["n0_pe"]); n0_err = float(fd["n0_error_pe"])
            lam_cm = float(fd["lambda_eff_cm"]); lam_err = float(fd["lambda_eff_error_cm"])
            d_fit = np.linspace(dist.min(), dist.max(), 200) / 10  # cm
            ax.plot(d_fit, _exp_model(d_fit, n0, lam_cm), "-", color=color, lw=1.5,
                    label=rf"Fit: $\lambda_\mathrm{{eff}}={lam_cm:.1f}\pm{lam_err:.1f}$ cm")
        # Show bulk attenuation length for reference
        ax.axvline(LAMBDA_BULK_CM, color="gray", ls="--", lw=1,
                   label=f"Bulk $\\lambda$={LAMBDA_BULK_CM:.0f} cm")
        ax.set_xlabel("Distance from end (cm)")
        ax.set_ylabel("Mean photon count per event (PE)")
        ax.set_title(f"Attenuation — {label}")
        ax.legend(fontsize=7)
    plt.tight_layout()
    lam_l = float(fits["left"]["lambda_eff_cm"])
    lam_r = float(fits["right"]["lambda_eff_cm"])
    lam_c = float(fits["combined"]["lambda_eff_cm"])
    _save(fig, fig_dir / "fig_attenuation.png",
          f"Photon count vs. distance from the respective end, with exponential fits. "
          f"$\\lambda_\\mathrm{{eff}}$: left = {lam_l:.1f} cm, right = {lam_r:.1f} cm, "
          f"combined = {lam_c:.1f} cm. Bulk EJ-204 $\\lambda$ = {LAMBDA_BULK_CM:.0f} cm "
          f"shown for reference (Mylar recirculation keeps photons in bar longer, "
          f"but absorption per bounce reduces $\\lambda_\\mathrm{{eff}}$).")

    # ─── FIG 3: Timing histograms at representative positions ─────────────────
    if repr_timing:
        n_repr = len(repr_timing)
        fig, axes_grid = plt.subplots(n_repr, 2, figsize=(10, 2.2 * n_repr),
                                      sharex=False, sharey=False)
        if n_repr == 1:
            axes_grid = axes_grid[np.newaxis, :]

        for row_idx, xv in enumerate(REPR_MM):
            if xv not in repr_timing:
                continue
            td = repr_timing[xv]
            for col_idx, (key, label, color) in enumerate([
                ("fpt_l",   "FPT End-Left (ns)",  "C0"),
                ("fpt_r",   "FPT End-Right (ns)", "C1"),
            ]):
                ax = axes_grid[row_idx, col_idx]
                data = td[key]
                data = data[np.isfinite(data)]
                if data.size == 0:
                    ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                    continue
                n_used = data.size
                lo, hi = np.percentile(data, 1), np.percentile(data, 99)
                bins = min(100, max(20, int(np.sqrt(n_used))))
                ax.hist(data, bins=bins, range=(lo, hi), color=color, alpha=0.75,
                        histtype="stepfilled", edgecolor="none")
                ax.set_xlabel("Arrival time (ns)")
                ax.set_ylabel("Events")
                ax.set_title(f"{_pos_tag(float(xv))}, {label}, N={n_used}")
        plt.tight_layout()
        _save(fig, fig_dir / "fig_timing_histograms.png",
              "Per-event first photon arrival time (FPT) at End-Left and End-Right "
              "for representative positions. Counts (not normalized). "
              "Each panel N = events with at least one photon on that end.")

        # Also FPT vs t50 side-by-side for v_eff anomaly diagnostic
        fig2, axes2 = plt.subplots(1, 2, figsize=(10, 4))
        for col_idx, (key_fpt, key_t50, end_label, color) in enumerate([
            ("fpt_l", "t50_l", "End-Left",  "C0"),
            ("fpt_r", "t50_r", "End-Right", "C1"),
        ]):
            ax = axes2[col_idx]
            for xv in REPR_MM:
                if xv not in repr_timing:
                    continue
                td = repr_timing[xv]
                for key, ls, lbl in [
                    (key_fpt, "-",  f"FPT"),
                    (key_t50, "--", f"$t_{{50}}$"),
                ]:
                    d = td[key]
                    d = d[np.isfinite(d)]
                    if d.size == 0:
                        continue
                    lo, hi = np.percentile(d, 1), np.percentile(d, 99)
                    bins = 80
                    counts, edges = np.histogram(d, bins=bins, range=(lo, hi))
                    cx = 0.5 * (edges[:-1] + edges[1:])
                    ax.plot(cx, counts / counts.max(), ls=ls, lw=1.2,
                            label=f"x={xv:+d} {lbl}")
            ax.set_xlabel("Arrival time (ns)")
            ax.set_ylabel("Normalized count")
            ax.set_title(f"FPT vs. $t_{{50}}$ — {end_label}")
            ax.legend(fontsize=6)
        plt.tight_layout()
        _save(fig2, fig_dir / "fig_fpt_vs_t50.png",
              "Comparison of FPT (first photon) and $t_{50}$ (median photon) "
              "arrival time distributions at representative positions. "
              "Normalized to peak. Shows systematic shift of $t_{50}$ relative to FPT.")

    # ─── FIG 4: Mean arrival time vs x and v_eff ─────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ok = np.isfinite(fpt_l_mean) & np.isfinite(fpt_r_mean)
    ok_l = np.isfinite(fpt_l_mean)
    ok_r = np.isfinite(fpt_r_mean)
    ok_t50_l = np.isfinite(t50_l_mean)
    ok_t50_r = np.isfinite(t50_r_mean)

    v_eff_fpt = v_eff_t50 = v_eff_sum4 = np.nan

    if ok.sum() >= 4:
        # v_eff from FPT: slope of ⟨FPT⟩ vs x
        # Left: ⟨FPT_L⟩ = t0_L + (x + L_half) / v  → slope = +1/v
        # Right: ⟨FPT_R⟩ = t0_R + (L_half - x) / v → slope = -1/v
        for ax, fpt_arr, label, color, sign in [
            (axes[0], fpt_l_mean, "End-Left FPT",  "C0", +1),
            (axes[0], fpt_r_mean, "End-Right FPT", "C1", -1),
        ]:
            ax.plot(x_mm[ok], fpt_arr[ok], "o-", ms=3, color=color, label=label)

        # Linear fits for v_eff
        v_eff_fpt_l_mm_ns = np.nan
        v_eff_fpt_r_mm_ns = np.nan

        # Left side
        ok_l = np.isfinite(fpt_l_mean)
        if ok_l.sum() >= 4:
            cl = np.polyfit(x_mm[ok_l], fpt_l_mean[ok_l], 1)
            slope_l = cl[0]  # ns/mm
            if abs(slope_l) > 1e-6:
                v_eff_fpt_l_mm_ns = abs(1.0 / slope_l)

        ok_r = np.isfinite(fpt_r_mean)
        if ok_r.sum() >= 4:
            cr = np.polyfit(x_mm[ok_r], fpt_r_mean[ok_r], 1)
            slope_r = cr[0]  # ns/mm (should be negative)
            if abs(slope_r) > 1e-6:
                v_eff_fpt_r_mm_ns = abs(1.0 / slope_r)

        if math.isfinite(v_eff_fpt_l_mm_ns) and math.isfinite(v_eff_fpt_r_mm_ns):
            v_eff_fpt = 0.5 * (v_eff_fpt_l_mm_ns + v_eff_fpt_r_mm_ns)

        # v_eff from SUM4 Δt = t_L - t_R ≈ 2x/v
        ok_sig = np.isfinite(mean_delta)
        if ok_sig.sum() >= 4:
            cs = np.polyfit(sig_x[ok_sig], mean_delta[ok_sig], 1)
            slope_s = cs[0]   # ns/mm; Δt = 2x/v → slope = 2/v
            if abs(slope_s) > 1e-6:
                v_eff_sum4 = abs(2.0 / slope_s)

        # v_eff from t50
        ok_t50_l = np.isfinite(t50_l_mean)
        ok_t50_r = np.isfinite(t50_r_mean)
        if ok_t50_l.sum() >= 4 and ok_t50_r.sum() >= 4:
            cl50 = np.polyfit(x_mm[ok_t50_l], t50_l_mean[ok_t50_l], 1)
            cr50 = np.polyfit(x_mm[ok_t50_r], t50_r_mean[ok_t50_r], 1)
            sl = cl50[0]; sr = cr50[0]
            if abs(sl) > 1e-6 and abs(sr) > 1e-6:
                v_eff_t50 = 0.5 * (abs(1/sl) + abs(1/sr))

        # Plot mean arrival time vs x
        axes[1].plot(x_mm[ok_l], fpt_l_mean[ok_l],  "o-", ms=3, color="C0",
                     label=r"$\langle t_\mathrm{FPT}\rangle$ Left")
        axes[1].plot(x_mm[ok_r], fpt_r_mean[ok_r],  "s-", ms=3, color="C1",
                     label=r"$\langle t_\mathrm{FPT}\rangle$ Right")
        axes[1].plot(x_mm[ok_t50_l], t50_l_mean[ok_t50_l], "^--", ms=3, color="C0",
                     alpha=0.6, label=r"$\langle t_{50}\rangle$ Left")
        axes[1].plot(x_mm[ok_t50_r], t50_r_mean[ok_t50_r], "v--", ms=3, color="C1",
                     alpha=0.6, label=r"$\langle t_{50}\rangle$ Right")

    axes[0].set_xlabel("Beam position $x$ (mm)"); axes[0].set_ylabel("Mean FPT (ns)")
    axes[0].set_title("Mean first photon time vs. position")
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("Beam position $x$ (mm)"); axes[1].set_ylabel("Mean arrival time (ns)")
    axes[1].set_title("Mean FPT and $t_{50}$ vs. position")
    axes[1].legend(fontsize=7)
    plt.tight_layout()
    v_str = (f"$v_\\mathrm{{eff}}$ (FPT) ≈ {v_eff_fpt:.1f} mm/ns"
             if math.isfinite(v_eff_fpt) else "v_eff not computed")
    v50_str = (f"$v_\\mathrm{{eff}}$ ($t_{{50}}$) ≈ {v_eff_t50:.1f} mm/ns"
               if math.isfinite(v_eff_t50) else "")
    vs4_str = (f"$v_\\mathrm{{eff}}$ (SUM4 $\\Delta t$) ≈ {v_eff_sum4:.1f} mm/ns"
               if math.isfinite(v_eff_sum4) else "")
    _save(fig, fig_dir / "fig_mean_arrival_time.png",
          f"Mean first-photon time (FPT) and median photon time ($t_{{50}}$) "
          f"vs. beam position, for both ends. "
          f"{v_str}. {v50_str}. {vs4_str}. "
          f"$v_\\mathrm{{eff}}$(FPT) > $v_\\mathrm{{eff}}$($t_{{50}}$) due to "
          f"distance-dependent tail of scattered (Mylar-reflected) photons.")

    # ─── FIG 5: σ_t(x) ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(sig_x, sig_s, yerr=sig_se, fmt="o-", ms=3, color="C3",
                label=r"$\sigma_\mathrm{single} = \sigma_{LR}/\sqrt{2}$ (SUM4 core fit)")
    ax.fill_between(sig_x, sig_s - sig_se, sig_s + sig_se, alpha=0.2, color="C3")
    ax.set_xlabel("Beam position $x$ (mm)")
    ax.set_ylabel(r"$\sigma_t$ (ps)")
    ax.set_title(r"Intrinsic timing resolution $\sigma_t(x)$ — End-only + Mylar, EJ-204")
    ax.legend()
    sigma_center = float(sig_s[np.argmin(np.abs(sig_x))])
    sigma_end690 = float(np.nanmean([sig_s[np.argmin(np.abs(sig_x - 690))],
                                     sig_s[np.argmin(np.abs(sig_x + 690))]]))
    _save(fig, fig_dir / "fig_sigma_t.png",
          rf"Intrinsic single-end timing resolution $\sigma_\mathrm{{single}} = "
          rf"\sigma_{{LR}}/\sqrt{{2}}$ (SUM4 leading-edge estimator, core Gaussian fit) "
          f"vs. beam position. "
          f"Center ($x=0$): {sigma_center:.0f} ps; "
          f"Near ends ($|x|=690$ mm): {sigma_end690:.0f} ps. "
          r"$\sigma_t$ shown is \textbf{intrinsic simulation prediction} "
          r"(optical propagation + SUM4 estimator only; "
          r"SiPM SPTR and electronics jitter not folded in).")

    # ─── FIG 6: Photon budget ─────────────────────────────────────────────────
    # Use aggregate census across all positions
    key_panel   = "Bar -> reflector panel"
    key_world   = "Bar -> World (escaped)"
    key_sipm    = "Bar -> SiPM (entering)"
    key_scint   = "Scint photons generated"

    n_gen  = total_summary.get(key_scint,  0.0)
    n_sipm = float(total_census.get(key_sipm, 0))
    n_esc  = float(total_census.get(key_world, 0))
    n_ref  = float(total_census.get(key_panel, 0))

    # Detected = N_sipm * PDE
    n_det  = n_sipm * PDE_PEAK   # approximate
    # Absorbed by Mylar = roughly what's not escaped and not SiPM
    # Note: n_ref is cumulative crossings (many bounces), not unique photons

    # Show photon fate as fraction of generated (unique photon fates)
    fig, ax = plt.subplots(figsize=(7, 4))
    if n_gen > 0:
        labels = [
            "Escaped\n(Bar→World)",
            "Enter SiPM\n(before PDE)",
            "Detected\n(after PDE)",
            "Mylar-\nabsorbed",
        ]
        vals = np.array([n_esc, n_sipm, n_det,
                         max(0, n_gen - n_esc - n_sipm)])
        pct  = 100.0 * vals / n_gen
        bars = ax.barh(labels[::-1], pct[::-1], color=["C1","C3","C2","C0"])
        ax.set_xlabel("Fraction of generated photons (%)")
        ax.set_title("Photon budget — End-only + Mylar (aggregated, 31 positions)")
        for bar, p in zip(bars, pct[::-1]):
            ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
                    f"{p:.2f}%", va="center", fontsize=8)
        ax.set_xlim(0, max(pct) * 1.2)
        # Annotate Mylar reflectivity
        ax.text(0.98, 0.02, f"Mylar R={mylar_r:.2f}, specular lobe={mylar_lobe}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                color="gray")
    else:
        ax.text(0.5, 0.5, "No log data available for budget", transform=ax.transAxes, ha="center")

    plt.tight_layout()
    _save(fig, fig_dir / "fig_photon_budget.png",
          f"Photon fate budget aggregated over all 31 positions × {n_events} events. "
          f"Generated={n_gen:.2e}; escaped (Bar→World)={n_esc/n_gen*100:.2f}%; "
          f"entering SiPM={n_sipm/n_gen*100:.2f}%; "
          f"detected (after PDE {PDE_PEAK*100:.0f}%)={n_det/n_gen*100:.2f}%. "
          f"Remainder absorbed by Mylar (R={mylar_r:.2f} per bounce). "
          f"Note: 'Bar→reflector panel' crossings ({n_ref:.2e}) "
          f"count all bounces (not unique photons).")

    # ─── TABLES ───────────────────────────────────────────────────────────────
    print("Writing tables...")

    # Table 1: simulation parameters (LaTeX strings only — no Unicode)
    bar_length_mm = 2 * BAR_HALF_MM
    params_latex = [
        # [Parameter, Value, Unit, Source]  — all ASCII/LaTeX, no Unicode
        ["Parameter", "Value", "Unit", "Source"],
        ["Bar material", "EJ-204 (OPSC-101)", "", "SSLG4"],
        ["Bar dimensions", f"{bar_length_mm:.0f}\\,$\\times$\\,60\\,$\\times$\\,10",
         "mm$^3$", "GDML"],
        ["Readout", "End-only ($\\pm X$ faces)", "", "geometry"],
        ["Wrapping ($\\pm Y,\\pm Z$)", "Mylar", "", "run\\_metadata"],
        ["Mylar reflectivity", f"{mylar_r:.2f}", "", "run\\_metadata"],
        ["Mylar specular lobe", f"{mylar_lobe:.1f}", "", "run\\_metadata"],
        ["Mylar $\\sigma_\\alpha$", f"{mylar_sigma:.1f}", "deg", "run\\_metadata"],
        ["Scint.~yield", f"{SCINT_YIELD}", "$\\gamma$/MeV", "SSLG4 (EJ-204)"],
        ["Rise time $\\tau_r$", f"{TAU_RISE_NS:.1f}", "ns", "SSLG4 (EJ-204)"],
        ["Decay time $\\tau_d$", f"{TAU_DECAY_NS:.1f}", "ns", "SSLG4 (EJ-204)"],
        ["Bulk att.~$\\lambda$", f"{LAMBDA_BULK_CM:.0f}", "cm", "EJ-204 datasheet"],
        ["SiPM model", SIPM_MODEL, "", "hardware"],
        ["SiPM PDE (peak)", f"{PDE_PEAK*100:.0f}", "\\% @ 420\\,nm", "datasheet"],
        ["Beam", "$\\mu^-$, 1\\,GeV", "", "macro"],
        ["Positions", f"{len(expected_pos)}", "", "run\\_metadata"],
        ["Events/position", f"{n_events}", "", "run\\_metadata"],
        ["Geant4 version", "11.4.0", "", "log"],
    ]
    # CSV with plain text (no LaTeX)
    params_csv = [
        ["Parameter", "Value", "Unit", "Source"],
        ["Bar material", "EJ-204 (OPSC-101)", "", "SSLG4"],
        ["Bar dimensions", f"{bar_length_mm:.0f}x60x10", "mm3", "GDML"],
        ["Readout", "End-only (+-X faces)", "", "geometry"],
        ["Wrapping (+-Y,+-Z)", "Mylar", "", "run_metadata"],
        ["Mylar reflectivity", f"{mylar_r:.2f}", "", "run_metadata"],
        ["Mylar specular lobe", f"{mylar_lobe:.1f}", "", "run_metadata"],
        ["Mylar sigma_alpha", f"{mylar_sigma:.1f}", "deg", "run_metadata"],
        ["Scint yield", f"{SCINT_YIELD}", "gamma/MeV", "SSLG4 (EJ-204)"],
        ["Rise time tau_r", f"{TAU_RISE_NS:.1f}", "ns", "SSLG4 (EJ-204)"],
        ["Decay time tau_d", f"{TAU_DECAY_NS:.1f}", "ns", "SSLG4 (EJ-204)"],
        ["Bulk attenuation lambda", f"{LAMBDA_BULK_CM:.0f}", "cm", "EJ-204 datasheet"],
        ["SiPM model", SIPM_MODEL, "", "hardware"],
        ["SiPM PDE (peak)", f"{PDE_PEAK*100:.0f}", "% @ 420 nm", "datasheet"],
        ["Beam", "mu-, 1 GeV", "", "macro"],
        ["Positions", f"{len(expected_pos)}", "", "run_metadata"],
        ["Events/position", f"{n_events}", "", "run_metadata"],
        ["Geant4 version", "11.4.0", "", "log"],
    ]
    with open(tab_dir / "tab_sim_params.csv", "w", newline="") as f:
        csv.writer(f).writerows(params_csv)

    # LaTeX version (use params_latex which has LaTeX markup)
    with open(tab_dir / "tab_sim_params.tex", "w") as f:
        f.write("\\begin{tabular}{llll}\\toprule\n")
        f.write(" & ".join(params_latex[0]) + " \\\\\n\\midrule\n")
        for row in params_latex[1:]:
            f.write(" & ".join(str(c) for c in row) + " \\\\\n")
        f.write("\\bottomrule\\end{tabular}\n")

    # Table 2: attenuation fits
    lam_data = []
    lam_data.append(["Side", "$N_0$ (PE)", "$\\sigma(N_0)$",
                      "$\\lambda_\\mathrm{eff}$ (cm)", "$\\sigma(\\lambda)$", "$N$", "Status"])
    for side in ["left", "right", "combined"]:
        fd = fits.get(side, {})
        lam_data.append([
            side.capitalize(),
            f"{float(fd.get('n0_pe', math.nan)):.1f}",
            f"{float(fd.get('n0_error_pe', math.nan)):.1f}",
            f"{float(fd.get('lambda_eff_cm', math.nan)):.2f}",
            f"{float(fd.get('lambda_eff_error_cm', math.nan)):.2f}",
            str(fd.get("n_points", "")),
            str(fd.get("status", "")),
        ])
    with open(tab_dir / "tab_attenuation_fits.csv", "w", newline="") as f:
        csv.writer(f).writerows(lam_data)

    with open(tab_dir / "tab_attenuation_fits.tex", "w") as f:
        f.write("\\begin{tabular}{lrrrrl}\\toprule\n")
        hdr = lam_data[0]
        f.write(f"{hdr[0]} & {hdr[1]} & {hdr[2]} & {hdr[3]} & {hdr[4]} & {hdr[5]} \\\\\n\\midrule\n")
        for row in lam_data[1:]:
            f.write(f"{row[0]} & {row[1]} & $\\pm${row[2]} & {row[3]} & $\\pm${row[4]} & {row[5]} \\\\\n")
        f.write("\\bottomrule\\end{tabular}\n")

    # Table 3: sigma_t at key positions (LaTeX notation only)
    key_xvals = [-690, -400, 0, 400, 690]
    sig_table_csv = [["x (mm)", "N triggered", "Trig. eff.", "sigma_single (ps)", "sigma_err (ps)"]]
    sig_table_tex_rows = [["$x$ (mm)", "$N_\\mathrm{trig}$",
                           "Trig.~eff.", "$\\sigma_\\mathrm{single}$ (ps)", "Err. (ps)"]]
    for kx in key_xvals:
        idx_arr = np.argmin(np.abs(sig_x - kx))
        row_csv = [
            f"{int(sig_x[idx_arr]):+d}",
            str(int(float(sig_rows[idx_arr]["n_triggered_lr"]))),
            f"{trig_eff[idx_arr]:.3f}",
            f"{sig_s[idx_arr]:.1f}",
            f"{sig_se[idx_arr]:.1f}",
        ]
        row_tex = [
            f"${int(sig_x[idx_arr]):+d}$",
            str(int(float(sig_rows[idx_arr]["n_triggered_lr"]))),
            f"{trig_eff[idx_arr]:.3f}",
            f"{sig_s[idx_arr]:.1f}",
            f"{sig_se[idx_arr]:.1f}",
        ]
        sig_table_csv.append(row_csv)
        sig_table_tex_rows.append(row_tex)
    with open(tab_dir / "tab_sigma_t.csv", "w", newline="") as f:
        csv.writer(f).writerows(sig_table_csv)

    with open(tab_dir / "tab_sigma_t.tex", "w") as f:
        f.write("\\begin{tabular}{rrrrl}\\toprule\n")
        hdr = sig_table_tex_rows[0]
        f.write(" & ".join(hdr) + " \\\\\n\\midrule\n")
        for row in sig_table_tex_rows[1:]:
            f.write(" & ".join(row) + " \\\\\n")
        f.write("\\bottomrule\\end{tabular}\n")

    # ─── deck_values.json ─────────────────────────────────────────────────────
    print("Writing deck_values.json...")

    def _safe(v) -> object:
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return None
        if isinstance(v, (np.floating, np.integer)):
            return float(v)
        return v

    # σ_t at representative positions
    sigma_at = {}
    for kx in key_xvals:
        idx_arr = int(np.argmin(np.abs(sig_x - kx)))
        sigma_at[str(kx)] = {
            "sigma_single_ps": _safe(float(sig_s[idx_arr])),
            "sigma_single_err_ps": _safe(float(sig_se[idx_arr])),
        }

    deck = {
        "source": "analysis_endonly_mylar.py",
        "data_dir": str(data_dir),
        "n_events_per_position": n_events,
        "n_positions": len(x_mm),
        "bar_half_length_mm": BAR_HALF_MM,
        "bar_length_mm": 2 * BAR_HALF_MM,
        "lambda_bulk_cm": LAMBDA_BULK_CM,
        "scint_yield_per_MeV": SCINT_YIELD,
        "tau_rise_ns": TAU_RISE_NS,
        "tau_decay_ns": TAU_DECAY_NS,
        "pde_peak_frac": PDE_PEAK,
        "sipm_model": SIPM_MODEL,
        "sipm_model_pde_curve": SIPM_MODEL_PDE,
        "mylar_reflectivity": mylar_r,
        "mylar_specular_lobe": mylar_lobe,
        "mylar_sigma_alpha_deg": mylar_sigma,
        "lambda_eff_left_cm": _safe(float(fits["left"]["lambda_eff_cm"])),
        "lambda_eff_left_err_cm": _safe(float(fits["left"]["lambda_eff_error_cm"])),
        "lambda_eff_right_cm": _safe(float(fits["right"]["lambda_eff_cm"])),
        "lambda_eff_right_err_cm": _safe(float(fits["right"]["lambda_eff_error_cm"])),
        "lambda_eff_combined_cm": _safe(float(fits["combined"]["lambda_eff_cm"])),
        "lambda_eff_combined_err_cm": _safe(float(fits["combined"]["lambda_eff_error_cm"])),
        "sigma_t_at_key_positions": sigma_at,
        "sigma_t_center_ps": _safe(sigma_center),
        "sigma_t_end690_ps": _safe(sigma_end690),
        "v_eff_fpt_mm_ns": _safe(v_eff_fpt),
        "v_eff_t50_mm_ns": _safe(v_eff_t50),
        "v_eff_sum4_mm_ns": _safe(v_eff_sum4),
        "photon_budget_generated": _safe(n_gen),
        "photon_budget_escaped": _safe(n_esc),
        "photon_budget_sipm_entering": _safe(n_sipm),
        "photon_budget_detected": _safe(n_det),
        "photon_budget_mylar_absorbed": _safe(max(0.0, n_gen - n_esc - n_sipm)),
        "photon_budget_reflector_crossings": _safe(float(n_ref)),
        "npe_left_at_x0": _safe(float(npe_l[np.argmin(np.abs(x_mm))])),
        "npe_right_at_x0": _safe(float(npe_r[np.argmin(np.abs(x_mm))])),
        "npe_left_at_x690": _safe(float(npe_l[np.argmin(np.abs(x_mm - 690))])),
        "npe_right_at_x690": _safe(float(npe_r[np.argmin(np.abs(x_mm - 690))])),
        "npe_left_at_xm690": _safe(float(npe_l[np.argmin(np.abs(x_mm + 690))])),
        "npe_right_at_xm690": _safe(float(npe_r[np.argmin(np.abs(x_mm + 690))])),
        "figures": {
            "fig_npe_vs_x":            "figures/fig_npe_vs_x.png",
            "fig_attenuation":         "figures/fig_attenuation.png",
            "fig_timing_histograms":   "figures/fig_timing_histograms.png",
            "fig_fpt_vs_t50":          "figures/fig_fpt_vs_t50.png",
            "fig_mean_arrival_time":   "figures/fig_mean_arrival_time.png",
            "fig_sigma_t":             "figures/fig_sigma_t.png",
            "fig_photon_budget":       "figures/fig_photon_budget.png",
        },
        "tables": {
            "tab_sim_params":       "tables/tab_sim_params",
            "tab_attenuation_fits": "tables/tab_attenuation_fits",
            "tab_sigma_t":          "tables/tab_sigma_t",
        },
    }

    json_path = out_dir / "deck_values.json"
    with open(json_path, "w") as f:
        json.dump(deck, f, indent=2)
    print(f"  [json] {json_path.name}")

    # Check for NaN/inf in fits
    print("\nSanity checks:")
    for side in ["left", "right", "combined"]:
        lam = float(fits[side].get("lambda_eff_cm", math.nan))
        status = fits[side].get("status", "?")
        ok_flag = "OK" if (math.isfinite(lam) and lam > 0) else "FAIL"
        print(f"  λ_eff {side}: {lam:.2f} cm [{ok_flag}] ({status})")
    nan_sigma = np.sum(~np.isfinite(sig_s))
    print(f"  σ_t NaN/inf: {nan_sigma}/{len(sig_s)}")
    if v_eff_fpt and math.isfinite(v_eff_fpt):
        print(f"  v_eff FPT: {v_eff_fpt:.1f} mm/ns = {v_eff_fpt/C_MM_NS:.3f}c")
    if v_eff_sum4 and math.isfinite(v_eff_sum4):
        print(f"  v_eff SUM4: {v_eff_sum4:.1f} mm/ns = {v_eff_sum4/C_MM_NS:.3f}c")
    print("\nDone.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614",
                    type=pathlib.Path)
    ap.add_argument("--out-dir", default=pathlib.Path(__file__).parent, type=pathlib.Path)
    args = ap.parse_args()
    run(args.data_dir, args.out_dir)
