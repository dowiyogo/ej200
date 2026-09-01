#!/usr/bin/env python3
"""
analyze_propagation_velocity.py

Extracts the effective light propagation velocity in the EJ-200/EJ-204 bar
from the EndTop scan.

Method:
  - For each gun position x, take the FIRST photon arrival time at END_L and END_R.
  - The median of that time distribution vs x is linear: t = x/v_eff + offset.
  - The slope gives 1/v_eff directly.
  - Also uses the half-difference (t_L - t_R)/2 = x/v_eff, which cancels the offset.

Outputs (in <out_dir>/):
  velocity/veff_curve.png          main figure (4 panels)
  velocity/veff_summary.csv        per-position table
  velocity/veff_results.txt        fitted v_eff, slope, intercept, uncertainty
"""

import argparse
import csv
import os
import sys
from pathlib import Path
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import uproot
from scipy.optimize import curve_fit

C_CM_NS = 29.9792458   # speed of light in cm/ns
BAR_HALF_LENGTH = 690.0  # mm


def first_photon_medians(root_path):
    """
    Load one ROOT file. Return per-event first-photon times for END_L and END_R.
    Also return (n_total_events, n_L, n_R).
    """
    with uproot.open(str(root_path)) as f:
        tree = f["sipm_hits"]
        arr = tree.arrays(["event_id", "face_type", "time_ns"], library="np")

    evt = arr["event_id"]
    ft  = arr["face_type"]
    t   = arr["time_ns"]
    n_total = len(np.unique(evt))

    def t_first(face_id):
        mask = ft == face_id
        e_ = evt[mask]; t_ = t[mask]
        if len(e_) == 0:
            return np.array([])
        order = np.argsort(e_, kind="stable")
        e_s = e_[order]; t_s = t_[order]
        splits = np.where(np.diff(e_s) != 0)[0] + 1
        return np.array([np.min(tg) for tg in np.split(t_s, splits)],
                        dtype=np.float64)

    tL = t_first(0)
    tR = t_first(1)
    return tL, tR, n_total


def robust_median_err(arr):
    """Return (median, MAD/sqrt(n)) as (value, err)."""
    if len(arr) == 0:
        return float("nan"), float("nan")
    med = np.median(arr)
    mad = 1.4826 * np.median(np.abs(arr - med))
    err = mad / np.sqrt(len(arr))
    return med, err


def linear(x, a, b):
    return a * x + b


def fit_line(xs, ys, errs):
    """Fit y = a*x + b with weighted least squares. Return popt, perr."""
    ok = np.isfinite(ys) & np.isfinite(errs) & (errs > 0)
    if ok.sum() < 3:
        return None, None
    try:
        popt, pcov = curve_fit(linear, xs[ok], ys[ok],
                               sigma=errs[ok], absolute_sigma=True,
                               p0=[0.013, 5.0])
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except Exception:
        return None, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    RUN_DIR = Path(args.run_dir)
    OUT_DIR = Path(args.out_dir) / "velocity"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Discover positions from summary.tsv
    tsv = RUN_DIR / "summary.tsv"
    scan_points = []
    with open(tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos_mm  = int(row["position_mm"])
            pos_tag = row["index"]
            rpath   = RUN_DIR / "outputs" / pos_tag / "photon_hits_run000.root"
            if rpath.exists():
                scan_points.append((pos_mm, rpath))
    scan_points.sort(key=lambda p: p[0])
    n_pos = len(scan_points)
    print(f"Found {n_pos} positions.")

    xs_all   = np.array([p[0] for p in scan_points], dtype=float)

    # ── Extract per-position medians ──────────────────────────
    medL = np.full(n_pos, np.nan)
    medR = np.full(n_pos, np.nan)
    errL = np.full(n_pos, np.nan)
    errR = np.full(n_pos, np.nan)
    nL   = np.zeros(n_pos, dtype=int)
    nR   = np.zeros(n_pos, dtype=int)
    n_ev = np.zeros(n_pos, dtype=int)
    tL_arrays = []
    tR_arrays = []

    for i, (pos_mm, rpath) in enumerate(scan_points):
        tL, tR, n_total = first_photon_medians(rpath)
        n_ev[i] = n_total
        nL[i] = len(tL); nR[i] = len(tR)
        tL_arrays.append(tL); tR_arrays.append(tR)
        if len(tL) > 5:
            medL[i], errL[i] = robust_median_err(tL)
        if len(tR) > 5:
            medR[i], errR[i] = robust_median_err(tR)
        print(f"  x={pos_mm:+5d} mm  nL={nL[i]:5d}  medL={medL[i]:7.4f} ns"
              f"  nR={nR[i]:5d}  medR={medR[i]:7.4f} ns")

    # ── Half-difference: (medL - medR)/2 = x/v_eff ───────────
    half_diff    = (medL - medR) / 2.0
    half_diff_err = np.sqrt(errL**2 + errR**2) / 2.0
    both_ok = np.isfinite(medL) & np.isfinite(medR)

    # ── Linear fits ───────────────────────────────────────────
    # Convert x from mm to cm for v in cm/ns
    xs_cm = xs_all / 10.0

    # Fit t_L vs x (all positions with END_L data)
    popt_L, perr_L = fit_line(xs_cm, medL, errL)
    popt_R, perr_R = fit_line(xs_cm, medR, errR)
    # Fit (t_L-t_R)/2 vs x — best for v_eff (offset-free)
    popt_half, perr_half = fit_line(xs_cm[both_ok], half_diff[both_ok],
                                    half_diff_err[both_ok])

    def veff_from_slope(slope, slope_err):
        """slope in ns/cm → v = 1/slope cm/ns."""
        if slope is None or slope == 0:
            return float("nan"), float("nan")
        v = 1.0 / abs(slope)
        v_err = v * (slope_err / abs(slope))
        return v, v_err

    vL, vL_err   = veff_from_slope(popt_L[0]    if popt_L    is not None else None,
                                   perr_L[0]    if perr_L    is not None else 0.0)
    vR, vR_err   = veff_from_slope(popt_R[0]    if popt_R    is not None else None,
                                   perr_R[0]    if perr_R    is not None else 0.0)
    # half_diff slope sign is positive (x increases → half_diff increases)
    vH, vH_err   = veff_from_slope(popt_half[0] if popt_half is not None else None,
                                   perr_half[0] if perr_half is not None else 0.0)

    # Best estimate: weighted average of the three
    vals = [v for v in [vL, vR, vH] if np.isfinite(v)]
    errs_v = [e for v, e in [(vL, vL_err), (vR, vR_err), (vH, vH_err)]
              if np.isfinite(v) and e > 0]
    if errs_v:
        weights = [1.0/e**2 for e in errs_v]
        v_best = sum(w*v for w, v in zip(weights, vals)) / sum(weights)
        v_best_err = 1.0 / np.sqrt(sum(weights))
    else:
        v_best = np.mean(vals) if vals else float("nan")
        v_best_err = float("nan")

    print(f"\n{'='*50}")
    print(f"v_eff from END_L slope   : {vL:.3f} ± {vL_err:.3f} cm/ns  = {vL/C_CM_NS:.4f} ± {vL_err/C_CM_NS:.4f} c")
    print(f"v_eff from END_R slope   : {vR:.3f} ± {vR_err:.3f} cm/ns  = {vR/C_CM_NS:.4f} ± {vR_err/C_CM_NS:.4f} c")
    print(f"v_eff from (t_L-t_R)/2  : {vH:.3f} ± {vH_err:.3f} cm/ns  = {vH/C_CM_NS:.4f} ± {vH_err/C_CM_NS:.4f} c")
    print(f"v_eff best estimate      : {v_best:.3f} ± {v_best_err:.3f} cm/ns  = {v_best/C_CM_NS:.4f} ± {v_best_err/C_CM_NS:.4f} c")
    print(f"{'='*50}")

    # ── Plot ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.30)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    # Color by N hits (visual weight)
    norm_nL = np.clip(nL / 5000, 0.02, 1.0)
    norm_nR = np.clip(nR / 5000, 0.02, 1.0)

    # ── Panel 1: t_EndL and t_EndR vs x ──────────────────────
    ok_L = np.isfinite(medL); ok_R = np.isfinite(medR)
    sc = ax1.scatter(xs_all[ok_L], medL[ok_L], c=norm_nL[ok_L], cmap="Blues",
                     vmin=0, vmax=1, s=40, zorder=3, label="END_L (−X)", marker="o")
    ax1.scatter(xs_all[ok_R], medR[ok_R], c=norm_nR[ok_R], cmap="Reds",
                vmin=0, vmax=1, s=40, zorder=3, label="END_R (+X)", marker="s")

    x_line = np.linspace(-700, 700, 300)
    if popt_L is not None:
        ax1.plot(x_line, linear(x_line/10, *popt_L),
                 "b--", lw=1.5,
                 label=f"fit L: v={vL:.2f} cm/ns = {vL/C_CM_NS:.3f}c")
    if popt_R is not None:
        ax1.plot(x_line, linear(x_line/10, *popt_R),
                 "r--", lw=1.5,
                 label=f"fit R: v={vR:.2f} cm/ns = {vR/C_CM_NS:.3f}c")

    ax1.set_xlabel("x_gun (mm)"); ax1.set_ylabel("median t₁ (ns)")
    ax1.set_title("First-photon median time at each end")
    ax1.legend(fontsize=7); ax1.grid(alpha=0.3)

    # ── Panel 2: half-difference (t_L - t_R)/2 vs x ─────────
    ax2.errorbar(xs_all[both_ok], half_diff[both_ok],
                 yerr=half_diff_err[both_ok],
                 fmt="ko", ms=5, lw=1, capsize=3, label="(t_L − t_R)/2")
    if popt_half is not None:
        ax2.plot(x_line, linear(x_line/10, *popt_half),
                 "g-", lw=2,
                 label=f"fit: v={vH:.2f}±{vH_err:.2f} cm/ns\n= {vH/C_CM_NS:.4f}±{vH_err/C_CM_NS:.4f} c")
    ax2.axhline(0, lw=0.7, color="gray", ls=":")
    ax2.set_xlabel("x_gun (mm)"); ax2.set_ylabel("(t_L − t_R)/2  (ns)")
    ax2.set_title("Half time-difference = x/v_eff  (offset-free)")
    ax2.legend(fontsize=7); ax2.grid(alpha=0.3)

    # ── Panel 3: hit counts vs x ──────────────────────────────
    ax3.semilogy(xs_all, nL, "b-o", ms=4, lw=1.2, label="END_L hits (N events with ≥1 hit)")
    ax3.semilogy(xs_all, nR, "r-s", ms=4, lw=1.2, label="END_R hits")
    ax3.axhline(5000, lw=0.8, color="gray", ls="--", label="5000 events")
    ax3.set_xlabel("x_gun (mm)"); ax3.set_ylabel("events with ≥1 hit")
    ax3.set_title("END hit coverage vs position")
    ax3.legend(fontsize=7); ax3.grid(alpha=0.3, which="both")

    # ── Panel 4: residuals of (t_L-t_R)/2 from linear fit ────
    if popt_half is not None:
        res = half_diff[both_ok] - linear(xs_all[both_ok]/10, *popt_half)
        res_ps = res * 1000
        ax4.scatter(xs_all[both_ok], res_ps, c="k", s=30, zorder=3)
        ax4.axhline(0, color="gray", lw=0.8)
        ax4.fill_between(xs_all[both_ok], -50, 50, alpha=0.1, color="green",
                         label="±50 ps band")
        rms_ps = float(np.std(res_ps))
        ax4.set_xlabel("x_gun (mm)"); ax4.set_ylabel("residual (ps)")
        ax4.set_title(f"Residuals from linear fit  (RMS = {rms_ps:.1f} ps)")
        ax4.legend(fontsize=7); ax4.grid(alpha=0.3)
        ax4.set_ylim(-max(150, rms_ps*3), max(150, rms_ps*3))

    # Global title
    fig.suptitle(
        f"Effective light propagation velocity in EJ-200 bar\n"
        f"v_eff = {v_best:.3f} ± {v_best_err:.3f} cm/ns  "
        f"= {v_best/C_CM_NS:.4f} ± {v_best_err/C_CM_NS:.4f} c\n"
        f"(first photon at each end, 5000 muons/position, {n_pos} positions)",
        fontsize=11
    )

    out_png = OUT_DIR / "veff_curve.png"
    fig.savefig(str(out_png), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"PNG: {out_png}")

    # ── Single clean summary figure ───────────────────────────
    fig2, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar(xs_all[both_ok], half_diff[both_ok] * 1000,
                yerr=half_diff_err[both_ok] * 1000,
                fmt="ko", ms=5, lw=1.2, capsize=3, label="(t_L − t_R)/2  (data)")
    if popt_half is not None:
        ax.plot(x_line,
                linear(x_line/10, *popt_half) * 1000,
                "g-", lw=2.5, zorder=5,
                label=f"Linear fit:  v_eff = {vH:.3f} ± {vH_err:.3f} cm/ns = {vH/C_CM_NS:.4f}c")
    ax.set_xlabel("Gun position  x  (mm)", fontsize=12)
    ax.set_ylabel("(t_EndL − t_EndR) / 2   (ps)", fontsize=12)
    ax.set_title("Effective light propagation velocity in scintillator bar\n"
                 "EJ-200 bar, EndTop geometry, 5000 μ⁻/position, 24 threads",
                 fontsize=11)
    ax.legend(fontsize=10); ax.grid(alpha=0.3)
    fig2.tight_layout()
    out_png2 = OUT_DIR / "veff_halfdiff_clean.png"
    fig2.savefig(str(out_png2), dpi=150)
    plt.close(fig2)
    print(f"PNG: {out_png2}")

    # ── CSV ───────────────────────────────────────────────────
    csv_path = OUT_DIR / "veff_summary.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pos_mm", "n_events_total",
                    "n_EndL", "med_tL_ns", "err_tL_ns",
                    "n_EndR", "med_tR_ns", "err_tR_ns",
                    "half_diff_ns", "half_diff_err_ns"])
        for i in range(n_pos):
            w.writerow([
                int(xs_all[i]), n_ev[i],
                nL[i], f"{medL[i]:.6f}", f"{errL[i]:.6f}",
                nR[i], f"{medR[i]:.6f}", f"{errR[i]:.6f}",
                f"{half_diff[i]:.6f}", f"{half_diff_err[i]:.6f}",
            ])
    print(f"CSV: {csv_path}")

    # ── Results text ──────────────────────────────────────────
    txt_path = OUT_DIR / "veff_results.txt"
    with open(txt_path, "w") as fh:
        fh.write(f"Effective light propagation velocity — EJ-200 bar\n")
        fh.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"Run: {args.run_dir}\n\n")
        fh.write(f"Method: first-photon median time (N=1) at END_L and END_R\n")
        fh.write(f"Fit: t = (1/v_eff) * x + offset,  x in cm, t in ns\n\n")

        def fmt(v, e):
            return f"{v:.4f} ± {e:.4f}" if np.isfinite(v) else "nan"

        fh.write(f"END_L slope     : {popt_L[0]:.6f} ± {perr_L[0]:.6f} ns/cm\n"
                 if popt_L is not None else "END_L fit: failed\n")
        fh.write(f"  v_eff (L)     : {vL:.4f} ± {vL_err:.4f} cm/ns = {vL/C_CM_NS:.5f} ± {vL_err/C_CM_NS:.5f} c\n\n")
        fh.write(f"END_R slope     : {-popt_R[0]:.6f} ± {perr_R[0]:.6f} ns/cm\n"
                 if popt_R is not None else "END_R fit: failed\n")
        fh.write(f"  v_eff (R)     : {vR:.4f} ± {vR_err:.4f} cm/ns = {vR/C_CM_NS:.5f} ± {vR_err/C_CM_NS:.5f} c\n\n")
        fh.write(f"Half-diff slope : {popt_half[0]:.6f} ± {perr_half[0]:.6f} ns/cm\n"
                 if popt_half is not None else "Half-diff fit: failed\n")
        fh.write(f"  v_eff (Δt/2)  : {vH:.4f} ± {vH_err:.4f} cm/ns = {vH/C_CM_NS:.5f} ± {vH_err/C_CM_NS:.5f} c\n\n")
        fh.write(f"BEST ESTIMATE\n")
        fh.write(f"  v_eff         : {v_best:.4f} ± {v_best_err:.4f} cm/ns\n")
        fh.write(f"               = {v_best/C_CM_NS:.5f} ± {v_best_err/C_CM_NS:.5f} c\n")
        fh.write(f"               = {v_best:.4f} cm/ns  (c = {C_CM_NS:.4f} cm/ns)\n")
        if np.isfinite(v_best):
            n_eff = C_CM_NS / v_best
            fh.write(f"  n_eff         : {n_eff:.4f}  (= c/v_eff)\n")
    print(f"TXT: {txt_path}")
    print(f"\nDone. All outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
