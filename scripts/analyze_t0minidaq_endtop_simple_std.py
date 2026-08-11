#!/usr/bin/env python3
"""
analyze_t0minidaq_endtop_simple_std.py

Simple exploratory analysis of the EndTop scan.
No Gaussian fits, no chi2. Just std, mean, median, MAD, efficiency.

Metrics computed:
  END: tL_N, tR_N, END_AVG=(tL+tR)/2, END_SUM=tL+tR, END_DIFF=tR-tL
  TOP individual: tTOP_N per channel (global_id 16..85)
  TOP aggregate:  all TOP hits merged per event, t_N

Usage:
  python3 scripts/analyze_t0minidaq_endtop_simple_std.py \
      --run-dir <path> --out-dir <path> [options]
"""

import argparse
import csv
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import uproot

os.environ.setdefault("ROOT_BATCH", "1")
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

# ── constants ─────────────────────────────────────────────────
END_L_IDS = set(range(0, 8))     # global_ids 0-7
END_R_IDS = set(range(8, 16))    # global_ids 8-15
TOP_IDS   = set(range(16, 86))   # global_ids 16-85
N_TOP     = 70

DIAG_POSITIONS = [-690, -400, 0, 400, 690]   # for histograms
TOP_PROFILE_IDS = [16, 31, 51, 70, 85]        # representative TOP IDs


# ── statistics helpers ────────────────────────────────────────
def _mad_sigma(a):
    if len(a) == 0:
        return float("nan")
    med = np.median(a)
    return 1.4826 * float(np.median(np.abs(a - med)))


def stats(a):
    """Return dict of basic statistics for array a (in ns → ps for std/mad)."""
    n = len(a)
    if n == 0:
        return dict(n=0, mean=np.nan, std_ps=np.nan,
                    median=np.nan, mad_ps=np.nan,
                    p5=np.nan, p95=np.nan, rmin=np.nan, rmax=np.nan)
    return dict(
        n=n,
        mean=float(np.mean(a)),
        std_ps=float(np.std(a, ddof=1)) * 1000.0 if n > 1 else np.nan,
        median=float(np.median(a)),
        mad_ps=_mad_sigma(a) * 1000.0,
        p5=float(np.percentile(a, 5)),
        p95=float(np.percentile(a, 95)),
        rmin=float(np.min(a)),
        rmax=float(np.max(a)),
    )


def tN_per_event(evt_arr, t_arr, N, face_mask):
    """
    For events with ≥N hits in face_mask, return t[N-1] per event.
    Returns (tN array, n_total_events).
    """
    mask = face_mask
    e_ = evt_arr[mask]
    t_ = t_arr[mask]
    if len(e_) == 0:
        return np.array([]), 0
    order = np.argsort(e_, kind="stable")
    e_s = e_[order]; t_s = t_[order]
    splits = np.where(np.diff(e_s) != 0)[0] + 1
    groups = np.split(t_s, splits)
    result = np.array([np.sort(g)[N-1] for g in groups if len(g) >= N],
                      dtype=np.float64)
    return result, len(groups)


def tN_per_event_by_gid(evt_arr, t_arr, gid_arr, gid, N):
    """Per-channel (gid) tN computation."""
    mask = gid_arr == gid
    return tN_per_event(evt_arr, t_arr, N, mask)


def tN_per_event_set(evt_arr, t_arr, face_mask, N):
    """tN for a set of channels (face_mask already applied)."""
    return tN_per_event(evt_arr, t_arr, N, face_mask)


# ── load one file ─────────────────────────────────────────────
def load_file(path):
    with uproot.open(str(path)) as f:
        tree = f["sipm_hits"]
        arr = tree.arrays(
            ["event_id", "face_type", "global_id", "time_ns"],
            library="np",
        )
    evt = arr["event_id"]
    ft  = arr["face_type"]
    gid = arr["global_id"]
    t   = arr["time_ns"]
    n_total = len(np.unique(evt))
    return evt, ft, gid, t, n_total


# ── END metrics for one position ──────────────────────────────
def compute_end_metrics(evt, ft, gid, t, n_total, end_ns):
    """
    Compute END_L, END_R, AVG, SUM, DIFF for each N in end_ns.
    Returns list of row dicts.
    """
    mask_L = (ft == 0)
    mask_R = (ft == 1)
    rows = []

    for N in end_ns:
        tL_all, nL_ev = tN_per_event_set(evt, t, mask_L, N)
        tR_all, nR_ev = tN_per_event_set(evt, t, mask_R, N)

        # Pair up events that have BOTH endpoints
        # Strategy: build event-indexed dicts
        e_L = evt[mask_L]; t_L = t[mask_L]
        e_R = evt[mask_R]; t_R = t[mask_R]

        # tN per event using dict approach
        def tN_dict(e_arr, t_arr):
            if len(e_arr) == 0:
                return {}
            o = np.argsort(e_arr, kind="stable")
            es = e_arr[o]; ts = t_arr[o]
            spl = np.where(np.diff(es) != 0)[0] + 1
            d = {}
            for eg, tg in zip(np.split(es, spl), np.split(ts, spl)):
                tg_sorted = np.sort(tg)
                if len(tg_sorted) >= N:
                    d[eg[0]] = tg_sorted[N-1]
            return d

        dL = tN_dict(e_L, t_L)
        dR = tN_dict(e_R, t_R)

        # Find events present in BOTH
        common = set(dL.keys()) & set(dR.keys())
        n_valid = len(common)
        eff = n_valid / n_total if n_total > 0 else 0.0

        if n_valid == 0:
            row = {"N": N, "n_events_total": n_total, "n_valid": 0,
                   "efficiency": 0.0}
            for k in ["tL","tR","avg","diff","dsum"]:
                for m in ["mean","std_ps","median","mad_ps","p5","p95","rmin","rmax"]:
                    row[f"{k}_{m}"] = np.nan
            rows.append(row)
            continue

        common_sorted = sorted(common)
        tL_v = np.array([dL[e] for e in common_sorted])
        tR_v = np.array([dR[e] for e in common_sorted])
        avg  = (tL_v + tR_v) / 2.0
        diff = tR_v - tL_v
        dsum = tL_v + tR_v

        s_L   = stats(tL_v)
        s_R   = stats(tR_v)
        s_avg = stats(avg)
        s_dif = stats(diff)
        s_sum = stats(dsum)

        row = {"N": N, "n_events_total": n_total, "n_valid": n_valid, "efficiency": eff}
        for label, s in [("tL", s_L), ("tR", s_R),
                          ("avg", s_avg), ("diff", s_dif), ("dsum", s_sum)]:
            for m in ["mean","std_ps","median","mad_ps","p5","p95","rmin","rmax"]:
                row[f"{label}_{m}"] = s[m]
        row["_tL"] = tL_v
        row["_tR"] = tR_v
        row["_avg"] = avg
        row["_diff"] = diff
        rows.append(row)
    return rows


# ── TOP aggregate for one position ────────────────────────────
def compute_top_aggregate(evt, ft, t, n_total, top_ns):
    mask_top = (ft == 2)
    rows = []
    for N in top_ns:
        tN_vals, _ = tN_per_event_set(evt, t, mask_top, N)
        n_valid = len(tN_vals)
        eff = n_valid / n_total if n_total > 0 else 0.0
        s = stats(tN_vals) if n_valid > 0 else stats(np.array([]))
        rows.append({"N": N, "observable": f"TOP_ALL_N{N}",
                     "n_events_total": n_total, "n_valid": n_valid, "efficiency": eff,
                     **{k: s[k] for k in ["mean","std_ps","median","mad_ps","p5","p95","rmin","rmax"]},
                     "_tN": tN_vals})
    return rows


# ── TOP individual for one position ───────────────────────────
def compute_top_individual(evt, gid, t, n_total, top_n_list):
    rows = []
    unique_gids = np.array(sorted(g for g in np.unique(gid) if g in TOP_IDS))
    for g in unique_gids:
        mask_g = gid == g
        for N in top_n_list:
            tN_vals, n_ev_g = tN_per_event_set(evt, t, mask_g, N)
            n_valid = len(tN_vals)
            eff = n_valid / n_total if n_total > 0 else 0.0
            s = stats(tN_vals) if n_valid > 0 else stats(np.array([]))
            rows.append({"top_id": int(g), "N": N,
                         "n_events_total": n_total, "n_valid": n_valid,
                         "efficiency": eff,
                         **{k: s[k] for k in
                            ["mean","std_ps","median","mad_ps","p5","p95","rmin","rmax"]}})
    return rows


# ── Histogram helper ──────────────────────────────────────────
def make_hist_png(ax, data_ns, title, xlabel, color, bin_ps=10):
    if len(data_ns) < 5:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title(title, fontsize=9)
        return
    bin_ns = bin_ps / 1000.0
    lo = np.percentile(data_ns, 0.5); hi = np.percentile(data_ns, 99.5)
    lo -= 2*bin_ns; hi += 2*bin_ns
    n_bins = max(10, int((hi-lo)/bin_ns))
    ax.hist(data_ns, bins=n_bins, range=(lo, hi), color=color, alpha=0.75,
            edgecolor="none")
    s = stats(data_ns)
    ax.axvline(s["mean"], color="k", lw=1.5, ls="--", label=f"mean={s['mean']:.3f} ns")
    ax.axvline(s["median"], color="gray", lw=1, ls=":", label=f"med={s['median']:.3f} ns")
    info = (f"n={s['n']}  std={s['std_ps']:.0f} ps\n"
            f"MAD={s['mad_ps']:.0f} ps")
    ax.text(0.97, 0.97, info, transform=ax.transAxes, ha="right", va="top",
            fontsize=7.5, family="monospace",
            bbox=dict(fc="white", alpha=0.8, boxstyle="round"))
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("counts", fontsize=9)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)


# ── Trend plot helper ─────────────────────────────────────────
def trend_plot(ax, xs, ys, label, color=None, marker="o", ms=5, lw=1.5,
               yerr=None, alpha=1.0):
    ok = np.isfinite(np.array(ys, dtype=float))
    xa = np.array(xs)[ok]; ya = np.array(ys, dtype=float)[ok]
    if len(xa) == 0:
        return
    kw = dict(marker=marker, ms=ms, lw=lw, label=label, alpha=alpha)
    if color: kw["color"] = color
    if yerr is not None:
        ye = np.array(yerr, dtype=float)[ok]
        ax.errorbar(xa, ya, yerr=ye, capsize=2, **kw)
    else:
        ax.plot(xa, ya, **kw)


# ── ROOT writers ──────────────────────────────────────────────
def write_root(out_root, all_end, all_top_agg, all_top_ind,
               diag_hists, bin_ps, positions):
    rfile = ROOT.TFile(str(out_root), "RECREATE")

    def safe_tge(name, xs_list, ys_list):
        """Write a TGraphErrors only if there are finite points."""
        xs = np.ascontiguousarray(xs_list, dtype=np.float64)
        ys = np.ascontiguousarray(ys_list, dtype=np.float64)
        ok = np.isfinite(ys)
        n = int(ok.sum())
        if n < 2:
            return None
        xf = np.ascontiguousarray(xs[ok])
        yf = np.ascontiguousarray(ys[ok])
        ez = np.zeros(n, dtype=np.float64)
        gr = ROOT.TGraphErrors(n, xf, yf, ez, ez)
        gr.SetName(name); gr.SetTitle(name)
        return gr

    # ── END graphs ────────────────────────────────────────────
    dir_end = rfile.mkdir("end_graphs")
    dir_end.cd()
    ns_end = sorted({r["N"] for r in all_end})
    for N in ns_end:
        rows = sorted([r for r in all_end if r["N"]==N], key=lambda r: r["x_mm"])
        xs_r = [r["x_mm"] for r in rows]
        for col, tag in [("avg_std_ps","END_AVG_std"),
                         ("diff_std_ps","END_DIFF_std"),
                         ("tL_std_ps","tL_std"),
                         ("tR_std_ps","tR_std"),
                         ("avg_mean","END_AVG_mean_ns"),
                         ("diff_mean","END_DIFF_mean_ns"),
                         ("efficiency","eff")]:
            gr = safe_tge(f"{tag}_N{N}", xs_r, [r.get(col, np.nan) for r in rows])
            if gr is not None:
                gr.Write()
    rfile.cd()

    # ── TOP aggregate graphs ──────────────────────────────────
    dir_top = rfile.mkdir("top_agg_graphs")
    dir_top.cd()
    obs_list = sorted({r["observable"] for r in all_top_agg})
    for obs in obs_list:
        rows = sorted([r for r in all_top_agg if r["observable"]==obs],
                      key=lambda r: r["x_mm"])
        xs_r = [r["x_mm"] for r in rows]
        for col, tag in [("std_ps","std"), ("efficiency","eff"),
                         ("mean","mean_ns"), ("median","median_ns")]:
            gr = safe_tge(f"{obs}_{tag}", xs_r, [r.get(col, np.nan) for r in rows])
            if gr is not None:
                gr.Write()
    rfile.cd()

    # ── Diagnostic histograms ─────────────────────────────────
    dir_h = rfile.mkdir("diag_histograms")
    dir_h.cd()
    bin_ns = bin_ps / 1000.0
    for key, data_ns in diag_hists.items():
        if len(data_ns) < 3:
            continue
        lo = float(np.percentile(data_ns, 0.5))
        hi = float(np.percentile(data_ns, 99.5))
        lo -= 3*bin_ns; hi += 3*bin_ns
        if hi <= lo:
            continue
        n_bins = max(5, int((hi - lo) / bin_ns))
        safe_key = key.replace("-","m").replace("+","p")
        h = ROOT.TH1D(safe_key, safe_key, n_bins, lo, hi)
        for v in data_ns:
            h.Fill(float(v))
        h.Write()
    rfile.cd()

    rfile.Close()
    print(f"  ROOT: {out_root}")


# ── PNG summary plots ─────────────────────────────────────────
def make_summary_pngs(out_png_dir, all_end, all_top_agg, all_top_ind,
                      diag_hists, bin_ps, positions):

    colors_N = {1:"steelblue", 2:"darkorange", 4:"seagreen", 8:"crimson", 20:"purple"}

    # ── end_std_avg_vs_x.png ────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    ns_end = sorted({r["N"] for r in all_end})
    for N in ns_end:
        rows = sorted([r for r in all_end if r["N"]==N], key=lambda r: r["x_mm"])
        trend_plot(ax, [r["x_mm"] for r in rows],
                   [r.get("avg_std_ps", np.nan) for r in rows],
                   f"N={N}", color=colors_N.get(N))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("std((t_L + t_R)/2)  (ps)")
    ax.set_title("END_AVG  std vs position  —  simple exploration")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(out_png_dir / "end_std_avg_vs_x.png"), dpi=120)
    plt.close(fig)

    # ── end_std_diff_vs_x.png ───────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    for N in ns_end:
        rows = sorted([r for r in all_end if r["N"]==N], key=lambda r: r["x_mm"])
        trend_plot(ax, [r["x_mm"] for r in rows],
                   [r.get("diff_std_ps", np.nan) for r in rows],
                   f"N={N}", color=colors_N.get(N))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("std(t_R − t_L)  (ps)")
    ax.set_title("END_DIFF  std vs position")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(out_png_dir / "end_std_diff_vs_x.png"), dpi=120)
    plt.close(fig)

    # ── end_efficiency_vs_x.png ──────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    for N in ns_end:
        rows = sorted([r for r in all_end if r["N"]==N], key=lambda r: r["x_mm"])
        trend_plot(ax, [r["x_mm"] for r in rows],
                   [r.get("efficiency", np.nan) for r in rows],
                   f"N={N}", color=colors_N.get(N))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("efficiency (events with both ends)")
    ax.set_title("END paired efficiency vs position")
    ax.set_ylim(-0.05, 1.1)
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(out_png_dir / "end_efficiency_vs_x.png"), dpi=120)
    plt.close(fig)

    # ── end_mean_diff_vs_x.png ───────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    for N in ns_end:
        rows = sorted([r for r in all_end if r["N"]==N], key=lambda r: r["x_mm"])
        trend_plot(ax, [r["x_mm"] for r in rows],
                   [r.get("diff_mean", np.nan) for r in rows],
                   f"N={N}", color=colors_N.get(N))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("mean(t_R − t_L)  (ns)")
    ax.set_title("END_DIFF  mean vs position  (position proxy)")
    ax.axhline(0, lw=0.8, ls=":", color="gray")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(out_png_dir / "end_mean_diff_vs_x.png"), dpi=120)
    plt.close(fig)

    # ── top_all_std_vs_x.png ─────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    obs_list = sorted({r["observable"] for r in all_top_agg})
    for obs in obs_list:
        rows = sorted([r for r in all_top_agg if r["observable"]==obs],
                      key=lambda r: r["x_mm"])
        N = int(obs.split("N")[-1])
        trend_plot(ax, [r["x_mm"] for r in rows],
                   [r.get("std_ps", np.nan) for r in rows],
                   obs, color=colors_N.get(N))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("std(t_N)  (ps)")
    ax.set_title("TOP_ALL  std vs position — simple exploration")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(out_png_dir / "top_all_std_vs_x.png"), dpi=120)
    plt.close(fig)

    # ── top_all_efficiency_vs_x.png ──────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    for obs in obs_list:
        rows = sorted([r for r in all_top_agg if r["observable"]==obs],
                      key=lambda r: r["x_mm"])
        N = int(obs.split("N")[-1])
        trend_plot(ax, [r["x_mm"] for r in rows],
                   [r.get("efficiency", np.nan) for r in rows],
                   obs, color=colors_N.get(N))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("efficiency")
    ax.set_title("TOP_ALL efficiency vs position")
    ax.set_ylim(-0.05, 1.1)
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(out_png_dir / "top_all_efficiency_vs_x.png"), dpi=120)
    plt.close(fig)

    # ── TOP individual heatmaps ───────────────────────────────
    if all_top_ind:
        pos_sorted = sorted(positions)
        rows_N1 = [r for r in all_top_ind if r["N"]==1]
        top_ids_present = sorted({r["top_id"] for r in rows_N1})

        # Build matrices [n_pos x n_ch]
        n_pos = len(pos_sorted); n_ch = len(top_ids_present)
        mat_std = np.full((n_ch, n_pos), np.nan)
        mat_eff = np.full((n_ch, n_pos), np.nan)
        pid_map = {gid: i for i, gid in enumerate(top_ids_present)}
        pos_map = {p: i for i, p in enumerate(pos_sorted)}
        for r in rows_N1:
            ix = pos_map.get(r["x_mm"])
            iy = pid_map.get(r["top_id"])
            if ix is not None and iy is not None:
                mat_std[iy, ix] = r.get("std_ps", np.nan)
                mat_eff[iy, ix] = r["efficiency"]

        # std heatmap
        fig, ax = plt.subplots(figsize=(14, 8))
        im = ax.imshow(mat_std, aspect="auto", origin="lower",
                       extent=[pos_sorted[0]-10, pos_sorted[-1]+10,
                                top_ids_present[0]-0.5, top_ids_present[-1]+0.5],
                       cmap="plasma_r", vmin=0,
                       vmax=np.nanpercentile(mat_std[np.isfinite(mat_std)], 95)
                            if np.any(np.isfinite(mat_std)) else 1000)
        plt.colorbar(im, ax=ax, label="std (ps)")
        ax.set_xlabel("x_gun (mm)"); ax.set_ylabel("TOP global_id")
        ax.set_title("TOP individual  std(t₁)  [ps]  N=1  — heatmap")
        fig.tight_layout()
        fig.savefig(str(out_png_dir / "top_individual_std_heatmap_N1.png"), dpi=120)
        plt.close(fig)

        # efficiency heatmap
        fig, ax = plt.subplots(figsize=(14, 8))
        im = ax.imshow(mat_eff, aspect="auto", origin="lower",
                       extent=[pos_sorted[0]-10, pos_sorted[-1]+10,
                                top_ids_present[0]-0.5, top_ids_present[-1]+0.5],
                       cmap="viridis", vmin=0, vmax=1)
        plt.colorbar(im, ax=ax, label="efficiency")
        ax.set_xlabel("x_gun (mm)"); ax.set_ylabel("TOP global_id")
        ax.set_title("TOP individual  efficiency  N=1  — heatmap")
        fig.tight_layout()
        fig.savefig(str(out_png_dir / "top_individual_efficiency_heatmap_N1.png"), dpi=120)
        plt.close(fig)

        # Profile: selected TOP IDs
        sel_ids = [g for g in TOP_PROFILE_IDS if g in pid_map]
        if not sel_ids:
            # fallback: pick 5 spread evenly
            step = max(1, n_ch//5)
            sel_ids = top_ids_present[::step][:5]
        fig, axes = plt.subplots(len(sel_ids), 2, figsize=(12, 3*len(sel_ids)),
                                  sharex=True)
        if len(sel_ids) == 1:
            axes = [axes]
        for i, gid in enumerate(sel_ids):
            rows_g = [r for r in rows_N1 if r["top_id"]==gid]
            rows_g.sort(key=lambda r: r["x_mm"])
            xs_g = [r["x_mm"] for r in rows_g]
            axes[i][0].plot(xs_g, [r.get("std_ps", np.nan) for r in rows_g],
                            "o-", ms=4, lw=1.3, color="steelblue")
            axes[i][0].set_ylabel(f"std (ps)\ngid={gid}", fontsize=8)
            axes[i][0].grid(alpha=0.3)
            axes[i][1].plot(xs_g, [r.get("efficiency", np.nan) for r in rows_g],
                            "s-", ms=4, lw=1.3, color="seagreen")
            axes[i][1].set_ylabel("efficiency", fontsize=8)
            axes[i][1].set_ylim(-0.05, 1.1); axes[i][1].grid(alpha=0.3)
        axes[-1][0].set_xlabel("x_gun (mm)")
        axes[-1][1].set_xlabel("x_gun (mm)")
        fig.suptitle("TOP individual std & efficiency profiles  (N=1)", fontsize=11)
        fig.tight_layout()
        fig.savefig(str(out_png_dir / "top_individual_std_profiles_selected_ids.png"),
                    dpi=120)
        plt.close(fig)

    # ── Diagnostic histograms ─────────────────────────────────
    hist_dir = out_png_dir / "hists"
    hist_dir.mkdir(exist_ok=True)

    # Group diag_hists by prefix
    from collections import defaultdict
    by_key = defaultdict(dict)
    for key, data in diag_hists.items():
        parts = key.rsplit("_x", 1)
        if len(parts) == 2:
            prefix = parts[0]
            try:
                xval = int(parts[1])
            except ValueError:
                continue
            by_key[prefix][xval] = data

    for prefix, pos_dict in by_key.items():
        pos_list = sorted(pos_dict.keys())
        n = len(pos_list)
        cols = min(3, n); rows_p = (n + cols - 1) // cols
        fig, axes = plt.subplots(rows_p, cols, figsize=(5*cols, 4*rows_p),
                                  squeeze=False)
        for idx, pos in enumerate(pos_list):
            r = idx // cols; c = idx % cols
            ax = axes[r][c]
            color = "steelblue" if "avg" in prefix else \
                    "darkorange" if "diff" in prefix else "seagreen"
            make_hist_png(ax, pos_dict[pos],
                          f"{prefix}  x={pos:+d} mm",
                          "time (ns)", color)
        # Hide unused
        for idx in range(n, rows_p*cols):
            axes[idx//cols][idx%cols].set_visible(False)
        fig.suptitle(f"{prefix}  — diagnostic histograms", fontsize=11)
        fig.tight_layout()
        safe = prefix.replace("/","_")
        fig.savefig(str(hist_dir / f"{safe}.png"), dpi=110)
        plt.close(fig)


# ── CSV writers ───────────────────────────────────────────────
def write_csvs(out_csv_dir, all_end, all_top_agg, all_top_ind):
    # END metrics
    end_cols = ["x_mm","N","n_events_total","n_valid","efficiency",
                "tL_mean","tL_std_ps","tL_median","tL_mad_ps",
                "tR_mean","tR_std_ps","tR_median","tR_mad_ps",
                "avg_mean","avg_std_ps","avg_median","avg_mad_ps",
                "diff_mean","diff_std_ps","diff_median","diff_mad_ps",
                "dsum_mean","dsum_std_ps","dsum_median","dsum_mad_ps"]
    with open(out_csv_dir / "end_simple_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=end_cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(all_end, key=lambda r: (r["x_mm"], r["N"])):
            w.writerow({k: (f"{r[k]:.6g}" if isinstance(r.get(k), float) else r.get(k,""))
                        for k in end_cols})

    # TOP aggregate
    top_agg_cols = ["x_mm","observable","N","n_events_total","n_valid","efficiency",
                    "mean","std_ps","median","mad_ps","p5","p95","rmin","rmax"]
    with open(out_csv_dir / "top_aggregate_simple_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=top_agg_cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(all_top_agg, key=lambda r: (r["x_mm"], r["N"])):
            w.writerow({k: (f"{r[k]:.6g}" if isinstance(r.get(k), float) else r.get(k,""))
                        for k in top_agg_cols})

    # TOP individual
    top_ind_cols = ["x_mm","top_id","N","n_events_total","n_valid","efficiency",
                    "mean","std_ps","median","mad_ps","p5","p95","rmin","rmax"]
    with open(out_csv_dir / "top_individual_simple_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=top_ind_cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(all_top_ind, key=lambda r: (r["x_mm"], r["top_id"], r["N"])):
            w.writerow({k: (f"{r[k]:.6g}" if isinstance(r.get(k), float) else r.get(k,""))
                        for k in top_ind_cols})

    # Long format all
    long_cols = ["type","x_mm","group","N","n_events_total","n_valid","efficiency",
                 "mean","std_ps","median","mad_ps"]
    with open(out_csv_dir / "all_simple_metrics_long.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=long_cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(all_end, key=lambda r: (r["x_mm"], r["N"])):
            for grp, mk, sk in [("END_L","tL_mean","tL_std_ps"),
                                  ("END_R","tR_mean","tR_std_ps"),
                                  ("END_AVG","avg_mean","avg_std_ps"),
                                  ("END_DIFF","diff_mean","diff_std_ps"),
                                  ("END_SUM","dsum_mean","dsum_std_ps")]:
                w.writerow({"type":"END","x_mm":r["x_mm"],"group":grp,"N":r["N"],
                             "n_events_total":r["n_events_total"],"n_valid":r["n_valid"],
                             "efficiency":f"{r['efficiency']:.4f}",
                             "mean":f"{r.get(mk,float('nan')):.6g}",
                             "std_ps":f"{r.get(sk,float('nan')):.4g}",
                             "median":f"{r.get(mk.replace('mean','median'),float('nan')):.6g}",
                             "mad_ps":f"{r.get(sk.replace('std_ps','mad_ps'),float('nan')):.4g}"})
        for r in sorted(all_top_agg, key=lambda r: (r["x_mm"], r["N"])):
            w.writerow({"type":"TOP_AGG","x_mm":r["x_mm"],"group":r["observable"],
                        "N":r["N"],"n_events_total":r["n_events_total"],
                        "n_valid":r["n_valid"],
                        "efficiency":f"{r['efficiency']:.4f}",
                        "mean":f"{r.get('mean',float('nan')):.6g}",
                        "std_ps":f"{r.get('std_ps',float('nan')):.4g}",
                        "median":f"{r.get('median',float('nan')):.6g}",
                        "mad_ps":f"{r.get('mad_ps',float('nan')):.4g}"})


# ── Markdown report ───────────────────────────────────────────
def write_report(out_summary, all_end, all_top_agg, all_top_ind, run_dir, n_pos, positions):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = []
    lines.append("# EndTop Simple STD Analysis — Exploratory Report\n")
    lines.append(f"**Generated:** {now}  ")
    lines.append(f"**Run:** `{run_dir}`  ")
    lines.append(f"**Positions:** {n_pos}  ")
    lines.append(f"**Events/position:** 5000\n")

    lines.append("## Channel mapping\n")
    lines.append("| Group | global_id range | face_type |")
    lines.append("|---|---|---|")
    lines.append("| END_L | 0–7   | 0 |")
    lines.append("| END_R | 8–15  | 1 |")
    lines.append("| TOP   | 16–85 | 2 |\n")

    lines.append("## Definitions\n")
    lines.append("- **END_AVG** = (t_L + t_R) / 2")
    lines.append("- **END_SUM** = t_L + t_R")
    lines.append("- **END_DIFF** = t_R − t_L")
    lines.append("- **std** = sample standard deviation")
    lines.append("- **MAD** = 1.4826 × median|t − median(t)|")
    lines.append("- All std/MAD reported in **ps**; times in **ns**\n")

    lines.append("## END metrics summary (N=1, 2, 4, 8)\n")
    lines.append("| N | x_range | eff_min | eff_max | std_AVG range (ps) | std_DIFF range (ps) |")
    lines.append("|---|---|---|---|---|---|")
    for N in sorted({r["N"] for r in all_end}):
        rows = [r for r in all_end if r["N"]==N and r["n_valid"]>0]
        if not rows: continue
        effs = [r["efficiency"] for r in rows]
        std_avg = [r.get("avg_std_ps", np.nan) for r in rows]
        std_dif = [r.get("diff_std_ps", np.nan) for r in rows]
        std_avg_ok = [v for v in std_avg if np.isfinite(v)]
        std_dif_ok = [v for v in std_dif if np.isfinite(v)]
        lines.append(
            f"| {N} | {min(r['x_mm'] for r in rows)}..{max(r['x_mm'] for r in rows)} "
            f"| {min(effs):.3f} | {max(effs):.3f} "
            f"| {min(std_avg_ok):.0f}–{max(std_avg_ok):.0f} "
            f"| {min(std_dif_ok):.0f}–{max(std_dif_ok):.0f} |"
            if std_avg_ok and std_dif_ok else
            f"| {N} | — | — | — | — | — |"
        )

    lines.append("\n## TOP aggregate summary\n")
    lines.append("| Observable | eff_min | eff_max | std_min (ps) | std_max (ps) |")
    lines.append("|---|---|---|---|---|")
    for obs in sorted({r["observable"] for r in all_top_agg}):
        rows = [r for r in all_top_agg if r["observable"]==obs and r["n_valid"]>0]
        if not rows: continue
        effs = [r["efficiency"] for r in rows]
        stds = [r.get("std_ps", np.nan) for r in rows]
        stds_ok = [v for v in stds if np.isfinite(v)]
        lines.append(
            f"| `{obs}` | {min(effs):.3f} | {max(effs):.3f} "
            f"| {min(stds_ok):.0f} | {max(stds_ok):.0f} |"
            if stds_ok else f"| `{obs}` | — | — | — | — |"
        )

    lines.append("\n## Physical observations\n")

    # END_DIFF ~ linear with x
    rows_N1 = [r for r in all_end if r["N"]==1 and np.isfinite(r.get("diff_mean", np.nan))]
    rows_N1.sort(key=lambda r: r["x_mm"])
    if rows_N1:
        xs_n = np.array([r["x_mm"] for r in rows_N1])/10.0
        ys_n = np.array([r["diff_mean"] for r in rows_N1])
        ok = np.isfinite(ys_n)
        if ok.sum() > 2:
            p = np.polyfit(xs_n[ok], ys_n[ok], 1)
            # END_DIFF = t_R - t_L = (L/2-x)/v - (L/2+x)/v = -2x/v
            # so d(DIFF)/dx = -2/v → v_eff = -2/slope
            v_eff = -2.0/p[0] if p[0]!=0 else np.nan
            lines.append(
                f"- **END_DIFF (N=1) vs x:** approximately linear, "
                f"slope ≈ {p[0]:.4f} ns/cm  "
                f"→ v_eff = −2/slope ≈ {v_eff:.1f} cm/ns "
                f"= {v_eff/29.98:.3f} c"
            )

    # TOP std trend
    rows_top1 = sorted([r for r in all_top_agg if r["observable"]=="TOP_ALL_N1"],
                        key=lambda r: r["x_mm"])
    if rows_top1:
        stds = [r.get("std_ps", np.nan) for r in rows_top1]
        ok_stds = [v for v in stds if np.isfinite(v)]
        if ok_stds:
            lines.append(
                f"- **TOP_ALL_N1 std:** ranges {min(ok_stds):.0f}–{max(ok_stds):.0f} ps "
                f"across all positions"
            )

    # END efficiency
    rows_eff_N1 = [r for r in all_end if r["N"]==1]
    if rows_eff_N1:
        effs = sorted(rows_eff_N1, key=lambda r: r["x_mm"])
        best_eff_pos = [r["x_mm"] for r in effs if r["efficiency"]>0.9]
        lines.append(
            f"- **END paired efficiency (N=1):** > 90% for positions {best_eff_pos[0]}..{best_eff_pos[-1]} mm"
            if best_eff_pos else "- **END efficiency:** < 90% at all positions"
        )

    lines.append("- **std_END_AVG vs std_END_DIFF:** the AVG std mixes t_0 offset "
                 "and position uncertainty; DIFF std is a position-independent spread "
                 "related to arrival time resolution.")
    lines.append("- This analysis uses raw std without Gaussian core fit. "
                 "For core-Gaussian resolution see `analysis_corefit/`.\n")

    md_path = out_summary / "SIMPLE_STD_ANALYSIS_REPORT.md"
    with open(md_path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"Report: {md_path}")


# ── main ─────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir",  required=True)
    parser.add_argument("--out-dir",  required=True)
    parser.add_argument("--max-files", type=int, default=999)
    parser.add_argument("--bin-ps",   type=float, default=10.0)
    parser.add_argument("--end-n",    default="1,2,4,8")
    parser.add_argument("--top-n",    default="1,4,8,20")
    parser.add_argument("--no-top-individual-hists", action="store_true")
    args = parser.parse_args()

    RUN_DIR = Path(args.run_dir)
    OUT_DIR = Path(args.out_dir)
    BIN_PS  = args.bin_ps
    END_NS  = [int(x) for x in args.end_n.split(",")]
    TOP_AGG_NS = [int(x) for x in args.top_n.split(",")]
    TOP_IND_NS = [1]  # per-channel: N=1 only (others can be sparse)

    for sub in ["csv","png","root","logs","summary","png/hists"]:
        (OUT_DIR / sub).mkdir(parents=True, exist_ok=True)

    # Discover positions
    tsv_path = RUN_DIR / "summary.tsv"
    scan_points = []
    with open(tsv_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pos_mm  = int(row["position_mm"])
            pos_tag = row["index"]
            rpath   = RUN_DIR / "outputs" / pos_tag / "photon_hits_run000.root"
            if rpath.exists():
                scan_points.append((pos_mm, rpath))
    scan_points.sort(key=lambda p: p[0])
    scan_points = scan_points[:args.max_files]
    n_pos = len(scan_points)
    positions = [p[0] for p in scan_points]
    print(f"Processing {n_pos} positions  (bin={BIN_PS}ps  END_N={END_NS}  TOP_AGG_N={TOP_AGG_NS})")

    all_end      = []
    all_top_agg  = []
    all_top_ind  = []
    diag_hists   = {}
    diag_pos_set = set(DIAG_POSITIONS)

    t0 = time.time()
    for idx, (pos_mm, rpath) in enumerate(scan_points):
        sys.stdout.write(f"\n[{idx+1:02d}/{n_pos}] x={pos_mm:+4d} mm  loading ... ")
        sys.stdout.flush()
        evt, ft, gid, t_arr, n_total = load_file(rpath)
        sys.stdout.write("done\n"); sys.stdout.flush()

        # END
        sys.stdout.write(f"  END metrics ... "); sys.stdout.flush()
        end_rows = compute_end_metrics(evt, ft, gid, t_arr, n_total, END_NS)
        for r in end_rows:
            r["x_mm"] = pos_mm
            # Collect diag histograms at selected positions
            if pos_mm in diag_pos_set:
                for tag, key in [("avg",f"end_avg_N{r['N']}_x{pos_mm}"),
                                  ("diff",f"end_diff_N{r['N']}_x{pos_mm}")]:
                    if f"_{tag}" in r and len(r.get(f"_{tag}",[]))>0:
                        diag_hists[key] = r[f"_{tag}"]
        all_end.extend(end_rows)
        sys.stdout.write("OK\n"); sys.stdout.flush()

        # TOP aggregate
        sys.stdout.write(f"  TOP aggregate ... "); sys.stdout.flush()
        top_agg_rows = compute_top_aggregate(evt, ft, t_arr, n_total, TOP_AGG_NS)
        for r in top_agg_rows:
            r["x_mm"] = pos_mm
            N = r["N"]
            if pos_mm in diag_pos_set and r["n_valid"]>0:
                diag_hists[f"top_all_N{N}_x{pos_mm}"] = r["_tN"]
        all_top_agg.extend(top_agg_rows)
        sys.stdout.write("OK\n"); sys.stdout.flush()

        # TOP individual
        if not args.no_top_individual_hists:
            sys.stdout.write(f"  TOP individual ({TOP_IND_NS}) ... "); sys.stdout.flush()
            top_ind_rows = compute_top_individual(evt, gid, t_arr, n_total, TOP_IND_NS)
            for r in top_ind_rows:
                r["x_mm"] = pos_mm
            all_top_ind.extend(top_ind_rows)
            sys.stdout.write(f"OK ({len(top_ind_rows)} rows)\n"); sys.stdout.flush()

    print(f"\nAll positions done in {time.time()-t0:.1f}s")

    # Write outputs
    print("\n── Writing CSV ──")
    write_csvs(OUT_DIR / "csv", all_end, all_top_agg, all_top_ind)
    print(f"  csv/end_simple_metrics.csv")
    print(f"  csv/top_aggregate_simple_metrics.csv")
    print(f"  csv/top_individual_simple_metrics.csv")
    print(f"  csv/all_simple_metrics_long.csv")

    print("\n── Writing PNGs ──")
    make_summary_pngs(OUT_DIR/"png", all_end, all_top_agg, all_top_ind,
                      diag_hists, BIN_PS, positions)
    for png in sorted((OUT_DIR/"png").glob("*.png")):
        print(f"  png/{png.name}")

    print("\n── Writing ROOT ──")
    out_root = OUT_DIR / "root" / "simple_std_histograms_and_graphs.root"
    write_root(out_root, all_end, all_top_agg, all_top_ind,
               diag_hists, BIN_PS, positions)

    print("\n── Writing report ──")
    write_report(OUT_DIR/"summary", all_end, all_top_agg, all_top_ind,
                 str(RUN_DIR), n_pos, positions)

    print(f"\n✓ Analysis complete. Outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
