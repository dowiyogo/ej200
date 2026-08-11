#!/usr/bin/env python3
"""
analyze_t0minidaq_endtop_corefit.py

Core-fit timing analysis for the EndTop scan.
Reads 31 ROOT files, computes t_N timing observables, fits Gaussian
to the core, runs bootstrap, and writes CSV + PNG + ROOT outputs.

Usage:
  python3 scripts/analyze_t0minidaq_endtop_corefit.py \
      --run-dir <path> --out-dir <path> [--bin-ps 10] [--n-bootstrap 200]
"""

import argparse
import csv
import os
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no X display
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import uproot
from scipy.optimize import curve_fit


def _mad_sigma(a):
    """MAD-based robust sigma: 1.4826 * median(|a - median(a)|)."""
    med = np.median(a)
    return 1.4826 * np.median(np.abs(a - med))

# ── suppress ROOT banner ──────────────────────────────────────
os.environ.setdefault("ROOT_BATCH", "1")
import ROOT  # noqa: E402
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

# ── constants ─────────────────────────────────────────────────
OBSERVABLES = [
    ("TOP_ALL", 4),
    ("TOP_ALL", 8),
    ("TOP_ALL", 20),
    ("END_ALL", 4),
    ("END_ALL", 8),
    ("END_L",   4),
    ("END_R",   4),
]
OBS_NAMES = [f"{g}_N{n}" for g, n in OBSERVABLES]

GROUP_FACE = {
    "TOP_ALL": {2},
    "END_ALL": {0, 1},
    "END_L":   {0},
    "END_R":   {1},
}

CSV_COLS = [
    "obs", "pos_mm", "n_events_total", "n_events_valid", "efficiency",
    "mean_ns", "std_ps", "median_ns", "mad_sigma_ps",
    "fit_mu_ns", "fit_sigma_ps", "fit_sigma_err_ps",
    "chi2_ndf", "chi2", "ndf",
    "bootstrap_sigma_ps", "bootstrap_sigma_err_ps", "n_bootstrap_success",
    "quality_flag",
]


# ── gaussian ─────────────────────────────────────────────────
def gaussian(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ── core fit ──────────────────────────────────────────────────
def fit_gaussian_core(times_ns, bin_width_ns, verbose=False):
    """
    Fit a Gaussian to the core of a t_N distribution.
    Returns dict with fit results and quality flag.
    """
    if len(times_ns) < 10:
        return _empty_fit("LOW_STATS")

    med = np.median(times_ns)
    mad = _mad_sigma(times_ns)
    if mad <= 0:
        mad = np.std(times_ns)
    if mad <= 0:
        return _empty_fit("FIT_FAILED")

    # Robust histogram range: 0.5%–99.5% percentiles
    p_lo = np.percentile(times_ns, 0.5)
    p_hi = np.percentile(times_ns, 99.5)
    lo = p_lo - 0.5 * mad
    hi = p_hi + 0.5 * mad
    if hi - lo < bin_width_ns:
        hi = lo + 20 * bin_width_ns

    n_bins = max(10, int(round((hi - lo) / bin_width_ns)))
    n_bins = min(n_bins, 2000)
    counts, edges = np.histogram(times_ns, bins=n_bins, range=(lo, hi))
    centers = 0.5 * (edges[:-1] + edges[1:])

    # peak bin
    peak_idx = np.argmax(counts)
    t_peak = centers[peak_idx]
    A0 = float(counts[peak_idx])

    # fit window: ±2σ_MAD around peak
    win_lo = t_peak - 2.0 * mad
    win_hi = t_peak + 2.0 * mad
    mask = (centers >= win_lo) & (centers <= win_hi) & (counts > 0)
    if mask.sum() < 3:
        mask = (centers >= t_peak - 5 * bin_width_ns) & (centers <= t_peak + 5 * bin_width_ns) & (counts > 0)
    if mask.sum() < 3:
        return _empty_fit("FIT_FAILED")

    x_fit = centers[mask]
    y_fit = counts[mask].astype(float)
    yerr = np.sqrt(np.maximum(y_fit, 1.0))

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(
                gaussian, x_fit, y_fit,
                p0=[A0, t_peak, mad],
                sigma=yerr,
                absolute_sigma=True,
                bounds=([0, lo, 1e-5], [np.inf, hi, 5 * mad + 1e-9]),
                maxfev=10000,
            )
        perr = np.sqrt(np.diag(pcov))
    except Exception:
        return _empty_fit("FIT_FAILED")

    A_fit, mu_fit, sigma_fit = popt
    sigma_err = perr[2]
    if sigma_fit <= 0 or not np.isfinite(sigma_fit) or not np.isfinite(mu_fit):
        return _empty_fit("FIT_FAILED")

    # chi2 in fit window
    y_pred = gaussian(x_fit, *popt)
    residuals = (y_fit - y_pred) / yerr
    chi2_val = float(np.sum(residuals ** 2))
    ndf = max(1, len(x_fit) - 3)
    chi2_ndf = chi2_val / ndf

    quality = "OK"
    if len(times_ns) < 100:
        quality = "LOW_STATS"
    elif chi2_ndf >= 3.0:
        quality = "NON_GAUSSIAN_CORE"

    return {
        "mean_ns": float(np.mean(times_ns)),
        "std_ps": float(np.std(times_ns)) * 1000.0,
        "median_ns": float(med),
        "mad_sigma_ps": float(mad) * 1000.0,
        "fit_mu_ns": float(mu_fit),
        "fit_sigma_ps": float(sigma_fit) * 1000.0,
        "fit_sigma_err_ps": float(sigma_err) * 1000.0,
        "chi2": chi2_val,
        "ndf": ndf,
        "chi2_ndf": chi2_ndf,
        "quality_flag": quality,
        # internals for plotting
        "_counts": counts, "_edges": edges, "_centers": centers,
        "_popt": popt, "_mask_fit": mask,
        "_win_lo": win_lo, "_win_hi": win_hi,
        "_mad": mad, "_t_peak": t_peak,
    }


def _empty_fit(flag):
    return {k: float("nan") for k in [
        "mean_ns", "std_ps", "median_ns", "mad_sigma_ps",
        "fit_mu_ns", "fit_sigma_ps", "fit_sigma_err_ps",
        "chi2", "ndf", "chi2_ndf",
    ]} | {
        "quality_flag": flag,
        "_counts": None, "_edges": None, "_centers": None,
        "_popt": None, "_mask_fit": None,
        "_win_lo": None, "_win_hi": None, "_mad": None, "_t_peak": None,
    }


# ── bootstrap ─────────────────────────────────────────────────
def bootstrap_sigma(times_ns, bin_width_ns, n_boot):
    """Return (mean_sigma_ps, std_sigma_ps, n_success)."""
    if n_boot == 0 or len(times_ns) < 10:
        return float("nan"), float("nan"), 0

    rng = np.random.default_rng(seed=42)
    sigmas = []
    for _ in range(n_boot):
        sample = rng.choice(times_ns, size=len(times_ns), replace=True)
        res = fit_gaussian_core(sample, bin_width_ns)
        if res["quality_flag"] not in ("FIT_FAILED", "LOW_STATS") \
                and np.isfinite(res["fit_sigma_ps"]) and res["fit_sigma_ps"] > 0:
            sigmas.append(res["fit_sigma_ps"])

    if not sigmas:
        return float("nan"), float("nan"), 0
    return float(np.mean(sigmas)), float(np.std(sigmas)), len(sigmas)


# ── load one ROOT file ────────────────────────────────────────
def load_times(root_path):
    """
    Returns dict: group_name -> np.array of (event_id, time_ns) per hit.
    Uses uproot for speed.
    """
    with uproot.open(root_path) as f:
        tree = f["sipm_hits"]
        arr = tree.arrays(["event_id", "face_type", "time_ns"], library="np")

    evt_ids   = arr["event_id"]
    face_types = arr["face_type"]
    times     = arr["time_ns"]
    return evt_ids, face_types, times


# ── compute t_N for all events ────────────────────────────────
def compute_tN(evt_ids, face_types, times, faces_set, N):
    """
    For each event, collect hits in `faces_set`, sort by time, take t[N-1].
    Returns np.array of t_N values (one per valid event), and n_total_events.
    """
    # group by event_id
    unique_events = np.unique(evt_ids)
    n_total = len(unique_events)

    mask_face = np.isin(face_types, list(faces_set))
    e_sel = evt_ids[mask_face]
    t_sel = times[mask_face]

    # build event -> sorted times mapping efficiently
    order = np.argsort(e_sel, kind="stable")
    e_sorted = e_sel[order]
    t_sorted = t_sel[order]

    # find boundaries of each event group
    splits = np.where(np.diff(e_sorted) != 0)[0] + 1
    e_groups = np.split(e_sorted, splits)
    t_groups = np.split(t_sorted, splits)

    tN_list = []
    for tg in t_groups:
        if len(tg) >= N:
            ts = np.sort(tg)
            tN_list.append(ts[N - 1])

    return np.array(tN_list, dtype=np.float64), n_total


# ── plot one observable at one position ───────────────────────
def make_plot(obs_name, pos_mm, times_ns, result, n_total, out_png, bin_width_ns):
    fig, ax = plt.subplots(figsize=(9, 6))

    counts  = result.get("_counts")
    edges   = result.get("_edges")
    centers = result.get("_centers")
    popt    = result.get("_popt")
    mask    = result.get("_mask_fit")

    if counts is not None:
        ax.bar(centers, counts, width=(edges[1] - edges[0]),
               color="#5599cc", alpha=0.75, label="histogram")

        if popt is not None:
            x_curve = np.linspace(edges[0], edges[-1], 2000)
            ax.plot(x_curve, gaussian(x_curve, *popt),
                    "r-", lw=2, label="Gaussian fit")

            win_lo = result["_win_lo"]
            win_hi = result["_win_hi"]
            ax.axvspan(win_lo, win_hi, alpha=0.10, color="red", label="fit window")
    else:
        ax.text(0.5, 0.5, result["quality_flag"],
                ha="center", va="center", transform=ax.transAxes, fontsize=14)

    n_valid = len(times_ns) if times_ns is not None else 0
    eff = n_valid / n_total if n_total > 0 else 0.0
    flag = result.get("quality_flag", "?")

    sig_std  = result.get("std_ps", float("nan"))
    sig_mad  = result.get("mad_sigma_ps", float("nan"))
    sig_fit  = result.get("fit_sigma_ps", float("nan"))
    sig_err  = result.get("fit_sigma_err_ps", float("nan"))
    chi2_ndf = result.get("chi2_ndf", float("nan"))
    boot_s   = result.get("bootstrap_sigma_ps", float("nan"))
    boot_e   = result.get("bootstrap_sigma_err_ps", float("nan"))

    def fmt(v, dec=1):
        return f"{v:.{dec}f}" if np.isfinite(v) else "—"

    info = (
        f"x = {pos_mm:+.0f} mm     {obs_name}\n"
        f"n_total={n_total}   n_valid={n_valid}   eff={eff:.3f}\n"
        f"σ_std = {fmt(sig_std)} ps\n"
        f"σ_MAD = {fmt(sig_mad)} ps\n"
        f"σ_fit = {fmt(sig_fit)} ± {fmt(sig_err)} ps\n"
        f"boot  = {fmt(boot_s)} ± {fmt(boot_e)} ps\n"
        f"χ²/ndf = {fmt(chi2_ndf, 2)}\n"
        f"flag: {flag}"
    )
    ax.text(0.97, 0.97, info,
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, family="monospace",
            bbox=dict(boxstyle="round", fc="white", alpha=0.85))

    ax.set_xlabel("t_N  (ns)", fontsize=11)
    ax.set_ylabel("Counts", fontsize=11)
    ax.set_title(f"{obs_name}  |  x = {pos_mm:+.0f} mm", fontsize=12)
    ax.legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


# ── summary plot ──────────────────────────────────────────────
def summary_plot(positions, values, errs, title, ylabel, out_png,
                 obs_list=None, colors=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    if obs_list is None:
        obs_list = [None]
        values = [values]
        errs   = [errs]
    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(obs_list)))

    for obs, vals, err_arr, col in zip(obs_list, values, errs, colors):
        pos_arr = np.array(positions)
        val_arr = np.array(vals, dtype=float)
        err_arr2 = np.array(err_arr, dtype=float)
        ok = np.isfinite(val_arr)
        if ok.any():
            ax.errorbar(pos_arr[ok], val_arr[ok],
                        yerr=err_arr2[ok] if np.isfinite(err_arr2[ok]).any() else None,
                        marker="o", ms=5, lw=1.4, capsize=3,
                        label=obs, color=col)

    ax.set_xlabel("x (mm)", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


# ── write ROOT output ─────────────────────────────────────────
def write_root_output(all_results, out_root):
    """Write histograms, canvases, and TGraphErrors to ROOT file."""
    rfile = ROOT.TFile(str(out_root), "RECREATE")

    # TGraphErrors: sigma_fit vs x for each observable
    gr_sigma = {}
    gr_mu    = {}
    gr_eff   = {}
    gr_chi2  = {}

    for obs_name in OBS_NAMES:
        rows = [r for r in all_results if r["obs"] == obs_name]
        rows.sort(key=lambda r: r["pos_mm"])
        n = len(rows)
        xs    = np.array([r["pos_mm"]         for r in rows], dtype=np.float64)
        sigs  = np.array([r["fit_sigma_ps"]   for r in rows], dtype=np.float64)
        serrs = np.array([r["fit_sigma_err_ps"] for r in rows], dtype=np.float64)
        mus   = np.array([r["fit_mu_ns"]      for r in rows], dtype=np.float64)
        effs  = np.array([r["efficiency"]     for r in rows], dtype=np.float64)
        chi2s = np.array([r["chi2_ndf"]       for r in rows], dtype=np.float64)
        zeros = np.zeros(n, dtype=np.float64)

        # replace nan with 0 for ROOT (ROOT doesn't handle nan in graphs well)
        sigs_r  = np.where(np.isfinite(sigs),  sigs,  0.0)
        serrs_r = np.where(np.isfinite(serrs), serrs, 0.0)
        mus_r   = np.where(np.isfinite(mus),   mus,   0.0)
        effs_r  = np.where(np.isfinite(effs),  effs,  0.0)
        chi2s_r = np.where(np.isfinite(chi2s), chi2s, 0.0)

        gr = ROOT.TGraphErrors(n, xs, sigs_r, zeros, serrs_r)
        gr.SetName(f"gr_sigma_{obs_name}")
        gr.SetTitle(f"#sigma_{{fit}} vs x;x (mm);#sigma_{{fit}} (ps)")
        gr_sigma[obs_name] = gr

        grm = ROOT.TGraphErrors(n, xs, mus_r, zeros, zeros)
        grm.SetName(f"gr_mu_{obs_name}")
        grm.SetTitle(f"#mu_{{fit}} vs x;x (mm);#mu_{{fit}} (ns)")
        gr_mu[obs_name] = grm

        gre = ROOT.TGraph(n, xs, effs_r)
        gre.SetName(f"gr_eff_{obs_name}")
        gre.SetTitle(f"Efficiency vs x;x (mm);efficiency")
        gr_eff[obs_name] = gre

        grc = ROOT.TGraph(n, xs, chi2s_r)
        grc.SetName(f"gr_chi2_{obs_name}")
        grc.SetTitle(f"#chi^{{2}}/ndf vs x;x (mm);#chi^{{2}}/ndf")
        gr_chi2[obs_name] = grc

    # Write graphs
    rfile.mkdir("graphs").cd()
    for d in [gr_sigma, gr_mu, gr_eff, gr_chi2]:
        for g in d.values():
            g.Write()
    rfile.cd()

    print("  ROOT: graphs written")
    rfile.Close()
    print(f"  ROOT file: {out_root}")


# ── main ─────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir",     required=True)
    parser.add_argument("--out-dir",     required=True)
    parser.add_argument("--bin-ps",      type=float, default=10.0)
    parser.add_argument("--n-bootstrap", type=int,   default=200)
    args = parser.parse_args()

    RUN_DIR  = Path(args.run_dir)
    OUT_DIR  = Path(args.out_dir)
    BIN_NS   = args.bin_ps / 1000.0   # convert ps → ns
    N_BOOT   = args.n_bootstrap

    # make output dirs
    for d in ["csv", "png", "root", "logs", "summary"]:
        (OUT_DIR / d).mkdir(parents=True, exist_ok=True)

    # parse summary.tsv
    positions = []  # list of (pos_mm, root_path)
    tsv_path = RUN_DIR / "summary.tsv"
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos_mm = int(row["position_mm"])
            pos_tag = row["index"]  # e.g. "x-690mm"
            root_path = RUN_DIR / "outputs" / pos_tag / "photon_hits_run000.root"
            if root_path.exists():
                positions.append((pos_mm, root_path))
            else:
                print(f"WARNING: ROOT not found for {pos_tag}: {root_path}")

    positions.sort(key=lambda p: p[0])
    n_pos = len(positions)
    print(f"Found {n_pos} positions with ROOT files")

    # Accumulate all results
    all_results = []

    # Per-observable containers for summary plots
    summary_data = {obs: {"pos": [], "sigma": [], "sigma_err": [], "eff": [],
                          "chi2": [], "mu": [], "std": []} for obs in OBS_NAMES}

    t0_total = time.time()

    for p_idx, (pos_mm, root_path) in enumerate(positions):
        t0_pos = time.time()
        print(f"\n[{p_idx+1:02d}/{n_pos}] x={pos_mm:+4d} mm  loading {root_path.name}...")

        try:
            evt_ids, face_types, times_all = load_times(str(root_path))
        except Exception as e:
            print(f"  ERROR loading {root_path}: {e}")
            continue

        for obs_name in OBS_NAMES:
            grp_name, N = next((g, n) for g, n in OBSERVABLES
                               if f"{g}_N{n}" == obs_name)
            faces_set = GROUP_FACE[grp_name]

            sys.stdout.write(f"  {obs_name} ... ")
            sys.stdout.flush()

            t_vals, n_total = compute_tN(evt_ids, face_types, times_all, faces_set, N)
            n_valid = len(t_vals)
            eff = n_valid / n_total if n_total > 0 else 0.0

            if n_valid < 5:
                result = _empty_fit("LOW_EFFICIENCY")
                result["quality_flag"] = "LOW_EFFICIENCY"
                boot_s, boot_e, n_bs = float("nan"), float("nan"), 0
            else:
                result = fit_gaussian_core(t_vals, BIN_NS)

                # Bootstrap
                boot_s, boot_e, n_bs = bootstrap_sigma(t_vals, BIN_NS, N_BOOT)
                result["bootstrap_sigma_ps"] = boot_s
                result["bootstrap_sigma_err_ps"] = boot_e
                result["n_bootstrap_success"] = n_bs

            result["bootstrap_sigma_ps"] = boot_s
            result["bootstrap_sigma_err_ps"] = boot_e
            result.setdefault("n_bootstrap_success", n_bs)

            flag = result.get("quality_flag", "?")
            sig  = result.get("fit_sigma_ps", float("nan"))
            sys.stdout.write(f"σ_fit={sig:.1f} ps  [{flag}]\n")
            sys.stdout.flush()

            # PNG per position/observable
            png_name = f"pos_{pos_mm:+04d}mm_{obs_name}.png".replace("+", "p").replace("-", "m")
            out_png  = OUT_DIR / "png" / png_name
            make_plot(obs_name, pos_mm,
                      t_vals if n_valid >= 5 else None,
                      result, n_total, str(out_png), BIN_NS)

            row = {
                "obs":                    obs_name,
                "pos_mm":                 pos_mm,
                "n_events_total":         n_total,
                "n_events_valid":         n_valid,
                "efficiency":             eff,
                "mean_ns":                result.get("mean_ns", float("nan")),
                "std_ps":                 result.get("std_ps", float("nan")),
                "median_ns":              result.get("median_ns", float("nan")),
                "mad_sigma_ps":           result.get("mad_sigma_ps", float("nan")),
                "fit_mu_ns":              result.get("fit_mu_ns", float("nan")),
                "fit_sigma_ps":           result.get("fit_sigma_ps", float("nan")),
                "fit_sigma_err_ps":       result.get("fit_sigma_err_ps", float("nan")),
                "chi2_ndf":               result.get("chi2_ndf", float("nan")),
                "chi2":                   result.get("chi2", float("nan")),
                "ndf":                    result.get("ndf", float("nan")),
                "bootstrap_sigma_ps":     boot_s,
                "bootstrap_sigma_err_ps": boot_e,
                "n_bootstrap_success":    n_bs,
                "quality_flag":           flag,
            }
            all_results.append(row)

            # accumulate for summary
            sd = summary_data[obs_name]
            sd["pos"].append(pos_mm)
            sd["sigma"].append(sig)
            sd["sigma_err"].append(result.get("fit_sigma_err_ps", float("nan")))
            sd["eff"].append(eff)
            sd["chi2"].append(result.get("chi2_ndf", float("nan")))
            sd["mu"].append(result.get("fit_mu_ns", float("nan")))
            sd["std"].append(result.get("std_ps", float("nan")))

        print(f"  position done in {time.time()-t0_pos:.1f}s")

    print(f"\nAll positions done in {time.time()-t0_total:.1f}s")

    # ── Write CSV ─────────────────────────────────────────────
    csv_all = OUT_DIR / "csv" / "corefit_results_all.csv"
    with open(csv_all, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS)
        w.writeheader()
        for row in all_results:
            w.writerow({k: (f"{row[k]:.6g}" if isinstance(row[k], float) else row[k])
                        for k in CSV_COLS})

    for obs_name in OBS_NAMES:
        rows_obs = [r for r in all_results if r["obs"] == obs_name]
        rows_obs.sort(key=lambda r: r["pos_mm"])
        csv_obs = OUT_DIR / "csv" / f"corefit_results_{obs_name}.csv"
        with open(csv_obs, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=CSV_COLS)
            w.writeheader()
            for row in rows_obs:
                w.writerow({k: (f"{row[k]:.6g}" if isinstance(row[k], float) else row[k])
                            for k in CSV_COLS})

    print(f"CSV written: {csv_all}")

    # ── Summary PNGs ─────────────────────────────────────────
    colors = plt.cm.tab10(np.linspace(0, 1, len(OBS_NAMES)))

    # sigma_fit vs x (one PNG per observable + combined)
    for i, obs_name in enumerate(OBS_NAMES):
        sd = summary_data[obs_name]
        summary_plot(
            sd["pos"], sd["sigma"], sd["sigma_err"],
            f"σ_fit vs x  [{obs_name}]",
            "σ_fit (ps)",
            str(OUT_DIR / "png" / f"sigma_fit_vs_x_{obs_name}.png"),
        )

    # Combined sigma
    summary_plot(
        summary_data[OBS_NAMES[0]]["pos"],
        [summary_data[o]["sigma"]     for o in OBS_NAMES],
        [summary_data[o]["sigma_err"] for o in OBS_NAMES],
        "σ_fit vs x — all observables",
        "σ_fit (ps)",
        str(OUT_DIR / "png" / "summary_sigma_fit_vs_x.png"),
        obs_list=OBS_NAMES, colors=colors,
    )

    # efficiency
    summary_plot(
        summary_data[OBS_NAMES[0]]["pos"],
        [summary_data[o]["eff"] for o in OBS_NAMES],
        [np.zeros(len(summary_data[o]["eff"])) for o in OBS_NAMES],
        "Efficiency vs x",
        "efficiency",
        str(OUT_DIR / "png" / "summary_efficiency_vs_x.png"),
        obs_list=OBS_NAMES, colors=colors,
    )

    # chi2/ndf
    summary_plot(
        summary_data[OBS_NAMES[0]]["pos"],
        [summary_data[o]["chi2"] for o in OBS_NAMES],
        [np.zeros(len(summary_data[o]["chi2"])) for o in OBS_NAMES],
        "χ²/ndf vs x",
        "χ²/ndf",
        str(OUT_DIR / "png" / "summary_chi2_ndf_vs_x.png"),
        obs_list=OBS_NAMES, colors=colors,
    )

    # mu_fit vs x
    summary_plot(
        summary_data[OBS_NAMES[0]]["pos"],
        [summary_data[o]["mu"] for o in OBS_NAMES],
        [np.zeros(len(summary_data[o]["mu"])) for o in OBS_NAMES],
        "μ_fit vs x",
        "μ_fit (ns)",
        str(OUT_DIR / "png" / "summary_mu_fit_vs_x.png"),
        obs_list=OBS_NAMES, colors=colors,
    )

    # std vs fit comparison (TOP_ALL_N4 only as example)
    obs_cmp = "TOP_ALL_N4"
    sd = summary_data[obs_cmp]
    fig, ax = plt.subplots(figsize=(10, 5))
    pos_arr = np.array(sd["pos"])
    ok = np.isfinite(np.array(sd["sigma"], dtype=float))
    ax.plot(pos_arr[ok], np.array(sd["std"], dtype=float)[ok],
            "s--", ms=5, lw=1, color="gray", label="σ_std")
    ax.errorbar(pos_arr[ok], np.array(sd["sigma"], dtype=float)[ok],
                yerr=np.array(sd["sigma_err"], dtype=float)[ok],
                fmt="o-", ms=5, lw=1.5, color="steelblue",
                capsize=3, label="σ_fit (core)")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("ps")
    ax.set_title(f"σ comparison: std vs Gaussian core fit  [{obs_cmp}]")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(OUT_DIR / "png" / "summary_std_vs_fit_comparison.png"), dpi=120)
    plt.close(fig)

    print("Summary PNGs written")

    # ── ROOT output ───────────────────────────────────────────
    out_root = OUT_DIR / "root" / "corefit_histograms_and_canvases.root"
    write_root_output(all_results, out_root)

    # ── Markdown report ───────────────────────────────────────
    write_report(all_results, summary_data, OUT_DIR, RUN_DIR,
                 args.bin_ps, N_BOOT, n_pos)

    print(f"\nAnalysis complete. Outputs in: {OUT_DIR}")


# ── markdown report ───────────────────────────────────────────
def write_report(all_results, summary_data, OUT_DIR, RUN_DIR, bin_ps, n_boot, n_pos):
    from collections import defaultdict

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # best/worst sigma per observable
    report_lines = []
    report_lines.append(f"# EndTop Core-Fit Analysis Report\n")
    report_lines.append(f"**Generated:** {now}  ")
    report_lines.append(f"**Run dir:** `{RUN_DIR}`  ")
    report_lines.append(f"**Positions:** {n_pos}  ")
    report_lines.append(f"**Events/position:** 5000  ")
    report_lines.append(f"**Binning:** {bin_ps} ps  ")
    report_lines.append(f"**Bootstrap:** {n_boot} iterations  \n")

    report_lines.append("## Observables\n")
    for obs in OBS_NAMES:
        report_lines.append(f"- `{obs}`")
    report_lines.append("")

    report_lines.append("## Summary Table\n")
    report_lines.append("| Observable | OK | LOW_EFF | LOW_STATS | NON_GAUSS | FAIL | min_σ_ps | max_σ_ps | mean_σ_ps |")
    report_lines.append("|---|---|---|---|---|---|---|---|---|")

    for obs_name in OBS_NAMES:
        rows = [r for r in all_results if r["obs"] == obs_name]
        flags = defaultdict(int)
        for r in rows:
            flags[r["quality_flag"]] += 1
        sigs = [r["fit_sigma_ps"] for r in rows if np.isfinite(r["fit_sigma_ps"])]
        if sigs:
            report_lines.append(
                f"| `{obs_name}` | {flags['OK']} | {flags['LOW_EFFICIENCY']} | "
                f"{flags['LOW_STATS']} | {flags['NON_GAUSSIAN_CORE']} | "
                f"{flags['FIT_FAILED']} | {min(sigs):.1f} | {max(sigs):.1f} | {np.mean(sigs):.1f} |"
            )
        else:
            report_lines.append(f"| `{obs_name}` | — | — | — | — | {len(rows)} | — | — | — |")

    report_lines.append("")
    report_lines.append("## Per-observable highlights\n")

    for obs_name in OBS_NAMES:
        rows = [r for r in all_results if r["obs"] == obs_name
                and np.isfinite(r["fit_sigma_ps"]) and r["quality_flag"] == "OK"]
        if not rows:
            report_lines.append(f"### {obs_name}\nNo OK fits.\n")
            continue
        rows_s = sorted(rows, key=lambda r: r["fit_sigma_ps"])
        best  = rows_s[0]
        worst = rows_s[-1]
        report_lines.append(f"### {obs_name}")
        report_lines.append(f"- **Best  σ_fit** = {best['fit_sigma_ps']:.1f} ps  at x={best['pos_mm']:+d} mm")
        report_lines.append(f"- **Worst σ_fit** = {worst['fit_sigma_ps']:.1f} ps  at x={worst['pos_mm']:+d} mm")
        effvals = [r["efficiency"] for r in [r2 for r2 in all_results if r2["obs"] == obs_name]]
        report_lines.append(f"- **Efficiency range:** {min(effvals):.3f} – {max(effvals):.3f}")
        report_lines.append("")

    report_lines.append("## Physical commentary\n")
    # TOP_ALL_N4 sigma trend
    sd4 = summary_data.get("TOP_ALL_N4", {})
    sigs4 = np.array(sd4.get("sigma", []), dtype=float)
    pos4  = np.array(sd4.get("pos",   []))
    ok4   = np.isfinite(sigs4)
    if ok4.any():
        sigs_ok = sigs4[ok4]; pos_ok = pos4[ok4]
        center_mask = np.abs(pos_ok) < 200
        edge_mask   = np.abs(pos_ok) > 500
        if center_mask.any() and edge_mask.any():
            sig_center = np.mean(sigs_ok[center_mask])
            sig_edge   = np.mean(sigs_ok[edge_mask])
            report_lines.append(
                f"- **TOP_ALL_N4 σ_fit:** center avg={sig_center:.1f} ps, "
                f"edge avg={sig_edge:.1f} ps"
            )
            if sig_center < sig_edge:
                report_lines.append("  → resolution degrades toward the ends (expected: fewer top photons at edges).")
            else:
                report_lines.append("  → no strong x-dependence in top readout resolution.")

    # END_ALL efficiency
    sd_end = summary_data.get("END_ALL_N4", {})
    eff_end = np.array(sd_end.get("eff", []), dtype=float)
    if len(eff_end) > 0:
        report_lines.append(
            f"- **END_ALL_N4 efficiency:** {np.nanmin(eff_end):.3f} – {np.nanmax(eff_end):.3f}  "
            f"(typical: {np.nanmedian(eff_end):.3f})"
        )

    # chi2 issues
    bad_chi2 = [(r["obs"], r["pos_mm"], r["chi2_ndf"])
                for r in all_results if r.get("chi2_ndf", 0) >= 3.0
                and np.isfinite(r.get("chi2_ndf", float("nan")))]
    if bad_chi2:
        report_lines.append(f"- **χ²/ndf ≥ 3:** {len(bad_chi2)} fits:")
        for obs, pos, chi2 in bad_chi2[:10]:
            report_lines.append(f"  - `{obs}` x={pos:+d} mm  χ²/ndf={chi2:.2f}")
    else:
        report_lines.append("- **All core fits have χ²/ndf < 3** — Gaussian core model is good.")

    report_lines.append("\n---\n*End of report*\n")

    md_path = OUT_DIR / "summary" / "ANALYSIS_CORE_FIT_REPORT.md"
    with open(md_path, "w") as fh:
        fh.write("\n".join(report_lines))
    print(f"Markdown report: {md_path}")


if __name__ == "__main__":
    main()
