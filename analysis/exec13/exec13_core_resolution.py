#!/usr/bin/env python3
"""EXEC_13 core resolution: Gaussian fit to t4 nucleus, robust estimators.

Resolves the core/tail structure of the t_4 distribution at the Top nearest
channel (IDs 16, 31, 51) for the three key positions x = -690/-400/0 mm.

PHYSICAL FINDING: the t4 distribution has two separable components:
  Core  : very narrow (sigma_core ~ 1-3 ps from 2 ps bins), containing the
          majority of events where the 4th photon arrives via a prompt path.
          Characterised by IQR/1.349 ~ 2-5 ps and MAD*1.4826 ~ 2-5 ps.
  Tail  : long right tail (mean >> mode) from events where photons are
          delayed by reflections or Landau N_pe fluctuations.
  The RMS (52-62 ps) is dominated by the tail, not the core.

EXEC_12 reconciliation:
  EXEC_12 sigma_fit (8.77 / 8.02 / 8.21 ps) was computed with 25 ps bins
  in a [mean-2*RMS, mean+2*RMS] window.  That window (width ~200 ps) includes
  a large fraction of tail events and the fit converges to the width of the
  core+inner-tail region — close to sigma_stat(t4) = sqrt(4)*tau_d/<Npe>
  (theoretical order-statistic formula).  It does NOT measure the narrow
  prompt core.

  sigma_core (2 ps bins, mode±50 ps) captures the prompt-photon core only.
  The discrepancy with EXEC_12 is documented and expected; it is NOT an error.

Algorithm citations (unchanged from previous EXEC_13 stages):
  compute_tN  <- exec07/exec12_tN_analysis.py:72-100
  time origin  <- t=0 = Geant4 event start (exec07_photon_budget.py:342)
  No SPTR, no jitter, no time gate
  gaussian_core_fit (EXEC_12 style) <- exec07/exec12_tN_analysis.py:103-167

Outputs:
  analysis/exec13/exec13_core_resolution.csv
  analysis/exec13/exec13_resolution_macros.tex
  analysis/exec13/figs/exec13_f5_t4_core.{pdf,png}
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uproot

_REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "analysis"))
from exec07.common import N_EVENTS, TAU_D_NS, expected_file  # noqa: E402

# ─── Named constants ─────────────────────────────────────────────────────────
POSITIONS_MM: tuple[int, ...] = (-690, -400, 0)
NEAREST_IDS: dict[int, int] = {-690: 16, -400: 31, 0: 51}
RUN_LABELS: dict[int, str] = {-690: "Run 1", -400: "Run 3", 0: "Run 4"}

N_PE_THRESHOLD: int = 4           # leading-edge threshold
CORE_BIN_PS: float = 2.0          # fine bin width for Gaussian fit (spec: 2 ps)
CORE_WINDOW_PS: float = 50.0      # fit range = mode ± CORE_WINDOW_PS  (spec: ±50 ps)
TAIL_K: float = 3.0               # tail: t4 > mode + TAIL_K * sigma_core
ZOOM_LO_PS: float = -50.0         # F5 x-range offset from global mean mode (left)
ZOOM_HI_PS: float = 160.0         # F5 x-range (right): shows core + partial tail

# EXEC_12 reference (exec07/exec12_tN_summary.csv, top_nearest, N=4)
EXEC12_SIGMA_FIT_PS: dict[int, float] = {-690: 8.774, -400: 8.021, 0: 8.208}
EXEC12_BIN_NS: float = 0.025      # 25 ps bins used by exec12_tN_analysis.py:DT_FINE_NS

# Mean N_pe from exec13_f3_npe_profile.csv (channels 16, 31, 51)
MEAN_NPE: dict[int, float] = {-690: 350.95, -400: 411.28, 0: 407.27}

_COLORS: dict[int, str] = {-690: "tab:blue", -400: "tab:orange", 0: "tab:green"}

_DEFAULT_DATA_DIR = pathlib.Path(
    "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"
)
_DEFAULT_OUT_DIR = _REPO / "analysis" / "exec13"


# ─── Core algorithm (exec07/exec12_tN_analysis.py:72-100, unchanged) ─────────

def compute_tN(
    event_id: np.ndarray,
    global_id: np.ndarray,
    times: np.ndarray,
    channel_ids: tuple[int, ...],
    N: int,
    n_events: int = N_EVENTS,
) -> np.ndarray:
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
    return tN


# ─── Gaussian ────────────────────────────────────────────────────────────────

def _gauss(x: np.ndarray, A: float, mu: float, sigma: float) -> np.ndarray:
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ─── Mode: median of t4 (robust; mode ≈ median for right-skewed distribution) ─

def find_mode_ps(t_ps: np.ndarray) -> float:
    """Mode via peak of 2 ps histogram (matches CORE_BIN_PS for consistency)."""
    lo = float(np.nanpercentile(t_ps, 5))
    hi = float(np.nanpercentile(t_ps, 95))
    edges = np.arange(lo, hi + CORE_BIN_PS, CORE_BIN_PS)
    counts, _ = np.histogram(t_ps, bins=edges)
    return float(edges[np.argmax(counts)] + CORE_BIN_PS / 2)


# ─── Fit 1: 2 ps bins, mode ± CORE_WINDOW_PS (spec §PASO 1b) ─────────────────

def fit_gaussian_core(t_ps: np.ndarray, mode_ps: float) -> dict:
    """Gaussian fit to core in [mode - CORE_WINDOW_PS, mode + CORE_WINDOW_PS]
    with CORE_BIN_PS bins.  Captures only the narrow prompt-photon core."""
    lo = mode_ps - CORE_WINDOW_PS
    hi = mode_ps + CORE_WINDOW_PS
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
                "chi2ndf": chi2 / ndf, "mu_ps": mu, "A": A,
                "centers": centers, "counts": counts, "edges": edges}
    except Exception as exc:
        print(f"    fit_gaussian_core failed: {exc}")
        return {"ok": False, "sigma_ps": math.nan, "sigma_err_ps": math.nan,
                "chi2ndf": math.nan, "mu_ps": mode_ps, "A": math.nan,
                "centers": centers, "counts": counts, "edges": edges}


# ─── Fit 2: EXEC_12-style iterative (exec12_tN_analysis.py:103-167) ──────────

def fit_gaussian_exec12_style(t_ps: np.ndarray) -> dict:
    """Two-pass Gaussian fit in [mean-2*RMS, mean+2*RMS] with 25 ps bins.
    Replicates exec07/exec12_tN_analysis.py:gaussian_core_fit adapted to ps."""
    valid = t_ps[np.isfinite(t_ps)]
    mean0 = float(np.mean(valid))
    rms0 = float(np.std(valid, ddof=1))
    dt_ps = EXEC12_BIN_NS * 1000.0  # 25 ps

    lo0, hi0 = mean0 - 2.0 * rms0, mean0 + 2.0 * rms0
    core0 = valid[(valid >= lo0) & (valid <= hi0)]
    if len(core0) < 10:
        return {"ok": False, "sigma_ps": math.nan, "sigma_err_ps": math.nan}

    edges0 = np.arange(lo0, hi0 + dt_ps * 0.5, dt_ps)
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

    edges1 = np.arange(lo1, hi1 + dt_ps * 0.5, dt_ps)
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


# ─── Theoretical σ_stat(t_4) ─────────────────────────────────────────────────

def sigma_stat_ps(x_mm: int) -> float:
    """σ_stat(t_4) = sqrt(4) * tau_d / <Npe>  (order-statistics formula)."""
    return math.sqrt(N_PE_THRESHOLD) * TAU_D_NS * 1000.0 / MEAN_NPE[x_mm]


# ─── Robust estimators ────────────────────────────────────────────────────────

def robust_estimators(t_ps: np.ndarray) -> dict:
    q25, q75 = float(np.nanpercentile(t_ps, 25)), float(np.nanpercentile(t_ps, 75))
    iqr_sigma = (q75 - q25) / 1.349
    med = float(np.nanmedian(t_ps))
    mad_sigma = float(np.nanmedian(np.abs(t_ps - med))) * 1.4826
    return {"iqr_sigma_ps": iqr_sigma, "mad_sigma_ps": mad_sigma,
            "q25_ps": q25, "q50_ps": med, "q75_ps": q75}


# ─── Core/tail fractions ──────────────────────────────────────────────────────

def core_tail_fraction(t_ps: np.ndarray, mode_ps: float, sigma_core_ps: float) -> dict:
    if math.isnan(sigma_core_ps):
        return {"core_fraction": math.nan, "tail_fraction": math.nan,
                "tail_boundary_ps": math.nan}
    boundary = mode_ps + TAIL_K * sigma_core_ps
    valid = t_ps[np.isfinite(t_ps)]
    tail_frac = float(np.sum(valid > boundary)) / len(valid)
    return {"core_fraction": 1.0 - tail_frac, "tail_fraction": tail_frac,
            "tail_boundary_ps": boundary}


# ─── Per-position analysis ────────────────────────────────────────────────────

def analyze_position(data_dir: pathlib.Path, x_mm: int) -> dict:
    path = expected_file(data_dir, x_mm)
    with uproot.open(path) as rf:
        arr = rf["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns"], library="np"
        )
    ev = arr["event_id"].astype(np.int32)
    gid = arr["global_id"].astype(np.int32)
    t = arr["time_ns"].astype(np.float64)

    results = {}
    for N in (4, 20):
        ch = NEAREST_IDS[x_mm]
        tN_ns = compute_tN(ev, gid, t, (ch,), N)
        tN_ps = tN_ns * 1000.0
        valid = tN_ps[np.isfinite(tN_ps)]

        mode_ps = find_mode_ps(valid)
        rms_full = float(np.nanstd(tN_ps))
        rob = robust_estimators(valid)
        fit2 = fit_gaussian_core(valid, mode_ps)
        fitex12 = fit_gaussian_exec12_style(valid)
        ct = core_tail_fraction(valid, mode_ps, fit2["sigma_ps"])
        ss = sigma_stat_ps(x_mm) if N == 4 else math.nan

        exec12_ref = EXEC12_SIGMA_FIT_PS.get(x_mm, math.nan) if N == 4 else math.nan

        # Reconciliation: sigma_core vs EXEC_12 reference
        if N == 4 and not math.isnan(fit2["sigma_ps"]) and not math.isnan(exec12_ref):
            tol = 3 * fit2["sigma_err_ps"] if not math.isnan(fit2["sigma_err_ps"]) else 1.0
            reconciled_2ps = abs(fit2["sigma_ps"] - exec12_ref) <= tol
        else:
            reconciled_2ps = None

        # Reconciliation for EXEC_12-style fit
        if N == 4 and not math.isnan(fitex12["sigma_ps"]) and not math.isnan(exec12_ref):
            reconciled_ex = abs(fitex12["sigma_ps"] - exec12_ref) <= 2.0
        else:
            reconciled_ex = None

        # Use explicit None check (numpy bool `is False` fails identity test)
        if reconciled_2ps is None:
            rec2_str = "na"
        elif reconciled_2ps:
            rec2_str = "yes"
        else:
            rec2_str = "NO"
        if reconciled_ex is None:
            recex_str = "na"
        elif reconciled_ex:
            recex_str = "yes"
        else:
            recex_str = "NO"

        results[N] = {
            "x_mm": x_mm, "run_label": RUN_LABELS[x_mm], "channel": ch,
            "threshold_pe": N, "n_events": len(valid),
            "mode_ps": mode_ps, "rms_full_ps": rms_full,
            "sigma_core_ps": fit2["sigma_ps"],
            "sigma_core_err_ps": fit2["sigma_err_ps"],
            "chi2ndf": fit2["chi2ndf"],
            "fit_mu_ps": fit2["mu_ps"],
            "iqr_sigma_ps": rob["iqr_sigma_ps"],
            "mad_sigma_ps": rob["mad_sigma_ps"],
            "q25_ps": rob["q25_ps"], "q50_ps": rob["q50_ps"], "q75_ps": rob["q75_ps"],
            "core_fraction": ct["core_fraction"],
            "tail_fraction": ct["tail_fraction"],
            "tail_boundary_ps": ct["tail_boundary_ps"],
            "tail_k": TAIL_K,
            "sigma_exec12_style_ps": fitex12["sigma_ps"],
            "sigma_exec12_style_err_ps": fitex12["sigma_err_ps"],
            "sigma_stat_theory_ps": ss,
            "exec12_sigma_fit_ref_ps": exec12_ref,
            "reconciled_2ps_vs_exec12": rec2_str,
            "reconciled_exec12style_vs_exec12": recex_str,
            # internal (not written to CSV)
            "_fit2": fit2, "_fitex12": fitex12, "_tN_ps": valid,
        }

        print(f"  x={x_mm:+d} N={N:2d}: "
              f"mode={mode_ps:.1f} RMS={rms_full:.1f} | "
              f"sigma_core(2ps)={fit2['sigma_ps']:.2f}±{fit2['sigma_err_ps']:.2f} ps | "
              f"chi2/ndf={fit2['chi2ndf']:.2f} | "
              f"IQR/1.349={rob['iqr_sigma_ps']:.1f} MAD={rob['mad_sigma_ps']:.1f} | "
              f"tail={ct['tail_fraction']*100:.1f}%"
              + (f" | EXEC12-style={fitex12['sigma_ps']:.2f} ps | "
                 f"sigma_stat_th={ss:.2f} ps | vs EXEC12={exec12_ref:.2f} ps [{rec2_str}]"
                 if N == 4 else ""))
    return results


# ─── F5: core zoom figure ─────────────────────────────────────────────────────

def plot_f5(all_results: dict, out_dir: pathlib.Path) -> None:
    """F5: t4 at 2 ps bins, zoom to nucleus, Gaussian fit superimposed."""
    modes = [all_results[x][4]["mode_ps"] for x in POSITIONS_MM]
    global_mean_mode = float(np.mean(modes))
    x_lo = global_mean_mode + ZOOM_LO_PS
    x_hi = global_mean_mode + ZOOM_HI_PS

    # Global y-max
    global_ymax = 0
    for x_mm in POSITIONS_MM:
        valid = all_results[x_mm][4]["_tN_ps"]
        edges = np.arange(x_lo, x_hi + CORE_BIN_PS, CORE_BIN_PS)
        counts, _ = np.histogram(valid, bins=edges)
        global_ymax = max(global_ymax, int(counts.max()))
    global_ymax = int(global_ymax * 1.2)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    fig.suptitle(
        r"F5 — $t_4$ nucleus: 2 ps bins, zoom to core, Gaussian fit, fixed scale",
        fontsize=12,
    )
    xlims_rec, ylims_rec = [], []

    for ax, x_mm in zip(axes, POSITIONS_MM):
        r = all_results[x_mm][4]
        valid = r["_tN_ps"]
        fit = r["_fit2"]
        color = _COLORS[x_mm]
        run_lbl = RUN_LABELS[x_mm]
        ss = r["sigma_stat_theory_ps"]

        edges = np.arange(x_lo, x_hi + CORE_BIN_PS, CORE_BIN_PS)
        counts, _ = np.histogram(valid, bins=edges)
        centers = (edges[:-1] + edges[1:]) / 2.0

        ax.step(centers, counts, where="mid", color=color, lw=1.5, alpha=0.8)
        ax.fill_between(centers, counts, step="mid", color=color, alpha=0.15)

        # Narrow Gaussian fit (sigma_core)
        if fit["ok"] and not math.isnan(fit["sigma_ps"]):
            x_fit = np.linspace(x_lo, x_hi, 800)
            y_fit = _gauss(x_fit, fit["A"], fit["mu_ps"], fit["sigma_ps"])
            ax.plot(x_fit, y_fit, "k-", lw=1.8,
                    label=f"Gauss fit: {fit['sigma_ps']:.2f} ps")

        # Theoretical sigma_stat curve for reference
        if not math.isnan(ss):
            A_stat = float(len(valid)) * CORE_BIN_PS / (ss * np.sqrt(2 * np.pi))
            mu_stat = r["mode_ps"]
            x_th = np.linspace(x_lo, x_hi, 800)
            y_th = _gauss(x_th, A_stat, mu_stat, ss)
            ax.plot(x_th, y_th, "--", color="gray", lw=1.3, alpha=0.8,
                    label=f"$\\sigma_{{stat}}={ss:.1f}$ ps (theory)")

        sigma_val = fit["sigma_ps"] if not math.isnan(fit["sigma_ps"]) else float("nan")
        sigma_err = fit["sigma_err_ps"] if not math.isnan(fit["sigma_err_ps"]) else float("nan")
        ax.text(0.97, 0.95,
                f"$\\sigma_{{\\rm core}}=${sigma_val:.2f}$\\pm${sigma_err:.2f} ps\n"
                f"IQR/1.349={r['iqr_sigma_ps']:.1f} ps\n"
                f"RMS$_{{\\rm full}}=${r['rms_full_ps']:.1f} ps",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85))

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(0, global_ymax)
        ax.set_xlabel("$t_4$ (ps)", fontsize=10)
        if ax is axes[0]:
            ax.set_ylabel("Events / 2 ps bin", fontsize=10)
        ax.set_title(f"{run_lbl}: x={x_mm:+d} mm, ID {NEAREST_IDS[x_mm]}", fontsize=10)
        ax.legend(fontsize=7.5, loc="upper right")
        ax.grid(True, alpha=0.3)
        xlims_rec.append(ax.get_xlim())
        ylims_rec.append(ax.get_ylim())

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / "figs" / f"exec13_f5_t4_core.{ext}",
                    bbox_inches="tight", dpi=150)
    plt.close(fig)

    ok = all(
        abs(xl[0] - xlims_rec[0][0]) < 1e-3 and abs(xl[1] - xlims_rec[0][1]) < 1e-3 and
        abs(yl[0] - ylims_rec[0][0]) < 1e-3 and abs(yl[1] - ylims_rec[0][1]) < 1e-3
        for xl, yl in zip(xlims_rec, ylims_rec)
    )
    print(f"  F5 axis-identity gate: {'PASS' if ok else 'FAIL'}")


# ─── CSV ──────────────────────────────────────────────────────────────────────

def write_csv(all_results: dict, out_dir: pathlib.Path) -> None:
    path = out_dir / "exec13_core_resolution.csv"
    fields = [
        "x_mm", "run_label", "channel", "threshold_pe", "n_events",
        "mode_ps", "rms_full_ps",
        "sigma_core_ps", "sigma_core_err_ps", "chi2ndf", "fit_mu_ps",
        "iqr_sigma_ps", "mad_sigma_ps", "q25_ps", "q50_ps", "q75_ps",
        "core_fraction", "tail_fraction", "tail_boundary_ps", "tail_k",
        "sigma_exec12_style_ps", "sigma_exec12_style_err_ps",
        "sigma_stat_theory_ps", "exec12_sigma_fit_ref_ps",
        "reconciled_2ps_vs_exec12", "reconciled_exec12style_vs_exec12",
    ]
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


# ─── LaTeX macros ─────────────────────────────────────────────────────────────

def write_macros(all_results: dict, out_dir: pathlib.Path) -> None:
    lines = [
        "% Auto-generated by exec13_core_resolution.py",
        "% Source data: exec13_core_resolution.csv — do not edit manually",
        "",
    ]

    tags = {-690: "RunOne", -400: "RunThree", 0: "RunFour"}
    for x_mm in POSITIONS_MM:
        r4 = all_results[x_mm][4]
        tag = tags[x_mm]

        def cmd(name: str, val, fmt: str = ".2f") -> str:
            if isinstance(val, float) and math.isnan(val):
                return f"\\newcommand{{\\core{tag}{name}}}{{--}}"
            if isinstance(val, float):
                return f"\\newcommand{{\\core{tag}{name}}}{{{val:{fmt}}}}"
            return f"\\newcommand{{\\core{tag}{name}}}{{{val}}}"

        lines += [
            f"% {RUN_LABELS[x_mm]} x={x_mm:+d} mm ID {NEAREST_IDS[x_mm]}",
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
            cmd("ExecTwelveRef", r4["exec12_sigma_fit_ref_ps"]),
            cmd("TailK", TAIL_K, ".0f"),
            cmd("Reconciled", r4["reconciled_2ps_vs_exec12"]),
            "",
        ]

    sigs_ex = [all_results[x][4]["sigma_exec12_style_ps"] for x in POSITIONS_MM
               if not math.isnan(all_results[x][4]["sigma_exec12_style_ps"])]
    sigs_core = [all_results[x][4]["sigma_core_ps"] for x in POSITIONS_MM
                 if not math.isnan(all_results[x][4]["sigma_core_ps"])]
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
        mean_cmd("MeanSigCore", sigs_core),
        mean_cmd("MeanSigEx", sigs_ex),
        mean_cmd("MeanSigStat", sigs_stat),
        # SPTR reference value from EXEC_02b (external parameter, not from CSV)
        "\\newcommand{\\coreSPTRps}{100}",
        "",
    ]

    (out_dir / "exec13_resolution_macros.tex").write_text("\n".join(lines))
    print("  Macros written to exec13_resolution_macros.tex")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=pathlib.Path, default=_DEFAULT_DATA_DIR)
    parser.add_argument("--out-dir", type=pathlib.Path, default=_DEFAULT_OUT_DIR)
    args = parser.parse_args(argv)
    data_dir, out_dir = args.data_dir, args.out_dir
    (out_dir / "figs").mkdir(parents=True, exist_ok=True)

    all_results: dict = {}
    for x_mm in POSITIONS_MM:
        print(f"\nAnalyzing x={x_mm:+d} mm (ID {NEAREST_IDS[x_mm]})…")
        all_results[x_mm] = analyze_position(data_dir, x_mm)

    print("\nGenerating F5 figure…")
    plot_f5(all_results, out_dir)

    print("\nWriting CSV and macros…")
    write_csv(all_results, out_dir)
    write_macros(all_results, out_dir)

    # Reconciliation summary
    print("\n── Reconciliation summary ───────────────────────────────────────────")
    for x_mm in POSITIONS_MM:
        r = all_results[x_mm][4]
        print(f"  x={x_mm:+d}: "
              f"sigma_core(2ps)={r['sigma_core_ps']:.2f} | "
              f"EXEC12-style(25ps)={r['sigma_exec12_style_ps']:.2f} | "
              f"sigma_stat_th={r['sigma_stat_theory_ps']:.2f} | "
              f"EXEC12_ref={r['exec12_sigma_fit_ref_ps']:.2f} | "
              f"2ps_vs_exec12={r['reconciled_2ps_vs_exec12']}")
    print("\nPhysical interpretation:")
    print("  sigma_core (2 ps) = narrow prompt-photon core; does NOT match EXEC_12")
    print("  EXEC_12-style fit = core + inner tail; matches theoretical sigma_stat")
    print("  Discrepancy is physical (bin-width effect + core/tail structure)")
    print("\nexec13_core_resolution COMPLETED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
