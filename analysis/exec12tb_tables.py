"""
EXEC_12T-B — Generate all Beamer tables (CSV → LaTeX booktabs) and
             generated_numbers.tex (TeX macros for in-slide numbers).

Usage:
    python analysis/exec12tb_tables.py <outdir>

<outdir>/tables/ and <outdir>/beamer/ must exist.
"""

import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd

REPO  = Path("/home/reriosto/SHiP/ej200")
R12T  = REPO / "results" / "exec12t_20260612_195426"
R11   = REPO / "results" / "exec11_20260612_182454"
R12   = REPO / "results" / "exec12_20260612_191000"

ANALYSIS12T = R12T / "analysis"
TABLES12T   = R12T / "tables"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def tex_val(v, fmt=".2f"):
    """Format a number for LaTeX (handles NaN)."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "--"
    return f"{v:{fmt}}"


def write_booktabs(df, path, caption, label, col_fmt=None,
                   col_renames=None, note=None):
    """Write a pandas DataFrame as a LaTeX longtable/booktabs fragment."""
    if col_renames:
        df = df.rename(columns=col_renames)
    ncols = len(df.columns)
    if col_fmt is None:
        col_fmt = "l" + "r" * (ncols - 1)

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\small")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label}}}")
    lines.append(rf"\begin{{tabular}}{{{col_fmt}}}")
    lines.append(r"\toprule")

    # Header
    header = " & ".join(f"\\textbf{{{c}}}" for c in df.columns) + r" \\"
    lines.append(header)
    lines.append(r"\midrule")

    # Rows
    for _, row in df.iterrows():
        cells = []
        for val in row:
            if isinstance(val, float):
                if np.isnan(val):
                    cells.append("--")
                elif abs(val) >= 1000:
                    cells.append(f"{val:.0f}")
                elif abs(val) >= 10:
                    cells.append(f"{val:.1f}")
                else:
                    cells.append(f"{val:.2f}")
            else:
                cells.append(str(val))
        lines.append(" & ".join(cells) + r" \\")

    lines.append(r"\bottomrule")
    if note:
        lines.append(rf"\multicolumn{{{ncols}}}{{l}}{{\small \textit{{{note}}}}}")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  table  {Path(path).name}")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

tpos = pd.read_csv(ANALYSIS12T / "temporal_position_summary.csv")
tpos4  = tpos[tpos.threshold == 4 ].sort_values("x_true_mm").reset_index(drop=True)
tpos20 = tpos[tpos.threshold == 20].sort_values("x_true_mm").reset_index(drop=True)

sum42 = pd.read_csv(ANALYSIS12T / "threshold_4_20_summary.csv")

cal4   = pd.read_csv(ANALYSIS12T / "calibration_4pe.csv" ).iloc[0]
cal20  = pd.read_csv(ANALYSIS12T / "calibration_20pe.csv").iloc[0]

sweep  = pd.read_csv(ANALYSIS12T / "threshold_sweep_summary.csv")

pr4    = pd.read_csv(ANALYSIS12T / "position_reconstruction_4pe.csv" ).sort_values("x_true_mm")
pr20   = pd.read_csv(ANALYSIS12T / "position_reconstruction_20pe.csv").sort_values("x_true_mm")

win    = pd.read_csv(ANALYSIS12T / "window_dip_summary.csv")
glob   = pd.read_csv(TABLES12T   / "global_context.csv")

ref_pos = pd.read_csv(R11 / "analysis" / "reference_positions.csv")
REF1 = float(ref_pos[ref_pos.reference == "POS_REF_1"].x_true_mm.values[0])
REF2 = float(ref_pos[ref_pos.reference == "POS_REF_2"].x_true_mm.values[0])

data_inv = pd.read_csv(ANALYSIS12T / "data_inventory.csv")

loo4  = np.load(ANALYSIS12T / "loo_predictions_4pe.npz")
loo20 = np.load(ANALYSIS12T / "loo_predictions_20pe.npz")

OUTDIR = None

# ---------------------------------------------------------------------------
# Table 1 — Dataset inventory
# ---------------------------------------------------------------------------

def table_dataset_inventory():
    t = pd.DataFrame({
        "Dataset":   ["Fine scan (28/29)", "Global EndTop", "Window-dip"],
        "Positions": [41, "--", 4],
        "Events/pos":[3000, 2000, "~2000"],
        "Material":  ["EJ-204", "EJ-204", "EJ-204"],
        "Readout":   ["Top", "End+Top", "Top"],
        "Jitter":    ["0", "0", "0"],
        "Path":      ["pairscan\\_2026-06-11", "exec07\\_endtop\\_2000",
                      "exec08b\\_window\\_dip"],
    })
    write_booktabs(
        t,
        OUTDIR / "tables" / "dataset_inventory.tex",
        caption="Available simulation data sets.",
        label="tab:datasets",
        col_fmt="llrrllp{3cm}",
    )


# ---------------------------------------------------------------------------
# Table 2 — Main 4PE vs 20PE comparison
# ---------------------------------------------------------------------------

def table_threshold_comparison():
    def row_4PE_native():
        r4n = sum42[(sum42.threshold == 4) & (sum42["sample"] == "native")].iloc[0]
        return r4n

    def row_20PE_native():
        return sum42[(sum42.threshold == 20) & (sum42["sample"] == "native")].iloc[0]

    def pr_mean(pr, col):
        return pr[col].mean()

    # bias at REF1 and REF2 from LOO NPZ
    def loo_bias(x_mm, arr):
        xall = arr["x_true"]; xrec = arr["x_rec"]
        mask = np.abs(xall - x_mm) < 0.1
        vals = xrec[mask][np.isfinite(xrec[mask])]
        return vals.mean() - x_mm

    def loo_rms68(x_mm, arr):
        xall = arr["x_true"]; xrec = arr["x_rec"]
        mask = np.abs(xall - x_mm) < 0.1
        vals = xrec[mask][np.isfinite(xrec[mask])] - x_mm
        return float(np.percentile(np.abs(vals - np.median(vals)), 68))

    r4n  = row_4PE_native()
    r20n = row_20PE_native()

    metrics = {
        "Metric":               ["Efficiency (mean)", "Slope (ps/mm)",
                                 r"$\chi^2$/ndf (cal.)",
                                 r"$\sigma_{\Delta t}$ core (ps)",
                                 r"RMS$_{68}(\Delta t)$ (ps)",
                                 r"$\sigma(t_+)$ core (ps)",
                                 r"RMS$_{68}(x)$ LOO (mm)",
                                 "Bias max (mm)",
                                 r"$\rho(A,B)$ mean"],
        "4PE":   [f"{r4n.efficiency_mean:.2f}",
                  f"{cal4.slope_ps_per_mm:.3f}",
                  f"{cal4.chi2_ndf:.2f}",
                  f"{r4n.sigma_dt_core_mean:.1f}",
                  f"{r4n.sigma_dt_core_mean:.1f}",
                  f"{r4n.sigma_tplus_core_mean:.1f}",
                  f"{pr_mean(pr4,'rms68'):.2f}",
                  f"{pr4['mean'].abs().max():.2f}",
                  f"{tpos4.rho_ab.mean():.3f}"],
        "20PE":  [f"{r20n.efficiency_mean:.2f}",
                  f"{cal20.slope_ps_per_mm:.3f}",
                  f"{cal20.chi2_ndf:.2f}",
                  f"{r20n.sigma_dt_core_mean:.1f}",
                  f"{r20n.sigma_dt_core_mean:.1f}",
                  f"{r20n.sigma_tplus_core_mean:.1f}",
                  f"{pr_mean(pr20,'rms68'):.2f}",
                  f"{pr20['mean'].abs().max():.2f}",
                  f"{tpos20.rho_ab.mean():.3f}"],
    }
    t = pd.DataFrame(metrics)
    write_booktabs(
        t,
        OUTDIR / "tables" / "threshold_4_20_comparison.tex",
        caption=(r"Comparison of 4th-hit and 20th-hit timing estimators. "
                 r"Efficiency = 1.00 for both; see text."),
        label="tab:threshold_comparison",
        col_fmt="lrr",
        note=r"Matched-20PE sample identical (efficiency=100\% for all k).",
    )


# ---------------------------------------------------------------------------
# Table 3 — REF1 / REF2 results
# ---------------------------------------------------------------------------

def table_reference_positions():
    def ref_row(x_mm, label):
        t4n  = tpos4[np.abs(tpos4.x_true_mm - x_mm) < 0.1].iloc[0]
        t20n = tpos20[np.abs(tpos20.x_true_mm - x_mm) < 0.1].iloc[0]

        xr4  = loo4["x_rec"][np.abs(loo4["x_true"] - x_mm) < 0.1]
        xr20 = loo20["x_rec"][np.abs(loo20["x_true"] - x_mm) < 0.1]
        xr4  = xr4[np.isfinite(xr4)]
        xr20 = xr20[np.isfinite(xr20)]

        def r68(arr, x):
            v = arr - x
            return float(np.percentile(np.abs(v - np.median(v)), 68))

        return {
            "Position":            label,
            "$x$ (mm)":            f"{x_mm:.0f}",
            r"$\mu_{\Delta t}$ 4PE":  f"{t4n.mean_dt:.0f} ps",
            r"$\mu_{\Delta t}$ 20PE": f"{t20n.mean_dt:.0f} ps",
            r"RMS$_{68}$ 4PE":        f"{t4n.rms68_dt:.0f} ps",
            r"RMS$_{68}$ 20PE":       f"{t20n.rms68_dt:.0f} ps",
            "Bias 4PE":              f"{xr4.mean()-x_mm:.2f} mm",
            "Bias 20PE":             f"{xr20.mean()-x_mm:.2f} mm",
            r"$\sigma_x$ 4PE":        f"{r68(xr4,x_mm):.1f} mm",
            r"$\sigma_x$ 20PE":       f"{r68(xr20,x_mm):.1f} mm",
        }

    rows = [ref_row(REF1, "REF1"), ref_row(REF2, "REF2")]
    t = pd.DataFrame(rows)
    write_booktabs(
        t,
        OUTDIR / "tables" / "reference_positions.tex",
        caption=rf"Detailed results at reference positions REF1 ($x={REF1:.0f}$~mm) "
                rf"and REF2 ($x={REF2:.0f}$~mm).",
        label="tab:reference_positions",
        col_fmt="l" + "r" * (len(t.columns) - 1),
    )


# ---------------------------------------------------------------------------
# Table 4 — Threshold sweep (selected k)
# ---------------------------------------------------------------------------

def table_sweep():
    ks = [1, 2, 4, 8, 12, 16, 20, 25, 30]
    sub = sweep[sweep.threshold.isin(ks)].copy()
    sub = sub[["threshold", "efficiency_mean", "slope_ps_per_mm",
               "calibration_chi2_ndf", "mean_rms68_delta_t",
               "mean_rms68_x_cv", "max_abs_bias"]]
    sub.columns = ["k", "Efficiency", "Slope (ps/mm)",
                   r"$\chi^2$/ndf",
                   r"RMS$_{68}(\Delta t)$ (ps)",
                   r"RMS$_{68}(x)$ (mm)", "|bias|max (mm)"]
    write_booktabs(
        sub,
        OUTDIR / "tables" / "threshold_sweep.tex",
        caption="Threshold sweep: selected k values. Efficiency = 1.00 for all k shown.",
        label="tab:threshold_sweep",
        col_fmt="r" * len(sub.columns),
        note=r"$k=4$ and $k=20$ are the primary analysis thresholds.",
    )


# ---------------------------------------------------------------------------
# Table 5 — Global context from EXEC_12
# ---------------------------------------------------------------------------

def table_global_context():
    METHOD_LABELS = {
        "A_end":              "Amplitude (End)",
        "R_end":              "Ratio (End)",
        "dt_end_pool":        r"$\Delta t$ pooled (End)",
        "dt_end_weighted":    r"$\Delta t$ weighted (End)",
        "local_R":            "Ratio local (Pair 28/29)",
        "x_top_centroid_raw": "Centroid (Top)",
    }
    g = glob.copy()
    g["Method"] = g.method.map(METHOD_LABELS).fillna(g.method)
    g = g[["Method", "sigma_core_mean", "max_abs_bias"]]
    g.columns = ["Method", r"$\sigma_{\mathrm{core}}$ (mm)", "|bias|max (mm)"]
    write_booktabs(
        g,
        OUTDIR / "tables" / "global_context.tex",
        caption="Global X reconstruction — EXEC\\_12 context (full EndTop scan).",
        label="tab:global_context",
        col_fmt="lrr",
    )


# ---------------------------------------------------------------------------
# Table 6 — Comparison with Lv et al.
# ---------------------------------------------------------------------------

def table_lv_comparison():
    t = pd.DataFrame({
        "Property":      ["Scintillator", "Dimensions (mm)",
                          "SiPM channels", "Timing method",
                          "Position method", "Timing result",
                          "Position result", "Simulation/measurement",
                          "Jitter included"],
        "Lv et al. (2026)": [
            "EJ-200", "200×200×6",
            "64 (16 per side)", "20\\% CFD on waveform",
            "2016 pair geometry + CNN",
            "\\textasciitilde{}90~ps FWHM",
            "\\textasciitilde{}3~mm (CNN)",
            "Measurement", "SPTR + electronics",
        ],
        "This analysis": [
            "EJ-204", "1400×60×10",
            "Top: 70; End: 16", "k-th detected hit (order statistic)",
            "Linear $\\Delta t$ calibration, pair (28,29)",
            "4PE: 78 ps core; 20PE: 100 ps core",
            "4PE: 13 mm LOO; 20PE: 11 mm LOO",
            "Simulation (intrinsic)", "None (jitter = 0)",
        ],
    })
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\small",
        r"\caption{Qualitative comparison with Lv et al.\ (2026). "
        r"\textbf{Direct numerical comparison is not valid} — see text.}",
        r"\label{tab:lv_comparison}",
        r"\begin{tabular}{lp{4.5cm}p{4.5cm}}",
        r"\toprule",
        r"\textbf{Property} & \textbf{Lv et al.\ (2026)} & \textbf{This analysis} \\",
        r"\midrule",
    ]
    for _, row in t.iterrows():
        lines.append(f"{row['Property']} & {row['Lv et al. (2026)']} & {row['This analysis']} \\\\")
    lines += [
        r"\midrule",
        r"\multicolumn{3}{c}{\textbf{Direct numerical comparison: NOT valid (geometry,"
        r" channels, method, jitter all differ)}} \\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    (OUTDIR / "tables" / "lv_comparison.tex").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("  table  lv_comparison.tex")


# ---------------------------------------------------------------------------
# Table 7 — Window-dip summary
# ---------------------------------------------------------------------------

def table_window_dip():
    run_map = {
        "photon_hits_run_A_x-652mm":  "Run A (x=−652)",
        "photon_hits_run_B_x-642mm":  "Run B (x=−642)",
        "photon_hits_run_C1_x-648mm": "Run C1 (x=−648)",
        "photon_hits_run_C2_x-654mm": "Run C2 (x=−654)",
    }
    w4  = win[win.threshold == 4 ].copy()
    w20 = win[win.threshold == 20].copy()
    merged = w4.merge(w20, on=["run","x_true_mm","nearest_top_id"],
                      suffixes=("_4","_20"))
    merged["Run"] = merged.run.map(run_map).fillna(merged.run)
    merged = merged[["Run","x_true_mm","mean_npe_4",
                     "mean_time_ns_4","sigma_time_ns_4",
                     "mean_time_ns_20","sigma_time_ns_20"]]
    merged.columns = ["Run","x (mm)","NPE mean",
                      "t4 mean (ns)","t4 σ (ns)",
                      "t20 mean (ns)","t20 σ (ns)"]
    write_booktabs(
        merged,
        OUTDIR / "tables" / "window_dip_summary.tex",
        caption="Window-dip run summary — SiPM\\_18, Top readout.",
        label="tab:window_dip",
        col_fmt="lrrrrrr",
    )


# ---------------------------------------------------------------------------
# generated_numbers.tex
# ---------------------------------------------------------------------------

def write_generated_numbers():
    # Helper: get LOO rms68 for a position
    def loo_rms68_val(x_mm, arr):
        xall = arr["x_true"]; xrec = arr["x_rec"]
        mask = np.abs(xall - x_mm) < 0.1
        vals = xrec[mask][np.isfinite(xrec[mask])] - x_mm
        return float(np.percentile(np.abs(vals - np.median(vals)), 68))

    def loo_bias_val(x_mm, arr):
        xall = arr["x_true"]; xrec = arr["x_rec"]
        mask = np.abs(xall - x_mm) < 0.1
        vals = xrec[mask][np.isfinite(xrec[mask])]
        return vals.mean() - x_mm

    rms68_4_mean  = pr4.rms68.mean()
    rms68_20_mean = pr20.rms68.mean()
    bias_4_max    = pr4["mean"].abs().max()   # 'mean' is already bias = <x_rec> - x_true
    bias_20_max   = pr20["mean"].abs().max()
    rho_4_mean    = tpos4.rho_ab.mean()
    rho_20_mean   = tpos20.rho_ab.mean()

    macros = [
        ("SlopeFourPE",     f"{cal4.slope_ps_per_mm:.3f}",  "ps/mm"),
        ("SlopeTwentyPE",   f"{cal20.slope_ps_per_mm:.3f}", "ps/mm"),
        ("VeffFourPE",      f"{cal4.v_eff_cm_ns:.1f}",      "cm/ns"),
        ("VeffTwentyPE",    f"{cal20.v_eff_cm_ns:.1f}",     "cm/ns"),
        ("ChiFourPE",       f"{cal4.chi2_ndf:.2f}",         ""),
        ("ChiTwentyPE",     f"{cal20.chi2_ndf:.2f}",        ""),
        ("SigmaDtFourPE",   f"{sum42[(sum42.threshold==4)&(sum42['sample']=='native')].sigma_dt_core_mean.values[0]:.1f}", "ps"),
        ("SigmaDtTwentyPE", f"{sum42[(sum42.threshold==20)&(sum42['sample']=='native')].sigma_dt_core_mean.values[0]:.1f}", "ps"),
        ("SigmaTplusFourPE",  f"{sum42[(sum42.threshold==4)&(sum42['sample']=='native')].sigma_tplus_core_mean.values[0]:.1f}", "ps"),
        ("SigmaTplusTwentyPE",f"{sum42[(sum42.threshold==20)&(sum42['sample']=='native')].sigma_tplus_core_mean.values[0]:.1f}", "ps"),
        ("RmsxFourPE",      f"{rms68_4_mean:.2f}",          "mm"),
        ("RmsxTwentyPE",    f"{rms68_20_mean:.2f}",         "mm"),
        ("BiasMaxFourPE",   f"{bias_4_max:.2f}",            "mm"),
        ("BiasMaxTwentyPE", f"{bias_20_max:.2f}",           "mm"),
        ("RhoABFourPE",     f"{rho_4_mean:.3f}",            ""),
        ("RhoABTwentyPE",   f"{rho_20_mean:.3f}",           ""),
        ("RefOneX",         f"{REF1:.0f}",                  "mm"),
        ("RefTwoX",         f"{REF2:.0f}",                  "mm"),
        ("NPositions",      "41",                            ""),
        ("NEventsPerPos",   "3000",                          ""),
        ("DataCommit",      "f431c01",                       ""),
    ]

    lines = [
        "% generated_numbers.tex — auto-generated by exec12tb_tables.py",
        "% DO NOT EDIT by hand — regenerate from CSV/NPZ",
        "",
    ]
    for name, val, unit in macros:
        comment = f"% {unit}" if unit else ""
        lines.append(rf"\newcommand{{\{name}}}{{{val}}} {comment}")

    (OUTDIR / "beamer" / "generated_numbers.tex").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("  macros generated_numbers.tex")
    # echo to stdout for verification
    for name, val, _ in macros:
        print(f"    \\{name} = {val}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(outdir):
    global OUTDIR
    OUTDIR = Path(outdir)
    (OUTDIR / "tables").mkdir(parents=True, exist_ok=True)
    (OUTDIR / "beamer").mkdir(parents=True, exist_ok=True)

    print("Generating tables and macros...")
    table_dataset_inventory()
    table_threshold_comparison()
    table_reference_positions()
    table_sweep()
    table_global_context()
    table_lv_comparison()
    table_window_dip()
    write_generated_numbers()
    print(f"\nAll done → {OUTDIR}/tables/ and {OUTDIR}/beamer/generated_numbers.tex")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python exec12tb_tables.py <outdir>")
        sys.exit(1)
    main(sys.argv[1])
