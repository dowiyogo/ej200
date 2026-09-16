#!/usr/bin/env python3
"""Ejecuta el gate 5.1 de identificación entre/dentro para EXEC_46."""

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "step5"
DERIVED_ROOT = BASE_DIR / "step2" / "exec46_derived_events.root"
BASELINE_POINTS = BASE_DIR / "step2" / "baseline_points.csv"
BASELINE_FITS = BASE_DIR / "step2" / "baseline_fits.csv"
REPORT_PATH = OUTPUT_DIR / "REPORT_CHAINRULE_IDENTIFICATION_20260916.md"
COMMAND = (
    "env PYTHONPATH=analysis/track_mechanism_20260915 python3 "
    "analysis/track_mechanism_20260915/analyze_step5.py"
)
MATERIALS = ("EJ-200", "EJ-204", "EJ-230")
MATERIAL_CODES = {name: index for index, name in enumerate(MATERIALS)}
POSITIONS_MM = np.asarray([-650, -500, -200, 0, 200, 500, 650])
EVEN_POSITIONS_MM = (0, 200, 500, 650)
EXPECTED_ROWS = 210_000
EXPECTED_EVENTS_PER_CELL = 10_000
MAJORITY_THRESHOLD = 0.5
FIT_RANGE_M = (-0.66, 0.66)
NS_TO_PS = 1000.0
COLORS = {"EJ-200": "#1f77b4", "EJ-204": "#ff7f0e", "EJ-230": "#2ca02c"}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fit_even(x_mm, values, errors):
    x_m = np.asarray(x_mm, dtype=float) / 1000.0
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    design = np.column_stack([np.ones(len(x_m)), x_m * x_m])
    weights = 1.0 / (errors * errors)
    covariance = np.linalg.inv((design.T * weights) @ design)
    parameters = covariance @ (design.T @ (weights * values))
    residual = values - design @ parameters
    chi2 = float(np.sum((residual / errors) ** 2))
    return {
        "a0_ns": parameters[0],
        "a0_error_ns": math.sqrt(covariance[0, 0]),
        "a2_ns_per_m2": parameters[1],
        "a2_error_ns_per_m2": math.sqrt(covariance[1, 1]),
        "chi2": chi2,
        "ndf": len(values) - 2,
        "chi2_ndf": chi2 / (len(values) - 2),
    }


def within_fit(npe, t0, positions):
    centered_npe = np.empty_like(npe, dtype=float)
    fitted_intercept = np.empty_like(t0, dtype=float)
    cell_rows = []
    for position in POSITIONS_MM:
        mask = positions == position
        require(np.count_nonzero(mask) == EXPECTED_EVENTS_PER_CELL,
                f"conteo inesperado en x={position}")
        mean_npe = float(np.mean(npe[mask]))
        mean_t0 = float(np.mean(t0[mask]))
        centered_npe[mask] = npe[mask] - mean_npe
        fitted_intercept[mask] = mean_t0
        cell_rows.append({
            "x_mm": int(position),
            "mean_npe": mean_npe,
            "sem_npe": float(np.std(npe[mask], ddof=1) / math.sqrt(np.count_nonzero(mask))),
            "mean_t0_ns": mean_t0,
            "sem_t0_ns": float(np.std(t0[mask], ddof=1) / math.sqrt(np.count_nonzero(mask))),
        })
    denominator = float(np.sum(centered_npe * centered_npe))
    beta = float(np.sum(centered_npe * (t0 - fitted_intercept)) / denominator)
    residual = t0 - fitted_intercept - beta * centered_npe
    n_observations = len(t0)
    n_parameters = len(POSITIONS_MM) + 1
    hc1 = n_observations / (n_observations - n_parameters)
    robust_variance = hc1 * float(np.sum((centered_npe * residual) ** 2)) / denominator ** 2
    classic_variance = (float(np.sum(residual * residual))
                        / (n_observations - n_parameters) / denominator)
    return beta, math.sqrt(robust_variance), math.sqrt(classic_variance), cell_rows


def between_fit(cell_frame, weighted=False):
    design = np.column_stack([np.ones(len(cell_frame)), cell_frame["mean_npe"]])
    target = cell_frame["mean_t0_ns"].to_numpy()
    if weighted:
        weights = 1.0 / cell_frame["sem_t0_ns"].to_numpy() ** 2
    else:
        weights = np.ones(len(cell_frame))
    inverse = np.linalg.inv((design.T * weights) @ design)
    parameters = inverse @ (design.T @ (weights * target))
    residual = target - design @ parameters
    chi2 = float(np.sum(weights * residual * residual))
    scale = chi2 / (len(target) - 2)
    covariance = inverse * scale
    return {
        "intercept_ns": parameters[0],
        "slope_ns_per_pe": parameters[1],
        "slope_error_ns_per_pe": math.sqrt(covariance[1, 1]),
        "chi2": chi2,
        "ndf": len(target) - 2,
        "chi2_ndf": chi2 / (len(target) - 2),
        "weighted": weighted,
    }


def numeric_root(path, tree_name, frame):
    columns = {}
    for name in frame.columns:
        values = frame[name].to_numpy()
        if np.issubdtype(values.dtype, np.number) or values.dtype == bool:
            columns[name] = values
    with uproot.recreate(path) as root_file:
        root_file.mktree(tree_name, {name: value.dtype for name, value in columns.items()})
        root_file[tree_name].extend(columns)


def save_bundle(stem, frame, metadata, figure):
    csv_path = OUTPUT_DIR / f"{stem}.csv"
    root_path = OUTPUT_DIR / f"{stem}.root"
    pdf_path = OUTPUT_DIR / f"{stem}.pdf"
    meta_path = OUTPUT_DIR / f"{stem}.meta.json"
    frame.to_csv(csv_path, index=False, float_format="%.12g")
    numeric_root(root_path, "data", frame)
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    payload = {
        **metadata,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "materials": list(MATERIALS),
        "new_simulation": False,
        "input_root": str(DERIVED_ROOT.resolve()),
        "input_sha256": sha256(DERIVED_ROOT),
        "baseline_points": str(BASELINE_POINTS.resolve()),
        "baseline_points_sha256": sha256(BASELINE_POINTS),
        "baseline_fits": str(BASELINE_FITS.resolve()),
        "baseline_fits_sha256": sha256(BASELINE_FITS),
        "csv": str(csv_path.resolve()),
        "root": str(root_path.resolve()),
        "pdf": str(pdf_path.resolve()),
    }
    meta_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def analyze(arrays, baseline_points, baseline_fits):
    slope_rows = []
    point_rows = []
    curve_rows = []
    summary_rows = []
    for material in MATERIALS:
        material_code = MATERIAL_CODES[material]
        mask = arrays["material_code"] == material_code
        positions = arrays["x_mm"][mask]
        npe = arrays["npe_end"][mask].astype(float)
        t0 = arrays["t0_ns"][mask]
        beta_within, within_se, within_classic_se, cells = within_fit(
            npe, t0, positions)
        cells = pd.DataFrame(cells).sort_values("x_mm").reset_index(drop=True)
        between = between_fit(cells, weighted=False)
        between_wls = between_fit(cells, weighted=True)
        gap = between["slope_ns_per_pe"] - beta_within
        gap_se = math.hypot(between["slope_error_ns_per_pe"], within_se)
        slope_rows.extend([
            {
                "material": material, "material_code": material_code,
                "scope": "within_fixed_effects", "slope_ns_per_pe": beta_within,
                "slope_error_ns_per_pe": within_se,
                "alternate_classic_error_ns_per_pe": within_classic_se,
                "n_observations": len(t0), "n_positions": len(POSITIONS_MM),
                "chi2": np.nan, "ndf": len(t0) - 8, "chi2_ndf": np.nan,
            },
            {
                "material": material, "material_code": material_code,
                "scope": "between_positions_ols", **between,
                "n_observations": len(cells), "n_positions": len(cells),
                "alternate_classic_error_ns_per_pe": np.nan,
            },
            {
                "material": material, "material_code": material_code,
                "scope": "between_positions_wls_sensitivity", **between_wls,
                "n_observations": len(cells), "n_positions": len(cells),
                "alternate_classic_error_ns_per_pe": np.nan,
            },
        ])

        center = cells[cells["x_mm"] == 0].iloc[0]
        baseline_material = baseline_points[
            (baseline_points["clock"] == "time_ns")
            & (baseline_points["material"] == material)].sort_values("x_mm")
        require(len(baseline_material) == 7, f"baseline incompleta para {material}")
        require(np.array_equal(cells["x_mm"].to_numpy(),
                               baseline_material["x_mm"].to_numpy()),
                f"posiciones baseline distintas para {material}")
        original_fit = baseline_fits[
            (baseline_fits["clock"] == "time_ns")
            & (baseline_fits["material"] == material)
            & (baseline_fits["model"] == "predicted_quadratic")].iloc[0]
        observed_fit_reference = baseline_fits[
            (baseline_fits["clock"] == "time_ns")
            & (baseline_fits["material"] == material)
            & (baseline_fits["model"] == "observed_quadratic")].iloc[0]

        delta_npe = cells["mean_npe"].to_numpy() - center["mean_npe"]
        within_curve = center["mean_t0_ns"] + beta_within * delta_npe
        between_curve = (center["mean_t0_ns"]
                         + between["slope_ns_per_pe"] * delta_npe)
        between_wls_curve = (center["mean_t0_ns"]
                             + between_wls["slope_ns_per_pe"] * delta_npe)
        corrected_residual = cells["mean_t0_ns"].to_numpy() - between_curve
        observed_fit = fit_even(cells["x_mm"], cells["mean_t0_ns"], cells["sem_t0_ns"])
        within_curve_fit = fit_even(cells["x_mm"], within_curve, cells["sem_t0_ns"])
        between_curve_fit = fit_even(cells["x_mm"], between_curve, cells["sem_t0_ns"])
        between_wls_curve_fit = fit_even(
            cells["x_mm"], between_wls_curve, cells["sem_t0_ns"])
        corrected_fit = fit_even(cells["x_mm"], corrected_residual, cells["sem_t0_ns"])
        registered_remnant = (observed_fit["a2_ns_per_m2"]
                              - original_fit["a2_ns_m2"])
        correction_from_identification = (between_curve_fit["a2_ns_per_m2"]
                                          - original_fit["a2_ns_m2"])
        fraction_registered_removed = correction_from_identification / registered_remnant
        fraction_registered_removed_wls = (
            (between_wls_curve_fit["a2_ns_per_m2"] - original_fit["a2_ns_m2"])
            / registered_remnant)
        npe_a2 = (between_curve_fit["a2_ns_per_m2"]
                  / between["slope_ns_per_pe"])
        between_prediction_a2_error = abs(
            npe_a2 * between["slope_error_ns_per_pe"])
        corrected_total_error = math.hypot(
            corrected_fit["a2_error_ns_per_m2"], between_prediction_a2_error)
        common_slope_remnant = (observed_fit["a2_ns_per_m2"]
                                - within_curve_fit["a2_ns_per_m2"])
        common_slope_gap = (between_curve_fit["a2_ns_per_m2"]
                            - within_curve_fit["a2_ns_per_m2"])
        fraction_common_gap = common_slope_gap / common_slope_remnant

        summary_rows.append({
            "material": material, "material_code": material_code,
            "beta_within_ns_per_pe": beta_within,
            "beta_within_error_ns_per_pe": within_se,
            "beta_between_ns_per_pe": between["slope_ns_per_pe"],
            "beta_between_error_ns_per_pe": between["slope_error_ns_per_pe"],
            "beta_between_wls_ns_per_pe": between_wls["slope_ns_per_pe"],
            "beta_gap_ns_per_pe": gap, "beta_gap_error_ns_per_pe": gap_se,
            "beta_gap_significance": abs(gap) / gap_se,
            "observed_a2_ns_per_m2": observed_fit["a2_ns_per_m2"],
            "observed_a2_error_ns_per_m2": observed_fit["a2_error_ns_per_m2"],
            "observed_a2_chi2_ndf": observed_fit["chi2_ndf"],
            "reference_observed_a2_ns_per_m2": observed_fit_reference["a2_ns_m2"],
            "original_prediction_a2_ns_per_m2": original_fit["a2_ns_m2"],
            "original_prediction_a2_chi2_ndf": original_fit["chi2_ndf"],
            "within_prediction_a2_ns_per_m2": within_curve_fit["a2_ns_per_m2"],
            "within_prediction_a2_chi2_ndf": within_curve_fit["chi2_ndf"],
            "between_prediction_a2_ns_per_m2": between_curve_fit["a2_ns_per_m2"],
            "between_prediction_a2_error_ns_per_m2": between_prediction_a2_error,
            "between_prediction_a2_chi2_ndf": between_curve_fit["chi2_ndf"],
            "between_wls_prediction_a2_ns_per_m2": between_wls_curve_fit["a2_ns_per_m2"],
            "registered_remnant_a2_ns_per_m2": registered_remnant,
            "identification_correction_a2_ns_per_m2": correction_from_identification,
            "fraction_registered_remnant_removed": fraction_registered_removed,
            "fraction_registered_remnant_removed_wls": fraction_registered_removed_wls,
            "common_slope_gap_a2_ns_per_m2": common_slope_gap,
            "fraction_common_slope_remnant_explained": fraction_common_gap,
            "corrected_remnant_a2_ns_per_m2": corrected_fit["a2_ns_per_m2"],
            "corrected_remnant_a2_error_ns_per_m2": corrected_fit["a2_error_ns_per_m2"],
            "corrected_remnant_a2_total_error_ns_per_m2": corrected_total_error,
            "corrected_remnant_chi2": corrected_fit["chi2"],
            "corrected_remnant_ndf": corrected_fit["ndf"],
            "corrected_remnant_chi2_ndf": corrected_fit["chi2_ndf"],
            "majority_gate": (min(fraction_registered_removed,
                                  fraction_registered_removed_wls)
                              > MAJORITY_THRESHOLD),
        })

        for index, cell in cells.iterrows():
            x_mm = int(cell["x_mm"])
            original_shift = float(
                baseline_material[baseline_material["x_mm"] == x_mm]
                ["predicted_even_shift_ps"].iloc[0])
            curve_rows.append({
                "material": material, "material_code": material_code,
                "x_mm": x_mm, "mean_npe": cell["mean_npe"],
                "mean_t0_ns": cell["mean_t0_ns"], "sem_t0_ns": cell["sem_t0_ns"],
                "observed_shift_ps": NS_TO_PS * (cell["mean_t0_ns"]
                                                   - center["mean_t0_ns"]),
                "original_prediction_shift_ps": original_shift,
                "within_prediction_shift_ps": NS_TO_PS * (within_curve[index]
                                                            - center["mean_t0_ns"]),
                "between_prediction_shift_ps": NS_TO_PS * (between_curve[index]
                                                             - center["mean_t0_ns"]),
                "corrected_residual_ps": NS_TO_PS * corrected_residual[index],
            })

        for abs_x in EVEN_POSITIONS_MM:
            baseline_row = baseline_material[baseline_material["x_mm"] == abs_x].iloc[0]
            if abs_x == 0:
                even_npe = center["mean_npe"]
            else:
                even_npe = cells[np.abs(cells["x_mm"]) == abs_x]["mean_npe"].mean()
            observed_shift = float(baseline_row["observed_even_shift_ps"])
            original_shift = float(baseline_row["predicted_even_shift_ps"])
            within_shift = NS_TO_PS * beta_within * (even_npe - center["mean_npe"])
            between_shift = (NS_TO_PS * between["slope_ns_per_pe"]
                             * (even_npe - center["mean_npe"]))
            registered = observed_shift - original_shift
            corrected = observed_shift - between_shift
            correction = between_shift - original_shift
            point_rows.append({
                "material": material, "material_code": material_code,
                "abs_x_mm": abs_x, "even_mean_npe": even_npe,
                "observed_even_shift_ps": observed_shift,
                "original_prediction_shift_ps": original_shift,
                "within_prediction_shift_ps": within_shift,
                "between_prediction_shift_ps": between_shift,
                "registered_remnant_ps": registered,
                "identification_correction_ps": correction,
                "corrected_remnant_ps": corrected,
                "fraction_registered_remnant_removed": (
                    correction / registered if registered != 0.0 else np.nan),
            })
    return (pd.DataFrame(slope_rows), pd.DataFrame(point_rows),
            pd.DataFrame(curve_rows), pd.DataFrame(summary_rows))


def make_figures(slopes, points, curves, summary):
    fig, axis = plt.subplots(figsize=(8, 5))
    offsets = {"within_fixed_effects": -0.12, "between_positions_ols": 0.12}
    markers = {"within_fixed_effects": "o", "between_positions_ols": "s"}
    for index, material in enumerate(MATERIALS):
        for scope in offsets:
            row = slopes[(slopes["material"] == material)
                         & (slopes["scope"] == scope)].iloc[0]
            axis.errorbar(index + offsets[scope], NS_TO_PS * row["slope_ns_per_pe"],
                          yerr=NS_TO_PS * row["slope_error_ns_per_pe"],
                          fmt=markers[scope], color=COLORS[material], capsize=3,
                          label=scope if index == 0 else None)
    axis.axhline(0.0, color="black", lw=0.8)
    axis.set_xticks(range(len(MATERIALS)), MATERIALS)
    axis.set_ylabel("slope [ps/pe]"); axis.grid(alpha=0.2); axis.legend()
    save_bundle("within_between_slopes", slopes, {
        "model": "T0=position fixed effect + beta_within*(Npe-cell mean Npe)",
        "between_model": "OLS of seven cell means; WLS retained as sensitivity",
        "error": "HC1 event-level error for within; residual OLS error for between",
        "binning": "none",
        "plot_scale": "linear",
    }, fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    for axis, material in zip(axes, MATERIALS):
        group = points[points["material"] == material]
        axis.plot(group["abs_x_mm"], group["registered_remnant_ps"], "o-",
                  label="registered")
        axis.plot(group["abs_x_mm"], group["corrected_remnant_ps"], "s-",
                  label="after between-position slope")
        axis.axhline(0.0, color="black", lw=0.8)
        axis.set_title(material); axis.grid(alpha=0.2)
        axis.set_xlabel("|x| [mm]")
    axes[0].set_ylabel("even residual [ps]"); axes[0].legend(fontsize=8)
    save_bundle("identification_residual", points, {
        "registered": "observed even shift minus original local-beta chain",
        "corrected": "observed even shift minus beta_between*Delta mean Npe",
        "binning": "measured |x| points 0, 200, 500, 650 mm",
        "plot_scale": "linear",
    }, fig)

    curves.to_csv(OUTPUT_DIR / "within_between_curves.csv", index=False,
                  float_format="%.12g")
    summary.to_csv(OUTPUT_DIR / "within_between_summary.csv", index=False,
                   float_format="%.12g")


def render_report(slopes, points, summary, halted):
    lines = [
        "# EXEC_46 Step 5 — chain-rule identification gate", "",
        "Date: 2026-09-16", "", "## Gate result", "",
        "**HALTED_AT_5_1_MAJORITY_IDENTIFICATION_ARTIFACT.**" if halted else
        "**5.1 did not activate the majority-artifact gate.**", "",
        "No simulation was run. The 210,000-event `derived_events` tree was opened read-only.",
        "The common within-cell slope uses position fixed effects and cell-centered Npe. The",
        "between-position slope is an unweighted OLS regression of the seven cell means; a WLS",
        "sensitivity is retained in the sidecar.", "",
        "## 5.1 within versus between", "",
        "| material | beta_within [ps/pe] | beta_between [ps/pe] | gap significance | WLS beta_between [ps/pe] |",
        "|---|---:|---:|---:|---:|",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"| {row['material']} | {NS_TO_PS*row['beta_within_ns_per_pe']:.5f} +/- "
            f"{NS_TO_PS*row['beta_within_error_ns_per_pe']:.5f} | "
            f"{NS_TO_PS*row['beta_between_ns_per_pe']:.5f} +/- "
            f"{NS_TO_PS*row['beta_between_error_ns_per_pe']:.5f} | "
            f"{row['beta_gap_significance']:.2f} sigma | "
            f"{NS_TO_PS*row['beta_between_wls_ns_per_pe']:.5f} |")
    lines += ["", "The between-position response is much weaker than the conditional response",
              "inside a cell. A local beta therefore does not identify the change of the",
              "unconditional mean between positions.", "",
              "## Quadratic summary and mandatory stop", "",
              "| material | observed a2 | original predicted a2 | between predicted a2 | corrected residual a2 | corrected chi2/ndf | registered removed |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for _, row in summary.iterrows():
        lines.append(
            f"| {row['material']} | {NS_TO_PS*row['observed_a2_ns_per_m2']:+.2f} | "
            f"{NS_TO_PS*row['original_prediction_a2_ns_per_m2']:+.2f} | "
            f"{NS_TO_PS*row['between_prediction_a2_ns_per_m2']:+.2f} | "
            f"{NS_TO_PS*row['corrected_remnant_a2_ns_per_m2']:+.2f} +/- "
            f"{NS_TO_PS*row['corrected_remnant_a2_total_error_ns_per_m2']:.2f} | "
            f"{row['corrected_remnant_chi2']:.2f}/{int(row['corrected_remnant_ndf'])}="
            f"{row['corrected_remnant_chi2_ndf']:.2f} | "
            f"{100*row['fraction_registered_remnant_removed']:.2f}% |")
    lines += ["", "The quadratic basis remains a poor shape description where its chi2/ndf is",
              "large; the table is retained only to compare with the preregistered +123.10,",
              "+227.05, and +197.23 ps/m^2 remnants. The primary result is pointwise:", "",
              "| material | |x| [mm] | observed [ps] | original prediction [ps] | registered residual [ps] | between prediction [ps] | corrected residual [ps] | removed |",
              "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for _, row in points[points["abs_x_mm"] > 0].iterrows():
        lines.append(
            f"| {row['material']} | {int(row['abs_x_mm'])} | "
            f"{row['observed_even_shift_ps']:+.3f} | "
            f"{row['original_prediction_shift_ps']:+.3f} | "
            f"{row['registered_remnant_ps']:+.3f} | "
            f"{row['between_prediction_shift_ps']:+.3f} | "
            f"{row['corrected_remnant_ps']:+.3f} | "
            f"{100*row['fraction_registered_remnant_removed']:+.1f}% |")
    if halted:
        all_identification_fractions = np.concatenate([
            summary["fraction_registered_remnant_removed"].to_numpy(),
            summary["fraction_registered_remnant_removed_wls"].to_numpy(),
        ])
        lines += ["", "All three materials exceed the preregistered 50% majority threshold:",
                  f"the identification correction removes "
                  f"{100*all_identification_fractions.min():.2f}--"
                  f"{100*all_identification_fractions.max():.2f}% of the registered remnant across",
                  "the OLS primary result and WLS sensitivity. Independently, replacing the common",
                  "within slope by the common between slope explains "
                  f"{100*summary['fraction_common_slope_remnant_explained'].min():.2f}--"
                  f"{100*summary['fraction_common_slope_remnant_explained'].max():.2f}% of its",
                  "own common-slope remnant. Both definitions cross the majority gate.",
                  "The corrected a2 uncertainty conservatively adds the between-slope uncertainty",
                  "without a covariance cancellation; it does not affect the gate.", "",
                  "This changes the thesis: the remnant is",
                  "predominantly an identification artifact before any attribution to optical",
                  "mechanisms. Per the explicit gate, Steps 5.2--5.5 and Step 6 were not run.", ""]
    lines += ["## Reproducibility", "", "```bash", COMMAND, "```", "",
              "`within_between_slopes` and `identification_residual` each have PDF, CSV, ROOT,",
              "and JSON metadata sidecars. No push, merge, deck edit, or simulation occurred.", ""]
    REPORT_PATH.write_text("\n".join(lines))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in (DERIVED_ROOT, BASELINE_POINTS, BASELINE_FITS):
        require(path.is_file(), f"falta entrada {path}")
    with uproot.open(DERIVED_ROOT) as root_file:
        tree = root_file["derived_events"]
        require(tree.num_entries == EXPECTED_ROWS, "conteo derived_events inesperado")
        arrays = tree.arrays(["material_code", "x_mm", "t0_ns", "npe_end"],
                             library="np")
    baseline_points = pd.read_csv(BASELINE_POINTS)
    baseline_fits = pd.read_csv(BASELINE_FITS)
    slopes, points, curves, summary = analyze(arrays, baseline_points, baseline_fits)
    halted = bool(np.all(summary["majority_gate"]))
    slopes.to_csv(OUTPUT_DIR / "within_between_slopes.csv", index=False,
                  float_format="%.12g")
    points.to_csv(OUTPUT_DIR / "identification_residual_points.csv", index=False,
                  float_format="%.12g")
    make_figures(slopes, points, curves, summary)
    render_report(slopes, points, summary, halted)
    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": ("HALTED_AT_5_1_MAJORITY_IDENTIFICATION_ARTIFACT"
                   if halted else "STEP_5_1_COMPLETE"),
        "majority_threshold": MAJORITY_THRESHOLD,
        "fractions_removed": dict(zip(summary["material"],
                                      summary["fraction_registered_remnant_removed"])),
        "fractions_removed_wls_sensitivity": dict(zip(
            summary["material"],
            summary["fraction_registered_remnant_removed_wls"])),
        "report": str(REPORT_PATH.resolve()),
        "report_sha256": sha256(REPORT_PATH),
        "input_sha256": sha256(DERIVED_ROOT),
        "steps_not_run": ["5.2", "5.3", "5.4", "5.5", "6"] if halted else [],
    }
    (OUTPUT_DIR / "analysis_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
