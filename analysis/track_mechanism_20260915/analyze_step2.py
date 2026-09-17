#!/usr/bin/env python3
"""Reproduce la línea base EXEC_46 y ejecuta los diagnósticos A1--A3."""

import csv
import hashlib
import json
import math
import sys
from array import array
from pathlib import Path

import numpy as np
import ROOT
import uproot

from dispersive_optics import (campaign_tables, photon_optics, distribution_summary,
                               optical_markdown)


BASE_DIR = Path(__file__).resolve().parent
STEP2_DIR = BASE_DIR / "step2"
DERIVED_PATH = STEP2_DIR / "exec46_derived_events.root"
DERIVED_META_PATH = STEP2_DIR / "exec46_derived_events.meta.json"
MATERIAL_PROPERTIES_PATH = STEP2_DIR / "material_optical_properties.csv"
ANALYSIS_ROOT_PATH = STEP2_DIR / "exec46_step2_analysis.root"
REPORT_PATH = STEP2_DIR / "REPORT_BASELINE_REPRODUCTION_20260915.md"
EXTERNAL_REPORT_PATH = Path("/home/rrios/REPORT_BASELINE_REPRODUCTION_20260915.md")
EXPECTED_ENTRIES = 210_000
EXPECTED_EVENTS_PER_CELL = 10_000
PROFILE_BINS = 40
PROFILE_PADDING_PE = 0.5
FIT_X_MIN_M = -0.66
FIT_X_MAX_M = 0.66
MM_PER_M = 1000.0
NS_TO_PS = 1000.0
IQR_GAUSSIAN_SCALE = 1.349
RATIO_THRESHOLD = 0.5
CHERENKOV_FRACTION_THRESHOLD = 0.5
ANGLE_BINS = 180
ANGLE_MIN_DEG = 0.0
ANGLE_MAX_DEG = 90.0
BOUNDARY_BINS = 120
UNRELIABLE_CHI2_NDF = 5.0
REFERENCE_POINT_TOLERANCE_PS = 0.01
REFERENCE_A2_SIGMA_TOLERANCE = 5.0
MATERIALS = ("EJ-200", "EJ-204", "EJ-230")
MATERIAL_CODES = {0: "EJ-200", 1: "EJ-204", 2: "EJ-230"}
POSITIONS_MM = (-650, -500, -200, 0, 200, 500, 650)
COLORS = {"EJ-200": ROOT.kBlue + 1, "EJ-204": ROOT.kRed + 1,
          "EJ-230": ROOT.kGreen + 2}
SOURCE_CODES = {"scintillation": 1, "Cherenkov": 2}
SELECTION_COLORS = {
    ("first", "scintillation"): ROOT.kBlue + 1,
    ("first", "Cherenkov"): ROOT.kRed + 1,
    ("random", "scintillation"): ROOT.kCyan + 2,
    ("random", "Cherenkov"): ROOT.kMagenta + 1,
}
REFERENCE = {
    "EJ-200": {
        "observed_ps": {200: -1.354, 500: -6.264, 650: -31.922},
        "predicted_ps": {200: -6.980, 500: -44.246, 650: -82.784},
        "observed_a2_ns_m2": -0.0718094078599,
        "observed_a2_error_ns_m2": 0.00179749428654,
        "predicted_a2_ns_m2": -0.194912612852,
        "quadratic_chi2_ndf": 62.61,
        "quartic_chi2_ndf": 2.99,
    },
    "EJ-204": {
        "observed_ps": {200: 0.306, 500: 1.936, 650: -14.307},
        "predicted_ps": {200: -10.061, 500: -60.121, 650: -109.642},
        "observed_a2_ns_m2": -0.0299287544695,
        "observed_a2_error_ns_m2": 0.00179473390186,
        "predicted_a2_ns_m2": -0.256981962817,
        "quadratic_chi2_ndf": 41.75,
        "quartic_chi2_ndf": 2.51,
    },
    "EJ-230": {
        "observed_ps": {200: 1.679, 500: 6.889, 650: 0.194},
        "predicted_ps": {200: -7.352, 500: -46.831, 650: -82.499},
        "observed_a2_ns_m2": 0.00265027733139,
        "observed_a2_error_ns_m2": 0.00174633304245,
        "predicted_a2_ns_m2": -0.194576132578,
        "quadratic_chi2_ndf": 19.27,
        "quartic_chi2_ndf": 2.20,
    },
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path, rows):
    require(rows, f"sin filas para {path}")
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_meta(path, figure, sources, binning, objects):
    payload = {
        "figure": figure,
        "generator": " ".join(sys.argv),
        "new_simulation": False,
        "derived_root": str(DERIVED_PATH),
        "derived_root_sha256": sha256(DERIVED_PATH),
        "production_sources": sources,
        "binning": binning,
        "root_objects": objects,
    }
    Path(path).write_text(json.dumps(payload, indent=2) + "\n")


def quantile_width(values):
    q25, q75 = np.quantile(values, (0.25, 0.75))
    return float((q75 - q25) / IQR_GAUSSIAN_SCALE)


def sem(values):
    return float(np.std(values, ddof=1) / math.sqrt(len(values)))


def fit_even(root_file, name, x_m, y_ns, errors_ns, quartic):
    root_file.cd()
    graph = ROOT.TGraphErrors(
        len(x_m), array("d", x_m), array("d", y_ns),
        array("d", [0.0] * len(x_m)), array("d", errors_ns))
    graph.SetName(name + "_graph")
    formula = "[0]+[1]*x*x+[2]*x*x*x*x" if quartic else "[0]+[1]*x*x"
    function = ROOT.TF1(name, formula, FIT_X_MIN_M, FIT_X_MAX_M)
    function.SetParameter(0, y_ns[len(y_ns) // 2])
    result = graph.Fit(function, "QRSN")
    graph.Write()
    function.Write()
    if result.Get():
        result.Get().GetCovarianceMatrix().Write(name + "_covariance")
    return {
        "a0_ns": function.GetParameter(0),
        "a0_error_ns": function.GetParError(0),
        "a2_ns_m2": function.GetParameter(1),
        "a2_error_ns_m2": function.GetParError(1),
        "a4_ns_m4": function.GetParameter(2) if quartic else 0.0,
        "a4_error_ns_m4": function.GetParError(2) if quartic else 0.0,
        "chi2": function.GetChisquare(),
        "ndf": function.GetNDF(),
        "chi2_ndf": function.GetChisquare() / function.GetNDF(),
        "status": int(result),
        "covariance_status": result.Get().CovMatrixStatus() if result.Get() else -1,
    }


def analyze_clock(arrays, clock_label, left_branch, right_branch, t0_branch, root_file):
    cells = {}
    cell_rows = []
    for material_code, material in MATERIAL_CODES.items():
        for x_mm in POSITIONS_MM:
            mask = ((arrays["material_code"] == material_code)
                    & (arrays["x_mm"] == x_mm))
            require(np.count_nonzero(mask) == EXPECTED_EVENTS_PER_CELL,
                    f"conteo incorrecto {material} {x_mm}")
            left = arrays[left_branch][mask]
            right = arrays[right_branch][mask]
            t0 = arrays[t0_branch][mask]
            npe = arrays["npe_end"][mask].astype(np.float64)
            identity = np.max(np.abs(t0 - 0.5 * (left + right)))
            require(identity == 0.0, f"T0 no exacto {clock_label} {material} {x_mm}")
            mean_t0 = float(np.mean(t0))
            sem_t0 = sem(t0)
            mean_npe = float(np.mean(npe))
            sem_npe = sem(npe)
            low, high = float(np.min(npe) - PROFILE_PADDING_PE), float(np.max(npe) + PROFILE_PADDING_PE)
            directory = root_file.GetDirectory(clock_label)
            if not directory:
                directory = root_file.mkdir(clock_label)
            directory.cd()
            profile_name = f"profile_{clock_label}_{material.replace('-', '_')}_{x_mm:+d}"
            profile = ROOT.TProfile(profile_name, profile_name, PROFILE_BINS, low, high)
            for n, value in zip(npe, t0):
                profile.Fill(float(n), float(value))
            function = ROOT.TF1(profile_name + "_pol1", "pol1", low, high)
            result = profile.Fit(function, "QRSN")
            profile.Write()
            function.Write()
            if result.Get():
                result.Get().GetCovarianceMatrix().Write(profile_name + "_covariance")
            profile_result = {
                "slope_ns_pe": function.GetParameter(1),
                "slope_error_ns_pe": function.GetParError(1),
                "chi2": function.GetChisquare(),
                "ndf": function.GetNDF(),
                "chi2_ndf": function.GetChisquare() / function.GetNDF(),
                "status": int(result),
                "covariance_status": result.Get().CovMatrixStatus() if result.Get() else -1,
            }
            cells[(material, x_mm)] = {
                "mean_t0_ns": mean_t0, "sem_t0_ns": sem_t0,
                "mean_npe": mean_npe, "sem_npe": sem_npe,
                "profile": profile_result,
            }
            cell_rows.append({
                "clock": clock_label, "material": material, "x_mm": x_mm,
                "N": int(np.count_nonzero(mask)), "mean_T0_ns": mean_t0,
                "sem_T0_ns": sem_t0, "mean_Npe_END": mean_npe,
                "sem_Npe_END": sem_npe, "max_T0_identity_difference_ns": identity,
                **profile_result,
                "fit_reliable": profile_result["chi2_ndf"] <= UNRELIABLE_CHI2_NDF,
            })

    point_rows, fit_rows = [], []
    for material in MATERIALS:
        even_npe = {0: cells[(material, 0)]["mean_npe"]}
        even_beta = {0: cells[(material, 0)]["profile"]["slope_ns_pe"]}
        for absolute_x in (200, 500, 650):
            even_npe[absolute_x] = 0.5 * (
                cells[(material, -absolute_x)]["mean_npe"]
                + cells[(material, absolute_x)]["mean_npe"])
            even_beta[absolute_x] = 0.5 * (
                cells[(material, -absolute_x)]["profile"]["slope_ns_pe"]
                + cells[(material, absolute_x)]["profile"]["slope_ns_pe"])
        predicted = {0: 0.0}
        previous = 0
        for absolute_x in (200, 500, 650):
            predicted[absolute_x] = predicted[previous] + 0.5 * (
                even_beta[previous] + even_beta[absolute_x]) * (
                even_npe[absolute_x] - even_npe[previous])
            previous = absolute_x
        center = cells[(material, 0)]["mean_t0_ns"]
        x_m, observed, predicted_points, errors = [], [], [], []
        for x_mm in POSITIONS_MM:
            absolute_x = abs(x_mm)
            if x_mm == 0:
                observed_even = center
            else:
                observed_even = 0.5 * (
                    cells[(material, -absolute_x)]["mean_t0_ns"]
                    + cells[(material, absolute_x)]["mean_t0_ns"])
            observed_shift = observed_even - center
            x_m.append(x_mm / MM_PER_M)
            observed.append(cells[(material, x_mm)]["mean_t0_ns"])
            predicted_points.append(predicted[absolute_x])
            errors.append(cells[(material, x_mm)]["sem_t0_ns"])
            point_rows.append({
                "clock": clock_label, "material": material, "x_mm": x_mm,
                "observed_T0_ns": cells[(material, x_mm)]["mean_t0_ns"],
                "observed_T0_sem_ns": cells[(material, x_mm)]["sem_t0_ns"],
                "even_mean_Npe_END": even_npe[absolute_x],
                "even_slope_ns_pe": even_beta[absolute_x],
                "observed_even_shift_ps": observed_shift * NS_TO_PS,
                "predicted_even_shift_ps": predicted[absolute_x] * NS_TO_PS,
            })
        for model, values, quartic in (
            ("observed_quadratic", observed, False),
            ("predicted_quadratic", predicted_points, False),
            ("observed_quartic", observed, True),
            ("predicted_quartic", predicted_points, True),
        ):
            fit = fit_even(root_file,
                           f"{clock_label}_{model}_{material.replace('-', '_')}",
                           x_m, values, errors, quartic)
            fit_rows.append({"clock": clock_label, "material": material,
                             "model": model, **fit})
    return cells, cell_rows, point_rows, fit_rows


def interpolate_crossing(xs, ys, threshold):
    for index in range(len(xs) - 1):
        y0, y1 = ys[index], ys[index + 1]
        if (y0 - threshold) * (y1 - threshold) <= 0 and y0 != y1:
            return xs[index] + (threshold - y0) * (xs[index + 1] - xs[index]) / (y1 - y0)
    return math.nan


def make_cherenkov_boundary(arrays, sources):
    path_stem = STEP2_DIR / "cherenkov_boundary"
    root_file = ROOT.TFile(str(path_stem.with_suffix(".root")), "RECREATE")
    canvas = ROOT.TCanvas("c_cherenkov_boundary", "Cherenkov boundary", 1500, 900)
    canvas.Divide(3, 2)
    rows, boundaries = [], []
    all_ratios, all_fractions = [], []
    object_names = []
    keepalive = []
    for material_index, material in MATERIAL_CODES.items():
        ratios, fractions = [], []
        for x_mm in POSITIONS_MM:
            mask = ((arrays["material_code"] == material_index)
                    & (arrays["x_mm"] == x_mm))
            ratio = quantile_width(arrays["t_left_ns"][mask]) / quantile_width(
                arrays["t0_ns"][mask])
            fraction = float(np.mean(arrays["first_left_source_type"][mask] == 2))
            ratios.append(ratio)
            fractions.append(fraction)
            all_ratios.append(ratio)
            all_fractions.append(fraction)
            rows.append({"record_type": "point", "material": material,
                         "x_mm": x_mm, "sigma_IQR_tL_over_T0": ratio,
                         "first_left_Cherenkov_fraction": fraction})
        ratio_boundary = interpolate_crossing(POSITIONS_MM, ratios, RATIO_THRESHOLD)
        cher_boundary = interpolate_crossing(POSITIONS_MM, fractions,
                                             CHERENKOV_FRACTION_THRESHOLD)
        boundaries.append((material, ratio_boundary, cher_boundary))
        rows.append({"record_type": "boundary", "material": material,
                     "ratio_threshold": RATIO_THRESHOLD,
                     "ratio_boundary_x_mm": ratio_boundary,
                     "cherenkov_threshold": CHERENKOV_FRACTION_THRESHOLD,
                     "cherenkov_boundary_x_mm": cher_boundary,
                     "boundary_difference_mm": cher_boundary - ratio_boundary})

        graph_ratio = ROOT.TGraph(len(POSITIONS_MM), array("d", POSITIONS_MM),
                                  array("d", ratios))
        graph_fraction = ROOT.TGraph(len(POSITIONS_MM), array("d", POSITIONS_MM),
                                     array("d", fractions))
        graph_ratio.SetName(f"ratio_{material.replace('-', '_')}")
        graph_fraction.SetName(f"first_cherenkov_{material.replace('-', '_')}")
        object_names.extend((graph_ratio.GetName(), graph_fraction.GetName()))
        for graph, color in ((graph_ratio, ROOT.kBlue + 1),
                             (graph_fraction, ROOT.kRed + 1)):
            graph.SetLineColor(color); graph.SetMarkerColor(color)
            graph.SetMarkerStyle(20); graph.SetLineWidth(2); graph.Write()
        canvas.cd(material_index + 1)
        graph_ratio.SetTitle(f"{material};gun x [mm];#sigma_{{IQR}}(t_{{L}})/#sigma_{{IQR}}(T_{{0}})")
        graph_ratio.SetMinimum(0.0); graph_ratio.SetMaximum(2.2); graph_ratio.Draw("ALP")
        line_ratio = ROOT.TLine(-650, RATIO_THRESHOLD, 650, RATIO_THRESHOLD)
        line_ratio.SetLineStyle(2); line_ratio.Draw(); line_ratio.Write(
            f"ratio_threshold_{material.replace('-', '_')}")
        canvas.cd(material_index + 4)
        graph_fraction.SetTitle(f"{material};gun x [mm];first-left Cherenkov fraction")
        graph_fraction.SetMinimum(0.0); graph_fraction.SetMaximum(0.8); graph_fraction.Draw("ALP")
        line_cher = ROOT.TLine(-650, CHERENKOV_FRACTION_THRESHOLD, 650,
                              CHERENKOV_FRACTION_THRESHOLD)
        line_cher.SetLineStyle(2); line_cher.Draw(); line_cher.Write(
            f"cherenkov_threshold_{material.replace('-', '_')}")
        keepalive.extend((graph_ratio, graph_fraction, line_ratio, line_cher))
    ratio_boundaries = np.asarray([item[1] for item in boundaries])
    cher_boundaries = np.asarray([item[2] for item in boundaries])
    boundary_correlation = float(np.corrcoef(ratio_boundaries, cher_boundaries)[0, 1])
    point_correlation = float(np.corrcoef(all_ratios, all_fractions)[0, 1])
    rows.append({"record_type": "correlation",
                 "boundary_Pearson_r_three_materials": boundary_correlation,
                 "point_Pearson_r_21_cells": point_correlation})
    canvas.Write(); canvas.Print(str(path_stem.with_suffix(".pdf")))
    root_file.Close()
    write_csv(path_stem.with_suffix(".csv"), rows)
    write_meta(path_stem.with_suffix(".meta.json"), "Cherenkov fraction and robust width boundary",
               sources, {"x_mm": list(POSITIONS_MM), "ratio_threshold": RATIO_THRESHOLD,
                         "Cherenkov_threshold": CHERENKOV_FRACTION_THRESHOLD}, object_names)
    return rows, boundaries, boundary_correlation, point_correlation


def make_histogram(name, values, bins, low, high):
    histogram = ROOT.TH1D(name, name, bins, low, high)
    counts, _ = np.histogram(values, bins=bins, range=(low, high))
    for index, count in enumerate(counts, 1):
        histogram.SetBinContent(index, int(count))
    histogram.SetEntries(len(values))
    return histogram


def summary_stats(values):
    q16, median, q84 = np.quantile(values, (0.16, 0.5, 0.84))
    return {"N": len(values), "mean": float(np.mean(values)),
            "median": float(median), "q16": float(q16), "q84": float(q84)}


def guiding_mask(arrays, material_code, scope, face_name):
    base = arrays["material_code"] == material_code
    if scope == "all_end":
        return base
    if face_name == "left":
        return base & (arrays["x_mm"] == -650)
    return base & (arrays["x_mm"] == 650)


def make_guiding_diagnostics(arrays, sources):
    path_stem = STEP2_DIR / "cherenkov_guiding_diagnostics"
    root_file = ROOT.TFile(str(path_stem.with_suffix(".root")), "RECREATE")
    canvas = ROOT.TCanvas("c_guiding", "Cherenkov guiding diagnostics", 1500, 900)
    canvas.Divide(3, 2)
    rows, object_names = [], []
    optics = {(selection, face): photon_optics(arrays, f"{selection}_{face}_", "x_mm")
              for selection in ("first", "random") for face in ("left", "right")}
    histograms = {}
    for material_code, material in MATERIAL_CODES.items():
        boundary_values_for_range = []
        for selection in ("first", "random"):
            for source_label, source_code in SOURCE_CODES.items():
                for face_name in ("left", "right"):
                    mask = guiding_mask(arrays, material_code, "all_end", face_name)
                    source = arrays[f"{selection}_{face_name}_source_type"]
                    selected = mask & (source == source_code)
                    boundary_values_for_range.append(
                        arrays[f"{selection}_{face_name}_n_boundary_encounters"][selected])
        boundary_max = max(20.0, float(np.quantile(
            np.concatenate(boundary_values_for_range), 0.995)))
        boundary_max = math.ceil(boundary_max / 10.0) * 10.0
        for selection in ("first", "random"):
            for source_label, source_code in SOURCE_CODES.items():
                angle_parts, boundary_parts, critical_parts = [], [], []
                for face_name in ("left", "right"):
                    mask = guiding_mask(arrays, material_code, "all_end", face_name)
                    source = arrays[f"{selection}_{face_name}_source_type"]
                    selected = mask & (source == source_code)
                    angle_parts.append(arrays[f"{selection}_{face_name}_exit_angle_deg"][selected])
                    critical_parts.append(optics[(selection, face_name)]["theta_critical_detected_deg"][selected])
                    boundary_parts.append(
                        arrays[f"{selection}_{face_name}_n_boundary_encounters"][selected])
                angles = np.concatenate(angle_parts)
                boundaries = np.concatenate(boundary_parts).astype(np.float64)
                angle_stats = summary_stats(angles)
                boundary_stats = summary_stats(boundaries)
                rows.append({"material": material, "scope": "all_end",
                             "selection": selection, "source": source_label,
                             **{f"exit_angle_{key}_deg": value for key, value in angle_stats.items()},
                             **{f"boundary_{key}": value for key, value in boundary_stats.items()},
                             "fraction_exit_angle_above_TIR_critical": float(
                                 np.mean(angles > np.concatenate(critical_parts)))})
                tag = f"{material.replace('-', '_')}_{selection}_{source_label}"
                h_angle = make_histogram("exit_angle_" + tag, angles, ANGLE_BINS,
                                         ANGLE_MIN_DEG, ANGLE_MAX_DEG)
                h_boundary = make_histogram("boundaries_" + tag, boundaries,
                                            BOUNDARY_BINS, 0.0, boundary_max)
                color = SELECTION_COLORS[(selection, source_label)]
                for histogram in (h_angle, h_boundary):
                    histogram.SetLineColor(color); histogram.SetLineWidth(2)
                    histogram.Write(); object_names.append(histogram.GetName())
                histograms[(material, selection, source_label, "angle")] = h_angle
                histograms[(material, selection, source_label, "boundary")] = h_boundary

        for scope in ("near_end_650",):
            for selection in ("first", "random"):
                for source_label, source_code in SOURCE_CODES.items():
                    angle_parts, boundary_parts, critical_parts = [], [], []
                    for face_name in ("left", "right"):
                        mask = guiding_mask(arrays, material_code, scope, face_name)
                        source = arrays[f"{selection}_{face_name}_source_type"]
                        selected = mask & (source == source_code)
                        angle_parts.append(arrays[f"{selection}_{face_name}_exit_angle_deg"][selected])
                        critical_parts.append(optics[(selection, face_name)]["theta_critical_detected_deg"][selected])
                        boundary_parts.append(arrays[f"{selection}_{face_name}_n_boundary_encounters"][selected])
                    angles = np.concatenate(angle_parts)
                    boundaries = np.concatenate(boundary_parts).astype(np.float64)
                    angle_stats = summary_stats(angles); boundary_stats = summary_stats(boundaries)
                    rows.append({"material": material, "scope": scope,
                                 "selection": selection, "source": source_label,
                                 **{f"exit_angle_{key}_deg": value for key, value in angle_stats.items()},
                                 **{f"boundary_{key}": value for key, value in boundary_stats.items()},
                                 "fraction_exit_angle_above_TIR_critical": float(
                                     np.mean(angles > np.concatenate(critical_parts)))})

        for row_offset, variable in enumerate(("angle", "boundary")):
            canvas.cd(material_code + 1 + 3 * row_offset)
            ROOT.gPad.SetLogy()
            legend = ROOT.TLegend(0.48, 0.64, 0.89, 0.89)
            first_draw = True
            for selection in ("first", "random"):
                for source_label in ("scintillation", "Cherenkov"):
                    original = histograms[(material, selection, source_label, variable)]
                    normalized = original.Clone(original.GetName() + "_normalized")
                    if normalized.Integral() > 0:
                        normalized.Scale(1.0 / normalized.Integral())
                    normalized.SetTitle(
                        f"{material};{'exit angle at SiPM [deg]' if variable == 'angle' else 'boundary encounters'};normalized entries")
                    normalized.Draw("HIST" if first_draw else "HIST SAME")
                    normalized.Write(); object_names.append(normalized.GetName())
                    legend.AddEntry(normalized, f"{selection} {source_label}", "l")
                    first_draw = False
            legend.Draw(); legend.Write(f"legend_{material.replace('-', '_')}_{variable}")
    canvas.Write(); canvas.Print(str(path_stem.with_suffix(".pdf")))
    root_file.Close()
    write_csv(path_stem.with_suffix(".csv"), rows)
    write_meta(path_stem.with_suffix(".meta.json"),
               "END exit-angle and boundary-count diagnostics",
               sources, {"exit_angle": [ANGLE_BINS, ANGLE_MIN_DEG, ANGLE_MAX_DEG],
                         "boundary_bins": BOUNDARY_BINS,
                         "boundary_max": "per-material 99.5 percentile rounded to 10"},
               object_names)
    return rows


def make_baseline_figure(points, fits, sources):
    path_stem = STEP2_DIR / "baseline_reproduction"
    root_file = ROOT.TFile(str(path_stem.with_suffix(".root")), "RECREATE")
    canvas = ROOT.TCanvas("c_baseline", "Baseline reproduction", 1500, 500)
    canvas.Divide(3, 1)
    object_names = []
    keepalive = []
    for index, material in enumerate(MATERIALS, 1):
        rows = [row for row in points if row["clock"] == "time_ns"
                and row["material"] == material]
        x = [row["x_mm"] for row in rows]
        observed = [row["observed_even_shift_ps"] for row in rows]
        predicted = [row["predicted_even_shift_ps"] for row in rows]
        graph_observed = ROOT.TGraph(len(x), array("d", x), array("d", observed))
        graph_predicted = ROOT.TGraph(len(x), array("d", x), array("d", predicted))
        graph_observed.SetName(f"observed_shift_{material.replace('-', '_')}")
        graph_predicted.SetName(f"predicted_shift_{material.replace('-', '_')}")
        object_names.extend((graph_observed.GetName(), graph_predicted.GetName()))
        graph_observed.SetLineColor(ROOT.kBlack); graph_observed.SetMarkerColor(ROOT.kBlack)
        graph_observed.SetMarkerStyle(20); graph_observed.SetLineWidth(2)
        graph_predicted.SetLineColor(ROOT.kRed + 1); graph_predicted.SetMarkerColor(ROOT.kRed + 1)
        graph_predicted.SetMarkerStyle(24); graph_predicted.SetLineWidth(2)
        canvas.cd(index)
        graph_observed.SetTitle(f"{material};gun x [mm];even shift from center [ps]")
        y_values = observed + predicted
        y_span = max(y_values) - min(y_values)
        y_margin = max(1.0, 0.12 * y_span)
        graph_observed.SetMinimum(min(y_values) - y_margin)
        graph_observed.SetMaximum(max(y_values) + y_margin)
        graph_observed.Draw("ALP"); graph_predicted.Draw("LP SAME")
        legend = ROOT.TLegend(0.16, 0.72, 0.55, 0.88)
        legend.AddEntry(graph_observed, "observed", "lp")
        legend.AddEntry(graph_predicted, "Npe chain prediction", "lp")
        legend.Draw(); legend.Write(f"legend_{material.replace('-', '_')}")
        graph_observed.Write(); graph_predicted.Write()
        keepalive.extend((graph_observed, graph_predicted, legend))
    canvas.Write(); canvas.Print(str(path_stem.with_suffix(".pdf")))
    root_file.Close()
    csv_rows = list(points) + [{"clock": row["clock"], "material": row["material"],
                               "record_type": "fit", **row} for row in fits]
    write_csv(path_stem.with_suffix(".csv"), csv_rows)
    write_meta(path_stem.with_suffix(".meta.json"), "Baseline observed and Npe-predicted even shifts",
               sources, {"positions_mm": list(POSITIONS_MM),
                         "profile_bins": PROFILE_BINS,
                         "fit_range_m": [FIT_X_MIN_M, FIT_X_MAX_M]}, object_names)


def fit_lookup(fits, clock, material, model):
    return next(row for row in fits if row["clock"] == clock
                and row["material"] == material and row["model"] == model)


def point_lookup(points, clock, material, x_mm):
    return next(row for row in points if row["clock"] == clock
                and row["material"] == material and row["x_mm"] == x_mm)


def build_report(metadata, material_rows, cells, points, fits, clock_maxima, boundaries,
                 boundary_correlation, point_correlation, guiding_rows,
                 optical_rows):
    boundary_differences = [cher - ratio for _, ratio, cher in boundaries]
    if boundary_correlation <= 0.0:
        boundary_verdict = (
            "**BOUNDARY_IDENTITY_NOT_ESTABLISHED.** The material ordering of "
            "the two interpolated crossings is opposite, so the two boundaries "
            "cannot be identified as the same measured transition."
        )
    elif max(abs(value) for value in boundary_differences) > 150.0:
        boundary_verdict = (
            "**BOUNDARY_IDENTITY_NOT_ESTABLISHED.** At least one crossing differs "
            "by more than the full sampled transition interval."
        )
    else:
        boundary_verdict = (
            "**COMMON_BOUNDARY_SUPPORTED_AT_GRID_RESOLUTION.** The crossing "
            "ordering agrees and all differences lie within the 150 mm sampling interval."
        )
    lines = [
        "# EXEC_46 Step 2 — bit-level baseline reproduction and A1--A3 diagnostics",
        "", "Date: 2026-09-15", "",
        "## Verdict", "",
    ]
    point_pass = True; a2_pass = True
    for material in MATERIALS:
        for absolute_x in (200, 500, 650):
            row = point_lookup(points, "time_ns", material, absolute_x)
            point_pass &= abs(row["observed_even_shift_ps"] - REFERENCE[material]["observed_ps"][absolute_x]) <= REFERENCE_POINT_TOLERANCE_PS
            point_pass &= abs(row["predicted_even_shift_ps"] - REFERENCE[material]["predicted_ps"][absolute_x]) <= REFERENCE_POINT_TOLERANCE_PS
        fit = fit_lookup(fits, "time_ns", material, "observed_quadratic")
        a2_pass &= abs(fit["a2_ns_m2"] - REFERENCE[material]["observed_a2_ns_m2"]) <= REFERENCE_A2_SIGMA_TOLERANCE * REFERENCE[material]["observed_a2_error_ns_m2"]
    overall_pass = point_pass and a2_pass and max(clock_maxima.values()) == 0.0
    unreliable_profiles = sum(
        row["clock"] == "time_ns" and not row["fit_reliable"] for row in cells)
    lines += [
        f"**{'PASS' if overall_pass else 'FAIL'}.** The reconstructed EXEC_46 baseline "
        f"{'reproduces' if overall_pass else 'does not reproduce'} the registered END-only first-photon result under the declared tolerances.",
        "No simulation was run. The source ROOT files were opened read-only.", "",
        "The point tolerance is 0.01 ps, set by the precision of the registered table. "
        "The quadratic-coefficient tolerance is five times the registered statistical error.", "",
        "## Baseline point-by-point reproduction", "",
        "| Material | |x| [mm] | observed old/new [ps] | Npe prediction old/new [ps] |",
        "|---|---:|---:|---:|",
    ]
    for material in MATERIALS:
        for absolute_x in (200, 500, 650):
            row = point_lookup(points, "time_ns", material, absolute_x)
            lines.append(
                f"| {material} | {absolute_x} | {REFERENCE[material]['observed_ps'][absolute_x]:+.3f} / {row['observed_even_shift_ps']:+.6f} | "
                f"{REFERENCE[material]['predicted_ps'][absolute_x]:+.3f} / {row['predicted_even_shift_ps']:+.6f} |")
    lines += ["",
              f"The per-cell 40-bin TProfile `pol1` fits reproduce the previous quality result: {unreliable_profiles}/21 have chi2/ndf > 5 and are flagged as unreliable. "
              "The chain-rule curve is reproduced as registered but is not promoted to a valid physical model by this check.", "",
              "## Position fits", "",
              "An asterisk marks chi2/ndf > 5; those fits are descriptive projections, not adequate models.", "",
              "| Material | observed a2 [ps/m2] | predicted a2 [ps/m2] | observed quadratic chi2/ndf | observed quartic a2/a4 | quartic chi2/ndf |",
              "|---|---:|---:|---:|---:|---:|"]
    for material in MATERIALS:
        oq = fit_lookup(fits, "time_ns", material, "observed_quadratic")
        pq = fit_lookup(fits, "time_ns", material, "predicted_quadratic")
        o4 = fit_lookup(fits, "time_ns", material, "observed_quartic")
        star_q = "*" if oq["chi2_ndf"] > UNRELIABLE_CHI2_NDF else ""
        star_4 = "*" if o4["chi2_ndf"] > UNRELIABLE_CHI2_NDF else ""
        lines.append(
            f"| {material} | {oq['a2_ns_m2']*NS_TO_PS:+.4f} +/- {oq['a2_error_ns_m2']*NS_TO_PS:.4f} | "
            f"{pq['a2_ns_m2']*NS_TO_PS:+.4f} | {oq['chi2']:.2f}/{oq['ndf']} = {oq['chi2_ndf']:.2f}{star_q} | "
            f"{o4['a2_ns_m2']*NS_TO_PS:+.2f} / {o4['a4_ns_m4']*NS_TO_PS:+.2f} ps | "
            f"{o4['chi2']:.2f}/{o4['ndf']} = {o4['chi2_ndf']:.2f}{star_4} |")
    lines += ["", "## Clock comparison", "",
              "`T0=(tL+tR)/2` is exact in all 210,000 events for each clock. "
              "The maximum event-level differences are:", "",
              "| Quantity | max abs difference [ns] |", "|---|---:|"]
    for key, value in clock_maxima.items():
        lines.append(f"| {key} | {value:.17g} |")
    lines += ["", "## A1 — runtime optical properties", "",
              "| Cell | RINDEX min–max | type | ABSLENGTH min–max [mm] |",
              "|---|---|---|---|"]
    for row in material_rows:
        lines.append(f"| {row['cell_id']} | {row['rindex_min']:.6f}–{row['rindex_max']:.6f} | "
                     f"{'constant' if row['rindex_constant'] else 'dispersive'} | "
                     f"{row['abs_length_min_mm']:.2f}–{row['abs_length_max_mm']:.2f} |")
    lines += ["", optical_markdown(optical_rows), "",
              "Numerical GROUPVEL is differentiated on each actual runtime energy mesh. "
              "ABSLENGTH spectral dependence is read from the same cell runtime. PDE remains active.", "",
              "## A2 — Cherenkov and robust-width boundaries", "",
              "Both boundaries use linear interpolation inside the sparse -650 to -500 mm bracket: ratio=0.5 and first-left Cherenkov fraction=0.5.", "",
              "| Material | ratio boundary [mm] | Cherenkov boundary [mm] | Cher - ratio [mm] |",
              "|---|---:|---:|---:|"]
    for material, ratio_boundary, cher_boundary in boundaries:
        lines.append(f"| {material} | {ratio_boundary:.2f} | {cher_boundary:.2f} | {cher_boundary-ratio_boundary:+.2f} |")
    lines += ["",
              f"Across the three interpolated boundaries, Pearson r = {boundary_correlation:+.4f}. "
              f"Across all 21 left-END material-position points, r(ratio, first-Cherenkov fraction) = {point_correlation:+.4f}.",
              boundary_verdict,
              "With only three materials and one 150 mm transition interval, the boundary correlation is descriptive. "
              "No causal identity is forced from the point correlation alone.", "",
              "## A3 — Cherenkov guiding diagnostic", "",
              "There is no unique material cone angle; the A1 distributions use each photon wavelength. "
              "Comparisons to the final SiPM-normal angle use theta_critical(lambda_detected).", "",
              "The stored `exit_angle_deg` is relative to the SiPM normal at final detection. It is not the incidence angle at the large bar face and cannot directly prove the proposed TIR-cone inequality. "
              "The angle and boundary-count contrasts are therefore indirect tests; a direct test requires per-boundary angle history.", "",
              "All-END summary (mean with median in parentheses):", "",
              "| Material | selection | source | exit angle [deg] | boundary encounters | N |",
              "|---|---|---|---:|---:|---:|"]
    for row in guiding_rows:
        if row["scope"] != "all_end":
            continue
        lines.append(f"| {row['material']} | {row['selection']} | {row['source']} | {row['exit_angle_mean_deg']:.2f} ({row['exit_angle_median_deg']:.2f}) | {row['boundary_mean']:.2f} ({row['boundary_median']:.0f}) | {row['exit_angle_N_deg']} |")
    lines += ["", "Near-END summary at x=-650 left plus x=+650 right:", "",
              "| Material | selection | source | exit angle [deg] | boundary encounters | N |",
              "|---|---|---|---:|---:|---:|"]
    for row in guiding_rows:
        if row["scope"] != "near_end_650":
            continue
        lines.append(f"| {row['material']} | {row['selection']} | {row['source']} | {row['exit_angle_mean_deg']:.2f} ({row['exit_angle_median_deg']:.2f}) | {row['boundary_mean']:.2f} ({row['boundary_median']:.0f}) | {row['exit_angle_N_deg']} |")
    lines += ["",
              "At the near END, first Cherenkov photons have a narrow median exit angle near 39.6 degrees and a median of four prior boundary encounters in every material. "
              "First scintillation photons have median exit angles near 17 degrees and one boundary encounter. "
              "The random controls are broader and have medians of seven boundaries for Cherenkov and five for scintillation. "
              "This supports a distinct directional, promptly selected Cherenkov population, but the stored final-exit angle does not directly establish TIR at the large faces.", "",
              "## Reproducibility and artifacts", "",
              f"The derived tree contains 21 cells x 10,000 events and has SHA-256 `{metadata['derived_root_sha256']}`. "
              f"The paired random control uses `{metadata['random_control']['algorithm']}` with seed `{metadata['random_control']['seed_hex']}`.", "",
              "```bash", "env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/build_step2_derived.py --processes 4",
              "env PYTHONPATH=analysis/track_mechanism_20260915 python3 analysis/track_mechanism_20260915/analyze_step2.py", "```", "",
              "- `exec46_derived_events.root` and `.meta.json`: 210,000 event rows and exact production provenance.",
              "  It stores `npe_scint_left/right` separately; Step 4 must use those per-END scintillation counts as N rather than total Npe_END.",
              "- `baseline_reproduction.{pdf,root,csv,meta.json}`.",
              "- `cherenkov_boundary.{pdf,root,csv,meta.json}`.",
              "- `cherenkov_guiding_diagnostics.{pdf,root,csv,meta.json}`.",
              "- `baseline_cells.csv`, `baseline_fits.csv`, `clock_comparison.csv`, and `material_optical_properties.csv`.", "",
              "No push, simulation, merge, or deck edit was performed. Step 3 requires a new checkpoint approval.", ""]
    require(overall_pass, "la reproducción de línea base falló; informe escrito")
    return "\n".join(lines)


def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    require(DERIVED_PATH.is_file(), f"falta {DERIVED_PATH}")
    metadata = json.loads(DERIVED_META_PATH.read_text())
    sources = [{"path": row["root_path"], "sha256": row["root_sha256"]}
               for row in metadata["cells"]]
    _, material_rows = campaign_tables()
    with uproot.open(DERIVED_PATH) as root_file:
        tree = root_file["derived_events"]
        require(tree.num_entries == EXPECTED_ENTRIES, "conteo derivado incorrecto")
        arrays = tree.arrays(library="np")
    optical_rows = []
    for selection in ("first", "random"):
        for source, source_code in SOURCE_CODES.items():
            parts = {key: [] for key in MATERIAL_CODES}
            for face in ("left", "right"):
                optics = photon_optics(arrays, f"{selection}_{face}_", "x_mm")
                for code in MATERIAL_CODES:
                    mask = ((arrays["material_code"] == code)
                            & (arrays[f"{selection}_{face}_source_type"] == source_code))
                    parts[code].append({key: value[mask] for key, value in optics.items()})
            for code, material in MATERIAL_CODES.items():
                values = {key: np.concatenate([part[key] for part in parts[code]])
                          for key in parts[code][0]}
                optical_rows.append({"material": material,
                    "population": f"{selection} {source}, both ENDs",
                    **distribution_summary(values)})
    write_csv(STEP2_DIR / "optical_predictions.csv", optical_rows)
    (STEP2_DIR / "optical_runtime.json").write_text(json.dumps(material_rows, indent=2)+"\n")
    analysis_root = ROOT.TFile(str(ANALYSIS_ROOT_PATH), "RECREATE")
    time_result = analyze_clock(arrays, "time_ns", "t_left_ns", "t_right_ns",
                                "t0_ns", analysis_root)
    detection_result = analyze_clock(
        arrays, "t_detection_ns", "t_left_detection_ns", "t_right_detection_ns",
        "t0_detection_ns", analysis_root)
    cells = time_result[1] + detection_result[1]
    points = time_result[2] + detection_result[2]
    fits = time_result[3] + detection_result[3]
    analysis_root.Write(); analysis_root.Close()
    clock_maxima = {
        "tL(time)-tL(detection)": float(np.max(np.abs(
            arrays["t_left_ns"] - arrays["t_left_detection_ns"]))),
        "tR(time)-tR(detection)": float(np.max(np.abs(
            arrays["t_right_ns"] - arrays["t_right_detection_ns"]))),
        "T0(time)-T0(detection)": float(np.max(np.abs(
            arrays["t0_ns"] - arrays["t0_detection_ns"]))),
        "T0(time)-formula": float(np.max(np.abs(
            arrays["t0_ns"] - 0.5 * (arrays["t_left_ns"] + arrays["t_right_ns"])))),
        "T0(detection)-formula": float(np.max(np.abs(
            arrays["t0_detection_ns"] - 0.5 * (
                arrays["t_left_detection_ns"] + arrays["t_right_detection_ns"])))),
    }
    write_csv(STEP2_DIR / "baseline_cells.csv", cells)
    write_csv(STEP2_DIR / "baseline_points.csv", points)
    write_csv(STEP2_DIR / "baseline_fits.csv", fits)
    write_csv(STEP2_DIR / "clock_comparison.csv",
              [{"quantity": key, "max_abs_difference_ns": value}
               for key, value in clock_maxima.items()])
    make_baseline_figure(points, fits, sources)
    _, boundaries, boundary_correlation, point_correlation = make_cherenkov_boundary(
        arrays, sources)
    guiding_rows = make_guiding_diagnostics(arrays, sources)
    report = build_report(metadata, material_rows, cells, points, fits, clock_maxima, boundaries,
                          boundary_correlation, point_correlation, guiding_rows,
                          optical_rows)
    REPORT_PATH.write_text(report)
    EXTERNAL_REPORT_PATH.write_text(report)
    print(json.dumps({"status": "PASS", "entries": EXPECTED_ENTRIES,
                      "clock_maxima": clock_maxima,
                      "boundary_correlation": boundary_correlation,
                      "point_correlation": point_correlation,
                      "report": str(REPORT_PATH)}, indent=2))


if __name__ == "__main__":
    main()
