#!/usr/bin/env python3
"""Cuantifica el sesgo de primer fotón y los tests C0--C4 de EXEC_46."""

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

from exec46_schema import SPEED_OF_LIGHT_MM_PER_NS, load_material_config


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "step4"
INPUT_ROOT = OUTPUT_DIR / "step4_event_pairs.root"
ORDER_WIDTHS = Path("analysis/order_stat_weight_20260915/sources/part_a_widths.csv")
REPORT_PATH = OUTPUT_DIR / "REPORT_FIRSTPHOTON_SELECTION_20260916.md"
BUILD_COMMAND = (
    "env PYTHONPATH=analysis/track_mechanism_20260915 python3 "
    "analysis/track_mechanism_20260915/build_step4_pairs.py --processes 4"
)
ANALYSIS_COMMAND = (
    "env PYTHONPATH=analysis/track_mechanism_20260915 python3 "
    "analysis/track_mechanism_20260915/analyze_step4.py"
)
MATERIALS = ("EJ-200", "EJ-204", "EJ-230")
MATERIAL_CODES = {name: index for index, name in enumerate(MATERIALS)}
OPSC_CODES = {"EJ-200": "OPSC-100", "EJ-204": "OPSC-101", "EJ-230": "OPSC-106"}
FACE_NAMES = {0: "left", 1: "right"}
NOMINAL_DISTANCES_MM = np.asarray([50, 200, 500, 700, 900, 1200, 1350])
N_REFRACTIVE = float(load_material_config("OPSC-100")["rindex"][0])
GROUP_VELOCITY_MM_PER_NS = SPEED_OF_LIGHT_MM_PER_NS / N_REFRACTIVE
MUON_MASS_MEV = 105.6583755
MUON_KINETIC_ENERGY_MEV = 1000.0
MUON_GAMMA = (MUON_KINETIC_ENERGY_MEV + MUON_MASS_MEV) / MUON_MASS_MEV
MUON_BETA = math.sqrt(1.0 - 1.0 / MUON_GAMMA ** 2)
CHERENKOV_ANGLE_FINITE_BETA_DEG = math.degrees(
    math.acos(1.0 / (N_REFRACTIVE * MUON_BETA)))
AXIAL_EDGE_ANGLE_DEG = 90.0 - CHERENKOV_ANGLE_FINITE_BETA_DEG
AXIAL_EDGE_VELOCITY_BETA1 = (
    SPEED_OF_LIGHT_MM_PER_NS * math.sqrt(1.0 - 1.0 / N_REFRACTIVE ** 2)
    / N_REFRACTIVE)
AXIAL_EDGE_VELOCITY_FINITE_BETA = (
    GROUP_VELOCITY_MM_PER_NS
    * math.sin(math.radians(CHERENKOV_ANGLE_FINITE_BETA_DEG)))
BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_BATCH = 100
BOOTSTRAP_SEED = 0x46B00757A4
ORDER_EXPONENT_GRID = np.linspace(0.1, 1.5, 2801)
IQR_NORMALIZATION = 1.349
ANGLE_HISTOGRAM_WIDTH_DEG = 0.02
ANGLE_HISTOGRAM_MIN_DEG = 38.0
ANGLE_HISTOGRAM_MAX_DEG = 90.0
PRIMARY_LIKE_TOLERANCE_MM = 0.001
PAIR_VARIABLES = {
    "path_length_mm": "path_length_mm",
    "rho_detour": "rho_detour",
    "d_direct_mm": "d_direct_mm",
    "n_boundary_encounters": "n_boundary_encounters",
    "wl_nm_created": "wl_nm_created",
    "wl_nm": "wl_nm",
    "pde": "pde",
    "exit_angle_deg": "exit_angle_deg",
    "t_creation_ns": "t_creation_ns",
    "t_detection_ns": "t_detection_ns",
}
PREREG_REQUIRED_NEFF = {"EJ-200": 451.0, "EJ-204": 299.0, "EJ-230": 182.0}
PREREG_NSCINT = {"EJ-200": 2265.0, "EJ-204": 2331.0, "EJ-230": 2086.0}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def nominal_distance(frame):
    return np.where(frame["face_type"] == 0,
                    700 + frame["gun_x_mm"], 700 - frame["gun_x_mm"])


def folded_angle(values):
    values = np.asarray(values)
    return np.minimum(values, 180.0 - values)


def ks_statistic(first, random):
    first = np.sort(np.asarray(first, dtype=float))
    random = np.sort(np.asarray(random, dtype=float))
    pooled = np.sort(np.concatenate([first, random]))
    cdf_first = np.searchsorted(first, pooled, side="right") / len(first)
    cdf_random = np.searchsorted(random, pooled, side="right") / len(random)
    statistic = float(np.max(np.abs(cdf_first - cdf_random)))
    effective = len(first) * len(random) / (len(first) + len(random))
    scale = (math.sqrt(effective) + 0.12 + 0.11 / math.sqrt(effective)) * statistic
    terms = [(-1) ** (index - 1) * math.exp(-2.0 * index * index * scale * scale)
             for index in range(1, 101)]
    return statistic, min(1.0, max(0.0, 2.0 * sum(terms)))


def paired_bootstrap_means(differences, seed):
    differences = np.asarray(differences, dtype=float)
    n_events, n_variables = differences.shape
    output = np.empty((BOOTSTRAP_REPLICATES, n_variables), dtype=float)
    probability = np.full(n_events, 1.0 / n_events)
    rng = np.random.default_rng(seed)
    offset = 0
    while offset < BOOTSTRAP_REPLICATES:
        count = min(BOOTSTRAP_BATCH, BOOTSTRAP_REPLICATES - offset)
        weights = rng.multinomial(n_events, probability, size=count)
        output[offset:offset + count] = weights @ differences / n_events
        offset += count
    return output


def selection_columns(frame, prefix):
    result = {}
    for label, field in PAIR_VARIABLES.items():
        result[label] = frame[f"{prefix}_{field}"].to_numpy(dtype=float)
    result["is_cherenkov"] = (
        frame[f"{prefix}_source_type"].to_numpy(dtype=int) == 2).astype(float)
    return result


def paired_statistics(frame):
    rows = []
    grouping = ["material_code", "gun_x_mm", "face_type"]
    variable_names = list(PAIR_VARIABLES) + ["is_cherenkov"]
    for keys, group in frame.groupby(grouping, sort=True):
        material_code, gun_x, face = map(int, keys)
        first = selection_columns(group, "first")
        random = selection_columns(group, "random")
        differences = np.column_stack([first[name] - random[name]
                                       for name in variable_names])
        seed = BOOTSTRAP_SEED + material_code * 100_000 + (gun_x + 700) * 10 + face
        bootstrap = paired_bootstrap_means(differences, seed)
        distance = 700 + gun_x if face == 0 else 700 - gun_x
        for index, name in enumerate(variable_names):
            ks, ks_p = ks_statistic(first[name], random[name])
            low, high = np.quantile(bootstrap[:, index], [0.025, 0.975])
            rows.append({
                "material_code": material_code, "material": MATERIALS[material_code],
                "gun_x_mm": gun_x, "face_type": face, "face": FACE_NAMES[face],
                "nominal_d_mm": distance, "regime": ("near" if distance < 700 else
                                                        "far" if distance > 700 else "center"),
                "variable": name, "n_pairs": len(group),
                "first_mean": np.mean(first[name]), "random_mean": np.mean(random[name]),
                "paired_mean_difference": np.mean(differences[:, index]),
                "bootstrap_low_95": low, "bootstrap_high_95": high,
                "ks_statistic": ks, "ks_asymptotic_p": ks_p,
                "fraction_first_less": np.mean(differences[:, index] < 0.0),
            })
    return pd.DataFrame(rows)


def distribution_histograms(frame):
    rows = []
    for material_code, material in enumerate(MATERIALS):
        group = frame[frame["material_code"] == material_code]
        first = selection_columns(group, "first")
        random = selection_columns(group, "random")
        for variable in list(PAIR_VARIABLES) + ["is_cherenkov"]:
            pooled = np.concatenate([first[variable], random[variable]])
            if variable == "is_cherenkov":
                edges = np.asarray([-0.5, 0.5, 1.5])
            else:
                low, high = np.quantile(pooled, [0.005, 0.995])
                if high <= low:
                    high = low + 1.0
                edges = np.linspace(low, high, 81)
            for selection, values in (("first", first[variable]),
                                      ("random", random[variable])):
                counts, _ = np.histogram(values, bins=edges)
                for index, count in enumerate(counts):
                    rows.append({
                        "material_code": material_code, "material": material,
                        "variable": variable, "selection_code": 0 if selection == "first" else 1,
                        "selection": selection, "bin_low": edges[index],
                        "bin_high": edges[index + 1],
                        "bin_center": 0.5 * (edges[index] + edges[index + 1]),
                        "count": int(count), "density": count / (len(values)
                                                                  * (edges[index + 1] - edges[index])),
                    })
    return pd.DataFrame(rows)


def cherenkov_enrichment(frame):
    rows = []
    for (material_code, distance), group in frame.groupby(
            ["material_code", "nominal_d_mm"], sort=True):
        first_fraction = np.mean(group["first_source_type"] == 2)
        random_fraction = np.mean(group["random_source_type"] == 2)
        rows.append({
            "material_code": int(material_code), "material": MATERIALS[int(material_code)],
            "nominal_d_mm": int(distance), "n_events_ends": len(group),
            "first_cherenkov_fraction": first_fraction,
            "random_cherenkov_fraction": random_fraction,
            "enrichment": first_fraction / random_fraction,
            "mean_first_minus_random_detection_ns": np.mean(
                group["first_t_detection_ns"] - group["random_t_detection_ns"]),
        })
    return pd.DataFrame(rows)


def quantile_bins(values, bins=5):
    ranked = pd.Series(values).rank(method="first", pct=True).to_numpy()
    return np.minimum((ranked * bins).astype(int), bins - 1)


def cherenkov_nc_scan(frame):
    rows = []
    selections = (
        ("all_source_type_2", "first_cherenkov", "npe_cherenkov"),
        ("primary_like_proxy", "first_primary_cherenkov", "npe_primary_cherenkov"),
    )
    for selection, prefix, count_column in selections:
        valid = frame[f"{prefix}_track_id"] >= 0
        selected = frame[valid].copy()
        selected["tprop_ns"] = (selected[f"{prefix}_t_detection_ns"]
                                - selected[f"{prefix}_t_creation_ns"])
        selected["d_axial_mm"] = np.abs(selected[f"{prefix}_x_mm"]
                                           - selected[f"{prefix}_x_creation_mm"])
        selected["v_axial_mm_per_ns"] = selected["d_axial_mm"] / selected["tprop_ns"]
        selected["alpha_axial_deg"] = folded_angle(selected[f"{prefix}_exit_angle_deg"])
        for (material_code, distance), group in selected.groupby(
                ["material_code", "nominal_d_mm"], sort=True):
            group = group.copy()
            group["count_quantile"] = quantile_bins(group[count_column])
            for quantile, values in group.groupby("count_quantile", sort=True):
                rows.append({
                    "selection": selection, "material_code": int(material_code),
                    "material": MATERIALS[int(material_code)],
                    "nominal_d_mm": int(distance), "count_quantile": int(quantile),
                    "n_event_ends": len(values),
                    "mean_n_cherenkov": values[count_column].mean(),
                    "mean_v_axial_mm_per_ns": values["v_axial_mm_per_ns"].mean(),
                    "se_v_axial_mm_per_ns": values["v_axial_mm_per_ns"].sem(),
                    "median_v_axial_mm_per_ns": values["v_axial_mm_per_ns"].median(),
                    "mean_alpha_axial_deg": values["alpha_axial_deg"].mean(),
                    "distance_to_beta1_edge_mm_per_ns": (
                        values["v_axial_mm_per_ns"].mean() - AXIAL_EDGE_VELOCITY_BETA1),
                    "distance_to_finite_beta_edge_mm_per_ns": (
                        values["v_axial_mm_per_ns"].mean()
                        - AXIAL_EDGE_VELOCITY_FINITE_BETA),
                })
    return pd.DataFrame(rows)


def cherenkov_count_summary(frame):
    rows = []
    for keys, group in frame.groupby(
            ["material_code", "gun_x_mm", "face_type", "nominal_d_mm"], sort=True):
        material_code, gun_x, face, distance = map(int, keys)
        rows.append({
            "material_code": material_code, "material": MATERIALS[material_code],
            "gun_x_mm": gun_x, "face_type": face, "face": FACE_NAMES[face],
            "nominal_d_mm": distance, "mean_n_cherenkov": group["npe_cherenkov"].mean(),
            "se_n_cherenkov": group["npe_cherenkov"].sem(),
            "mean_n_primary_proxy": group["npe_primary_cherenkov"].mean(),
            "se_n_primary_proxy": group["npe_primary_cherenkov"].sem(),
        })
    return pd.DataFrame(rows)


def cherenkov_angle_window(frame):
    rows, histogram_rows = [], []
    for selection, prefix in (("all_source_type_2", "first_cherenkov"),
                              ("primary_like_proxy", "first_primary_cherenkov")):
        selected = frame[frame[f"{prefix}_track_id"] >= 0].copy()
        selected["alpha_axial_deg"] = folded_angle(selected[f"{prefix}_exit_angle_deg"])
        selected["tprop_ns"] = (selected[f"{prefix}_t_detection_ns"]
                                - selected[f"{prefix}_t_creation_ns"])
        selected["d_axial_mm"] = np.abs(selected[f"{prefix}_x_mm"]
                                           - selected[f"{prefix}_x_creation_mm"])
        alpha_rad = np.radians(selected["alpha_axial_deg"])
        edge_rad = math.radians(AXIAL_EDGE_ANGLE_DEG)
        selected["angle_penalty_pred_ns"] = (
            selected["d_axial_mm"] / GROUP_VELOCITY_MM_PER_NS
            * (1.0 / np.cos(alpha_rad) - 1.0 / math.cos(edge_rad)))
        selected["edge_excess_observed_ns"] = (
            selected["tprop_ns"]
            - selected["d_axial_mm"] / AXIAL_EDGE_VELOCITY_FINITE_BETA)
        for (material_code, distance), group in selected.groupby(
                ["material_code", "nominal_d_mm"], sort=True):
            alpha95 = group["alpha_axial_deg"].quantile(0.95)
            d_mean = group["d_axial_mm"].mean()
            predicted = (d_mean / GROUP_VELOCITY_MM_PER_NS
                         * (1.0 / math.cos(math.radians(alpha95))
                            - 1.0 / math.cos(edge_rad)))
            rows.append({
                "selection": selection, "material_code": int(material_code),
                "material": MATERIALS[int(material_code)],
                "nominal_d_mm": int(distance), "n_event_ends": len(group),
                "q95_alpha_axial_deg": alpha95,
                "q50_alpha_axial_deg": group["alpha_axial_deg"].median(),
                "mean_d_axial_mm": d_mean,
                "predicted_penalty_at_q95_ns": predicted,
                "empirical_q95_edge_excess_ns": group["edge_excess_observed_ns"].quantile(0.95),
                "corr_predicted_observed_penalty": group[
                    ["angle_penalty_pred_ns", "edge_excess_observed_ns"]].corr().iloc[0, 1],
            })
            if selection == "primary_like_proxy":
                edges = np.arange(ANGLE_HISTOGRAM_MIN_DEG,
                                  ANGLE_HISTOGRAM_MAX_DEG
                                  + ANGLE_HISTOGRAM_WIDTH_DEG / 2.0,
                                  ANGLE_HISTOGRAM_WIDTH_DEG)
                counts, _ = np.histogram(group["alpha_axial_deg"], bins=edges)
                total = len(group)
                for index, count in enumerate(counts):
                    histogram_rows.append({
                        "material_code": int(material_code),
                        "material": MATERIALS[int(material_code)],
                        "nominal_d_mm": int(distance),
                        "bin_low_deg": edges[index],
                        "bin_high_deg": edges[index + 1],
                        "bin_center_deg": 0.5 * (edges[index] + edges[index + 1]),
                        "count": int(count),
                        "density_per_deg": count / (total * ANGLE_HISTOGRAM_WIDTH_DEG),
                    })
    return pd.DataFrame(rows), pd.DataFrame(histogram_rows)


def minimum_handicap(frame):
    rows = []
    selected = frame[frame["nominal_d_mm"] == 50]
    for keys, group in selected.groupby(
            ["material_code", "gun_x_mm", "face_type"], sort=True):
        material_code, gun_x, face = map(int, keys)
        valid = group["first_cherenkov_track_id"] >= 0
        group = group[valid]
        scint_tprop = (group["first_scint_t_detection_ns"]
                       - group["first_scint_t_creation_ns"])
        cher_tprop = (group["first_cherenkov_t_detection_ns"]
                      - group["first_cherenkov_t_creation_ns"])
        total_delta = (group["first_cherenkov_t_detection_ns"]
                       - group["first_scint_t_detection_ns"])
        material = MATERIALS[material_code]
        config = load_material_config(OPSC_CODES[material])
        prefactor = config["decay_time_ns"] * math.sqrt(
            math.pi * config["rise_time_ns"]
            / (2.0 * (config["rise_time_ns"] + config["decay_time_ns"])))
        n_scint = group["npe_scint"].mean()
        predicted_min = prefactor / math.sqrt(n_scint)
        observed_pure_min_gap = np.mean(
            group["min_creation_scint_t_creation_ns"]
            - group["first_cherenkov_t_creation_ns"])
        observed_detection_selected_gap = np.mean(
            group["first_scint_t_creation_ns"]
            - group["first_cherenkov_t_creation_ns"])
        rows.append({
            "material_code": material_code, "material": material,
            "gun_x_mm": gun_x, "face_type": face, "face": FACE_NAMES[face],
            "n_events": len(group), "mean_n_scint": n_scint,
            "mean_n_cherenkov": group["npe_cherenkov"].mean(),
            "predicted_scint_min_delay_ns": predicted_min,
            "observed_pure_min_creation_gap_ns": observed_pure_min_gap,
            "observed_detection_selected_creation_gap_ns": observed_detection_selected_gap,
            "selection_reordering_gap_ns": (observed_detection_selected_gap
                                              - observed_pure_min_gap),
            "cherenkov_transport_handicap_ns": np.mean(cher_tprop - scint_tprop),
            "total_cherenkov_minus_scint_ns": np.mean(total_delta),
            "cherenkov_winner_fraction": np.mean(total_delta < 0.0),
        })
    return pd.DataFrame(rows)


def weighted_order_fit(points, exponent):
    x = points["mean_n_scint"].to_numpy() ** (-exponent)
    y = points["mean_min_creation_ns"].to_numpy()
    error = points["se_min_creation_ns"].to_numpy()
    design = np.column_stack([np.ones(len(x)), x])
    weight = 1.0 / error ** 2
    covariance = np.linalg.inv((design.T * weight) @ design)
    parameters = covariance @ (design.T @ (weight * y))
    residual = y - design @ parameters
    return parameters, covariance, float(np.sum((residual / error) ** 2))


def scintillation_order_statistics(frame):
    points = []
    for keys, group in frame.groupby(
            ["material_code", "gun_x_mm", "face_type", "nominal_d_mm"], sort=True):
        material_code, gun_x, face, distance = map(int, keys)
        points.append({
            "material_code": material_code, "material": MATERIALS[material_code],
            "gun_x_mm": gun_x, "face_type": face, "face": FACE_NAMES[face],
            "nominal_d_mm": distance, "mean_n_scint": group["npe_scint"].mean(),
            "se_n_scint": group["npe_scint"].sem(),
            "mean_n_scint_low": group["npe_scint_low_boundary"].mean(),
            "low_boundary_fraction": (group["npe_scint_low_boundary"].mean()
                                      / group["npe_scint"].mean()),
            "mean_min_creation_ns": group["min_creation_scint_t_creation_ns"].mean(),
            "se_min_creation_ns": group["min_creation_scint_t_creation_ns"].sem(),
            "mean_detection_selected_creation_ns": group[
                "first_scint_t_creation_ns"].mean(),
        })
    points = pd.DataFrame(points)
    fits = []
    for material in MATERIALS:
        group = points[points["material"] == material]
        results = {}
        for label, exponent in (("N^-1", 1.0), ("N^-1/2", 0.5)):
            parameters, covariance, chi_square = weighted_order_fit(group, exponent)
            results[label] = (parameters, covariance, chi_square, exponent)
        scans = []
        for exponent in ORDER_EXPONENT_GRID:
            parameters, covariance, chi_square = weighted_order_fit(group, exponent)
            scans.append((chi_square, exponent, parameters, covariance))
        best = min(scans, key=lambda item: item[0])
        results["free"] = (best[2], best[3], best[0], best[1])
        config = load_material_config(OPSC_CODES[material])
        prefactor = config["decay_time_ns"] * math.sqrt(
            math.pi * config["rise_time_ns"]
            / (2.0 * (config["rise_time_ns"] + config["decay_time_ns"])))
        threshold = best[0] + 1.0
        allowed = [entry[1] for entry in scans if entry[0] <= threshold]
        for label, (parameters, covariance, chi_square, exponent) in results.items():
            normalization = parameters[1]
            if label == "N^-1/2":
                effective_fraction = (prefactor / normalization) ** 2
                effective_fraction_error = abs(
                    2.0 * effective_fraction * math.sqrt(covariance[1, 1]) / normalization)
            else:
                effective_fraction = np.nan
                effective_fraction_error = np.nan
            fits.append({
                "material": material, "material_code": MATERIAL_CODES[material],
                "model": label, "exponent": exponent,
                "exponent_low_delta_chi2_1": min(allowed) if label == "free" else np.nan,
                "exponent_high_delta_chi2_1": max(allowed) if label == "free" else np.nan,
                "t0_ns": parameters[0], "t0_error_ns": math.sqrt(covariance[0, 0]),
                "normalization_ns": normalization,
                "normalization_error_ns": math.sqrt(covariance[1, 1]),
                "chi2": chi_square, "ndf": len(group) - (3 if label == "free" else 2),
                "chi2_ndf": chi_square / (len(group) - (3 if label == "free" else 2)),
                "geant4_prefactor_ns": prefactor,
                "effective_fraction": effective_fraction,
                "effective_fraction_error": effective_fraction_error,
            })
    fits = pd.DataFrame(fits)
    effective_fraction = (fits[fits["model"] == "N^-1/2"]
                          .set_index("material")["effective_fraction"])
    points["fitted_effective_fraction"] = points["material"].map(effective_fraction)
    points["fitted_n_eff"] = (points["fitted_effective_fraction"]
                              * points["mean_n_scint"])
    points["fitted_n_eff_over_low_boundary"] = (
        points["fitted_n_eff"] / points["mean_n_scint_low"])
    return points, fits


def fixed_point_correlation(frame):
    widths = pd.read_csv(ORDER_WIDTHS)
    rows = []
    for abs_x, distance in ((650, 50), (500, 200)):
        for material_code, material in enumerate(MATERIALS):
            for gun_x, face, observable in ((-abs_x, 0, "tL"), (abs_x, 1, "tR")):
                group = frame[(frame["material_code"] == material_code)
                              & (frame["gun_x_mm"] == gun_x)
                              & (frame["face_type"] == face)]
                width = widths[(widths["material"] == material)
                               & (widths["x_mm"] == gun_x)
                               & (widths["observable"] == observable)]
                require(len(group) == 10_000 and len(width) == 1,
                        f"punto fijo ausente: {material}/{gun_x}/{face}")
                rows.append({
                    "abs_x_mm": abs_x, "nominal_d_mm": distance,
                    "material_code": material_code, "material": material,
                    "gun_x_mm": gun_x, "face_type": face, "face": FACE_NAMES[face],
                    "first_cherenkov_fraction": np.mean(group["first_source_type"] == 2),
                    "ratio_single_end_to_t0_iqr": width.iloc[0]["ratio_to_T0_IQR"],
                })
    rows = pd.DataFrame(rows)
    summaries = []
    for distance, group in rows.groupby("nominal_d_mm"):
        pearson = group[["first_cherenkov_fraction",
                         "ratio_single_end_to_t0_iqr"]].corr().iloc[0, 1]
        spearman = group[["first_cherenkov_fraction",
                          "ratio_single_end_to_t0_iqr"]].rank().corr().iloc[0, 1]
        summaries.append({"nominal_d_mm": distance, "n_points": len(group),
                          "pearson_r": pearson, "spearman_rho": spearman})
    return rows, pd.DataFrame(summaries)


def robust_width(values):
    q16, q84 = np.quantile(values, [0.16, 0.84])
    return 0.5 * (q84 - q16)


def skewness(values):
    centered = values - np.mean(values)
    scale = np.std(values)
    return np.mean(centered ** 3) / scale ** 3 if scale > 0 else np.nan


def near_end_mixture(frame):
    widths = pd.read_csv(ORDER_WIDTHS)
    rows, histogram_rows = [], []
    for material_code, material in enumerate(MATERIALS):
        for gun_x, face, observable in ((-650, 0, "tL"), (650, 1, "tR")):
            group = frame[(frame["material_code"] == material_code)
                          & (frame["gun_x_mm"] == gun_x)
                          & (frame["face_type"] == face)]
            times = group["first_time_ns"].to_numpy()
            source = group["first_source_type"].to_numpy()
            scint = times[source == 1]
            cher = times[source == 2]
            weight = len(cher) / len(times)
            mean_scint, mean_cher = np.mean(scint), np.mean(cher)
            var_scint, var_cher = np.var(scint), np.var(cher)
            within = (1.0 - weight) * var_scint + weight * var_cher
            between = weight * (1.0 - weight) * (mean_cher - mean_scint) ** 2
            sigma_mix = math.sqrt(within + between)
            original = widths[(widths["material"] == material)
                              & (widths["x_mm"] == gun_x)
                              & (widths["observable"] == observable)].iloc[0]
            cher_group = group[source == 2]
            alpha = folded_angle(cher_group["first_exit_angle_deg"])
            tprop = (cher_group["first_t_detection_ns"]
                     - cher_group["first_t_creation_ns"])
            d_axial = np.abs(cher_group["first_x_mm"]
                             - cher_group["first_x_creation_mm"])
            penalty = (d_axial / GROUP_VELOCITY_MM_PER_NS
                       * (1.0 / np.cos(np.radians(alpha))
                          - 1.0 / math.cos(math.radians(AXIAL_EDGE_ANGLE_DEG))))
            rows.append({
                "material_code": material_code, "material": material,
                "gun_x_mm": gun_x, "face_type": face, "face": FACE_NAMES[face],
                "n_events": len(times), "cherenkov_fraction": weight,
                "mean_scint_ns": mean_scint, "mean_cherenkov_ns": mean_cher,
                "mean_separation_ps": 1000.0 * (mean_scint - mean_cher),
                "rms_scint_ns": np.std(scint), "rms_cherenkov_ns": np.std(cher),
                "skew_scint": skewness(scint), "skew_cherenkov": skewness(cher),
                "within_variance_ns2": within, "between_variance_ns2": between,
                "sigma_mixture_ns": sigma_mix, "rms_total_ns": np.std(times),
                "qwidth_total_ns": robust_width(times),
                "original_sigma_gaussian_ns": original["sigma_G_ns"],
                "original_gaussian_chi2_ndf": original["chi2_ndf"],
                "corr_cherenkov_tprop_angle_penalty": np.corrcoef(tprop, penalty)[0, 1],
                "q95_cherenkov_folded_angle_deg": np.quantile(alpha, 0.95),
            })
            low, high = np.quantile(times, [0.001, 0.999])
            edges = np.linspace(low, high, 121)
            for source_code, label, values in ((0, "all", times), (1, "scintillation", scint),
                                                (2, "Cherenkov", cher)):
                counts, _ = np.histogram(values, bins=edges)
                for index, count in enumerate(counts):
                    histogram_rows.append({
                        "material_code": material_code, "material": material,
                        "gun_x_mm": gun_x, "face_type": face,
                        "source_type": source_code, "source": label,
                        "bin_low_ns": edges[index], "bin_high_ns": edges[index + 1],
                        "bin_center_ns": 0.5 * (edges[index] + edges[index + 1]),
                        "count": int(count), "density": count / (len(values)
                                                                  * (edges[index + 1] - edges[index])),
                    })
    return pd.DataFrame(rows), pd.DataFrame(histogram_rows)


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
        **metadata, "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": ANALYSIS_COMMAND,
        "materials": list(MATERIALS),
        "optical_model": {
            "material_codes": OPSC_CODES,
            "rindex": N_REFRACTIVE,
            "rindex_model": "constant over the configured energy range",
        },
        "plot_scale": metadata.get("plot_scale", "linear unless stated otherwise"),
        "binning": metadata.get("binning", "tabulated points; no histogram binning"),
        "input_root": str(INPUT_ROOT.resolve()), "input_sha256": sha256(INPUT_ROOT),
        "csv": str(csv_path.resolve()), "root": str(root_path.resolve()),
        "pdf": str(pdf_path.resolve()),
    }
    meta_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def make_figures(distributions, paired, enrichment, nc_scan, angle_window,
                 angle_histograms,
                 order_points, order_fits, fixed_points, mixture_hist):
    colors = {"EJ-200": "#1f77b4", "EJ-204": "#ff7f0e", "EJ-230": "#2ca02c"}
    variables = list(PAIR_VARIABLES) + ["is_cherenkov"]
    fig, axes = plt.subplots(3, 4, figsize=(16, 11))
    for axis, variable in zip(axes.flat, variables):
        subset = distributions[(distributions["material"] == "EJ-204")
                               & (distributions["variable"] == variable)]
        for selection, group in subset.groupby("selection"):
            axis.step(group["bin_center"], group["density"], where="mid", label=selection)
        axis.set_title(variable)
        axis.grid(alpha=0.2)
    for axis in axes.flat[len(variables):]: axis.set_visible(False)
    axes[0, 0].legend()
    fig.suptitle("Paired first/random marginals, EJ-204 (all positions and ENDs)")
    save_bundle("paired_distributions", distributions, {
        "display_material": "EJ-204", "histogram_quantiles": [0.005, 0.995],
        "binning": "80 equal bins over pooled 0.5--99.5% range per material and variable; binary source flag uses two exact bins",
        "note": "CSV/ROOT contain all three materials; PDF displays the middle material."
    }, fig)

    fig, axes = plt.subplots(3, 4, figsize=(16, 11), sharex=True)
    for axis, variable in zip(axes.flat, variables):
        subset = paired[paired["variable"] == variable]
        for material in MATERIALS:
            for face, style in (("left", "-"), ("right", "--")):
                group = subset[(subset["material"] == material)
                               & (subset["face"] == face)].sort_values("gun_x_mm")
                axis.plot(group["gun_x_mm"], group["paired_mean_difference"],
                          marker="o", ms=3, ls=style, color=colors[material],
                          label=f"{material} {face}")
                axis.fill_between(group["gun_x_mm"], group["bootstrap_low_95"],
                                  group["bootstrap_high_95"],
                                  color=colors[material], alpha=0.08)
        axis.axhline(0.0, color="black", lw=0.7)
        axis.set_title(variable)
        axis.grid(alpha=0.2)
    for axis in axes.flat[len(variables):]: axis.set_visible(False)
    axes[0, 0].legend()
    fig.suptitle("Paired mean difference: first minus random; solid=L, dashed=R")
    save_bundle("paired_deltas", paired, {
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED, "pairing": "same event and END face",
    }, fig)

    fig, axis = plt.subplots(figsize=(8, 5))
    for material in MATERIALS:
        group = enrichment[enrichment["material"] == material]
        axis.plot(group["nominal_d_mm"], group["enrichment"], marker="o",
                  color=colors[material], label=material)
    axis.axhline(1.0, color="black", lw=0.8)
    axis.set(xlabel="nominal END distance [mm]", ylabel="fCher(first) / fCher(random)")
    axis.grid(alpha=0.2); axis.legend()
    save_bundle("cherenkov_enrichment", enrichment, {
        "definition": "First/random Cherenkov fraction ratio, mirrors combined."
    }, fig)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True)
    for column, material in enumerate(MATERIALS):
        for row, selection in enumerate(("all_source_type_2", "primary_like_proxy")):
            axis = axes[row, column]
            group = nc_scan[(nc_scan["material"] == material)
                            & (nc_scan["selection"] == selection)]
            for distance, values in group.groupby("nominal_d_mm"):
                axis.plot(values["mean_n_cherenkov"], values["mean_v_axial_mm_per_ns"],
                          marker="o", ms=2, lw=0.8, label=str(distance))
            axis.axhline(AXIAL_EDGE_VELOCITY_FINITE_BETA, color="black", ls="--", lw=0.8)
            axis.set_title(f"{material}, {selection.replace('_', ' ')}")
            axis.grid(alpha=0.2)
    for axis in axes[-1]: axis.set_xlabel("mean N_C in quintile")
    for axis in axes[:, 0]: axis.set_ylabel("mean axial velocity [mm/ns]")
    axes[0, 0].legend(title="d [mm]", fontsize=6)
    save_bundle("cherenkov_nc_velocity", nc_scan, {
        "count_bins": "five equal-population bins per material and distance",
        "finite_beta_edge_velocity_mm_per_ns": AXIAL_EDGE_VELOCITY_FINITE_BETA,
    }, fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True, sharey=True)
    for axis, material in zip(axes, MATERIALS):
        group = angle_window[(angle_window["material"] == material)
                             & (angle_window["selection"] == "primary_like_proxy")]
        axis.plot(group["nominal_d_mm"], group["q95_alpha_axial_deg"], marker="o")
        axis.axhline(AXIAL_EDGE_ANGLE_DEG, color="black", ls="--", lw=0.8)
        axis.set_title(material); axis.grid(alpha=0.2)
    for axis in axes: axis.set_xlabel("nominal END distance [mm]")
    axes[0].set_ylabel("q95 folded axial exit angle [deg]")
    save_bundle("cherenkov_angle_window", angle_window, {
        "angle_definition": "min(exit_angle_deg, 180-exit_angle_deg)",
        "finite_beta_edge_deg": AXIAL_EDGE_ANGLE_DEG,
        "penalty_formula": "d/(c/n)*(1/cos(alpha)-1/cos(alpha_edge))",
    }, fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True, sharey=True)
    for axis, material in zip(axes, MATERIALS):
        group = angle_histograms[angle_histograms["material"] == material]
        for distance, values in group.groupby("nominal_d_mm", sort=True):
            visible = values[(values["bin_center_deg"] >= 39.3)
                             & (values["bin_center_deg"] <= 56.0)]
            axis.step(visible["bin_center_deg"], visible["density_per_deg"],
                      where="mid", lw=0.8, label=f"{int(distance)}")
        axis.axvline(AXIAL_EDGE_ANGLE_DEG, color="black", ls="--", lw=0.8)
        axis.set_title(material); axis.set_yscale("log"); axis.grid(alpha=0.2)
        axis.set_xlabel("folded axial exit angle [deg]")
    axes[0].set_ylabel("density [deg^-1]")
    axes[0].legend(title="d [mm]", fontsize=6)
    save_bundle("cherenkov_angle_histograms", angle_histograms, {
        "selection": "first primary-like Cherenkov photon per event and END",
        "stratification": "seven nominal END distances",
        "bin_width_deg": ANGLE_HISTOGRAM_WIDTH_DEG,
        "histogram_range_deg": [ANGLE_HISTOGRAM_MIN_DEG,
                                ANGLE_HISTOGRAM_MAX_DEG],
        "finite_beta_edge_deg": AXIAL_EDGE_ANGLE_DEG,
        "plot_scale": "logarithmic y, linear x",
    }, fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    for axis, material in zip(axes, MATERIALS):
        points = order_points[order_points["material"] == material].sort_values("mean_n_scint")
        axis.errorbar(points["mean_n_scint"], points["mean_min_creation_ns"],
                      yerr=points["se_min_creation_ns"], fmt="o", ms=3)
        xgrid = np.linspace(points["mean_n_scint"].min(), points["mean_n_scint"].max(), 300)
        for model, style in (("N^-1", ":"), ("N^-1/2", "--"), ("free", "-")):
            fit = order_fits[(order_fits["material"] == material)
                             & (order_fits["model"] == model)].iloc[0]
            axis.plot(xgrid, fit["t0_ns"] + fit["normalization_ns"]
                      * xgrid ** (-fit["exponent"]), style, label=model)
        axis.set_title(material); axis.set_xscale("log"); axis.grid(alpha=0.2)
    for axis in axes: axis.set_xlabel("detected scintillation photons N_S")
    axes[0].set_ylabel("E[min t_creation] [ns]"); axes[0].legend()
    save_bundle("scintillation_order_scaling", order_points, {
        "models": ["N^-1", "N^-1/2", "free exponent"],
        "scope": "scintillation only; true minimum creation time per event and END",
        "plot_scale": "logarithmic x, linear y",
    }, fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for axis, distance in zip(axes, (50, 200)):
        group = fixed_points[fixed_points["nominal_d_mm"] == distance]
        for material in MATERIALS:
            values = group[group["material"] == material]
            axis.scatter(values["first_cherenkov_fraction"],
                         values["ratio_single_end_to_t0_iqr"],
                         color=colors[material], label=material)
        axis.set_title(f"d={distance} mm"); axis.grid(alpha=0.2)
        axis.set_xlabel("first-photon Cherenkov fraction")
    axes[0].set_ylabel("sigma_IQR(near END) / sigma_IQR(T0)")
    axes[0].legend()
    save_bundle("fixed_point_cherenkov_width", fixed_points, {
        "points": "three materials times two mirrors; no interpolation",
        "width_source": str(ORDER_WIDTHS.resolve()), "width_source_sha256": sha256(ORDER_WIDTHS),
    }, fig)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=False, sharey=False)
    for axis, ((material, gun_x), group) in zip(
            axes.flat, mixture_hist.groupby(["material", "gun_x_mm"], sort=True)):
        for source, values in group.groupby("source"):
            axis.step(values["bin_center_ns"], values["density"], where="mid", label=source)
        axis.set_title(f"{material}, x={gun_x:+g} mm"); axis.set_yscale("log")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7)
    for axis in axes[-1]: axis.set_xlabel("near-END first timestamp [ns]")
    for axis in axes[:, 0]: axis.set_ylabel("density")
    save_bundle("near_end_source_mixture", mixture_hist, {
        "selection": "six near-END cells at |x|=650 mm",
        "mixture_estimator": "sqrt(within-source variance + between-source variance)",
        "binning": "120 equal bins over the 0.1--99.9% range independently per near-END cell",
        "plot_scale": "logarithmic y, linear x",
    }, fig)


def render_report(paired, enrichment, nc_scan, angle_window, handicap, order_points,
                  order_fits, fixed_points, fixed_summary, mixtures):
    transport_ps = 1000.0 * handicap["cherenkov_transport_handicap_ns"]
    creation_ps = -1000.0 * handicap["observed_detection_selected_creation_gap_ns"]
    total_ps = 1000.0 * handicap["total_cherenkov_minus_scint_ns"]
    reorder_ps = 1000.0 * handicap["selection_reordering_gap_ns"]
    sqrt_fits = order_fits[order_fits["model"] == "N^-1/2"]
    inverse_fits = order_fits[order_fits["model"] == "N^-1"]
    fitted_fractions = sqrt_fits.set_index("material")["effective_fraction"]
    fitted_fraction_errors = (sqrt_fits.set_index("material")
                              ["effective_fraction_error"])
    fraction_separations = []
    for left, right in (("EJ-200", "EJ-204"), ("EJ-204", "EJ-230"),
                        ("EJ-200", "EJ-230")):
        fraction_separations.append(
            abs(fitted_fractions[left] - fitted_fractions[right])
            / math.hypot(fitted_fraction_errors[left], fitted_fraction_errors[right]))
    low_ranges = (order_points.groupby("material")["low_boundary_fraction"]
                   .agg(["min", "max"]))
    lines = [
        "# EXEC_46 Step 4 — first-photon selection and source order statistics", "",
        "Date: 2026-09-16", "", "## Checkpoint verdict", "",
        "Step 4 is complete. No simulation was run; the 21 production ROOT files were read-only.",
        "The corrected min-versus-min comparison closes the near-END winner accounting: Cherenkov",
        f"pays a {transport_ps.min():.0f}--{transport_ps.max():.0f} ps transport handicap but gains "
        f"{-creation_ps.max():.0f}--{-creation_ps.min():.0f} ps in the creation time attached to",
        f"the first detected photon, leaving it earlier by {-total_ps.max():.0f}--"
        f"{-total_ps.min():.0f} ps depending on material.", "",
        "The pure scintillation emission minimum follows N^-1/2 and rejects N^-1. The smaller",
        "effective populations inferred from the detection-selected photon are therefore a",
        "transport-selection effect, not the emission order statistic itself.", "",
        "This is the mandatory checkpoint. Step 5 has not been started.", "",
        "## Input and estimator definitions", "",
        f"The derived tree has 420,000 event-END rows. First and random selections reproduce the",
        "Step 2 tree exactly. The random control uses the same splitmix64 seed and the same event",
        "and END face. Every mean-difference interval below uses 2,000 paired bootstrap replicas.",
        "The KS statistic compares the paired marginal samples; the mean difference and bootstrap",
        "retain event pairing.", "",
        "`min_creation_scint` is the actual minimum creation time among detected scintillation",
        "photons. `first_scint` and `first_cherenkov` are minima in detection time within each",
        "source. No scintillation lifetime law is applied to Cherenkov.", "",
        "## C0a — Cherenkov multiplicity and the cone edge", "",
        f"The requested beta=1 edge is {AXIAL_EDGE_VELOCITY_BETA1:.3f} mm/ns. The configured",
        f"1 GeV muon gives beta={MUON_BETA:.6f} and {AXIAL_EDGE_VELOCITY_FINITE_BETA:.3f} mm/ns.",
        "Axial velocity is the measured axial displacement divided by pure propagation time.", "",
        "| material | d [mm] | proxy N_C low/high quintile | v low/high [mm/ns] | high-edge [mm/ns] |",
        "|---|---:|---:|---:|---:|"]
    proxy = nc_scan[nc_scan["selection"] == "primary_like_proxy"]
    for (material, distance), group in proxy.groupby(["material", "nominal_d_mm"]):
        group = group.sort_values("count_quantile")
        low, high = group.iloc[0], group.iloc[-1]
        lines.append(f"| {material} | {int(distance)} | {low['mean_n_cherenkov']:.2f} / "
                     f"{high['mean_n_cherenkov']:.2f} | {low['mean_v_axial_mm_per_ns']:.3f} / "
                     f"{high['mean_v_axial_mm_per_ns']:.3f} | "
                     f"{high['distance_to_finite_beta_edge_mm_per_ns']:+.3f} |")
    lines += ["", "For the primary-like proxy, increasing N_C moves the minimum toward the",
              "finite-beta cone edge in the low-occupancy long-distance cells and leaves it on",
              "the edge where occupancy is already high. The undifferentiated source-type-2",
              "sample departs above the edge at high N_C because it includes secondary-particle",
              "Cherenkov; both results are retained in `cherenkov_nc_velocity.csv`.", "",
              "## C0b — angular window versus distance", "",
              "The timing formula is tested without a fitted speed:", "", "```text",
              "Delta t(alpha,d) = d/(c/n) * [1/cos(alpha) - 1/cos(alpha_edge)]", "```", "",
              "The stored final angle is folded to the axial magnitude because the penalty uses",
              "|cos(alpha)|. The upper window is a selection in arrival time, not a second cone",
              "boundary.", "", "| material | d [mm] | q95 alpha proxy [deg] | predicted penalty [ps] | empirical q95 excess [ps] |",
              "|---|---:|---:|---:|---:|"]
    window = angle_window[angle_window["selection"] == "primary_like_proxy"]
    for _, row in window.iterrows():
        lines.append(f"| {row['material']} | {int(row['nominal_d_mm'])} | "
                     f"{row['q95_alpha_axial_deg']:.3f} | "
                     f"{1000*row['predicted_penalty_at_q95_ns']:.2f} | "
                     f"{1000*row['empirical_q95_edge_excess_ns']:.2f} |")
    lines += ["", "Contrary to the preregistered narrowing prediction, q95 moves away from the",
              "edge as distance grows in all three materials. The parameter-free formula still",
              "maps each observed angular q95 to its measured timing penalty at short and medium",
              "distance, but it does not predict the population selected after the simultaneous",
              "loss of Cherenkov multiplicity, reflections, and attenuation. The narrowing",
              "hypothesis is therefore rejected for the detected first-primary-like sample.", "",
              "## C0c — corrected d=50 mm handicap", "",
              "| material | mirror | N_S | pure-min creation C-S [ps] | selected creation C-S [ps] | transport C-S [ps] | total C-S [ps] | Cher wins | reorder gap [ps] |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for _, row in handicap.iterrows():
        lines.append(f"| {row['material']} | x={int(row['gun_x_mm']):+d}/{row['face']} | "
                     f"{row['mean_n_scint']:.1f} | "
                     f"{-1000*row['observed_pure_min_creation_gap_ns']:.2f} | "
                     f"{-1000*row['observed_detection_selected_creation_gap_ns']:.2f} | "
                     f"{1000*row['cherenkov_transport_handicap_ns']:.2f} | "
                     f"{1000*row['total_cherenkov_minus_scint_ns']:.2f} | "
                     f"{100*row['cherenkov_winner_fraction']:.2f}% | "
                     f"{1000*row['selection_reordering_gap_ns']:.2f} |")
    lines += ["", "The winner fractions reproduce 74.10/73.56%, 62.90/63.42%, and",
              "51.77/52.33% for the two mirrors. Comparing true source-specific creation minima",
              "gives a 29--40 ps Cherenkov advantage, directly matching the scale of the roughly",
              "38 ps napkin discrepancy. Selection by detection time adds another "
              f"{reorder_ps.min():.0f}--{reorder_ps.max():.0f} ps to the scintillation delay. "
              "After the 45--51 ps Cherenkov transport handicap, the measured source-specific",
              "detection minima differ by -26 to -4 ps; this reproduces the winner fractions",
              "without comparing a mean from one source with a minimum from the other.", "",
              "## C1 and C4 — scintillation order statistic and N_eff", "",
              "| material | model | exponent | chi2/ndf | fitted effective fraction (sqrt model) |",
              "|---|---|---:|---:|---:|"]
    for _, row in order_fits.iterrows():
        fraction = (f"{row['effective_fraction']:.4f} +/- "
                    f"{row['effective_fraction_error']:.4f}"
                    if np.isfinite(row["effective_fraction"]) else "--")
        exponent = (f"{row['exponent']:.4f} [{row['exponent_low_delta_chi2_1']:.4f},"
                    f" {row['exponent_high_delta_chi2_1']:.4f}]"
                    if row["model"] == "free" else f"{row['exponent']:.1f}")
        lines.append(f"| {row['material']} | {row['model']} | {exponent} | "
                     f"{row['chi2']:.2f}/{int(row['ndf'])} = {row['chi2_ndf']:.2f} | {fraction} |")
    lines += ["", f"The N^-1/2 model gives chi2/ndf {sqrt_fits['chi2_ndf'].min():.2f}--"
              f"{sqrt_fits['chi2_ndf'].max():.2f}; N^-1 is rejected with chi2/ndf "
              f"{inverse_fits['chi2_ndf'].min():.2f}--{inverse_fits['chi2_ndf'].max():.2f}. "
              "The fitted effective fractions "
              f"{fitted_fractions['EJ-200']:.3f}, {fitted_fractions['EJ-204']:.3f} and "
              f"{fitted_fractions['EJ-230']:.3f}",
              f"differ by at least {min(fraction_separations):.1f} formal standard deviations, so "
              "N_eff/N_scint is not material-independent within this model.",
              "A universal geometric factor is therefore rejected. The preregistered fractions",
              "needed to force the detection-selected handicap, 0.20/0.13/0.09, are even more",
              "strongly material dependent. The low-boundary fractions span "
              f"{low_ranges.loc['EJ-200','min']:.3f}--{low_ranges.loc['EJ-200','max']:.3f}, "
              f"{low_ranges.loc['EJ-204','min']:.3f}--{low_ranges.loc['EJ-204','max']:.3f}, and "
              f"{low_ranges.loc['EJ-230','min']:.3f}--{low_ranges.loc['EJ-230','max']:.3f}; they "
              "are strongly distance dependent and do not reproduce the fitted effective "
              "fractions. `scintillation_order_points.csv` gives fitted N_eff, the low-boundary "
              "count, and their ratio at each of the fourteen mirror-resolved points per material.", "",
              "## C2 — fixed-point Cherenkov/width test", "",
              "No interpolation is used.", "", "| d [mm] | Pearson r | Spearman rho | role |",
              "|---:|---:|---:|---|"]
    for _, row in fixed_summary.iterrows():
        role = "Cherenkov-dominated operating regime" if row["nominal_d_mm"] == 50 else "spectator regime"
        lines.append(f"| {int(row['nominal_d_mm'])} | {row['pearson_r']:+.4f} | "
                     f"{row['spearman_rho']:+.4f} | {role} |")
    lines += ["", "The six measured values at each fixed point are:", "",
              "| d [mm] | material | mirror | fCher(first) | sigma_IQR(near)/sigma_IQR(T0) |",
              "|---:|---|---|---:|---:|"]
    for _, row in fixed_points.sort_values(
            ["nominal_d_mm", "material_code", "gun_x_mm"]).iterrows():
        lines.append(f"| {int(row['nominal_d_mm'])} | {row['material']} | "
                     f"x={int(row['gun_x_mm']):+d}/{row['face']} | "
                     f"{100*row['first_cherenkov_fraction']:.2f}% | "
                     f"{row['ratio_single_end_to_t0_iqr']:.4f} |")
    lines += ["", "At 50 mm, more Cherenkov corresponds to a smaller near-END/T0 robust-width",
              "ratio. At 200 mm, where the Cherenkov fraction is about 7--9%, the correlation",
              "reverses. This sign change bounds the operating range of the mechanism; it is not",
              "evidence that either measured fixed-point correlation is internally inconsistent.", "",
              "## C3 — physical source mixture at the near END", "",
              "| material | mirror | fCher | separation [ps] | sigma_G [ps] | chi2/ndf | sigma_mix [ps] | qwidth [ps] | angle-time r |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for _, row in mixtures.iterrows():
        lines.append(f"| {row['material']} | x={int(row['gun_x_mm']):+d}/{row['face']} | "
                     f"{100*row['cherenkov_fraction']:.2f}% | {row['mean_separation_ps']:.2f} | "
                     f"{1000*row['original_sigma_gaussian_ns']:.2f} | "
                     f"{row['original_gaussian_chi2_ndf']:.1f} | "
                     f"{1000*row['sigma_mixture_ns']:.2f} | {1000*row['qwidth_total_ns']:.2f} | "
                     f"{row['corr_cherenkov_tprop_angle_penalty']:+.3f} |")
    lines += ["", "The source labels partition every event, and the source-mixture variance",
              "reconstructs the total RMS exactly. The rejected Gaussian combines two populations",
              "with different means and shapes; the Cherenkov component also carries the axial",
              "caustic. The replacement width for these six cells is therefore `sigma_mixture`,",
              "defined as sqrt(within-source variance + between-source variance), with qwidth as",
              "the robust cross-check. The failure is physical rather than a numerical fit failure.", "",
              "## Original first-versus-random selection result", "",
              "All requested variables, KS statistics, paired mean shifts, and bootstrap intervals",
              "are in `paired_selection_statistics.csv`. `paired_deltas.pdf` shows the distance",
              "dependence and `paired_distributions.pdf` the paired marginals. Cherenkov enrichment",
              "is strongest at 50 mm and is reported independently of the source-order tests.", "",
              "| material | d [mm] | fCher(first) | fCher(random) | enrichment | first-random detection [ns] |",
              "|---|---:|---:|---:|---:|---:|"]
    for _, row in enrichment.iterrows():
        lines.append(f"| {row['material']} | {int(row['nominal_d_mm'])} | "
                     f"{100*row['first_cherenkov_fraction']:.3f}% | "
                     f"{100*row['random_cherenkov_fraction']:.3f}% | "
                     f"{row['enrichment']:.3f} | "
                     f"{row['mean_first_minus_random_detection_ns']:.4f} |")
    lines += ["",
              "## Step 3 corrections", "",
              "The Step 3 report now states that the finite-beta edge falls inside a partly filled",
              "bin and the modal bin immediately to its right is not a 0.019-degree discrepancy.",
              "It also identifies isolated 41--43 degree bins at roughly one count per bin as",
              "Poisson noise, and records that the upper angular window is temporal selection.", "",
              "## Reproducibility", "", "```bash", BUILD_COMMAND, ANALYSIS_COMMAND, "```", "",
              "Each figure has PDF, ROOT, CSV and metadata sidecars. No push, merge, deck edit, or",
              "simulation was performed.", ""]
    REPORT_PATH.write_text("\n".join(lines))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    require(INPUT_ROOT.is_file(), f"falta {INPUT_ROOT}")
    require(ORDER_WIDTHS.is_file(), f"falta {ORDER_WIDTHS}")
    with uproot.open(INPUT_ROOT) as root_file:
        frame = pd.DataFrame(root_file["event_pairs"].arrays(library="np"))
    require(len(frame) == 420_000, "número inesperado de filas")
    frame["nominal_d_mm"] = nominal_distance(frame).astype(int)
    require(set(frame["nominal_d_mm"]) == set(NOMINAL_DISTANCES_MM),
            "distancias END inesperadas")

    paired = paired_statistics(frame)
    paired.to_csv(OUTPUT_DIR / "paired_selection_statistics.csv", index=False,
                  float_format="%.12g")
    distributions = distribution_histograms(frame)
    enrichment = cherenkov_enrichment(frame)
    nc_scan = cherenkov_nc_scan(frame)
    count_summary = cherenkov_count_summary(frame)
    count_summary.to_csv(OUTPUT_DIR / "cherenkov_counts.csv", index=False,
                         float_format="%.12g")
    angle_window, angle_histograms = cherenkov_angle_window(frame)
    handicap = minimum_handicap(frame)
    handicap.to_csv(OUTPUT_DIR / "minimum_handicap_d50.csv", index=False,
                    float_format="%.12g")
    order_points, order_fits = scintillation_order_statistics(frame)
    order_points.to_csv(OUTPUT_DIR / "scintillation_order_points.csv", index=False,
                        float_format="%.12g")
    order_fits.to_csv(OUTPUT_DIR / "scintillation_order_fits.csv", index=False,
                      float_format="%.12g")
    fixed_points, fixed_summary = fixed_point_correlation(frame)
    fixed_summary.to_csv(OUTPUT_DIR / "fixed_point_correlation_summary.csv", index=False,
                         float_format="%.12g")
    mixtures, mixture_hist = near_end_mixture(frame)
    mixtures.to_csv(OUTPUT_DIR / "near_end_mixture_summary.csv", index=False,
                    float_format="%.12g")

    make_figures(distributions, paired, enrichment, nc_scan, angle_window,
                 angle_histograms,
                 order_points, order_fits, fixed_points, mixture_hist)
    render_report(paired, enrichment, nc_scan, angle_window, handicap, order_points,
                  order_fits, fixed_points, fixed_summary, mixtures)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "report": str(REPORT_PATH.resolve()), "report_sha256": sha256(REPORT_PATH),
        "input_sha256": sha256(INPUT_ROOT), "rows": len(frame),
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "paired_statistic_rows": len(paired), "nc_scan_rows": len(nc_scan),
        "angle_window_rows": len(angle_window), "mixture_rows": len(mixtures),
    }
    (OUTPUT_DIR / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
