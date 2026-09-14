#!/usr/bin/env python3
"""Diagnose material-dependent first-encounter escape in existing EXEC42 ROOTs."""
import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import uproot
from scipy.stats import chi2


N = 10000
NBOOT = 500
SEED = 43091401
NBINS = 2000
FACE_IDS = np.array([3, 4, 2, 5, 6, 7, 8], dtype=np.int32)
FACE_NAMES = ["+Z", "-Z", "-Y", "-X_SENSOR", "+X_SENSOR", "+Y_SENSOR",
              "OPEN_+Y_OR_+/-X_UNRESOLVED"]
NONSENSOR = {"+Z", "-Z", "-Y", "OPEN_+Y_OR_+/-X_UNRESOLVED"}
MATERIALS = ["EJ-200", "EJ-204", "EJ-230"]
OPSC = {"EJ-200": "opsc-100", "EJ-204": "opsc-101", "EJ-230": "opsc-106"}
ATTENUATION_MM = {"EJ-200": 3800.0, "EJ-204": 1600.0, "EJ-230": 1200.0}
POSITIONS = [-650, -500, -200, 0, 200, 500, 650]
SCENARIOS = ["world_y_30mm", "world_near_x", "world_far_x"]


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sidecars(directory, stem, columns, metadata):
    data = {name: np.asarray(values) for name, values in columns.items()}
    lengths = {len(values) for values in data.values()}
    if len(lengths) != 1:
        raise RuntimeError(f"{stem}: unequal column lengths {lengths}")
    csv_path = directory / f"{stem}.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(data)
        writer.writerows(zip(*data.values()))
    root_path = directory / f"{stem}.root"
    with uproot.recreate(root_path) as output:
        output["data"] = data
    dump(directory / f"{stem}.meta.json", {
        **metadata,
        "columns": list(data),
        "rows": next(iter(lengths)),
        "csv_sha256": sha(csv_path),
        "root_sha256": sha(root_path),
    })


def bootstrap_ratio(num, den, weights):
    boot = (weights @ num) / (weights @ den)
    return float(num.sum() / den.sum()), float(boot.std(ddof=1)), boot


def secant_quantile(hist, probability):
    """Approximate a secant quantile from uniform cosine bins."""
    total = hist.sum()
    if total == 0:
        return None
    target_cos_quantile = 1.0 - probability
    cumulative = np.cumsum(hist)
    index = int(np.searchsorted(cumulative, target_cos_quantile * total, side="left"))
    index = min(max(index, 0), NBINS - 1)
    cosine = (index + 0.5) / NBINS
    return 1.0 / cosine


def uniform_distance_survival(cosine, attenuation, width=10.0):
    cosine = np.maximum(cosine, 1.e-15)
    argument = width / (cosine * attenuation)
    return -np.expm1(-argument) / argument


def reweight_nonworld(face, cosine, x_mm, target_attenuation):
    baseline = ATTENUATION_MM["EJ-200"]
    cosine = np.maximum(cosine, 1.e-15)
    weight = np.empty(len(cosine), dtype=np.float64)
    z = (face == 0) | (face == 1)
    weight[z] = (uniform_distance_survival(cosine[z], target_attenuation) /
                 uniform_distance_survival(cosine[z], baseline))
    distances = {2: 30.0, 3: 700.0 + x_mm, 4: 700.0 - x_mm, 5: 30.0}
    delta = 1.0 / target_attenuation - 1.0 / baseline
    for face_index, distance in distances.items():
        selected = face == face_index
        weight[selected] = np.exp(-distance * delta / cosine[selected])
    return weight


def read_optical_properties(sslg4, exec42_results):
    rows = []
    for material in MATERIALS:
        directory = sslg4 / "data/oscnt" / OPSC[material]
        spectrum_path = directory / "scntComp1.txt"
        rindex_path = directory / "rIndex.txt"
        absorption_path = directory / "absLength.txt"
        spectrum = np.loadtxt(spectrum_path)
        rindex = np.loadtxt(rindex_path)
        absorption = np.loadtxt(absorption_path)
        wavelength, intensity = spectrum[:, 0], spectrum[:, 1]
        interpolated_n = np.interp(wavelength, rindex[:, 0], rindex[:, 1])
        norm = np.trapz(intensity, wavelength)
        mean_wavelength = np.trapz(wavelength * intensity, wavelength) / norm
        n_eff = np.trapz(interpolated_n * intensity, wavelength) / norm
        geometric_escape = 1.0 - math.sqrt(1.0 - 1.0 / n_eff**2)
        cells = [cell for cell in exec42_results["cells"].values()
                 if cell["configuration"]["material"] == material]
        aggregate = [cell["V2_H1"]["value"] for cell in cells]
        z_num = z_den = 0.0
        for cell in cells:
            for face in cell["V2_faces"]:
                if face["face"] in ("+Z", "-Z"):
                    z_num += face["numerator"]
                    z_den += face["denominator"]
        rows.append({
            "material": material,
            "opsc_code": OPSC[material],
            "emission_peak_nm": float(wavelength[np.argmax(intensity)]),
            "emission_mean_nm": float(mean_wavelength),
            "n_effective": float(n_eff),
            "geometric_escape": geometric_escape,
            "measured_escape_mean_x": float(np.mean(aggregate)),
            "measured_escape_min_x": float(np.min(aggregate)),
            "measured_escape_max_x": float(np.max(aggregate)),
            "measured_Z_conditional_escape": float(z_num / z_den),
            "attenuation_mm": float(absorption[0, 1] * 10.0),
            "rindex_path": str(rindex_path),
            "rindex_sha256": sha(rindex_path),
            "spectrum_path": str(spectrum_path),
            "spectrum_sha256": sha(spectrum_path),
            "absorption_path": str(absorption_path),
            "absorption_sha256": sha(absorption_path),
        })
    return rows


def analyze_cell(row, output, common_metadata):
    cell_id = row["cell_id"]
    material = row["material"]
    x_mm = int(row["x_mm"])
    cell_dir = output / "cells" / cell_id
    cell_dir.mkdir(parents=True, exist_ok=False)
    hist = np.zeros((7, 2, NBINS), dtype=np.int64)
    encounter_by_event = np.zeros(N, dtype=np.int64)
    closure = {}
    if material == "EJ-200":
        for target in ("EJ-204", "EJ-230"):
            closure[target] = {
                scenario: {"num": np.zeros(N), "den": np.zeros(N)}
                for scenario in SCENARIOS
            }
    face_map = np.full(9, -1, dtype=np.int8)
    face_map[FACE_IDS] = np.arange(7, dtype=np.int8)
    source = uproot.open(row["root_path"])
    print(f"START {cell_id} rows={row['encounter_rows']}", flush=True)
    selected_rows = 0
    for chunk in source["first_bar_encounters"].iterate(
            ["event_id", "source", "post_volume_id", "outcome", "exiting_bar",
             "cos_incidence"], step_size="256 MB", library="np"):
        selected = (chunk["source"] == 1) & (chunk["exiting_bar"] == 1)
        event_id = chunk["event_id"][selected].astype(np.int64)
        post = chunk["post_volume_id"][selected]
        outcome = chunk["outcome"][selected]
        cosine = chunk["cos_incidence"][selected].astype(np.float64)
        if np.any(~np.isfinite(cosine)) or np.any((cosine < 0) | (cosine > 1)):
            raise RuntimeError(cell_id + ": invalid signed cosine")
        face = face_map[post]
        if np.any(face < 0):
            raise RuntimeError(cell_id + ": unmapped face")
        escaped = ((outcome == 2) &
                   (((post >= 2) & (post <= 4)) | (post == 8)))
        bins = np.minimum((cosine * NBINS).astype(np.int64), NBINS - 1)
        keys = (face.astype(np.int64) * 2 + escaped.astype(np.int64)) * NBINS + bins
        hist += np.bincount(keys, minlength=14 * NBINS).reshape(7, 2, NBINS)
        encounter_by_event += np.bincount(event_id, minlength=N)
        selected_rows += len(event_id)

        if material == "EJ-200":
            world = face == 6
            nonworld = ~world
            for target in ("EJ-204", "EJ-230"):
                common_num = np.zeros(N)
                common_den = np.zeros(N)
                nw_weight = reweight_nonworld(face[nonworld], cosine[nonworld], x_mm,
                                               ATTENUATION_MM[target])
                common_den += np.bincount(event_id[nonworld], weights=nw_weight, minlength=N)
                common_num += np.bincount(event_id[nonworld & escaped],
                                          weights=nw_weight[escaped[nonworld]], minlength=N)
                world_distances = {
                    "world_y_30mm": 30.0,
                    "world_near_x": 700.0 - abs(x_mm),
                    "world_far_x": 700.0 + abs(x_mm),
                }
                delta = 1.0 / ATTENUATION_MM[target] - 1.0 / ATTENUATION_MM["EJ-200"]
                for scenario, distance in world_distances.items():
                    world_weight = np.exp(-distance * delta /
                                          np.maximum(cosine[world], 1.e-15))
                    closure[target][scenario]["den"] += common_den
                    closure[target][scenario]["num"] += common_num
                    closure[target][scenario]["den"] += np.bincount(
                        event_id[world], weights=world_weight, minlength=N)
                    closure[target][scenario]["num"] += np.bincount(
                        event_id[world & escaped], weights=world_weight[escaped[world]],
                        minlength=N)
    if selected_rows != int(hist.sum()) or selected_rows != int(encounter_by_event.sum()):
        raise RuntimeError(cell_id + ": encounter counts do not close")

    summary = []
    hist_columns = {"face": [], "outcome": [], "cos_bin_low": [],
                    "cos_bin_high": [], "count": []}
    for face_index, face_name in enumerate(FACE_NAMES):
        for outcome_index, outcome_name in enumerate(("NON_ESCAPE", "ESCAPE")):
            counts = hist[face_index, outcome_index]
            q10 = secant_quantile(counts, 0.10)
            q50 = secant_quantile(counts, 0.50)
            q90 = secant_quantile(counts, 0.90)
            fixed_distance = None
            if face_name in ("-Y", "+Y_SENSOR"):
                fixed_distance = 30.0
            elif face_name == "-X_SENSOR":
                fixed_distance = 700.0 + x_mm
            elif face_name == "+X_SENSOR":
                fixed_distance = 700.0 - x_mm
            summary.append({
                "cell_id": cell_id, "material": material, "x_mm": x_mm,
                "face": face_name, "outcome": outcome_name, "count": int(counts.sum()),
                "secant_q10": q10, "secant_median": q50, "secant_q90": q90,
                "fixed_distance_mm": fixed_distance,
                "path_q10_mm": None if fixed_distance is None or q10 is None else fixed_distance*q10,
                "path_median_mm": None if fixed_distance is None or q50 is None else fixed_distance*q50,
                "path_q90_mm": None if fixed_distance is None or q90 is None else fixed_distance*q90,
                "absolute_path_status": ("EVALUABLE_FIXED_SOURCE_PLANE" if fixed_distance is not None
                                         else "NOT_EVALUABLE_MISSING_CREATION_POINT_OR_FACE"),
            })
            nonzero = np.nonzero(counts)[0]
            for index in nonzero:
                hist_columns["face"].append(face_name)
                hist_columns["outcome"].append(outcome_name)
                hist_columns["cos_bin_low"].append(index / NBINS)
                hist_columns["cos_bin_high"].append((index + 1) / NBINS)
                hist_columns["count"].append(int(counts[index]))
    numeric_summary = {
        key: [np.nan if item[key] is None else item[key] for item in summary]
        for key in summary[0]
    }
    metadata = {**common_metadata, "cell_id": cell_id, "material": material,
                "x_mm": x_mm, "source_root": row["root_path"],
                "source_root_sha256": row["root_sha256"],
                "cosine_bins": NBINS, "selected_first_encounters": selected_rows}
    sidecars(cell_dir, "path_proxy_summary", numeric_summary, metadata)
    sidecars(cell_dir, "path_proxy_histogram", hist_columns, metadata)
    if material == "EJ-200":
        columns = {"event_id": np.arange(N)}
        for target in ("EJ-204", "EJ-230"):
            for scenario in SCENARIOS:
                columns[f"{OPSC[target]}_{scenario}_numerator"] = closure[target][scenario]["num"]
                columns[f"{OPSC[target]}_{scenario}_denominator"] = closure[target][scenario]["den"]
        sidecars(cell_dir, "attenuation_reweight_events", columns, metadata)
    dump(cell_dir / "path_results.json", {"metadata": metadata, "summary": summary})
    print(f"DONE {cell_id} selected={selected_rows}", flush=True)
    return summary, encounter_by_event, closure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--exec42-results", type=Path, required=True)
    parser.add_argument("--sslg4", type=Path, required=True)
    parser.add_argument("--preregistration", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise RuntimeError("output directory exists; refusing to overwrite")
    args.out.mkdir(parents=True)
    (args.out / "cells").mkdir()
    script_hash = sha(Path(__file__))
    prereg_hash = sha(args.preregistration)
    with args.inventory.open() as stream:
        inventory = list(csv.DictReader(stream))
    for row in inventory:
        for key in ("x_mm", "events", "encounter_rows"):
            row[key] = int(row[key])
    inventory.sort(key=lambda row: (row["material"], row["x_mm"]))
    if len(inventory) != 21:
        raise RuntimeError("input inventory does not contain 21 cells")
    exec42 = json.loads(args.exec42_results.read_text())
    rng = np.random.default_rng(SEED)
    weights = np.stack([np.bincount(rng.integers(0, N, N), minlength=N)
                        for _ in range(NBOOT)]).astype(np.float64)
    created = dt.datetime.now(dt.timezone.utc).isoformat()
    metadata = {
        "created_utc": created, "N_cells": 21, "N_per_cell": N,
        "seeds": [26092601, 8349041], "bootstrap_seed": SEED,
        "bootstrap_replicates": NBOOT,
        "bootstrap_unit": "generated event; common weights across all cells",
        "inventory": str(args.inventory), "inventory_sha256": sha(args.inventory),
        "exec42_results": str(args.exec42_results),
        "exec42_results_sha256": sha(args.exec42_results),
        "preregistration": str(args.preregistration),
        "preregistration_sha256": prereg_hash, "script_sha256": script_hash,
    }
    schema_audit = {
        "first_bar_encounters_branches": ["event_id", "track_id", "source", "pre_copy",
            "post_copy", "pre_volume_id", "post_volume_id", "outcome", "exiting_bar",
            "cos_incidence"],
        "required_but_absent_for_exact_path": ["creation_position", "encounter_position",
                                                "global_time", "track_length"],
        "absolute_path_plus_minus_Z": "NOT_EVALUABLE",
        "reason_plus_minus_Z": "production z spans -5..+5 mm but photon creation z is absent",
        "absolute_path_WorldPV": "NOT_EVALUABLE",
        "reason_WorldPV": "volume ID does not distinguish open +Y from uncovered +/-X",
        "path_proxy": "sec(theta)=1/cos_incidence",
    }
    dump(args.out / "schema_audit.json", {"metadata": metadata, **schema_audit})

    optical = read_optical_properties(args.sslg4, exec42)
    n_spread = max(row["n_effective"] for row in optical) - min(row["n_effective"] for row in optical)
    m2_status = "EXCLUDED" if n_spread < 1.e-4 else "NOT_EXCLUDED"
    sidecars(args.out, "optical_properties", {key: [row[key] for row in optical]
        for key in optical[0]}, metadata)

    all_path_summary = []
    event_data = {}
    closure_internal = {}
    for row in inventory:
        summary, encounters, closure = analyze_cell(row, args.out, metadata)
        all_path_summary.extend(summary)
        metrics_path = Path("/home/rrios/exec42_20260913/analysis/cells") / row["cell_id"] / "event_metrics.root"
        with uproot.open(metrics_path) as source:
            arrays = source["data"].arrays(["produced_scint", "V2_encounters", "V2_escaped"],
                                             library="np")
        if not np.array_equal(encounters, arrays["V2_encounters"]):
            raise RuntimeError(row["cell_id"] + ": streamed encounters differ from EXEC42")
        produced = arrays["produced_scint"].astype(np.float64)
        met_enc = arrays["V2_encounters"].astype(np.float64)
        escaped = arrays["V2_escaped"].astype(np.float64)
        absent = produced - met_enc
        absent_value, absent_se, absent_boot = bootstrap_ratio(absent, produced, weights)
        escape_value, escape_se, escape_boot = bootstrap_ratio(escaped, met_enc, weights)
        event_data[row["cell_id"]] = {
            "material": row["material"], "x_mm": row["x_mm"], "produced": produced,
            "encounters": met_enc, "escaped": escaped, "absent_value": absent_value,
            "absent_se": absent_se, "absent_boot": absent_boot,
            "escape_value": escape_value, "escape_se": escape_se, "escape_boot": escape_boot,
        }
        closure_internal[row["cell_id"]] = closure

    path_by = {(row["cell_id"], row["face"], row["outcome"]): row
               for row in all_path_summary}
    directional_rows = []
    for row in inventory:
        for face in FACE_NAMES:
            nonescape = path_by[row["cell_id"], face, "NON_ESCAPE"]
            escape = path_by[row["cell_id"], face, "ESCAPE"]
            evaluable = nonescape["count"] > 0 and escape["count"] > 0
            status = ("PASS" if evaluable and nonescape["secant_median"] > escape["secant_median"]
                      else "FAIL" if evaluable else "NOT_EVALUABLE")
            directional_rows.append({
                "cell_id": row["cell_id"], "material": row["material"], "x_mm": row["x_mm"],
                "face": face, "n_non_escape": nonescape["count"], "n_escape": escape["count"],
                "non_escape_secant_median": nonescape["secant_median"],
                "escape_secant_median": escape["secant_median"], "status": status,
                "enters_M1_decision": face in NONSENSOR and evaluable,
            })
    sidecars(args.out, "directional_path_test", {key: [np.nan if row[key] is None else row[key]
        for row in directional_rows] for key in directional_rows[0]}, metadata)

    absent_rows = []
    for cell_id, data in sorted(event_data.items(), key=lambda item: (item[1]["material"], item[1]["x_mm"])):
        absent_rows.append({"cell_id": cell_id, "material": data["material"],
            "x_mm": data["x_mm"], "attenuation_mm": ATTENUATION_MM[data["material"]],
            "produced": int(data["produced"].sum()), "first_encounters": int(data["encounters"].sum()),
            "absent": int((data["produced"]-data["encounters"]).sum()),
            "absent_fraction": data["absent_value"], "absent_fraction_se": data["absent_se"]})
    sidecars(args.out, "absent_photons", {key: [row[key] for row in absent_rows]
        for key in absent_rows[0]}, metadata)

    by_material_x = {(data["material"], data["x_mm"]): (cell_id, data)
                     for cell_id, data in event_data.items()}
    absent_pairs = []
    for x in POSITIONS:
        for longer, shorter in [("EJ-200", "EJ-204"), ("EJ-200", "EJ-230"),
                                ("EJ-204", "EJ-230")]:
            cell_l, data_l = by_material_x[longer, x]
            cell_s, data_s = by_material_x[shorter, x]
            difference = data_s["absent_value"] - data_l["absent_value"]
            boot_difference = data_s["absent_boot"] - data_l["absent_boot"]
            se = float(boot_difference.std(ddof=1))
            adjacent = (longer, shorter) in [("EJ-200", "EJ-204"), ("EJ-204", "EJ-230")]
            absent_pairs.append({"x_mm": x, "longer_material": longer,
                "shorter_material": shorter, "longer_cell": cell_l, "shorter_cell": cell_s,
                "absent_longer": data_l["absent_value"], "absent_shorter": data_s["absent_value"],
                "difference_shorter_minus_longer": difference, "paired_se": se,
                "difference_over_se": difference/se, "status": "PASS" if difference > 3*se else "FAIL",
                "adjacent_pair_for_M1_decision": adjacent})
    sidecars(args.out, "absent_pairwise", {key: [row[key] for row in absent_pairs]
        for key in absent_pairs[0]}, metadata)

    closure_rows = []
    for x in POSITIONS:
        baseline_cell, baseline_data = by_material_x["EJ-200", x]
        for target in ("EJ-204", "EJ-230"):
            target_cell, target_data = by_material_x[target, x]
            predictions = {}
            for scenario in SCENARIOS:
                component = closure_internal[baseline_cell][target][scenario]
                value, se, boot = bootstrap_ratio(component["num"], component["den"], weights)
                predictions[scenario] = {"value": value, "se": se, "boot": boot}
            low_name = min(SCENARIOS, key=lambda name: predictions[name]["value"])
            high_name = max(SCENARIOS, key=lambda name: predictions[name]["value"])
            observed = target_data["escape_value"]
            low = predictions[low_name]
            high = predictions[high_name]
            low_residual = observed - low["value"]
            high_residual = observed - high["value"]
            low_se = float((target_data["escape_boot"] - low["boot"]).std(ddof=1))
            high_se = float((target_data["escape_boot"] - high["boot"]).std(ddof=1))
            compatible = low_residual >= -3*low_se and high_residual <= 3*high_se
            central = predictions["world_y_30mm"]
            central_residual = observed - central["value"]
            central_se = float((target_data["escape_boot"] - central["boot"]).std(ddof=1))
            closure_rows.append({
                "x_mm": x, "baseline_cell": baseline_cell, "target_cell": target_cell,
                "target_material": target, "observed_escape": observed,
                "observed_escape_se": target_data["escape_se"],
                "predicted_world_y": central["value"], "predicted_world_y_se": central["se"],
                "predicted_low": low["value"], "predicted_low_scenario": low_name,
                "predicted_high": high["value"], "predicted_high_scenario": high_name,
                "envelope_width": high["value"]-low["value"],
                "central_residual": central_residual, "central_residual_paired_se": central_se,
                "central_abs_residual_over_se": abs(central_residual)/central_se,
                "low_residual": low_residual, "low_residual_paired_se": low_se,
                "high_residual": high_residual, "high_residual_paired_se": high_se,
                "status": "PASS" if compatible else "FAIL",
            })
    sidecars(args.out, "attenuation_closure", {key: [row[key] for row in closure_rows]
        for key in closure_rows[0]}, metadata)

    position_rows = []
    for first, second in [("EJ-200", "EJ-204"), ("EJ-200", "EJ-230"),
                          ("EJ-204", "EJ-230")]:
        point = np.array([by_material_x[first, x][1]["escape_value"] -
                          by_material_x[second, x][1]["escape_value"] for x in POSITIONS])
        boot = np.column_stack([by_material_x[first, x][1]["escape_boot"] -
                                by_material_x[second, x][1]["escape_boot"] for x in POSITIONS])
        covariance = np.cov(boot, rowvar=False, ddof=1)
        inverse = np.linalg.pinv(covariance)
        ones = np.ones(len(POSITIONS))
        mean = float((ones @ inverse @ point) / (ones @ inverse @ ones))
        residual = point - mean
        statistic = float(residual @ inverse @ residual)
        rank = int(np.linalg.matrix_rank(covariance))
        dof = max(rank - 1, 1)
        p_value = float(chi2.sf(statistic, dof))
        position_rows.append({"material_a": first, "material_b": second,
            "arithmetic_mean_difference": float(point.mean()), "GLS_constant": mean,
            "minimum_difference": float(point.min()), "maximum_difference": float(point.max()),
            "peak_to_peak": float(point.max()-point.min()),
            "fractional_peak_to_peak": float((point.max()-point.min())/abs(point.mean())),
            "chi_square": statistic, "covariance_rank": rank, "dof": dof, "p_value": p_value,
            "position_dependence": "RESOLVED" if p_value < 0.01 else "NOT_RESOLVED"})
    sidecars(args.out, "position_dependence", {key: [row[key] for row in position_rows]
        for key in position_rows[0]}, metadata)

    path_decision_rows = [row for row in directional_rows
                          if row["enters_M1_decision"]]
    directional_pass = bool(path_decision_rows) and all(row["status"] == "PASS"
                                                        for row in path_decision_rows)
    adjacent_absent = [row for row in absent_pairs if row["adjacent_pair_for_M1_decision"]]
    absent_pass = len(adjacent_absent) == 14 and all(row["status"] == "PASS"
                                                     for row in adjacent_absent)
    closure_pass = len(closure_rows) == 14 and all(row["status"] == "PASS"
                                                   for row in closure_rows)
    m1_confirmed = m2_status == "EXCLUDED" and directional_pass and absent_pass and closure_pass
    if m1_confirmed:
        m1_status = "CONFIRMED_WITHIN_PERSISTED_DATA_LIMITS"
    elif directional_pass and absent_pass:
        m1_status = "SUPPORTED_BUT_NOT_FULLY_CONFIRMED"
    else:
        m1_status = "NOT_CONFIRMED"
    results = {
        "metadata": metadata,
        "schema_audit": schema_audit,
        "M2": {"status": m2_status, "n_effective_spread": n_spread,
               "threshold": 1.e-4, "optical_properties": optical},
        "M1": {"status": m1_status, "directional_path_status": "PASS" if directional_pass else "FAIL",
               "directional_comparisons_pass": sum(row["status"] == "PASS" for row in path_decision_rows),
               "directional_comparisons_total": len(path_decision_rows),
               "survival_ordering_status": "PASS" if absent_pass else "FAIL",
               "adjacent_absent_pairs_pass": sum(row["status"] == "PASS" for row in adjacent_absent),
               "adjacent_absent_pairs_total": len(adjacent_absent),
               "quantitative_closure_status": "PASS" if closure_pass else "FAIL",
               "closure_cells_pass": sum(row["status"] == "PASS" for row in closure_rows),
               "closure_cells_total": len(closure_rows),
               "exact_full_path_distribution": "NOT_EVALUABLE_FROM_PRODUCTION_SCHEMA"},
        "absent_photons": absent_rows,
        "absent_pairwise": absent_pairs,
        "directional_path_test": directional_rows,
        "attenuation_closure": closure_rows,
        "position_dependence": position_rows,
        "proposed_H3_prime": {
            "applied": False,
            "text": "First-encounter escape may depend on material only through pre-encounter bulk-absorption filtering. Each material-pair difference must agree with a preregistered reweighting based on the attenuation lengths and an independently specified path distribution within a declared tolerance.",
            "production_requirement": "Persist exact path length or creation and first-encounter coordinates before considering golden-contract adoption.",
        },
    }
    dump(args.out / "results.json", results)
    print(json.dumps({"M2": results["M2"], "M1": results["M1"],
                      "position_dependence": position_rows}, indent=2))


if __name__ == "__main__":
    main()
