#!/usr/bin/env python3
"""Evaluate the preregistered EXEC40/41/42 physical contract on 21 ROOTs."""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import uproot
from scipy.stats import chi2


N = 10000
NBOOT = 500
SEED = 42091401
LOWER_BOUND = 0.226172
FACE_IDS = np.array([3, 4, 2, 5, 6, 7, 8], dtype=np.int32)
FACE_NAMES = ["+Z", "-Z", "-Y", "-X_SENSOR", "+X_SENSOR", "+Y_SENSOR",
              "OPEN_+Y_OR_+/-X_UNRESOLVED"]
YIELDS = {"EJ-200": 10000.0, "EJ-204": 10400.0, "EJ-230": 9700.0}


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(16*1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def sidecars(directory, stem, columns, metadata):
    data = {name: np.asarray(values) for name, values in columns.items()}
    csv_path = directory/(stem+".csv")
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(data); writer.writerows(zip(*data.values()))
    root_path = directory/(stem+".root")
    with uproot.recreate(root_path) as output:
        output["data"] = data
    dump(directory/(stem+".meta.json"), {
        **metadata, "columns": list(data), "rows": len(next(iter(data.values()))),
        "csv_sha256": sha(csv_path), "root_sha256": sha(root_path),
    })


def bootstrap_ratio(num, den, weights):
    boot = (weights @ num)/(weights @ den)
    return {"value": float(num.sum()/den.sum()), "se": float(boot.std(ddof=1)),
            "numerator": float(num.sum()), "denominator": float(den.sum())}, boot


def analyze_cell(inventory, weights, output_root, prereg, script_hash):
    cell_id = inventory["cell_id"]
    cell_dir = output_root/"cells"/cell_id
    cell_dir.mkdir(parents=True, exist_ok=False)
    source_path = Path(inventory["root_path"])
    source = uproot.open(source_path)
    event = source["event_observables"].arrays(library="np")
    order = np.argsort(event["event_id"])
    event = {name: values[order] for name, values in event.items()}
    if not np.array_equal(event["event_id"], np.arange(N)):
        raise RuntimeError(cell_id+": event IDs do not close")
    sensor = source["sipm_event_counts"].arrays(library="np")
    order = np.lexsort((sensor["global_id"], sensor["event_id"]))
    sensor = {name: values[order].reshape(N, 86) for name, values in sensor.items()}
    if not np.array_equal(sensor["event_id"][:, 0], np.arange(N)):
        raise RuntimeError(cell_id+": sensor event IDs do not close")

    produced = event["produced_scint"].astype(np.float64)
    total_edep = event["edep_total_MeV"].astype(np.float64)
    optical_edep = event["edep_optical_MeV"].astype(np.float64)
    nonoptical_edep = total_edep-optical_edep
    nominal_yield = YIELDS[inventory["material"]]
    v1, v1_boot = bootstrap_ratio(produced, nominal_yield*total_edep, weights)
    v1["lower_3se"], v1["upper_3se"] = v1["value"]-3*v1["se"], v1["value"]+3*v1["se"]
    v1["nominal_yield_ph_per_MeV"] = nominal_yield
    v1["status"] = "PASS" if v1["lower_3se"] <= 1 <= v1["upper_3se"] else "FAIL"
    v1_nonopt, v1_nonopt_boot = bootstrap_ratio(produced, nominal_yield*nonoptical_edep, weights)
    v1_nonopt["status"] = "DESCRIPTIVE_SECONDARY"

    incident = sensor["incident"].sum(axis=1).astype(np.float64)
    matched = sensor["matched_detected"].sum(axis=1).astype(np.float64)
    expected = sensor["expected_surface_pde_sum"].sum(axis=1).astype(np.float64)
    detected = sensor["detected"].sum(axis=1).astype(np.float64)
    orphan = sensor["unmatched_detected"].sum(axis=1).astype(np.float64)
    v5, v5_boot = bootstrap_ratio(matched, incident, weights)
    v5_expected, expected_boot = bootstrap_ratio(expected, incident, weights)
    v5_diff_boot = v5_boot-expected_boot
    v5_diff = v5["value"]-v5_expected["value"]
    v5_diff_se = float(v5_diff_boot.std(ddof=1))
    v5.update({"expected_incident_PDE": v5_expected, "difference": v5_diff,
               "difference_se": v5_diff_se,
               "status": "PASS" if abs(v5_diff) <= 3*v5_diff_se else "FAIL"})
    orphan_fraction, orphan_boot = bootstrap_ratio(orphan, detected, weights)
    orphan_fraction["status"] = "DESCRIPTIVE_ACCOUNTING_PATH"

    encounter = np.zeros((N, 7, 2), dtype=np.int64)
    hist = np.zeros((N, 20), dtype=np.int64)
    face_map = np.full(9, -1, dtype=np.int8)
    face_map[FACE_IDS] = np.arange(7, dtype=np.int8)
    selected_rows = 0
    print(f"START {cell_id} encounter_rows={inventory['encounter_rows']}", flush=True)
    for chunk in source["first_bar_encounters"].iterate(
            ["event_id", "source", "post_volume_id", "outcome", "exiting_bar", "cos_incidence"],
            step_size="256 MB", library="np"):
        selected = (chunk["source"] == 1) & (chunk["exiting_bar"] == 1)
        event_id = chunk["event_id"][selected].astype(np.int64)
        post = chunk["post_volume_id"][selected]
        outcome = chunk["outcome"][selected]
        cosine = chunk["cos_incidence"][selected]
        selected_rows += len(event_id)
        face = face_map[post]
        if np.any(face < 0):
            raise RuntimeError(cell_id+": unmapped production volume ID")
        encounter[:, :, 0] += np.bincount(event_id*7+face, minlength=N*7).reshape(N, 7)
        escaped = (outcome == 2) & (((post >= 2) & (post <= 4)) | (post == 8))
        encounter[:, :, 1] += np.bincount(event_id[escaped]*7+face[escaped],
                                          minlength=N*7).reshape(N, 7)
        valid = np.isfinite(cosine) & (cosine >= 0) & (cosine <= 1)
        bins = np.minimum((cosine[valid]*20).astype(np.int64), 19)
        hist += np.bincount(event_id[valid]*20+bins, minlength=N*20).reshape(N, 20)
    if selected_rows != int(encounter[:, :, 0].sum()):
        raise RuntimeError(cell_id+": face counts do not close")

    v2, v2_boot = bootstrap_ratio(encounter[:, :, 1].sum(axis=1),
                                   encounter[:, :, 0].sum(axis=1), weights)
    v2["lower_3se"] = v2["value"]-3*v2["se"]
    v2["lower_bound"] = LOWER_BOUND
    v2["status"] = "PASS" if v2["lower_3se"] >= LOWER_BOUND else "FAIL"
    face_results = []
    total_encounters = encounter[:, :, 0].sum(axis=1)
    total_escapes = encounter[:, :, 1].sum()
    for index, name in enumerate(FACE_NAMES):
        conditional, _ = bootstrap_ratio(encounter[:, index, 1], encounter[:, index, 0], weights)
        share_boot = (weights @ encounter[:, index, 0])/(weights @ total_encounters)
        conditional.update({"face": name,
            "encounter_share": float(encounter[:, index, 0].sum()/total_encounters.sum()),
            "encounter_share_se": float(share_boot.std(ddof=1)),
            "escape_share": float(encounter[:, index, 1].sum()/total_escapes) if total_escapes else 0.0})
        face_results.append(conditional)

    edges = np.linspace(0, 1, 21)
    observed = hist.sum(axis=0)/hist.sum()
    boot_hist = weights @ hist
    boot_hist = boot_hist/boot_hist.sum(axis=1)[:, None]
    covariance = np.cov(boot_hist, rowvar=False, ddof=1)
    rank = int(np.linalg.matrix_rank(covariance))
    def gof(reference):
        delta = observed-reference
        statistic = float(delta @ np.linalg.pinv(covariance) @ delta)
        return {"chi_square": statistic, "rank": rank, "p_value": float(chi2.sf(statistic, rank))}
    uniform = np.full(20, .05)
    flux = np.diff(edges**2)
    angular = {"status": "DESCRIPTIVE_NO_PASS_FAIL", "uniform": gof(uniform),
               "old_p_2mu": gof(flux), "photons": int(hist.sum()),
               "low_mu_probability": float(observed[:4].sum()),
               "low_mu_uniform": .2,
               "low_mu_difference": float(observed[:4].sum()-.2),
               "low_mu_se": float(boot_hist[:, :4].sum(axis=1).std(ddof=1))}

    npe_left = sensor["detected"][:, :8].sum(axis=1)
    npe_right = sensor["detected"][:, 8:16].sum(axis=1)
    metadata = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "cell_id": cell_id, "source_root": str(source_path),
        "source_root_sha256": inventory["root_sha256"], "N": N,
        "seeds": [26092601, 8349041], "workers_from_manifest": int(inventory["workers_manifest"]),
        "bootstrap_seed": SEED, "bootstrap_replicates": NBOOT,
        "bootstrap_unit": "generated event; common weights across all cells",
        "preregistration": str(prereg), "preregistration_sha256": sha(prereg),
        "script_sha256": script_hash}
    sidecars(cell_dir, "event_metrics", {
        "event_id": np.arange(N), "produced_scint": produced, "edep_total_MeV": total_edep,
        "edep_optical_MeV": optical_edep, "edep_nonoptical_MeV": nonoptical_edep,
        "V2_encounters": encounter[:, :, 0].sum(axis=1),
        "V2_escaped": encounter[:, :, 1].sum(axis=1), "V5_incident": incident,
        "V5_matched": matched, "V5_expected_PDE_sum": expected,
        "SD_detected": detected, "SD_orphan": orphan,
        "npe_left": npe_left, "npe_right": npe_right,
    }, metadata)
    sidecars(cell_dir, "face_decomposition", {
        "face": FACE_NAMES, "encounters": [row["denominator"] for row in face_results],
        "encounter_share": [row["encounter_share"] for row in face_results],
        "encounter_share_se": [row["encounter_share_se"] for row in face_results],
        "escaped": [row["numerator"] for row in face_results],
        "conditional_escape": [row["value"] for row in face_results],
        "conditional_escape_se": [row["se"] for row in face_results],
        "escape_share": [row["escape_share"] for row in face_results],
    }, metadata)
    sidecars(cell_dir, "angular_distribution", {
        "bin_low": edges[:-1], "bin_high": edges[1:], "photons": hist.sum(axis=0),
        "probability": observed, "probability_se": np.sqrt(np.diag(covariance)),
        "uniform_probability": uniform, "old_p_2mu_probability": flux,
    }, metadata)
    result = {"metadata": metadata,
        "configuration": {key: inventory[key] for key in
            ["cell_id", "material", "opsc_code", "x_mm", "events", "workers_manifest",
             "eventModulo_manifest", "seed1", "seed2", "PDE_path", "PDE_sha256", "root_path",
             "root_sha256", "root_bytes", "wall_s"]},
        "Npe": {"left_mean": float(npe_left.mean()), "right_mean": float(npe_right.mean()),
                "per_end_mean": float((npe_left.sum()+npe_right.sum())/(2*N))},
        "V1_primary": v1, "V1_nonoptical_secondary": v1_nonopt,
        "V2_H1": v2, "V2_faces": face_results, "V2_H2": angular,
        "V5_matched": v5, "SD_orphan_fraction": orphan_fraction}
    dump(cell_dir/"cell_results.json", result)
    print(f"DONE {cell_id} V1={v1['status']} V2={v2['status']} V5={v5['status']}", flush=True)
    return result, {"V2_boot": v2_boot}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--preregistration", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if (args.out/"cells").exists():
        raise RuntimeError("analysis/cells already exists; refusing to overwrite partial results")
    script_hash = sha(Path(__file__))
    with args.inventory.open() as stream:
        inventory = list(csv.DictReader(stream))
    for row in inventory:
        for name in ["x_mm", "events", "workers_manifest", "eventModulo_manifest", "seed1", "seed2",
                     "root_bytes", "event_rows", "encounter_rows", "sipm_event_rows", "hit_rows"]:
            row[name] = int(row[name])
        row["wall_s"] = float(row["wall_s"])
    inventory.sort(key=lambda row: (row["material"], row["x_mm"]))
    if len(inventory) != 21:
        raise RuntimeError("inventory does not have 21 cells")
    rng = np.random.default_rng(SEED)
    weights = np.stack([np.bincount(rng.integers(0, N, N), minlength=N)
                        for _ in range(NBOOT)]).astype(np.float64)
    results, internal = {}, {}
    for row in inventory:
        results[row["cell_id"]], internal[row["cell_id"]] = analyze_cell(
            row, weights, args.out, args.preregistration, script_hash)

    h3_rows = []
    materials = ["EJ-200", "EJ-204", "EJ-230"]
    cell_by = {(result["configuration"]["material"], result["configuration"]["x_mm"]): cell_id
               for cell_id, result in results.items()}
    for x in [-650, -500, -200, 0, 200, 500, 650]:
        for first, second in [(materials[0], materials[1]), (materials[0], materials[2]),
                              (materials[1], materials[2])]:
            a, b = cell_by[first, x], cell_by[second, x]
            value_a = results[a]["V2_H1"]["value"]
            value_b = results[b]["V2_H1"]["value"]
            diff = value_a-value_b
            boot_diff = internal[a]["V2_boot"]-internal[b]["V2_boot"]
            se = float(boot_diff.std(ddof=1))
            status = "PASS" if abs(diff) <= 3*se else "FAIL"
            h3_rows.append({"x_mm": x, "material_a": first, "material_b": second,
                "cell_a": a, "cell_b": b, "escape_a": value_a, "escape_b": value_b,
                "difference": diff, "paired_se": se, "abs_difference_over_se": abs(diff)/se,
                "status": status})
    h3_status = "PASS" if all(row["status"] == "PASS" for row in h3_rows) else "FAIL"
    metadata = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inventory": str(args.inventory), "inventory_sha256": sha(args.inventory),
        "N_cells": 21, "N_per_cell": N, "seeds": [26092601, 8349041],
        "bootstrap_seed": SEED, "bootstrap_replicates": NBOOT,
        "preregistration": str(args.preregistration),
        "preregistration_sha256": sha(args.preregistration), "script_sha256": script_hash}
    summary_rows = []
    for cell_id, result in results.items():
        cfg = result["configuration"]
        summary_rows.append({"cell_id": cell_id, "material": cfg["material"], "x_mm": cfg["x_mm"],
            "Npe_end": result["Npe"]["per_end_mean"],
            "V1": result["V1_primary"]["value"], "V1_se": result["V1_primary"]["se"],
            "V1_status": result["V1_primary"]["status"],
            "V1_nonoptical": result["V1_nonoptical_secondary"]["value"],
            "V1_nonoptical_se": result["V1_nonoptical_secondary"]["se"],
            "V2_escape": result["V2_H1"]["value"], "V2_se": result["V2_H1"]["se"],
            "V2_lower_3se": result["V2_H1"]["lower_3se"],
            "V2_status": result["V2_H1"]["status"],
            "V5_matched": result["V5_matched"]["value"],
            "V5_expected": result["V5_matched"]["expected_incident_PDE"]["value"],
            "V5_difference": result["V5_matched"]["difference"],
            "V5_difference_se": result["V5_matched"]["difference_se"],
            "V5_status": result["V5_matched"]["status"],
            "orphan_fraction": result["SD_orphan_fraction"]["value"],
            "orphan_fraction_se": result["SD_orphan_fraction"]["se"]})
    summary_rows.sort(key=lambda row: (row["material"], row["x_mm"]))
    sidecars(args.out, "contract_cells", {key: [row[key] for row in summary_rows]
        for key in summary_rows[0]}, metadata)
    sidecars(args.out, "H3_material_invariance", {key: [row[key] for row in h3_rows]
        for key in h3_rows[0]}, metadata)
    sidecars(args.out, "escape_profiles", {
        "cell_id": [row["cell_id"] for row in summary_rows],
        "material": [row["material"] for row in summary_rows],
        "x_mm": [row["x_mm"] for row in summary_rows],
        "escape": [row["V2_escape"] for row in summary_rows],
        "escape_se": [row["V2_se"] for row in summary_rows],
    }, metadata)
    required_statuses = ([row["V1_status"] for row in summary_rows] +
                         [row["V2_status"] for row in summary_rows] +
                         [row["V5_status"] for row in summary_rows] + [h3_status])
    ready = all(status == "PASS" for status in required_statuses)
    final = {"metadata": metadata,
        "status_counts": {
            "V1_PASS": sum(row["V1_status"] == "PASS" for row in summary_rows),
            "V2_H1_PASS": sum(row["V2_status"] == "PASS" for row in summary_rows),
            "V5_PASS": sum(row["V5_status"] == "PASS" for row in summary_rows),
            "H3_pair_PASS": sum(row["status"] == "PASS" for row in h3_rows),
            "H3_pair_total": len(h3_rows)},
        "H3_status": h3_status, "H3": h3_rows,
        "ready_for_acceptance": ready,
        "failure_reasons": [] if ready else [
            name for name, status in [("V1", all(row["V1_status"] == "PASS" for row in summary_rows)),
                                      ("V2_H1", all(row["V2_status"] == "PASS" for row in summary_rows)),
                                      ("V5_matched", all(row["V5_status"] == "PASS" for row in summary_rows)),
                                      ("H3", h3_status == "PASS")] if not status],
        "cells": results,
        "paired_design": "All cells share seeds; positions within a material are correlated."}
    dump(args.out/"contract_results.json", final)
    print(json.dumps({key: final[key] for key in
                      ["status_counts", "H3_status", "ready_for_acceptance", "failure_reasons"]}, indent=2))


if __name__ == "__main__":
    main()
