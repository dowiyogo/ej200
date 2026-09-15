#!/usr/bin/env python3
"""Ejecuta el inventario y las compuertas de integridad del Paso 1 de EXEC_46."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import sys

import numpy as np
import uproot

from exec46_schema import (
    BRANCH_TYPES,
    CAMPAIGN_DIR,
    LEFT_FACE,
    MATERIAL_BY_OPSC,
    RIGHT_FACE,
    SENSOR_MAP,
    SOURCE_LABELS,
    TOP_FACE,
    TREE_NAME,
    group_velocity_mm_per_ns,
    load_material_config,
)


OUTPUT_DIR = Path(__file__).resolve().parent
INVENTORY_PATH = OUTPUT_DIR / "exec46_inventory.csv"
SOURCE_CENSUS_PATH = OUTPUT_DIR / "exec46_source_census.csv"
SCHEMA_DUMP_PATH = OUTPUT_DIR / "exec46_schema_dump.json"
AUDIT_PATH = OUTPUT_DIR / "exec46_step1_audit.json"
STEP_SIZE = "256 MB"
PROCESS_COUNT = 4
EXPECTED_EVENTS = 10_000
PATH_TOLERANCE_MM = 1.0e-6
TIME_TOLERANCE_NS = 1.0e-12
SPEED_TOLERANCE_MM = 1.0e-6
SENSOR_POSITION_TOLERANCE_MM = 1.0e-9
NPE_REFERENCE_SIGMA_MULTIPLIER = 5.0
NPE_REFERENCE_ROUNDING_TOLERANCE = 0.005
TAU_TAIL_MULTIPLIER = 5.0
TAU_COMPATIBILITY_RELATIVE = 0.03
SCINTILLATION_SOURCE = 1
ALLOWED_SOURCE_TYPES = tuple(sorted(SOURCE_LABELS))

NPE_REFERENCE = {
    ("EJ-200", 0): 1056.81,
    ("EJ-200", -650): 2589.01,
    ("EJ-200", 650): 2602.54,
    ("EJ-204", 0): 796.28,
    ("EJ-204", -650): 2546.45,
    ("EJ-204", 650): 2535.60,
    ("EJ-230", 0): 608.50,
    ("EJ-230", -650): 2249.60,
    ("EJ-230", 650): 2252.03,
}

READ_BRANCHES = [
    "event_id", "track_id", "face_type", "global_id", "local_id",
    "time_ns", "t_detection_ns", "t_creation_ns",
    "x_creation_mm", "y_creation_mm", "z_creation_mm",
    "x_mm", "y_mm", "z_mm", "wl_nm_created", "path_length_mm",
    "n_boundary_encounters", "source_type",
]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_sha_file(path):
    if path is None:
        return {}
    values = {}
    for line in Path(path).read_text().splitlines():
        digest, filename = line.split(maxsplit=1)
        values[str(Path(filename).resolve())] = digest
    return values


def discover_cells(campaign_dir):
    campaign = json.loads((campaign_dir / "campaign.json").read_text())
    require(campaign["N_generated"] == EXPECTED_EVENTS, "N_generated no es 10000")
    require(len(campaign["cells"]) == 21, "la campaña no contiene 21 celdas")
    cells = []
    for cell in campaign["cells"]:
        done_path = Path(cell["output"]) / ".DONE"
        require(done_path.is_file(), f"falta {done_path}")
        done = json.loads(done_path.read_text())
        require(done["status"] == "SIMULATION_COMPLETE", f"celda incompleta: {cell['cell_id']}")
        root_path = Path(done["root_path"]).resolve()
        require(root_path.is_file(), f"falta ROOT: {root_path}")
        cells.append({**cell, "done": done, "root_path": str(root_path)})
    return sorted(cells, key=lambda item: (item["material"], item["x_mm"]))


def validate_schema(cells):
    dumps = {}
    expected = dict(BRANCH_TYPES)
    for cell in cells:
        with uproot.open(cell["root_path"]) as root_file:
            require(TREE_NAME in root_file, f"{cell['cell_id']}: falta {TREE_NAME}")
            actual = root_file[TREE_NAME].typenames()
        require(actual == expected, f"{cell['cell_id']}: esquema distinto: {actual}")
        dumps[cell["cell_id"]] = actual
    return dumps


def update_first(first_time, first_source, event_id, detection_time, source_type):
    if event_id.size == 0:
        return
    order = np.lexsort((detection_time, event_id))
    sorted_event = event_id[order]
    first = np.r_[True, sorted_event[1:] != sorted_event[:-1]]
    selected = order[first]
    events = event_id[selected]
    times = detection_time[selected]
    sources = source_type[selected]
    better = times < first_time[events]
    first_time[events[better]] = times[better]
    first_source[events[better]] = sources[better]


def analyze_cell(payload):
    cell, known_sha = payload
    cell_id = cell["cell_id"]
    root_path = Path(cell["root_path"])
    material_config = load_material_config(cell["opsc"])
    event_hits = np.zeros(EXPECTED_EVENTS, dtype=np.int64)
    face_hits = np.zeros((3, EXPECTED_EVENTS), dtype=np.int64)
    first_time = np.full((2, EXPECTED_EVENTS), np.inf)
    first_source = np.full((2, EXPECTED_EVENTS), -1, dtype=np.int8)
    source_counts = np.zeros(len(SOURCE_LABELS), dtype=np.int64)
    sensor_face = np.asarray([SENSOR_MAP[index][0] for index in range(len(SENSOR_MAP))])
    sensor_local = np.asarray([SENSOR_MAP[index][1] for index in range(len(SENSOR_MAP))])
    keys = []

    photons = 0
    delta_sum_sq = 0.0
    delta_nonzero = 0
    delta_max_abs = 0.0
    time_order_violations = 0
    path_violations = 0
    boundary_violations = 0
    speed_violations = 0
    source_violations = 0
    sensor_map_violations = 0
    minimum_time_margin = np.inf
    minimum_path_margin = np.inf
    minimum_boundary = np.inf
    maximum_speed_excess = -np.inf
    tau_count = 0
    tau_sum = 0.0
    tau_sum_sq = 0.0
    tau_cut = TAU_TAIL_MULTIPLIER * material_config["decay_time_ns"]

    with uproot.open(root_path) as root_file:
        tree = root_file[TREE_NAME]
        for arrays in tree.iterate(READ_BRANCHES, step_size=STEP_SIZE, library="np"):
            count = len(arrays["event_id"])
            photons += count
            event_id = arrays["event_id"].astype(np.int64, copy=False)
            require(np.all((0 <= event_id) & (event_id < EXPECTED_EVENTS)),
                    f"{cell_id}: event_id fuera de rango")
            event_hits += np.bincount(event_id, minlength=EXPECTED_EVENTS)

            face_type = arrays["face_type"].astype(np.int64, copy=False)
            require(np.all((0 <= face_type) & (face_type <= TOP_FACE)),
                    f"{cell_id}: face_type fuera de rango")
            for face in (LEFT_FACE, RIGHT_FACE, TOP_FACE):
                selected = face_type == face
                face_hits[face] += np.bincount(event_id[selected], minlength=EXPECTED_EVENTS)
                if face in (LEFT_FACE, RIGHT_FACE):
                    update_first(first_time[face], first_source[face], event_id[selected],
                                 arrays["t_detection_ns"][selected], arrays["source_type"][selected])

            source_type = arrays["source_type"].astype(np.int64, copy=False)
            valid_source = np.isin(source_type, ALLOWED_SOURCE_TYPES)
            source_violations += int(np.count_nonzero(~valid_source))
            if np.any(valid_source):
                source_counts += np.bincount(source_type[valid_source], minlength=len(SOURCE_LABELS))

            delta = arrays["time_ns"] - arrays["t_detection_ns"]
            finite_delta = np.isfinite(delta)
            delta_nonzero += int(np.count_nonzero((delta != 0.0) | ~finite_delta))
            delta_sum_sq += float(np.sum(delta[finite_delta] ** 2))
            if np.any(finite_delta):
                delta_max_abs = max(delta_max_abs, float(np.max(np.abs(delta[finite_delta]))))

            propagation = arrays["t_detection_ns"] - arrays["t_creation_ns"]
            chord = np.sqrt(
                (arrays["x_mm"] - arrays["x_creation_mm"]) ** 2
                + (arrays["y_mm"] - arrays["y_creation_mm"]) ** 2
                + (arrays["z_mm"] - arrays["z_creation_mm"]) ** 2
            )
            path_margin = arrays["path_length_mm"] - chord
            finite_time = np.isfinite(propagation)
            finite_path = np.isfinite(path_margin)
            time_order_violations += int(np.count_nonzero(~finite_time | (propagation < -TIME_TOLERANCE_NS)))
            path_violations += int(np.count_nonzero(~finite_path | (path_margin < -PATH_TOLERANCE_MM)))
            if np.any(finite_time):
                minimum_time_margin = min(minimum_time_margin, float(np.min(propagation[finite_time])))
            if np.any(finite_path):
                minimum_path_margin = min(minimum_path_margin, float(np.min(path_margin[finite_path])))

            boundaries = arrays["n_boundary_encounters"]
            boundary_violations += int(np.count_nonzero(boundaries < 0))
            minimum_boundary = min(minimum_boundary, int(np.min(boundaries)))

            group_speed = group_velocity_mm_per_ns(arrays["wl_nm_created"], material_config)
            speed_excess = chord - group_speed * propagation
            finite_speed = np.isfinite(speed_excess)
            speed_violations += int(np.count_nonzero(~finite_speed | (speed_excess > SPEED_TOLERANCE_MM)))
            if np.any(finite_speed):
                maximum_speed_excess = max(maximum_speed_excess, float(np.max(speed_excess[finite_speed])))

            global_id = arrays["global_id"].astype(np.int64, copy=False)
            local_id = arrays["local_id"].astype(np.int64, copy=False)
            valid_id = (0 <= global_id) & (global_id < len(SENSOR_MAP))
            sensor_map_violations += int(np.count_nonzero(~valid_id))
            if np.any(valid_id):
                expected_face = sensor_face[global_id[valid_id]]
                expected_local = sensor_local[global_id[valid_id]]
                sensor_map_violations += int(np.count_nonzero(
                    (face_type[valid_id] != expected_face) | (local_id[valid_id] != expected_local)))

            track_id = arrays["track_id"].astype(np.uint32, copy=False)
            keys.append((event_id.astype(np.uint64) << np.uint64(32)) | track_id.astype(np.uint64))

            tail = (source_type == SCINTILLATION_SOURCE) & (arrays["t_creation_ns"] >= tau_cut)
            excess = arrays["t_creation_ns"][tail] - tau_cut
            tau_count += excess.size
            tau_sum += float(np.sum(excess))
            tau_sum_sq += float(np.sum(excess ** 2))

        require(photons == tree.num_entries, f"{cell_id}: lectura incompleta")

    all_keys = np.concatenate(keys)
    del keys
    all_keys.sort(kind="quicksort")
    duplicate_count = int(np.count_nonzero(all_keys[1:] == all_keys[:-1]))
    del all_keys

    events_observed = int(np.count_nonzero(event_hits))
    event_id_complete = bool(events_observed == EXPECTED_EVENTS and np.all(event_hits > 0))
    npe_end = face_hits[LEFT_FACE] + face_hits[RIGHT_FACE]
    npe_end_mean = float(np.mean(npe_end))
    npe_end_sem = float(np.std(npe_end, ddof=1) / np.sqrt(EXPECTED_EVENTS))
    reference = NPE_REFERENCE.get((cell["material"], cell["x_mm"]))
    reference_difference = ""
    reference_tolerance = ""
    reference_pass = True
    if reference is not None:
        reference_difference = npe_end_mean - reference
        reference_tolerance = (NPE_REFERENCE_SIGMA_MULTIPLIER * npe_end_sem
                               + NPE_REFERENCE_ROUNDING_TOLERANCE)
        reference_pass = abs(reference_difference) <= reference_tolerance

    require(tau_count > 1, f"{cell_id}: cola de centelleo vacía")
    tau_fit = tau_sum / tau_count
    tau_variance = max(0.0, (tau_sum_sq - tau_count * tau_fit ** 2) / (tau_count - 1))
    tau_error = np.sqrt(tau_variance / tau_count)
    tau_relative_difference = (tau_fit - material_config["decay_time_ns"]) / material_config["decay_time_ns"]
    tau_pass = abs(tau_relative_difference) <= TAU_COMPATIBILITY_RELATIVE

    stat = root_path.stat()
    root_key = str(root_path.resolve())
    digest = known_sha.get(root_key) or sha256(root_path)
    hash_matches_done = digest == cell["done"]["root_sha256"]
    violations = {
        "time_order": time_order_violations,
        "path": path_violations,
        "boundary": boundary_violations,
        "apparent_speed": speed_violations,
        "duplicate_event_track": duplicate_count,
        "source_type": source_violations,
        "sensor_map": sensor_map_violations,
    }
    all_gates_pass = bool(
        event_id_complete and hash_matches_done and reference_pass and delta_nonzero == 0
        and all(value == 0 for value in violations.values()) and tau_pass
        and all(np.all(np.isfinite(first_time[face])) for face in (LEFT_FACE, RIGHT_FACE))
    )

    inventory = {
        "analysis_status": "FULLY_ANALYZED",
        "cell_id": cell_id,
        "material": cell["material"],
        "opsc_code": cell["opsc"],
        "x_mm": cell["x_mm"],
        "root_path": root_key,
        "root_size_bytes": stat.st_size,
        "root_mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "root_sha256": digest,
        "hash_matches_done": hash_matches_done,
        "events": events_observed,
        "photons": photons,
        "mean_npe_end": npe_end_mean,
        "sem_npe_end": npe_end_sem,
        "mean_npe_left": float(np.mean(face_hits[LEFT_FACE])),
        "mean_npe_right": float(np.mean(face_hits[RIGHT_FACE])),
        "mean_npe_top": float(np.mean(face_hits[TOP_FACE])),
        "fraction_events_left_ge1": float(np.mean(face_hits[LEFT_FACE] > 0)),
        "fraction_events_right_ge1": float(np.mean(face_hits[RIGHT_FACE] > 0)),
        "fraction_events_both_ends_ge1": float(np.mean(
            (face_hits[LEFT_FACE] > 0) & (face_hits[RIGHT_FACE] > 0))),
        "npe_reference": reference if reference is not None else "",
        "npe_reference_difference": reference_difference,
        "npe_reference_tolerance": reference_tolerance,
        "npe_reference_pass": reference_pass,
        "time_delta_max_abs_ns": delta_max_abs,
        "time_delta_rms_ns": float(np.sqrt(delta_sum_sq / photons)),
        "time_delta_nonzero_fraction": delta_nonzero / photons,
        "time_order_violations": time_order_violations,
        "minimum_propagation_time_ns": minimum_time_margin,
        "path_violations": path_violations,
        "minimum_path_minus_chord_mm": minimum_path_margin,
        "boundary_violations": boundary_violations,
        "minimum_boundary_encounters": minimum_boundary,
        "apparent_speed_violations": speed_violations,
        "maximum_chord_minus_vgroup_time_mm": maximum_speed_excess,
        "duplicate_event_track": duplicate_count,
        "source_type_violations": source_violations,
        "sensor_map_violations": sensor_map_violations,
        "tau_config_ns": material_config["decay_time_ns"],
        "tau_tail_cut_ns": tau_cut,
        "tau_tail_count": tau_count,
        "tau_fit_ns": tau_fit,
        "tau_fit_error_ns": float(tau_error),
        "tau_relative_difference": tau_relative_difference,
        "tau_pass": tau_pass,
        "all_gates_pass": all_gates_pass,
    }

    census = []
    for source_code, source_label in SOURCE_LABELS.items():
        census.append({
            "cell_id": cell_id, "material": cell["material"], "x_mm": cell["x_mm"],
            "scope": "all_hits", "face_type": "all", "source_type": source_code,
            "source_label": source_label, "count": int(source_counts[source_code]),
            "denominator": photons, "fraction": float(source_counts[source_code] / photons),
        })
        for face in (LEFT_FACE, RIGHT_FACE):
            denominator = int(np.count_nonzero(np.isfinite(first_time[face])))
            count = int(np.count_nonzero(first_source[face] == source_code))
            census.append({
                "cell_id": cell_id, "material": cell["material"], "x_mm": cell["x_mm"],
                "scope": "first_photon_end", "face_type": face, "source_type": source_code,
                "source_label": source_label, "count": count, "denominator": denominator,
                "fraction": float(count / denominator),
            })
    return inventory, census


def write_csv(path, rows):
    require(rows, f"sin filas para {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--sha-file", type=Path,
                        help="salida sha256sum ya calculada; si falta, calcula cada hash")
    parser.add_argument("--processes", type=int, default=PROCESS_COUNT)
    parser.add_argument("--cell", action="append",
                        help="limita una sonda de hard-abort; repetible")
    args = parser.parse_args()
    start = datetime.now(timezone.utc)
    cells = discover_cells(args.campaign.resolve())
    schema_dump = validate_schema(cells)
    campaign_cells = cells
    campaign_cell_count = len(campaign_cells)
    if args.cell:
        selected = set(args.cell)
        cells = [cell for cell in cells if cell["cell_id"] in selected]
        require(len(cells) == len(selected), "alguna celda solicitada no existe")
    known_sha = load_sha_file(args.sha_file)
    require(args.processes >= 1, "processes debe ser positivo")

    with mp.get_context("spawn").Pool(args.processes) as pool:
        results = pool.map(analyze_cell, [(cell, known_sha) for cell in cells])
    analyzed_inventory = [result[0] for result in results]
    inventory = list(analyzed_inventory)
    if args.cell:
        analyzed_ids = {row["cell_id"] for row in analyzed_inventory}
        for cell in campaign_cells:
            if cell["cell_id"] in analyzed_ids:
                continue
            root_path = Path(cell["root_path"])
            placeholder = {key: "" for key in analyzed_inventory[0]}
            placeholder.update({
                "analysis_status": "NOT_ANALYZED_AFTER_HARD_ABORT",
                "cell_id": cell["cell_id"],
                "material": cell["material"],
                "opsc_code": cell["opsc"],
                "x_mm": cell["x_mm"],
                "root_path": str(root_path),
                "root_size_bytes": root_path.stat().st_size,
                "root_mtime_utc": datetime.fromtimestamp(
                    root_path.stat().st_mtime, timezone.utc).isoformat(),
                "root_sha256": known_sha.get(str(root_path), ""),
                "hash_matches_done": known_sha.get(str(root_path)) == cell["done"]["root_sha256"],
                "events": cell["done"]["events_run"],
                "photons": cell["done"]["root_entries"],
                "mean_npe_end": (
                    cell["done"]["left_total"] + cell["done"]["right_total"]
                ) / EXPECTED_EVENTS,
                "mean_npe_left": cell["done"]["left_total"] / EXPECTED_EVENTS,
                "mean_npe_right": cell["done"]["right_total"] / EXPECTED_EVENTS,
                "mean_npe_top": cell["done"]["top_total"] / EXPECTED_EVENTS,
                "all_gates_pass": "",
            })
            inventory.append(placeholder)
    inventory.sort(key=lambda row: (row["material"], int(row["x_mm"])))
    census = [row for result in results for row in result[1]]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / INVENTORY_PATH.name, inventory)
    write_csv(args.output_dir / SOURCE_CENSUS_PATH.name, census)
    schema_payload = {
        "tree": TREE_NAME,
        "expected": BRANCH_TYPES,
        "verified_cells": schema_dump,
        "sensor_map": {str(key): value for key, value in SENSOR_MAP.items()},
    }
    (args.output_dir / SCHEMA_DUMP_PATH.name).write_text(json.dumps(schema_payload, indent=2) + "\n")
    audit = {
        "start_utc": start.isoformat(),
        "end_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "constants": {
            "step_size": STEP_SIZE,
            "processes": args.processes,
            "expected_events": EXPECTED_EVENTS,
            "path_tolerance_mm": PATH_TOLERANCE_MM,
            "time_tolerance_ns": TIME_TOLERANCE_NS,
            "speed_tolerance_mm": SPEED_TOLERANCE_MM,
            "npe_reference_sigma_multiplier": NPE_REFERENCE_SIGMA_MULTIPLIER,
            "npe_reference_rounding_tolerance": NPE_REFERENCE_ROUNDING_TOLERANCE,
            "tau_tail_multiplier": TAU_TAIL_MULTIPLIER,
            "tau_compatibility_relative": TAU_COMPATIBILITY_RELATIVE,
        },
        "cells": len(inventory),
        "analyzed_cells": len(analyzed_inventory),
        "campaign_cells": campaign_cell_count,
        "campaign_complete": len(analyzed_inventory) == campaign_cell_count,
        "all_gates_pass": (len(analyzed_inventory) == campaign_cell_count
                           and all(row["all_gates_pass"] for row in analyzed_inventory)),
    }
    (args.output_dir / AUDIT_PATH.name).write_text(json.dumps(audit, indent=2) + "\n")
    require(audit["all_gates_pass"], "falló al menos una compuerta; ver exec46_inventory.csv")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
