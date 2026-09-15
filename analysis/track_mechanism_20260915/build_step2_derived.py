#!/usr/bin/env python3
"""Construye el árbol por evento de EXEC_46 y controles pareados de fotones."""

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot

from analyze_step1 import discover_cells
from exec46_schema import (
    CAMPAIGN_DIR,
    LEFT_FACE,
    MATERIAL_BY_OPSC,
    RIGHT_FACE,
    TOP_FACE,
    TREE_NAME,
    load_material_config,
)


OUTPUT_DIR = Path(__file__).resolve().parent / "step2"
DERIVED_PATH = OUTPUT_DIR / "exec46_derived_events.root"
META_PATH = OUTPUT_DIR / "exec46_derived_events.meta.json"
MATERIAL_PATH = OUTPUT_DIR / "material_optical_properties.csv"
EXPECTED_EVENTS = 10_000
EXPECTED_CELLS = 21
PROCESS_COUNT = 4
STEP_SIZE = "256 MB"
SCINTILLATION_SOURCE = 1
CHERENKOV_SOURCE = 2
RANDOM_SELECTION_SEED = np.uint64(0x46A2C0DE5EED1234)
UINT64_MASK = np.uint64(0xFFFFFFFFFFFFFFFF)
READ_BRANCHES = [
    "event_id", "face_type", "track_id", "time_ns", "t_detection_ns",
    "t_creation_ns", "source_type", "exit_angle_deg",
    "n_boundary_encounters", "path_length_mm", "wl_nm", "wl_nm_created",
    "pde", "x_creation_mm", "y_creation_mm", "z_creation_mm",
]
PHOTON_FIELDS = [
    "track_id", "time_ns", "t_detection_ns", "t_creation_ns", "source_type",
    "exit_angle_deg", "n_boundary_encounters", "path_length_mm", "wl_nm",
    "wl_nm_created", "pde", "x_creation_mm", "y_creation_mm", "z_creation_mm",
]
MATERIAL_CODES = {"EJ-200": 0, "EJ-204": 1, "EJ-230": 2}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def splitmix64(values):
    """Hash determinista uniforme usado para el control aleatorio pareado."""
    with np.errstate(over="ignore"):
        values = (values + np.uint64(0x9E3779B97F4A7C15)) & UINT64_MASK
        values = ((values ^ (values >> np.uint64(30)))
                  * np.uint64(0xBF58476D1CE4E5B9)) & UINT64_MASK
        values = ((values ^ (values >> np.uint64(27)))
                  * np.uint64(0x94D049BB133111EB)) & UINT64_MASK
    return values ^ (values >> np.uint64(31))


def empty_selection():
    selection = {
        "priority": np.full(2 * EXPECTED_EVENTS, np.iinfo(np.uint64).max,
                            dtype=np.uint64),
        "track_id": np.full(2 * EXPECTED_EVENTS, np.iinfo(np.int32).max,
                            dtype=np.int32),
    }
    for field in PHOTON_FIELDS[1:]:
        dtype = np.int32 if field in ("source_type", "n_boundary_encounters") else np.float64
        fill = -1 if np.issubdtype(dtype, np.integer) else np.nan
        selection[field] = np.full(2 * EXPECTED_EVENTS, fill, dtype=dtype)
    return selection


def update_selection(selection, arrays, indices, priorities, first_by_time):
    if indices.size == 0:
        return
    event_id = arrays["event_id"][indices].astype(np.int64, copy=False)
    face = arrays["face_type"][indices].astype(np.int64, copy=False)
    group = 2 * event_id + face
    track = arrays["track_id"][indices].astype(np.int32, copy=False)
    if first_by_time:
        time = arrays["time_ns"][indices]
        order = np.lexsort((track, time, group))
    else:
        order = np.lexsort((track, priorities, group))
    sorted_group = group[order]
    keep = np.r_[True, sorted_group[1:] != sorted_group[:-1]]
    chosen_local = order[keep]
    chosen = indices[chosen_local]
    chosen_group = group[chosen_local]
    chosen_track = track[chosen_local]
    chosen_priority = priorities[chosen_local]
    if first_by_time:
        old_time = selection["time_ns"][chosen_group]
        new_time = arrays["time_ns"][chosen]
        better = ((~np.isfinite(old_time)) | (new_time < old_time)
                  | ((new_time == old_time)
                     & (chosen_track < selection["track_id"][chosen_group])))
    else:
        better = ((chosen_priority < selection["priority"][chosen_group])
                  | ((chosen_priority == selection["priority"][chosen_group])
                     & (chosen_track < selection["track_id"][chosen_group])))
    destination = chosen_group[better]
    source = chosen[better]
    selection["priority"][destination] = chosen_priority[better]
    selection["track_id"][destination] = chosen_track[better]
    for field in PHOTON_FIELDS[1:]:
        selection[field][destination] = arrays[field][source]


def reshape_face(values, face):
    return values.reshape(EXPECTED_EVENTS, 2)[:, face]


def analyze_cell(payload):
    cell, cell_index = payload
    root_path = Path(cell["root_path"])
    face_counts = np.zeros((3, EXPECTED_EVENTS), dtype=np.int32)
    scint_counts = np.zeros((2, EXPECTED_EVENTS), dtype=np.int32)
    cherenkov_counts = np.zeros((2, EXPECTED_EVENTS), dtype=np.int32)
    event_counts = np.zeros(EXPECTED_EVENTS, dtype=np.int64)
    first = empty_selection()
    random = empty_selection()
    clock_nonzero = 0

    with uproot.open(root_path) as root_file:
        tree = root_file[TREE_NAME]
        for arrays in tree.iterate(READ_BRANCHES, step_size=STEP_SIZE, library="np"):
            event_id = arrays["event_id"].astype(np.int64, copy=False)
            face = arrays["face_type"].astype(np.int64, copy=False)
            source = arrays["source_type"].astype(np.int64, copy=False)
            require(np.all((0 <= event_id) & (event_id < EXPECTED_EVENTS)),
                    f"{cell['cell_id']}: event_id fuera de rango")
            event_counts += np.bincount(event_id, minlength=EXPECTED_EVENTS)
            for face_code in (LEFT_FACE, RIGHT_FACE, TOP_FACE):
                mask = face == face_code
                face_counts[face_code] += np.bincount(
                    event_id[mask], minlength=EXPECTED_EVENTS).astype(np.int32)
            for face_code in (LEFT_FACE, RIGHT_FACE):
                for source_code, target in (
                    (SCINTILLATION_SOURCE, scint_counts),
                    (CHERENKOV_SOURCE, cherenkov_counts),
                ):
                    mask = (face == face_code) & (source == source_code)
                    target[face_code] += np.bincount(
                        event_id[mask], minlength=EXPECTED_EVENTS).astype(np.int32)

            end_indices = np.flatnonzero((face == LEFT_FACE) | (face == RIGHT_FACE))
            track_u64 = arrays["track_id"][end_indices].astype(np.uint64, copy=False)
            event_u64 = event_id[end_indices].astype(np.uint64, copy=False)
            face_u64 = face[end_indices].astype(np.uint64, copy=False)
            key = (track_u64 ^ (event_u64 << np.uint64(32))
                   ^ (face_u64 << np.uint64(60))
                   ^ (np.uint64(cell_index) << np.uint64(48))
                   ^ RANDOM_SELECTION_SEED)
            priority = splitmix64(key)
            update_selection(first, arrays, end_indices, arrays["time_ns"][end_indices], True)
            update_selection(random, arrays, end_indices, priority, False)
            clock_nonzero += int(np.count_nonzero(
                arrays["time_ns"] != arrays["t_detection_ns"]))

    require(np.all(event_counts > 0), f"{cell['cell_id']}: eventos ausentes")
    require(np.all(face_counts[LEFT_FACE] > 0) and np.all(face_counts[RIGHT_FACE] > 0),
            f"{cell['cell_id']}: evento sin END")
    require(clock_nonzero == 0, f"{cell['cell_id']}: relojes distintos")
    require(np.all(np.isfinite(first["time_ns"])),
            f"{cell['cell_id']}: primer fotón END ausente")
    require(np.all(np.isfinite(random["time_ns"])),
            f"{cell['cell_id']}: control aleatorio END ausente")

    output = {
        "cell_code": np.full(EXPECTED_EVENTS, cell_index, dtype=np.int32),
        "material_code": np.full(EXPECTED_EVENTS, MATERIAL_CODES[cell["material"]],
                                 dtype=np.int32),
        "x_mm": np.full(EXPECTED_EVENTS, cell["x_mm"], dtype=np.int32),
        "event_id": np.arange(EXPECTED_EVENTS, dtype=np.int32),
        "npe_left": face_counts[LEFT_FACE],
        "npe_right": face_counts[RIGHT_FACE],
        "npe_end": face_counts[LEFT_FACE] + face_counts[RIGHT_FACE],
        "npe_top": face_counts[TOP_FACE],
        "npe_scint_left": scint_counts[LEFT_FACE],
        "npe_scint_right": scint_counts[RIGHT_FACE],
        "npe_cherenkov_left": cherenkov_counts[LEFT_FACE],
        "npe_cherenkov_right": cherenkov_counts[RIGHT_FACE],
    }
    for label, selection in (("first", first), ("random", random)):
        for face_name, face_code in (("left", LEFT_FACE), ("right", RIGHT_FACE)):
            for field in PHOTON_FIELDS:
                output[f"{label}_{face_name}_{field}"] = reshape_face(
                    selection[field], face_code)
    output["t_left_ns"] = output["first_left_time_ns"]
    output["t_right_ns"] = output["first_right_time_ns"]
    output["t0_ns"] = 0.5 * (output["t_left_ns"] + output["t_right_ns"])
    output["t_left_detection_ns"] = output["first_left_t_detection_ns"]
    output["t_right_detection_ns"] = output["first_right_t_detection_ns"]
    output["t0_detection_ns"] = 0.5 * (
        output["t_left_detection_ns"] + output["t_right_detection_ns"])
    return output


def write_material_properties():
    rows = []
    for opsc_code, material in MATERIAL_BY_OPSC.items():
        config = load_material_config(opsc_code)
        rindex_constant = bool(np.all(config["rindex"] == config["rindex"][0]))
        absorption = config["absorption_length_mm"]
        absorption_constant = bool(np.all(absorption == absorption[0]))
        rows.append({
            "material": material,
            "opsc_code": opsc_code,
            "rindex_wavelength_min_nm": float(np.min(config["rindex_wavelength_nm"])),
            "rindex_wavelength_max_nm": float(np.max(config["rindex_wavelength_nm"])),
            "rindex_min": float(np.min(config["rindex"])),
            "rindex_max": float(np.max(config["rindex"])),
            "rindex_constant": rindex_constant,
            "abs_wavelength_min_nm": float(np.min(config["absorption_wavelength_nm"])),
            "abs_wavelength_max_nm": float(np.max(config["absorption_wavelength_nm"])),
            "abs_length_min_mm": float(np.min(absorption)),
            "abs_length_max_mm": float(np.max(absorption)),
            "abs_length_constant": absorption_constant,
            "rindex_path": str(config["paths"]["rindex"]),
            "absorption_path": str(config["paths"]["absorption"]),
        })
    with MATERIAL_PATH.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN_DIR)
    parser.add_argument("--processes", type=int, default=PROCESS_COUNT)
    args = parser.parse_args()
    require(args.processes >= 1, "processes debe ser positivo")
    start = datetime.now(timezone.utc)
    cells = discover_cells(args.campaign.resolve())
    require(len(cells) == EXPECTED_CELLS, "se esperaban 21 celdas")
    cells.sort(key=lambda cell: (cell["material"], int(cell["x_mm"])))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with mp.get_context("spawn").Pool(args.processes) as pool:
        results = pool.map(analyze_cell, [(cell, index) for index, cell in enumerate(cells)])
    branch_names = list(results[0])
    require(all(list(result) == branch_names for result in results),
            "ramas derivadas inconsistentes")
    combined = {name: np.concatenate([result[name] for result in results])
                for name in branch_names}
    require(len(combined["event_id"]) == EXPECTED_CELLS * EXPECTED_EVENTS,
            "conteo derivado incorrecto")
    with uproot.recreate(DERIVED_PATH) as output:
        output["derived_events"] = combined
    material_rows = write_material_properties()
    payload = {
        "start_utc": start.isoformat(),
        "end_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "input_mode": "READ",
        "new_simulation": False,
        "tree": "derived_events",
        "entries": len(combined["event_id"]),
        "cells": [{"cell_code": index, "cell_id": cell["cell_id"],
                   "material": cell["material"], "x_mm": cell["x_mm"],
                   "root_path": cell["root_path"],
                   "root_sha256": cell["done"]["root_sha256"]}
                  for index, cell in enumerate(cells)],
        "random_control": {
            "algorithm": "minimum splitmix64 priority per event and END face",
            "seed_hex": hex(int(RANDOM_SELECTION_SEED)),
            "unit": "one photon uniformly selected from each event and END face",
        },
        "material_properties": material_rows,
        "derived_root_sha256": sha256(DERIVED_PATH),
    }
    META_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in
                      ("start_utc", "end_utc", "entries", "derived_root_sha256")},
                     indent=2))


if __name__ == "__main__":
    main()
