#!/usr/bin/env python3
"""Construye pares primero/aleatorio y mínimos por fuente para el Paso 4."""

import argparse
import hashlib
import json
import multiprocessing as mp
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot

from analyze_step1 import discover_cells
from exec46_schema import (
    BAR_HALF_Z_MM,
    CAMPAIGN_DIR,
    LEFT_FACE,
    RIGHT_FACE,
    SPEED_OF_LIGHT_MM_PER_NS,
    TREE_NAME,
)


OUTPUT_DIR = Path(__file__).resolve().parent / "step4"
OUTPUT_ROOT = OUTPUT_DIR / "step4_event_pairs.root"
OUTPUT_META = OUTPUT_DIR / "step4_event_pairs.meta.json"
EXPECTED_EVENTS = 10_000
EXPECTED_CELLS = 21
PROCESS_COUNT = 4
STEP_SIZE = "192 MB"
SCINTILLATION_SOURCE = 1
CHERENKOV_SOURCE = 2
PRIMARY_LIKE_TOLERANCE_MM = 0.001
LOW_BOUNDARY_MAX = 2
RANDOM_SELECTION_SEED = np.uint64(0x46A2C0DE5EED1234)
UINT64_MASK = np.uint64(0xFFFFFFFFFFFFFFFF)
READ_BRANCHES = [
    "event_id", "track_id", "face_type", "source_type", "time_ns",
    "t_detection_ns", "t_creation_ns", "x_mm", "y_mm", "z_mm",
    "x_creation_mm", "y_creation_mm", "z_creation_mm", "path_length_mm",
    "exit_angle_deg", "n_boundary_encounters", "wl_nm_created", "wl_nm", "pde",
]
SELECTION_FIELDS = [
    "track_id", "source_type", "time_ns", "t_detection_ns", "t_creation_ns",
    "x_mm", "y_mm", "z_mm", "x_creation_mm", "y_creation_mm", "z_creation_mm",
    "path_length_mm", "exit_angle_deg", "n_boundary_encounters",
    "wl_nm_created", "wl_nm", "pde", "d_direct_mm", "rho_detour",
    "t_creation_corrected_ns",
]
SELECTION_NAMES = (
    "first", "random", "first_scint", "first_cherenkov",
    "min_creation_scint", "min_corrected_creation_scint",
    "first_primary_cherenkov",
)
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
    """Hash determinista idéntico al control aleatorio del Paso 2."""
    with np.errstate(over="ignore"):
        values = (values + np.uint64(0x9E3779B97F4A7C15)) & UINT64_MASK
        values = ((values ^ (values >> np.uint64(30)))
                  * np.uint64(0xBF58476D1CE4E5B9)) & UINT64_MASK
        values = ((values ^ (values >> np.uint64(27)))
                  * np.uint64(0x94D049BB133111EB)) & UINT64_MASK
    return values ^ (values >> np.uint64(31))


def empty_selection():
    output = {"priority": np.full(2 * EXPECTED_EVENTS,
                                  np.iinfo(np.uint64).max, dtype=np.uint64)}
    for field in SELECTION_FIELDS:
        integer = field in ("track_id", "source_type", "n_boundary_encounters")
        output[field] = np.full(2 * EXPECTED_EVENTS, -1 if integer else np.nan,
                                dtype=np.int32 if integer else np.float64)
    return output


def update_selection(selection, arrays, selected, priorities, quantities,
                     use_priority=False):
    if selected.size == 0:
        return
    event = arrays["event_id"][selected].astype(np.int64, copy=False)
    face = arrays["face_type"][selected].astype(np.int64, copy=False)
    group = 2 * event + face
    track = arrays["track_id"][selected].astype(np.int32, copy=False)
    order = np.lexsort((track, priorities, group))
    sorted_group = group[order]
    keep = np.r_[True, sorted_group[1:] != sorted_group[:-1]]
    local = order[keep]
    chosen = selected[local]
    destination = group[local]
    candidate_priority = priorities[local]
    if use_priority:
        better = ((candidate_priority < selection["priority"][destination])
                  | ((candidate_priority == selection["priority"][destination])
                     & (track[local] < selection["track_id"][destination])))
    else:
        old = selection["priority"][destination].view(np.float64)
        better = ((~np.isfinite(old)) | (candidate_priority < old)
                  | ((candidate_priority == old)
                     & (track[local] < selection["track_id"][destination])))
    destination = destination[better]
    chosen = chosen[better]
    candidate_priority = candidate_priority[better]
    if use_priority:
        selection["priority"][destination] = candidate_priority.astype(np.uint64)
    else:
        selection["priority"][destination] = candidate_priority.astype(np.float64).view(np.uint64)
    for field in SELECTION_FIELDS:
        source = quantities[field] if field in quantities else arrays[field]
        selection[field][destination] = source[chosen]


def analyze_cell(payload):
    cell, cell_index, scratch_dir = payload
    selections = {name: empty_selection() for name in SELECTION_NAMES}
    counts_scint = np.zeros((2, EXPECTED_EVENTS), dtype=np.int32)
    counts_cherenkov = np.zeros((2, EXPECTED_EVENTS), dtype=np.int32)
    counts_scint_low = np.zeros((2, EXPECTED_EVENTS), dtype=np.int32)
    counts_primary_cherenkov = np.zeros((2, EXPECTED_EVENTS), dtype=np.int32)
    event_seen = np.zeros(EXPECTED_EVENTS, dtype=np.int64)
    gun_x = float(cell["x_mm"])

    with uproot.open(cell["root_path"]) as root_file:
        for arrays in root_file[TREE_NAME].iterate(
                READ_BRANCHES, step_size=STEP_SIZE, library="np"):
            event = arrays["event_id"].astype(np.int64, copy=False)
            face = arrays["face_type"].astype(np.int64, copy=False)
            source = arrays["source_type"].astype(np.int64, copy=False)
            require(np.all((event >= 0) & (event < EXPECTED_EVENTS)),
                    f"{cell['cell_id']}: event_id fuera de rango")
            event_seen += np.bincount(event, minlength=EXPECTED_EVENTS)
            end = (face == LEFT_FACE) | (face == RIGHT_FACE)
            selected = np.flatnonzero(end)
            dx = arrays["x_mm"] - arrays["x_creation_mm"]
            dy = arrays["y_mm"] - arrays["y_creation_mm"]
            dz = arrays["z_mm"] - arrays["z_creation_mm"]
            direct = np.sqrt(dx * dx + dy * dy + dz * dz)
            rho = arrays["path_length_mm"] / direct
            deposit_time = ((BAR_HALF_Z_MM - arrays["z_creation_mm"])
                            / SPEED_OF_LIGHT_MM_PER_NS)
            corrected_creation = arrays["t_creation_ns"] - deposit_time
            quantities = {
                "d_direct_mm": direct,
                "rho_detour": rho,
                "t_creation_corrected_ns": corrected_creation,
            }

            detection_priority = arrays["t_detection_ns"][selected]
            update_selection(selections["first"], arrays, selected,
                             detection_priority, quantities)
            track_u64 = arrays["track_id"][selected].astype(np.uint64, copy=False)
            event_u64 = event[selected].astype(np.uint64, copy=False)
            face_u64 = face[selected].astype(np.uint64, copy=False)
            key = (track_u64 ^ (event_u64 << np.uint64(32))
                   ^ (face_u64 << np.uint64(60))
                   ^ (np.uint64(cell_index) << np.uint64(48))
                   ^ RANDOM_SELECTION_SEED)
            update_selection(selections["random"], arrays, selected,
                             splitmix64(key), quantities, use_priority=True)

            for source_code, name, target in (
                    (SCINTILLATION_SOURCE, "first_scint", counts_scint),
                    (CHERENKOV_SOURCE, "first_cherenkov", counts_cherenkov)):
                source_selected = np.flatnonzero(end & (source == source_code))
                update_selection(selections[name], arrays, source_selected,
                                 arrays["t_detection_ns"][source_selected], quantities)
                for face_code in (LEFT_FACE, RIGHT_FACE):
                    mask = (face == face_code) & (source == source_code)
                    target[face_code] += np.bincount(
                        event[mask], minlength=EXPECTED_EVENTS).astype(np.int32)

            scint_selected = np.flatnonzero(end & (source == SCINTILLATION_SOURCE))
            update_selection(selections["min_creation_scint"], arrays, scint_selected,
                             arrays["t_creation_ns"][scint_selected], quantities)
            update_selection(selections["min_corrected_creation_scint"], arrays,
                             scint_selected, corrected_creation[scint_selected], quantities)
            primary_like = (end & (source == CHERENKOV_SOURCE)
                            & (np.abs(arrays["x_creation_mm"] - gun_x)
                               < PRIMARY_LIKE_TOLERANCE_MM)
                            & (np.abs(arrays["y_creation_mm"])
                               < PRIMARY_LIKE_TOLERANCE_MM))
            primary_selected = np.flatnonzero(primary_like)
            update_selection(selections["first_primary_cherenkov"], arrays,
                             primary_selected,
                             arrays["t_detection_ns"][primary_selected], quantities)
            for face_code in (LEFT_FACE, RIGHT_FACE):
                low = ((face == face_code) & (source == SCINTILLATION_SOURCE)
                       & (arrays["n_boundary_encounters"] <= LOW_BOUNDARY_MAX))
                counts_scint_low[face_code] += np.bincount(
                    event[low], minlength=EXPECTED_EVENTS).astype(np.int32)
                primary_face = primary_like & (face == face_code)
                counts_primary_cherenkov[face_code] += np.bincount(
                    event[primary_face], minlength=EXPECTED_EVENTS).astype(np.int32)

    require(np.all(event_seen > 0), f"{cell['cell_id']}: eventos ausentes")
    require(np.all(selections["first"]["track_id"] >= 0),
            f"{cell['cell_id']}: first ausente")
    require(np.all(selections["random"]["track_id"] >= 0),
            f"{cell['cell_id']}: random ausente")
    require(np.all(selections["first_scint"]["track_id"] >= 0),
            f"{cell['cell_id']}: first scint ausente")

    output = {
        "cell_code": np.repeat(np.int32(cell_index), 2 * EXPECTED_EVENTS),
        "material_code": np.repeat(np.int32(MATERIAL_CODES[cell["material"]]),
                                   2 * EXPECTED_EVENTS),
        "gun_x_mm": np.repeat(np.int32(cell["x_mm"]), 2 * EXPECTED_EVENTS),
        "event_id": np.repeat(np.arange(EXPECTED_EVENTS, dtype=np.int32), 2),
        "face_type": np.tile(np.asarray([LEFT_FACE, RIGHT_FACE], np.int32),
                             EXPECTED_EVENTS),
        "npe_scint": counts_scint.T.reshape(-1),
        "npe_cherenkov": counts_cherenkov.T.reshape(-1),
        "npe_scint_low_boundary": counts_scint_low.T.reshape(-1),
        "npe_primary_cherenkov": counts_primary_cherenkov.T.reshape(-1),
    }
    for name, selection in selections.items():
        for field in SELECTION_FIELDS:
            output[f"{name}_{field}"] = selection[field]
    scratch_path = Path(scratch_dir) / f"{cell['cell_id']}.npz"
    np.savez_compressed(scratch_path, **output)
    return str(scratch_path)


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
    scratch = OUTPUT_DIR / "scratch"
    scratch.mkdir(exist_ok=True)
    payloads = [(cell, index, scratch) for index, cell in enumerate(cells)]
    with mp.get_context("spawn").Pool(args.processes) as pool:
        paths = pool.map(analyze_cell, payloads)
    chunks = [np.load(path) for path in paths]
    columns = {key: np.concatenate([chunk[key] for chunk in chunks])
               for key in chunks[0].files}
    require(len(columns["event_id"]) == EXPECTED_CELLS * EXPECTED_EVENTS * 2,
            "número inesperado de filas")
    with uproot.recreate(OUTPUT_ROOT) as root_file:
        root_file.mktree("event_pairs", {key: value.dtype for key, value in columns.items()})
        root_file["event_pairs"].extend(columns)
    for chunk in chunks:
        chunk.close()
    for path in paths:
        Path(path).unlink()
    scratch.rmdir()
    finish = datetime.now(timezone.utc)
    metadata = {
        "created_utc": finish.isoformat(), "started_utc": start.isoformat(),
        "wall_seconds": (finish - start).total_seconds(),
        "command": ("env PYTHONPATH=analysis/track_mechanism_20260915 python3 "
                    "analysis/track_mechanism_20260915/build_step4_pairs.py --processes "
                    f"{args.processes}"),
        "campaign": str(args.campaign.resolve()), "cells": EXPECTED_CELLS,
        "events_per_cell": EXPECTED_EVENTS, "rows": len(columns["event_id"]),
        "processes": args.processes, "step_size": STEP_SIZE,
        "random_seed_hex": hex(int(RANDOM_SELECTION_SEED)),
        "primary_like_tolerance_mm": PRIMARY_LIKE_TOLERANCE_MM,
        "low_boundary_max": LOW_BOUNDARY_MAX,
        "muon_transit_correction": {
            "gun_direction": [0.0, 0.0, -1.0],
            "bar_entry_z_mm": BAR_HALF_Z_MM,
            "formula": "t_creation-(bar_entry_z-z_creation)/c",
            "speed_mm_per_ns": SPEED_OF_LIGHT_MM_PER_NS,
        },
        "output": str(OUTPUT_ROOT.resolve()), "output_sha256": sha256(OUTPUT_ROOT),
    }
    OUTPUT_META.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
