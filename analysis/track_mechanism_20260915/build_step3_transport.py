#!/usr/bin/env python3
"""Resume el transporte fotón-a-fotón de EXEC_46 sin modificar los ROOT fuente."""

import argparse
import csv
import json
import multiprocessing as mp
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot

from analyze_step1 import discover_cells
from exec46_schema import CAMPAIGN_DIR, LEFT_FACE, RIGHT_FACE, TREE_NAME


OUTPUT_DIR = Path(__file__).resolve().parent / "step3"
EXPECTED_EVENTS = 10_000
EXPECTED_CELLS = 21
PROCESS_COUNT = 4
STEP_SIZE = "192 MB"
SOURCE_CODES = (1, 2)
SOURCE_NAMES = {1: "scintillation", 2: "Cherenkov"}
FACE_NAMES = {LEFT_FACE: "left", RIGHT_FACE: "right"}
MATERIAL_CODES = {"EJ-200": 0, "EJ-204": 1, "EJ-230": 2}
D_EDGES_MM = np.arange(0.0, 1450.0 + 10.0, 10.0)
TPROP_EDGES_NS = np.arange(0.0, 60.0 + 0.2, 0.2)
BOUNDARY_LABELS = ("0", "1-2", "3-5", "6-10", ">10")
ANGLE_EDGES_DEG = np.asarray([0.0, 15.0, 30.0, 39.2652479139,
                              45.0, 60.0, 90.0, 120.0, 180.0000001])
N_ANGLE_BINS = len(ANGLE_EDGES_DEG) - 1
READ_BRANCHES = [
    "event_id", "track_id", "face_type", "source_type",
    "t_detection_ns", "t_creation_ns", "x_mm", "y_mm", "z_mm",
    "x_creation_mm", "y_creation_mm", "z_creation_mm", "path_length_mm",
    "exit_angle_deg", "n_boundary_encounters",
]
FIRST_FIELDS = (
    "t_detection_ns", "t_creation_ns", "tprop_ns", "d_direct_mm",
    "path_length_mm", "rho_detour", "v_apparent_mm_per_ns",
    "exit_angle_deg", "n_boundary_encounters", "x_creation_mm",
    "y_creation_mm", "z_creation_mm", "x_mm", "y_mm", "z_mm",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def boundary_category(values):
    output = np.full(values.shape, 4, dtype=np.int8)
    output[values == 0] = 0
    output[(values >= 1) & (values <= 2)] = 1
    output[(values >= 3) & (values <= 5)] = 2
    output[(values >= 6) & (values <= 10)] = 3
    return output


def cluster_mean_se(event_count, event_sum):
    count = int(np.sum(event_count))
    if count == 0:
        return np.nan, np.nan
    mean = float(np.sum(event_sum) / count)
    influence = event_sum - mean * event_count
    se = float(np.sqrt(EXPECTED_EVENTS / (EXPECTED_EVENTS - 1)
                       * np.sum(influence * influence)) / count)
    return mean, se


def update_first(first, arrays, selected, group, event_id, quantities):
    if selected.size == 0:
        return
    combined = group * EXPECTED_EVENTS + event_id
    detection = arrays["t_detection_ns"][selected]
    track = arrays["track_id"][selected]
    order = np.lexsort((track, detection, combined))
    sorted_combined = combined[order]
    keep = np.r_[True, sorted_combined[1:] != sorted_combined[:-1]]
    chosen_local = order[keep]
    chosen = selected[chosen_local]
    destination = combined[chosen_local]
    chosen_detection = detection[chosen_local]
    chosen_track = track[chosen_local]
    old_detection = first["t_detection_ns"].reshape(-1)[destination]
    old_track = first["track_id"].reshape(-1)[destination]
    better = ((~np.isfinite(old_detection))
              | (chosen_detection < old_detection)
              | ((chosen_detection == old_detection) & (chosen_track < old_track)))
    destination = destination[better]
    chosen = chosen[better]
    first["track_id"].reshape(-1)[destination] = arrays["track_id"][chosen]
    for field in FIRST_FIELDS:
        source_values = quantities[field] if field in quantities else arrays[field]
        first[field].reshape(-1)[destination] = source_values[chosen]


def make_first_storage():
    output = {"track_id": np.full((4, EXPECTED_EVENTS), -1, dtype=np.int32)}
    for field in FIRST_FIELDS:
        dtype = np.int32 if field == "n_boundary_encounters" else np.float64
        fill = -1 if dtype == np.int32 else np.nan
        output[field] = np.full((4, EXPECTED_EVENTS), fill, dtype=dtype)
    return output


def analyze_cell(payload):
    cell, cell_index, scratch_dir = payload
    n_d = len(D_EDGES_MM) - 1
    n_t = len(TPROP_EDGES_NS) - 1
    count = np.zeros((4, n_d), dtype=np.int64)
    sum_t = np.zeros((4, n_d), dtype=np.float64)
    sum_t2 = np.zeros((4, n_d), dtype=np.float64)
    hist2d = np.zeros((4, n_d, n_t), dtype=np.int64)
    under_over = np.zeros((4, 4), dtype=np.int64)

    event_count = np.zeros((4, EXPECTED_EVENTS), dtype=np.int32)
    event_sum_t = np.zeros((4, EXPECTED_EVENTS), dtype=np.float64)
    totals = {name: np.zeros(4, dtype=np.float64) for name in
              ("count", "sum_t", "sum_t2", "sum_d", "sum_d2", "sum_path",
               "sum_rho", "sum_v")}
    stratum_event_count = np.zeros((4, 5, EXPECTED_EVENTS), dtype=np.int32)
    stratum_event_sum = np.zeros((4, 5, EXPECTED_EVENTS), dtype=np.float64)
    stratum_sum2 = np.zeros((4, 5), dtype=np.float64)
    angle_event_count = np.zeros((4, N_ANGLE_BINS, EXPECTED_EVENTS), dtype=np.int32)
    angle_event_sum = np.zeros((4, N_ANGLE_BINS, EXPECTED_EVENTS), dtype=np.float64)
    angle_sum2 = np.zeros((4, N_ANGLE_BINS), dtype=np.float64)
    first = make_first_storage()

    with uproot.open(cell["root_path"]) as root_file:
        tree = root_file[TREE_NAME]
        for arrays in tree.iterate(READ_BRANCHES, step_size=STEP_SIZE, library="np"):
            face = arrays["face_type"].astype(np.int64, copy=False)
            source = arrays["source_type"].astype(np.int64, copy=False)
            selected = np.flatnonzero(
                ((face == LEFT_FACE) | (face == RIGHT_FACE))
                & ((source == SOURCE_CODES[0]) | (source == SOURCE_CODES[1])))
            if selected.size == 0:
                continue
            event_id = arrays["event_id"][selected].astype(np.int64, copy=False)
            require(np.all((event_id >= 0) & (event_id < EXPECTED_EVENTS)),
                    f"{cell['cell_id']}: event_id fuera de rango")
            group = ((source[selected] - 1) * 2 + face[selected]).astype(np.int64)
            dx = arrays["x_mm"] - arrays["x_creation_mm"]
            dy = arrays["y_mm"] - arrays["y_creation_mm"]
            dz = arrays["z_mm"] - arrays["z_creation_mm"]
            d_direct = np.sqrt(dx * dx + dy * dy + dz * dz)
            tprop = arrays["t_detection_ns"] - arrays["t_creation_ns"]
            path = arrays["path_length_mm"]
            require(np.all(tprop[selected] > 0.0), f"{cell['cell_id']}: tprop no positivo")
            require(np.all(d_direct[selected] > 0.0), f"{cell['cell_id']}: d_direct no positivo")
            rho = path / d_direct
            vapp = d_direct / tprop
            quantities = {
                "tprop_ns": tprop, "d_direct_mm": d_direct,
                "rho_detour": rho, "v_apparent_mm_per_ns": vapp,
            }
            update_first(first, arrays, selected, group, event_id, quantities)

            dsel = d_direct[selected]
            tsel = tprop[selected]
            psel = path[selected]
            rsel = rho[selected]
            vsel = vapp[selected]
            flat_event_group = group * EXPECTED_EVENTS + event_id
            event_count.reshape(-1)[:] += np.bincount(
                flat_event_group, minlength=4 * EXPECTED_EVENTS).astype(np.int32)
            event_sum_t.reshape(-1)[:] += np.bincount(
                flat_event_group, weights=tsel, minlength=4 * EXPECTED_EVENTS)
            for group_index in range(4):
                mask = group == group_index
                if not np.any(mask):
                    continue
                totals["count"][group_index] += np.count_nonzero(mask)
                totals["sum_t"][group_index] += np.sum(tsel[mask])
                totals["sum_t2"][group_index] += np.sum(tsel[mask] ** 2)
                totals["sum_d"][group_index] += np.sum(dsel[mask])
                totals["sum_d2"][group_index] += np.sum(dsel[mask] ** 2)
                totals["sum_path"][group_index] += np.sum(psel[mask])
                totals["sum_rho"][group_index] += np.sum(rsel[mask])
                totals["sum_v"][group_index] += np.sum(vsel[mask])

            d_bin = np.searchsorted(D_EDGES_MM, dsel, side="right") - 1
            t_bin = np.searchsorted(TPROP_EDGES_NS, tsel, side="right") - 1
            valid_d = (d_bin >= 0) & (d_bin < n_d)
            flat_d = group[valid_d] * n_d + d_bin[valid_d]
            count.reshape(-1)[:] += np.bincount(flat_d, minlength=4 * n_d)
            sum_t.reshape(-1)[:] += np.bincount(
                flat_d, weights=tsel[valid_d], minlength=4 * n_d)
            sum_t2.reshape(-1)[:] += np.bincount(
                flat_d, weights=tsel[valid_d] ** 2, minlength=4 * n_d)
            valid_2d = valid_d & (t_bin >= 0) & (t_bin < n_t)
            flat_2d = ((group[valid_2d] * n_d + d_bin[valid_2d]) * n_t
                       + t_bin[valid_2d])
            hist2d.reshape(-1)[:] += np.bincount(
                flat_2d, minlength=4 * n_d * n_t)
            for group_index in range(4):
                gm = group == group_index
                under_over[group_index] += [
                    np.count_nonzero(gm & (d_bin < 0)),
                    np.count_nonzero(gm & (d_bin >= n_d)),
                    np.count_nonzero(gm & (t_bin < 0)),
                    np.count_nonzero(gm & (t_bin >= n_t)),
                ]

            bcat = boundary_category(arrays["n_boundary_encounters"][selected])
            acat = np.searchsorted(ANGLE_EDGES_DEG,
                                   arrays["exit_angle_deg"][selected], side="right") - 1
            require(np.all((acat >= 0) & (acat < N_ANGLE_BINS)),
                    f"{cell['cell_id']}: ángulo fuera de rango")
            flat_b = (group * 5 + bcat) * EXPECTED_EVENTS + event_id
            stratum_event_count.reshape(-1)[:] += np.bincount(
                flat_b, minlength=4 * 5 * EXPECTED_EVENTS).astype(np.int32)
            stratum_event_sum.reshape(-1)[:] += np.bincount(
                flat_b, weights=tsel, minlength=4 * 5 * EXPECTED_EVENTS)
            stratum_sum2.reshape(-1)[:] += np.bincount(
                group * 5 + bcat, weights=tsel ** 2, minlength=4 * 5)
            flat_a = (group * N_ANGLE_BINS + acat) * EXPECTED_EVENTS + event_id
            angle_event_count.reshape(-1)[:] += np.bincount(
                flat_a, minlength=4 * N_ANGLE_BINS * EXPECTED_EVENTS).astype(np.int32)
            angle_event_sum.reshape(-1)[:] += np.bincount(
                flat_a, weights=tsel, minlength=4 * N_ANGLE_BINS * EXPECTED_EVENTS)
            angle_sum2.reshape(-1)[:] += np.bincount(
                group * N_ANGLE_BINS + acat, weights=tsel ** 2,
                minlength=4 * N_ANGLE_BINS)

    rows = []
    boundary_rows = []
    angle_rows = []
    for source_index, source_code in enumerate(SOURCE_CODES):
        for face_code in (LEFT_FACE, RIGHT_FACE):
            group_index = 2 * source_index + face_code
            n = int(totals["count"][group_index])
            require(n > 0, f"{cell['cell_id']}: grupo vacío {source_code}/{face_code}")
            mean_t, se_t = cluster_mean_se(event_count[group_index],
                                           event_sum_t[group_index])
            variance = max(0.0, totals["sum_t2"][group_index] / n - mean_t ** 2)
            nominal_d = 700.0 + int(cell["x_mm"]) if face_code == LEFT_FACE \
                else 700.0 - int(cell["x_mm"])
            rows.append({
                "cell_id": cell["cell_id"], "material": cell["material"],
                "x_mm": int(cell["x_mm"]), "face": FACE_NAMES[face_code],
                "source_type": source_code, "source": SOURCE_NAMES[source_code],
                "nominal_d_mm": nominal_d, "n_photons": n,
                "mean_d_direct_mm": totals["sum_d"][group_index] / n,
                "sd_d_direct_mm": np.sqrt(max(0.0, totals["sum_d2"][group_index] / n
                                                  - (totals["sum_d"][group_index] / n) ** 2)),
                "mean_tprop_ns": mean_t, "se_tprop_cluster_ns": se_t,
                "var_tprop_ns2": variance,
                "mean_path_length_mm": totals["sum_path"][group_index] / n,
                "mean_rho_detour": totals["sum_rho"][group_index] / n,
                "mean_v_apparent_mm_per_ns": totals["sum_v"][group_index] / n,
            })
            for category, label in enumerate(BOUNDARY_LABELS):
                cat_n = int(np.sum(stratum_event_count[group_index, category]))
                cat_mean, cat_se = cluster_mean_se(
                    stratum_event_count[group_index, category],
                    stratum_event_sum[group_index, category])
                cat_var = (max(0.0, stratum_sum2[group_index, category] / cat_n
                               - cat_mean ** 2) if cat_n else np.nan)
                boundary_rows.append({**rows[-1], "boundary_category": label,
                                      "stratum_n": cat_n,
                                      "stratum_mean_tprop_ns": cat_mean,
                                      "stratum_se_cluster_ns": cat_se,
                                      "stratum_var_tprop_ns2": cat_var})
            for category in range(N_ANGLE_BINS):
                cat_n = int(np.sum(angle_event_count[group_index, category]))
                cat_mean, cat_se = cluster_mean_se(
                    angle_event_count[group_index, category],
                    angle_event_sum[group_index, category])
                cat_var = (max(0.0, angle_sum2[group_index, category] / cat_n
                               - cat_mean ** 2) if cat_n else np.nan)
                angle_rows.append({**rows[-1],
                                   "angle_low_deg": ANGLE_EDGES_DEG[category],
                                   "angle_high_deg": ANGLE_EDGES_DEG[category + 1],
                                   "stratum_n": cat_n,
                                   "stratum_mean_tprop_ns": cat_mean,
                                   "stratum_se_cluster_ns": cat_se,
                                   "stratum_var_tprop_ns2": cat_var})

    scratch_path = Path(scratch_dir) / f"{cell['cell_id']}.npz"
    np.savez_compressed(scratch_path, count=count, sum_t=sum_t, sum_t2=sum_t2,
                        hist2d=hist2d, under_over=under_over,
                        **{f"first_{key}": value for key, value in first.items()})
    return rows, boundary_rows, angle_rows, str(scratch_path), cell_index


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN_DIR)
    parser.add_argument("--processes", type=int, default=PROCESS_COUNT)
    args = parser.parse_args()
    require(args.processes >= 1, "processes debe ser positivo")
    start = datetime.now(timezone.utc)
    cells = discover_cells(args.campaign.resolve())
    require(len(cells) == EXPECTED_CELLS, "se esperaban 21 celdas")
    cells.sort(key=lambda item: (item["material"], int(item["x_mm"])))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    scratch = OUTPUT_DIR / "scratch"
    scratch.mkdir(exist_ok=True)
    payloads = [(cell, index, scratch) for index, cell in enumerate(cells)]
    with mp.get_context("spawn").Pool(args.processes) as pool:
        results = pool.map(analyze_cell, payloads)

    cell_rows, boundary_rows, angle_rows = [], [], []
    n_d = len(D_EDGES_MM) - 1
    n_t = len(TPROP_EDGES_NS) - 1
    aggregate_count = np.zeros((3, 4, n_d), dtype=np.int64)
    aggregate_sum = np.zeros((3, 4, n_d), dtype=np.float64)
    aggregate_sum2 = np.zeros((3, 4, n_d), dtype=np.float64)
    aggregate_hist = np.zeros((3, 4, n_d, n_t), dtype=np.int64)
    aggregate_under_over = np.zeros((3, 4, 4), dtype=np.int64)
    first_columns = {"cell_code": [], "material_code": [], "gun_x_mm": [],
                     "event_id": [], "face_type": [], "source_type": []}
    for field in ("track_id",) + FIRST_FIELDS:
        first_columns[field] = []
    for rows, brows, arows, scratch_path, cell_index in results:
        cell_rows.extend(rows)
        boundary_rows.extend(brows)
        angle_rows.extend(arows)
        cell = cells[cell_index]
        material_code = MATERIAL_CODES[cell["material"]]
        data = np.load(scratch_path)
        aggregate_count[material_code] += data["count"]
        aggregate_sum[material_code] += data["sum_t"]
        aggregate_sum2[material_code] += data["sum_t2"]
        aggregate_hist[material_code] += data["hist2d"]
        aggregate_under_over[material_code] += data["under_over"]
        for group_index in range(4):
            source_index, face_code = divmod(group_index, 2)
            valid = data["first_track_id"][group_index] >= 0
            count_valid = int(np.count_nonzero(valid))
            first_columns["cell_code"].append(np.full(count_valid, cell_index, np.int32))
            first_columns["material_code"].append(np.full(count_valid, material_code, np.int32))
            first_columns["gun_x_mm"].append(
                np.full(count_valid, int(cell["x_mm"]), np.int32))
            first_columns["event_id"].append(np.flatnonzero(valid).astype(np.int32))
            first_columns["face_type"].append(np.full(count_valid, face_code, np.int32))
            first_columns["source_type"].append(np.full(count_valid, source_index + 1, np.int32))
            for field in ("track_id",) + FIRST_FIELDS:
                first_columns[field].append(data[f"first_{field}"][group_index, valid])

    write_csv(OUTPUT_DIR / "all_photon_cell_face.csv", cell_rows)
    write_csv(OUTPUT_DIR / "boundary_strata.csv", boundary_rows)
    write_csv(OUTPUT_DIR / "angle_strata.csv", angle_rows)
    np.savez_compressed(OUTPUT_DIR / "microscopic_transport_bins.npz",
                        d_edges_mm=D_EDGES_MM, tprop_edges_ns=TPROP_EDGES_NS,
                        count=aggregate_count, sum_tprop_ns=aggregate_sum,
                        sum_tprop2_ns2=aggregate_sum2, hist2d=aggregate_hist,
                        under_over=aggregate_under_over)
    first_columns = {key: np.concatenate(values) for key, values in first_columns.items()}
    with uproot.recreate(OUTPUT_DIR / "first_by_source.root") as root_file:
        root_file["first_by_source"] = first_columns
    for _, _, _, scratch_path, _ in results:
        Path(scratch_path).unlink()
    scratch.rmdir()
    finish = datetime.now(timezone.utc)
    metadata = {
        "created_utc": finish.isoformat(), "started_utc": start.isoformat(),
        "wall_seconds": (finish - start).total_seconds(),
        "command": ("env PYTHONPATH=analysis/track_mechanism_20260915 python3 "
                    "analysis/track_mechanism_20260915/build_step3_transport.py --processes "
                    f"{args.processes}"),
        "campaign": str(args.campaign.resolve()), "cells": EXPECTED_CELLS,
        "events_per_cell": EXPECTED_EVENTS, "processes": args.processes,
        "step_size": STEP_SIZE, "sources": SOURCE_NAMES,
        "d_edges_mm": D_EDGES_MM.tolist(),
        "tprop_edges_ns": TPROP_EDGES_NS.tolist(),
        "boundary_categories": BOUNDARY_LABELS,
        "angle_edges_deg": ANGLE_EDGES_DEG.tolist(),
        "first_rows": len(first_columns["event_id"]),
    }
    (OUTPUT_DIR / "transport_build.meta.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
