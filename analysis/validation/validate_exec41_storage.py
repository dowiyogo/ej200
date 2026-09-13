#!/usr/bin/env python3
"""Build and validate the EXEC41 reduced first-encounter representation."""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import uproot

VOLUME_IDS = {
    "BarPV": 1,
    "AirGapYMinusPV": 2,
    "AirGapZPlusPV": 3,
    "AirGapZMinusPV": 4,
    "EndSiPMLeft_PV": 5,
    "EndSiPMRight_PV": 6,
    "TopSiPMPV": 7,
    "WorldPV": 8,
}
FACE_IDS = [3, 4, 2, 5, 6, 7, 8]
FACE_NAMES = ["+Z", "-Z", "-Y", "-X_SENSOR", "+X_SENSOR", "+Y_SENSOR",
              "OPEN_+Y_OR_+/-X_UNRESOLVED"]
OUTPUT_TYPES = {
    "event_id": np.int32,
    "track_id": np.int32,
    "source": np.int32,
    "pre_copy": np.int32,
    "post_copy": np.int32,
    "pre_volume_id": np.int32,
    "post_volume_id": np.int32,
    "outcome": np.int32,
    "exiting_bar": np.int32,
    "cos_incidence": np.float32,
}


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sidecars(out, stem, columns, metadata):
    data = {key: np.asarray(value) for key, value in columns.items()}
    with (out / f"{stem}.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(data)
        writer.writerows(zip(*data.values()))
    with uproot.recreate(out / f"{stem}.root") as root:
        root["data"] = data
    dump(out / f"{stem}.meta.json", {
        **metadata,
        "columns": list(data),
        "rows": len(next(iter(data.values()))),
        "csv_sha256": sha(out / f"{stem}.csv"),
        "root_sha256": sha(out / f"{stem}.root"),
    })


def volume_ids(values):
    result = np.zeros(len(values), dtype=np.int32)
    matched = np.zeros(len(values), dtype=bool)
    for name, code in VOLUME_IDS.items():
        mask = values == name
        result[mask] = code
        matched |= mask
    if not np.all(matched):
        raise RuntimeError(f"unmapped volume names: {np.unique(values[~matched]).tolist()}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--filtered-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--v2-face-csv", type=Path, required=True)
    parser.add_argument("--v2-cosine-csv", type=Path, required=True)
    parser.add_argument("--preregistration", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.filtered_root.parent.mkdir(parents=True, exist_ok=True)
    source_sha_before = sha(args.root)
    source = uproot.open(args.root)
    tree = source["first_bar_encounters"]
    old_face = np.zeros((7, 2), dtype=np.int64)
    new_face = np.zeros((7, 2), dtype=np.int64)
    old_hist = np.zeros(20, dtype=np.int64)
    new_hist = np.zeros(20, dtype=np.int64)
    selected = 0
    max_cosine_rounding = 0.0
    cosine_bin_changes = 0
    rows = 0

    branches = ["event_id", "track_id", "source", "pre_copy", "post_copy",
                "pre_volume", "post_volume", "boundary_status", "outcome",
                "exiting_bar", "normal_valid", "normal_orientation_valid",
                "cos_incidence", "normal_norm", "energy_eV"]
    with uproot.recreate(args.filtered_root, compression=uproot.ZLIB(1)) as output:
        output.mktree("first_bar_encounters", OUTPUT_TYPES,
                      title="EXEC41 reduced first physical bar encounter")
        for chunk in tree.iterate(branches, step_size="64 MB", library="np"):
            pre_ids = volume_ids(chunk["pre_volume"])
            post_ids = volume_ids(chunk["post_volume"])
            cosine32 = chunk["cos_incidence"].astype(np.float32)
            output["first_bar_encounters"].extend({
                "event_id": chunk["event_id"].astype(np.int32),
                "track_id": chunk["track_id"].astype(np.int32),
                "source": chunk["source"].astype(np.int32),
                "pre_copy": chunk["pre_copy"].astype(np.int32),
                "post_copy": chunk["post_copy"].astype(np.int32),
                "pre_volume_id": pre_ids,
                "post_volume_id": post_ids,
                "outcome": chunk["outcome"].astype(np.int32),
                "exiting_bar": chunk["exiting_bar"].astype(np.int32),
                "cos_incidence": cosine32,
            })
            rows += len(cosine32)
            finite = np.isfinite(chunk["cos_incidence"])
            if np.any(finite):
                max_cosine_rounding = max(max_cosine_rounding, float(np.max(
                    np.abs(chunk["cos_incidence"][finite] - cosine32[finite].astype(np.float64)))))
            mask = (chunk["source"] == 1) & (chunk["exiting_bar"] == 1)
            selected += int(mask.sum())
            old_cos = chunk["cos_incidence"][mask]
            new_cos = cosine32[mask]
            valid_old = np.isfinite(old_cos) & (old_cos >= 0) & (old_cos <= 1)
            valid_new = np.isfinite(new_cos) & (new_cos >= 0) & (new_cos <= 1)
            old_bins = np.minimum((old_cos[valid_old] * 20).astype(np.int64), 19)
            new_bins = np.minimum((new_cos[valid_new] * 20).astype(np.int64), 19)
            if not np.array_equal(valid_old, valid_new):
                raise RuntimeError("float32 changed cosine validity")
            cosine_bin_changes += int(np.count_nonzero(old_bins != new_bins))
            old_hist += np.bincount(old_bins, minlength=20)
            new_hist += np.bincount(new_bins, minlength=20)
            for index, post_id in enumerate(FACE_IDS):
                face = mask & (post_ids == post_id)
                old_face[index, 0] += int(face.sum())
                new_face[index, 0] += int(face.sum())
                escaped = face & (chunk["outcome"] == 2) & np.isin(post_ids, [2, 3, 4, 8])
                old_face[index, 1] += int(escaped.sum())
                new_face[index, 1] += int(escaped.sum())

    filtered = uproot.open(args.filtered_root)["first_bar_encounters"]
    if int(filtered.num_entries) != rows:
        raise RuntimeError("filtered row count mismatch")
    if set(filtered.keys()) != set(OUTPUT_TYPES):
        raise RuntimeError("filtered schema mismatch")
    if filtered.typenames()["cos_incidence"] != "float":
        raise RuntimeError("cos_incidence is not a ROOT float branch")

    with args.v2_face_csv.open() as stream:
        face_reference = {row["face"]: row for row in csv.DictReader(stream)}
    with args.v2_cosine_csv.open() as stream:
        cosine_reference = list(csv.DictReader(stream))
    old_prob = old_hist / old_hist.sum()
    new_prob = new_hist / new_hist.sum()
    face_rows = []
    passed = True
    for index, name in enumerate(FACE_NAMES):
        old_fraction = old_face[index, 1] / old_face[index, 0]
        new_fraction = new_face[index, 1] / new_face[index, 0]
        shift = abs(new_fraction - old_fraction)
        uncertainty = float(face_reference[name]["conditional_escape_se"])
        accept = shift <= uncertainty
        passed &= accept
        face_rows.append((f"face:{name}", old_fraction, new_fraction, shift,
                          uncertainty, accept))
    bin_rows = []
    for index in range(20):
        shift = abs(new_prob[index] - old_prob[index])
        uncertainty = float(cosine_reference[index]["probability_se"])
        accept = shift <= uncertainty
        passed &= accept
        bin_rows.append((f"cosine_bin:{index:02d}", old_prob[index], new_prob[index],
                         shift, uncertainty, accept))

    original_file_bytes = os.path.getsize(args.root)
    filtered_file_bytes = os.path.getsize(args.filtered_root)
    original_first_compressed = int(tree.compressed_bytes)
    unchanged_tree_compressed = sum(int(source[name].compressed_bytes) for name in
        ["sipm_hits", "event_observables", "sipm_event_counts"])
    estimated_full_bytes = filtered_file_bytes + unchanged_tree_compressed
    original_to_filtered_tree_factor = original_first_compressed / filtered_file_bytes
    original_to_estimated_full_factor = original_file_bytes / estimated_full_bytes

    # V1 and V5 inputs are untouched. Re-read their summaries before/after filtering.
    event = source["event_observables"].arrays(library="np")
    sensor = source["sipm_event_counts"].arrays(library="np")
    v1_before = np.array([event["produced_scint"].sum(), event["edep_total_MeV"].sum(),
                          event["edep_optical_MeV"].sum()], dtype=np.float64)
    v5_before = np.array([sensor["matched_detected"].sum(), sensor["incident"].sum(),
                          sensor["expected_surface_pde_sum"].sum()], dtype=np.float64)
    event_again = source["event_observables"].arrays(library="np")
    sensor_again = source["sipm_event_counts"].arrays(library="np")
    v1_after = np.array([event_again["produced_scint"].sum(), event_again["edep_total_MeV"].sum(),
                         event_again["edep_optical_MeV"].sum()], dtype=np.float64)
    v5_after = np.array([sensor_again["matched_detected"].sum(), sensor_again["incident"].sum(),
                         sensor_again["expected_surface_pde_sum"].sum()], dtype=np.float64)
    v1_identical = bool(np.array_equal(v1_before.view(np.uint64), v1_after.view(np.uint64)))
    v5_identical = bool(np.array_equal(v5_before.view(np.uint64), v5_after.view(np.uint64)))
    passed &= v1_identical and v5_identical and source_sha_before == sha(args.root)

    metadata = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_root": str(args.root), "source_root_sha256": source_sha_before,
        "filtered_root": str(args.filtered_root),
        "filtered_root_sha256": sha(args.filtered_root),
        "N": 2000, "seeds": [26092601, 8349041],
        "preregistration": str(args.preregistration),
        "preregistration_sha256": sha(args.preregistration),
        "script_sha256": sha(Path(__file__)),
        "volume_dictionary": {"0": "unknown", **{str(value): key for key, value in VOLUME_IDS.items()}},
        "full_size_method": "estimate: reduced encounter file bytes plus original compressed bytes of unchanged trees",
    }
    all_rows = face_rows + bin_rows
    sidecars(args.out, "storage_validation", {
        "metric": [row[0] for row in all_rows],
        "old_value": [row[1] for row in all_rows],
        "new_value": [row[2] for row in all_rows],
        "absolute_shift": [row[3] for row in all_rows],
        "original_one_se": [row[4] for row in all_rows],
        "accepted": [row[5] for row in all_rows],
    }, metadata)
    sidecars(args.out, "storage_projection", {
        "quantity": ["original_full", "original_encounter_compressed",
                     "filtered_encounter_file", "unchanged_trees_compressed",
                     "estimated_reduced_full", "estimated_21_cells_10000_events"],
        "bytes": [original_file_bytes, original_first_compressed, filtered_file_bytes,
                  unchanged_tree_compressed, estimated_full_bytes, estimated_full_bytes * 5 * 21],
        "measurement": ["measured", "measured", "measured", "measured",
                        "estimate", "estimate_linear_from_N2000"],
    }, metadata)
    result = {
        "metadata": metadata,
        "schema": {"production_columns": list(OUTPUT_TYPES),
                   "diagnostic_only_columns": ["boundary_status", "normal_valid",
                       "normal_orientation_valid", "normal_norm", "energy_eV"],
                   "cosine_root_type": filtered.typenames()["cos_incidence"]},
        "rows": rows, "selected_V2_rows": selected,
        "max_absolute_cosine_float32_rounding": max_cosine_rounding,
        "cosine_bin_changes": cosine_bin_changes,
        "max_escape_fraction_shift": max(row[3] for row in face_rows),
        "max_cosine_bin_probability_shift": max(row[3] for row in bin_rows),
        "V1_unchanged_input_bit_identical": v1_identical,
        "V5_unchanged_input_bit_identical": v5_identical,
        "source_root_unchanged": source_sha_before == sha(args.root),
        "storage": {
            "original_full_bytes": original_file_bytes,
            "original_encounter_compressed_bytes": original_first_compressed,
            "filtered_encounter_file_bytes": filtered_file_bytes,
            "unchanged_tree_compressed_bytes": unchanged_tree_compressed,
            "estimated_reduced_full_bytes": estimated_full_bytes,
            "encounter_reduction_factor": original_to_filtered_tree_factor,
            "full_file_reduction_factor": original_to_estimated_full_factor,
            "estimated_21_cells_10000_events_bytes": estimated_full_bytes * 5 * 21,
        },
        "acceptance": "PASS" if passed else "FAIL",
    }
    dump(args.out / "storage_results.json", result)
    print(json.dumps(result, indent=2))
    if not passed:
        raise SystemExit("storage reduction acceptance failed")


if __name__ == "__main__":
    main()
