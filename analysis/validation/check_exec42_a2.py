#!/usr/bin/env python3
"""Validate the EXEC42 production-schema pilot against the EXEC40/41 cell."""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import uproot


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ordered(tree):
    arrays = tree.arrays(library="np")
    order = np.argsort(arrays["event_id"], kind="stable")
    return {name: values[order] for name, values in arrays.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--exec41", type=Path, required=True)
    parser.add_argument("--d1", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    candidate = uproot.open(args.candidate)
    reference = uproot.open(args.reference)
    old_v2 = json.loads((args.exec41 / "v2_results.json").read_text())
    old_diagnosis = json.loads((args.exec41 / "diagnosis_results.json").read_text())
    old_storage = json.loads((args.exec41 / "storage_results.json").read_text())

    cand_event = ordered(candidate["event_observables"])
    ref_event = ordered(reference["event_observables"])
    count_fields = ["event_id", "produced_scint", "produced_optical", "detected_left",
                    "detected_right", "detected_top"]
    count_differences = {name: int(np.count_nonzero(cand_event[name] != ref_event[name]))
                         for name in count_fields}
    totals = {name: int(cand_event[name].sum()) for name in count_fields[1:]}
    npe_end = (totals["detected_left"] + totals["detected_right"]) / (2*len(cand_event["event_id"]))
    gp1 = {
        "status": "PASS" if max(count_differences.values()) == 0 and npe_end == 397.1525 else "FAIL",
        "npe_end": npe_end, "reference_npe_end": 397.1525,
        "exact_difference": npe_end - 397.1525,
        "events_with_different_counts": count_differences,
        "totals": totals,
    }

    produced = cand_event["produced_scint"].astype(np.float64).sum()
    total_edep = cand_event["edep_total_MeV"].sum()
    nonoptical_edep = total_edep - cand_event["edep_optical_MeV"].sum()
    v1_total = produced/(10400*total_edep)
    v1_nonoptical = produced/(10400*nonoptical_edep)

    face_counts = np.zeros((7, 2), dtype=np.int64)
    cosine_hist = np.zeros(20, dtype=np.int64)
    face_ids = [3, 4, 2, 5, 6, 7, 8]
    for chunk in candidate["first_bar_encounters"].iterate(
            ["source", "post_volume_id", "outcome", "exiting_bar", "cos_incidence"],
            step_size="64 MB", library="np"):
        selected = (chunk["source"] == 1) & (chunk["exiting_bar"] == 1)
        post = chunk["post_volume_id"]
        for index, volume_id in enumerate(face_ids):
            face = selected & (post == volume_id)
            face_counts[index, 0] += int(face.sum())
            face_counts[index, 1] += int(np.count_nonzero(
                face & (chunk["outcome"] == 2) & np.isin(post, [2, 3, 4, 8])))
        cosine = chunk["cos_incidence"][selected]
        valid = np.isfinite(cosine) & (cosine >= 0) & (cosine <= 1)
        bins = np.minimum((cosine[valid]*20).astype(np.int64), 19)
        cosine_hist += np.bincount(bins, minlength=20)
    v2_h1 = face_counts[:, 1].sum()/face_counts[:, 0].sum()

    sensor = candidate["sipm_event_counts"].arrays(library="np")
    v5_matched = sensor["matched_detected"].sum()/sensor["incident"].sum()
    v5_expected = sensor["expected_surface_pde_sum"].sum()/sensor["incident"].sum()
    old_all = next(row for row in old_diagnosis["arms"] if row["arm"] == "all")
    comparisons = {
        "V1_total": {"candidate": v1_total,
                     "reference": old_diagnosis["V1"]["total_denominator_ratio"]["value"],
                     "tolerance": old_diagnosis["V1"]["total_denominator_ratio"]["se"]},
        "V1_nonoptical": {"candidate": v1_nonoptical,
                          "reference": old_diagnosis["V1"]["nonoptical_denominator_ratio"]["value"],
                          "tolerance": old_diagnosis["V1"]["nonoptical_denominator_ratio"]["se"]},
        "V2_H1": {"candidate": v2_h1, "reference": old_v2["H1"]["value"],
                  "tolerance": old_v2["H1"]["se"]},
        "V5_matched": {"candidate": v5_matched,
                       "reference": old_all["matched_over_incident"]["value"],
                       "tolerance": old_all["matched_over_incident"]["se"]},
        "V5_expected": {"candidate": v5_expected,
                        "reference": old_all["incident_spectrum_pde"]["value"],
                        "tolerance": old_all["incident_spectrum_pde"]["se"]},
    }
    for value in comparisons.values():
        value["absolute_difference"] = abs(value["candidate"] - value["reference"])
        value["status"] = "PASS" if value["absolute_difference"] <= value["tolerance"] else "FAIL"
    gp2 = {"status": "PASS" if all(value["status"] == "PASS" for value in comparisons.values()) else "FAIL",
           "comparisons": comparisons,
           "face_counts": face_counts.tolist(), "cosine_hist": cosine_hist.tolist()}

    actual_bytes = os.path.getsize(args.candidate)
    estimate_bytes = old_storage["storage"]["estimated_reduced_full_bytes"]
    d1_bytes = os.path.getsize(args.d1)
    gp3 = {"status": "MEASURED_ESTIMATE_UNDERSHOT",
           "actual_root_bytes": actual_bytes, "EXEC41_estimated_bytes": estimate_bytes,
           "actual_over_estimate": actual_bytes/estimate_bytes,
           "actual_over_D1": actual_bytes/d1_bytes, "D1_bytes": d1_bytes,
           "tree_compressed_bytes": {name: int(candidate[name].compressed_bytes)
              for name in ["sipm_hits", "event_observables", "first_bar_encounters", "sipm_event_counts"]}}

    metadata = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "candidate": str(args.candidate), "candidate_sha256": sha(args.candidate),
        "reference": str(args.reference), "reference_sha256": sha(args.reference),
        "N": 2000, "seeds": [26092601, 8349041], "workers": 4, "eventModulo": 1,
        "script_sha256": sha(Path(__file__)),
    }
    result = {"metadata": metadata, "G_P_1": gp1, "G_P_2": gp2, "G_P_3": gp3,
              "grid_preparation_allowed": gp1["status"] == "PASS"}
    (args.out/"a2_results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    rows = [("G-P.1 Npe/end", npe_end, 397.1525, abs(npe_end-397.1525), 0.0, gp1["status"])]
    rows += [(name, value["candidate"], value["reference"], value["absolute_difference"],
              value["tolerance"], value["status"]) for name, value in comparisons.items()]
    with (args.out/"a2_validation.csv").open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(
            ["metric", "candidate", "reference", "absolute_difference", "tolerance", "status"])
        writer.writerows(rows)
    arrays = {"metric": np.asarray([row[0] for row in rows]),
              "candidate": np.asarray([row[1] for row in rows]),
              "reference": np.asarray([row[2] for row in rows]),
              "absolute_difference": np.asarray([row[3] for row in rows]),
              "tolerance": np.asarray([row[4] for row in rows]),
              "status": np.asarray([row[5] for row in rows])}
    with uproot.recreate(args.out/"a2_validation.root") as output:
        output["data"] = arrays
    (args.out/"a2_validation.meta.json").write_text(json.dumps({
        **metadata, "csv_sha256": sha(args.out/"a2_validation.csv"),
        "root_sha256": sha(args.out/"a2_validation.root"), "rows": len(rows),
    }, indent=2)+"\n")
    print(json.dumps(result, indent=2))
    if gp1["status"] != "PASS":
        raise SystemExit("G-P.1 failed: grid preparation is prohibited")


if __name__ == "__main__":
    main()
