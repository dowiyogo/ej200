#!/usr/bin/env python3
"""Fail-closed input verification for the completed EXEC42 grid."""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import uproot


REQUIRED_TREES = {
    "sipm_hits": {"event_id", "global_id", "time_ns"},
    "event_observables": {"event_id", "produced_scint", "edep_total_MeV", "edep_optical_MeV"},
    "first_bar_encounters": {"event_id", "track_id", "source", "post_volume_id", "outcome",
                             "exiting_bar", "cos_incidence"},
    "sipm_event_counts": {"event_id", "global_id", "incident", "matched_detected",
                          "unmatched_detected", "expected_surface_pde_sum"},
}


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(16*1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--a2-meta", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    campaign = json.loads(args.campaign.read_text())
    a2 = json.loads(args.a2_meta.read_text())
    manifest = [json.loads(line) for line in args.manifest.read_text().splitlines() if line.strip()]
    completed = [row for row in manifest if row["status"] == "SIMULATION_COMPLETE"]
    failed = [row for row in manifest if row["status"] == "FAILED"]
    errors = []
    if len(completed) != 21 or failed:
        errors.append(f"manifest completion/failure count is {len(completed)}/{len(failed)}")
    if len({row["cell_id"] for row in completed}) != 21:
        errors.append("manifest does not have one completed row per cell")
    roots = [row["root_path"] for row in completed]
    if len(set(roots)) != 21:
        errors.append("completed ROOT paths are not unique")
    rows = []
    for index, row in enumerate(sorted(completed, key=lambda item: item["cell_id"]), 1):
        cell = next(cell for cell in campaign["cells"] if cell["cell_id"] == row["cell_id"])
        done_path = Path(cell["output"])/".DONE"
        if not done_path.is_file():
            errors.append(row["cell_id"]+": missing .DONE")
            continue
        done = json.loads(done_path.read_text())
        for key in ("root_path", "root_sha256", "root_size_bytes", "workers", "eventModulo",
                    "seeds", "PDE_path", "PDE_sha256", "exit_code", "events_run"):
            if done.get(key) != row.get(key):
                errors.append(f"{row['cell_id']}: .DONE/manifest mismatch in {key}")
        root_path = Path(row["root_path"])
        if not root_path.is_file() or root_path.stat().st_size != row["root_size_bytes"]:
            errors.append(row["cell_id"]+": ROOT missing or wrong size")
            continue
        actual_hash = sha(root_path)
        print(f"HASH {index:02d}/21 {row['cell_id']} {actual_hash}", flush=True)
        if actual_hash != row["root_sha256"]:
            errors.append(row["cell_id"]+": ROOT SHA256 mismatch")
        try:
            root = uproot.open(root_path)
            tree_entries = {}
            for name, branches in REQUIRED_TREES.items():
                tree = root[name]
                if not branches <= set(tree.keys()):
                    errors.append(f"{row['cell_id']}: missing {name} branches")
                tree_entries[name] = int(tree.num_entries)
            event_ids = root["event_observables"]["event_id"].array(library="np")
            if len(event_ids) != 10000 or not np.array_equal(np.sort(event_ids), np.arange(10000)):
                errors.append(row["cell_id"]+": event_observables is not exactly events 0..9999")
            if tree_entries["sipm_event_counts"] != 860000:
                errors.append(row["cell_id"]+": sipm_event_counts is not 10000*86")
        except Exception as error:
            errors.append(f"{row['cell_id']}: ROOT read failed: {error}")
            tree_entries = {}
        expected = {
            "events_run": 10000, "N_generated": 10000, "workers": 4, "eventModulo": 1,
            "seeds": [26092601, 8349041], "PDE_path": a2["PDE_path"],
            "PDE_sha256": a2["PDE_sha256"], "binary_sha256": a2["binary_sha256"],
            "simulation_commit": "b35ee84acadef12c93506e0720580c91f901fcbf",
            "diagnostics": False, "exit_code": 0,
        }
        for key, value in expected.items():
            if row.get(key) != value:
                errors.append(f"{row['cell_id']}: {key}={row.get(key)!r}, expected {value!r}")
        rows.append({
            "cell_id": row["cell_id"], "material": row["material"], "opsc_code": row["opsc_code"],
            "x_mm": row["x_mm"], "events": row["events_run"], "workers_manifest": row["workers"],
            "eventModulo_manifest": row["eventModulo"], "seed1": row["seeds"][0],
            "seed2": row["seeds"][1], "PDE_path": row["PDE_path"],
            "PDE_sha256": row["PDE_sha256"], "root_path": row["root_path"],
            "root_sha256": actual_hash, "root_bytes": row["root_size_bytes"],
            "event_rows": tree_entries.get("event_observables", -1),
            "encounter_rows": tree_entries.get("first_bar_encounters", -1),
            "sipm_event_rows": tree_entries.get("sipm_event_counts", -1),
            "hit_rows": tree_entries.get("sipm_hits", -1), "wall_s": row["wall_s"],
        })
    rows.sort(key=lambda item: item["cell_id"])
    csv_path = args.out/"input_inventory.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    arrays = {key: np.asarray([row[key] for row in rows]) for key in rows[0]}
    with uproot.recreate(args.out/"input_inventory.root") as output:
        output["data"] = arrays
    metadata = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "campaign": str(args.campaign), "campaign_sha256": sha(args.campaign),
        "manifest": str(args.manifest), "manifest_sha256": sha(args.manifest),
        "A2_meta": str(args.a2_meta), "A2_meta_sha256": sha(args.a2_meta),
        "rows": len(rows), "csv_sha256": sha(csv_path),
        "root_sha256": sha(args.out/"input_inventory.root"), "script_sha256": sha(Path(__file__)),
        "hash_method": "SHA256 recomputed over every complete ROOT",
    }
    dump(args.out/"input_inventory.meta.json", metadata)
    alternatives = [
        {"path": "/home/rrios/exec34r_20260912", "selection": "NOT_SELECTED",
         "reason": "older EXEC34R manifest, simulation commit 420addf, before production observations"},
        {"path": str(args.campaign.parent), "selection": "SELECTED",
         "reason": "EXEC42 manifest, production commit b35ee84, completed 2026-09-14"},
    ]
    result = {"metadata": metadata, "status": "PASS" if not errors else "FAIL",
              "completed": len(rows), "failed_manifest_rows": len(failed),
              "unique_root_paths": len(set(item["root_path"] for item in rows)),
              "PDE_paths": sorted(set(item["PDE_path"] for item in rows)),
              "PDE_hashes": sorted(set(item["PDE_sha256"] for item in rows)),
              "workers_from_manifest": sorted(set(item["workers_manifest"] for item in rows)),
              "seeds": sorted(set((item["seed1"], item["seed2"]) for item in rows)),
              "campaign_selection": alternatives, "errors": errors}
    dump(args.out/"input_verification.json", result)
    print(json.dumps(result, indent=2), flush=True)
    if errors:
        raise SystemExit("EXEC42 input verification failed")


if __name__ == "__main__":
    main()
