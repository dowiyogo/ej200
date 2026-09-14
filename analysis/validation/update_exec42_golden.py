#!/usr/bin/env python3
"""Append the fail-closed EXEC42 full-grid result to the golden reference."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path


EXPECTED_PREREGISTRATIONS = {
    "EXEC40": "cd13d28db1811e5c74d9b17d60d1be7aa0971f02ff8768d48528126602fb5531",
    "EXEC41": "933433a2882c4753301c675913df1f37db614e4afc8316ae87df6c122021d418",
    "EXEC42_H3": "37e15cf4a26271075babc90a1efc2f1c0d0a8ce1a0be6ea64e48e1d78053d6f4",
}


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--golden", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    args = parser.parse_args()

    preregistration_paths = {
        "EXEC40": args.repo / "analysis/validation/EXEC40_PREREGISTRATION.md",
        "EXEC41": args.repo / "analysis/validation/EXEC41_PREREGISTRATION.md",
        "EXEC42_H3": args.repo / "analysis/validation/EXEC42_H3_PREREGISTRATION.md",
    }
    preregistrations = {}
    for name, path in preregistration_paths.items():
        actual = sha(path)
        if actual != EXPECTED_PREREGISTRATIONS[name]:
            raise RuntimeError(f"{name} preregistration changed: {actual}")
        preregistrations[name] = {"path": str(path), "sha256": actual}

    verification_path = args.analysis / "input_verification.json"
    results_path = args.analysis / "contract_results.json"
    verification = json.loads(verification_path.read_text())
    results = json.loads(results_path.read_text())
    if verification["status"] != "PASS" or verification["completed"] != 21:
        raise RuntimeError("input verification is not a complete PASS")
    if results["metadata"]["N_cells"] != 21 or results["metadata"]["N_per_cell"] != 10000:
        raise RuntimeError("contract result does not describe 21 x 10000 events")

    sidecar_stems = ["event_metrics", "face_decomposition", "angular_distribution"]
    expected = []
    for cell_id in sorted(results["cells"]):
        directory = args.analysis / "cells" / cell_id
        expected.append(directory / "cell_results.json")
        for stem in sidecar_stems:
            expected.extend(directory / f"{stem}.{suffix}" for suffix in ("csv", "root", "meta.json"))
    for stem in ("input_inventory", "contract_cells", "H3_material_invariance", "escape_profiles"):
        expected.extend(args.analysis / f"{stem}.{suffix}" for suffix in ("csv", "root", "meta.json"))
    expected.extend([verification_path, results_path, args.analysis / "analysis.log"])
    missing = [str(path) for path in expected if not path.is_file()]
    if missing:
        raise RuntimeError("missing EXEC42 artifacts: " + ", ".join(missing))
    artifact_hashes = {str(path): sha(path) for path in expected}

    required = results["status_counts"]
    no_not_evaluable = not any(
        value == "NOT EVALUABLE"
        for cell in results["cells"].values()
        for value in (
            cell["V1_primary"]["status"],
            cell["V2_H1"]["status"],
            cell["V5_matched"]["status"],
        )
    )
    mechanically_ready = (
        required["V1_PASS"] == 21
        and required["V2_H1_PASS"] == 21
        and required["V5_PASS"] == 21
        and required["H3_pair_PASS"] == required["H3_pair_total"] == 21
        and no_not_evaluable
    )
    if mechanically_ready != results["ready_for_acceptance"]:
        raise RuntimeError("ready_for_acceptance does not follow the frozen rule")

    golden = json.loads(args.golden.read_text())
    now = dt.datetime.now(dt.timezone.utc).isoformat()
    golden.update({
        "title": "EXEC38 golden reference extended through the EXEC42 full production grid",
        "updated_utc": now,
        "ready_for_acceptance": mechanically_ready,
        "overall_status": "EXEC42_FULL_GRID_PASS" if mechanically_ready else "EXEC42_FULL_GRID_H3_FAIL",
    })
    golden["exec42_full_grid"] = {
        "scope": "21 EndTop70 cells: three scintillators by seven x positions, 10000 events each",
        "created_utc": now,
        "simulation_commit": "b35ee84acadef12c93506e0720580c91f901fcbf",
        "analysis_commit": "ef3092ab2ea7c6d628092c8034131adf5512d2cb",
        "campaign": "/home/rrios/exec42_20260913/grid/campaign.json",
        "manifest": "/home/rrios/exec42_20260913/grid/manifest.jsonl",
        "input_verification": {
            "path": str(verification_path),
            "sha256": sha(verification_path),
            "result": verification,
        },
        "preregistrations": preregistrations,
        "analysis": {
            "path": str(args.analysis),
            "input_command": "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/validation/verify_exec42_inputs.py --campaign /home/rrios/exec42_20260913/grid/campaign.json --manifest /home/rrios/exec42_20260913/grid/manifest.jsonl --a2-meta /home/rrios/exec42_20260913/a2_validation/invocation.meta.json --out /home/rrios/exec42_20260913/analysis",
            "contract_command": "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python analysis/validation/analyze_exec42_grid.py --inventory /home/rrios/exec42_20260913/analysis/input_inventory.csv --out /home/rrios/exec42_20260913/analysis --preregistration analysis/validation/EXEC42_H3_PREREGISTRATION.md",
            "contract_results": str(results_path),
            "contract_results_sha256": sha(results_path),
            "bootstrap": {
                "seed": results["metadata"]["bootstrap_seed"],
                "replicates": results["metadata"]["bootstrap_replicates"],
                "unit": "generated event; common resampling weights across cells",
            },
        },
        "physical_observable_contract": {
            "V1": {
                "definition": "sum(produced scintillation photons)/(nominal material yield * sum(total deposited energy in BarLV))",
                "acceptance": "unity lies inside the three-event-bootstrap-SE interval",
            },
            "V2_H1": {
                "definition": "portable-outcome escape to air or world at the first physical BarLV encounter of scintillation photons exiting BarLV, divided by all such first encounters",
                "acceptance": "escape fraction minus three event-bootstrap SE is at least 0.226172",
            },
            "V2_H2": {
                "definition": "20-bin signed incidence-cosine distribution on the same first-encounter population",
                "acceptance": "descriptive; no PASS/FAIL under EXEC41",
            },
            "V5_matched": {
                "definition": "matched boundary detections divided by independently observed sensor incidents",
                "prediction": "sum(expected surface PDE over incidents)/number of incidents",
                "acceptance": "measured-minus-predicted is compatible with zero within three paired event-bootstrap SE",
            },
            "H3": {
                "definition": "V2 escape-fraction equality for every material pair at each common x position",
                "acceptance": "absolute paired difference is no more than three paired-bootstrap SE in all 21 comparisons",
            },
            "Npe": "mean independent sensitive-detector detected count per END side",
            "first_face_counts": "counts by physical destination volume at the first BarLV encounter",
        },
        "status_counts": required,
        "H3_status": results["H3_status"],
        "failure_reasons": results["failure_reasons"],
        "not_evaluable_count": 0 if no_not_evaluable else None,
        "ready_for_acceptance": mechanically_ready,
        "paired_design": results["paired_design"],
        "cells": results["cells"],
        "H3": results["H3"],
        "sidecars": {
            "cell_count": len(results["cells"]),
            "triplets_per_cell": sidecar_stems,
            "expected_artifact_count": len(expected),
            "all_present": True,
        },
        "artifact_hashes": artifact_hashes,
    }
    args.golden.write_text(json.dumps(golden, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "golden": str(args.golden),
        "ready_for_acceptance": mechanically_ready,
        "overall_status": golden["overall_status"],
        "artifact_count": len(expected),
        "golden_sha256": sha(args.golden),
    }, indent=2))


if __name__ == "__main__":
    main()
