#!/usr/bin/env python3
"""Append the preregistered EXEC41 single-cell results to the golden reference."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--golden", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--preregistration", type=Path, required=True)
    parser.add_argument("--source-parent-commit", required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    golden = json.loads(args.golden.read_text())
    v2 = json.loads((args.results_dir / "v2_results.json").read_text())
    diagnosis = json.loads((args.results_dir / "diagnosis_results.json").read_text())
    storage = json.loads((args.results_dir / "storage_results.json").read_text())
    all_arm = next(row for row in diagnosis["arms"] if row["arm"] == "all")
    artifact_names = [
        "v2_results.json", "v2_face_decomposition.csv", "v2_face_decomposition.root",
        "v2_face_decomposition.meta.json", "v2_uniform_cosine.csv",
        "v2_uniform_cosine.root", "v2_uniform_cosine.meta.json",
        "v2_face_events.csv", "v2_face_events.root", "v2_face_events.meta.json",
        "diagnosis_results.json", "detection_gap_by_sensor.csv",
        "detection_gap_by_sensor.root", "detection_gap_by_sensor.meta.json",
        "detection_gap_by_event.csv", "detection_gap_by_event.root",
        "detection_gap_by_event.meta.json", "detection_gap_by_arm.csv",
        "detection_gap_by_arm.root", "detection_gap_by_arm.meta.json",
        "v1_event_accounting.csv", "v1_event_accounting.root",
        "v1_event_accounting.meta.json", "code_path_audit.csv",
        "code_path_audit.root", "code_path_audit.meta.json",
        "storage_results.json", "storage_validation.csv", "storage_validation.root",
        "storage_validation.meta.json", "storage_projection.csv",
        "storage_projection.root", "storage_projection.meta.json",
        "first_bar_encounters_production.root",
    ]
    artifacts = {str(args.results_dir / name): sha(args.results_dir / name)
                 for name in artifact_names}
    golden["title"] = "EXEC38 golden reference with EXEC40 and corrected EXEC41 single-cell observations"
    golden["ready_for_acceptance"] = False
    golden["overall_status"] = "INCOMPLETE_SINGLE_CELL_V2_CORRECTED_V5_MATCHED_V1_EXPLAINED"
    golden["exec41_single_cell"] = {
        "scope": "EJ-204 x=0 EndTop70, N=2000 existing EXEC40 ROOT only; no simulation and no grid extension",
        "source_parent_commit": args.source_parent_commit,
        "preregistration_path": str(args.preregistration.resolve()),
        "preregistration_sha256": sha(args.preregistration),
        "superseded_interpretation_only": "exec40_single_cell.V2 p(mu)=2mu and 0.36-0.40 expectation",
        "report_path": str(args.report),
        "ready_for_acceptance": False,
        "V1": {
            "status": "EXPLAINED_BY_OPTICAL_PHOTON_ENERGY_RESCINTILLATION",
            **diagnosis["V1"],
            "EXEC40_reported_nonoptical_ratio_se": 0.0001603379137171514,
        },
        "V2": {
            "status": "CORRECTED_H1_PASS_H2_DESCRIPTIVE",
            "H1": v2["H1"],
            "faces": v2["faces"],
            "mechanism": v2["mechanism"],
            "H2": v2["H2"],
            "old_EXEC40_hypothesis_retained": v2["old_EXEC40_hypothesis"],
            "face_resolution_limitation": v2["face_resolution_limitation"],
        },
        "V5": {
            "status": "MATCHED_BOUNDARY_POPULATION_PASS_WITH_SEPARATE_SD_ACCOUNTING_GAP",
            "matched_detection_over_incident": all_arm["matched_over_incident"],
            "incident_spectrum_expected_PDE": all_arm["incident_spectrum_pde"],
            "difference": all_arm["matched_over_incident"]["value"] - all_arm["incident_spectrum_pde"]["value"],
            "legacy_SD_gap": diagnosis["detection_gap"],
            "cause": diagnosis["cause"],
            "spectral_filtering": diagnosis["spectral_filtering"],
        },
        "storage": storage,
        "analysis_commands": {
            "V2": "analysis/validation/analyze_exec41_v2.py --root /home/rrios/exec40_20260913/cell_validated/photon_hits_run000.root --out /home/rrios/exec41_20260913 --preregistration analysis/validation/EXEC41_PREREGISTRATION.md",
            "diagnosis": "analysis/validation/diagnose_exec41.py --root /home/rrios/exec40_20260913/cell_validated/photon_hits_run000.root --out /home/rrios/exec41_20260913 --repo /home/rrios/ej200_exec40_20260913 --exec27-summary /home/rrios/ej200_exec26_20260909/build_exec27_20260910/audit/analysis_summary.json --exec27-top-states /home/rrios/ej200_exec26_20260909/build_exec27_20260910/run500/tables/top_boundary_states.csv --preregistration analysis/validation/EXEC41_PREREGISTRATION.md",
            "storage": "analysis/validation/validate_exec41_storage.py --root /home/rrios/exec40_20260913/cell_validated/photon_hits_run000.root --filtered-root /home/rrios/exec41_20260913/first_bar_encounters_production.root --out /home/rrios/exec41_20260913 --v2-face-csv /home/rrios/exec41_20260913/v2_face_decomposition.csv --v2-cosine-csv /home/rrios/exec41_20260913/v2_uniform_cosine.csv --preregistration analysis/validation/EXEC41_PREREGISTRATION.md",
        },
        "artifact_hashes": artifacts,
    }
    golden["updated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    args.golden.write_text(json.dumps(golden, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
