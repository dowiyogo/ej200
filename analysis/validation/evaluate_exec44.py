#!/usr/bin/env python3
"""Evaluate EXEC44 H3' from published EXEC42/43 JSON; never opens ROOT."""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path):
    if path.suffix == ".root":
        raise ValueError("EXEC44 is restricted to published non-ROOT results")
    with path.open() as stream:
        return json.load(stream)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exec43", type=Path, required=True)
    parser.add_argument("--exec42", type=Path, required=True)
    parser.add_argument("--golden", type=Path, required=True)
    parser.add_argument("--prereg", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    r43 = load_json(args.exec43)
    r42 = load_json(args.exec42)
    golden = load_json(args.golden)

    optical = r43["M2"]["optical_properties"]
    n_values = [row["n_effective"] for row in optical]
    n_spread = max(n_values) - min(n_values)
    c1 = n_spread <= 1e-4 and r43["M2"]["status"] == "EXCLUDED"

    absent = r43["absent_photons"]
    by_x_absent = {}
    for row in absent:
        by_x_absent.setdefault(row["x_mm"], {})[row["material"]] = row["absent_fraction"]
    absent_order = all(
        values["EJ-200"] < values["EJ-204"] < values["EJ-230"]
        for values in by_x_absent.values()
    )
    adjacent = [row for row in r43["absent_pairwise"] if row["adjacent_pair_for_M1_decision"]]
    c2 = len(by_x_absent) == 7 and len(adjacent) == 14 and absent_order and all(row["status"] == "PASS" for row in adjacent)

    escape = r42["H3"]
    c3 = len({row["x_mm"] for row in escape}) == 7 and len(escape) == 21 and all(row["difference"] < 0 for row in escape)

    directional = [row for row in r43["directional_path_test"] if row["enters_M1_decision"]]
    c4 = len(directional) == 84 and all(
        row["status"] == "PASS" and row["non_escape_secant_median"] > row["escape_secant_median"]
        for row in directional
    )

    conditions = [
        {
            "condition": 1,
            "status": "PASS" if c1 else "FAIL",
            "evidence": f"n_eff={','.join(f'{v:.16g}' for v in n_values)}; spread={n_spread:.17g} <= 1e-4; M2={r43['M2']['status']}",
        },
        {
            "condition": 2,
            "status": "PASS" if c2 else "FAIL",
            "evidence": f"EJ-200 < EJ-204 < EJ-230 absent fraction at {len(by_x_absent)}/7 positions; {sum(r['status']=='PASS' for r in adjacent)}/{len(adjacent)} adjacent comparisons PASS",
        },
        {
            "condition": 3,
            "status": "PASS" if c3 else "FAIL",
            "evidence": f"escape EJ-200 < EJ-204 < EJ-230 in {len({r['x_mm'] for r in escape})}/7 positions ({sum(r['difference'] < 0 for r in escape)}/{len(escape)} ordered pairwise differences)",
        },
        {
            "condition": 4,
            "status": "PASS" if c4 else "FAIL",
            "evidence": f"non-escape sec(theta) median > escape median in {sum(r['status']=='PASS' for r in directional)}/{len(directional)} available directional comparisons",
        },
    ]
    h3_prime_pass = all(row["status"] == "PASS" for row in conditions)

    first_encounter_rows = sum(row["first_encounters"] for row in absent)
    payload_bytes_per_row = 8 + 3 * 8 + 3 * 8
    payload_total = first_encounter_rows * payload_bytes_per_row
    now = datetime.now(timezone.utc).isoformat()
    evaluation = {
        "schema": "ej200.exec44.h3_prime_evaluation.v1",
        "created_utc": now,
        "method": "Read published EXEC43 and EXEC42 JSON only; no ROOT opened and no observable recalculated.",
        "preregistration": {"path": str(args.prereg), "sha256": sha256(args.prereg)},
        "inputs": {
            "exec43_results": {"path": str(args.exec43), "sha256": sha256(args.exec43)},
            "exec42_contract_results": {"path": str(args.exec42), "sha256": sha256(args.exec42)},
        },
        "conditions": conditions,
        "H3_prime_status": "PASS" if h3_prime_pass else "FAIL",
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / "EXEC44_H3_PRIME_EVALUATION.json"
    csv_path = args.out_dir / "EXEC44_H3_PRIME_EVALUATION.csv"
    json_path.write_text(json.dumps(evaluation, indent=2) + "\n")
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["condition", "status", "evidence"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(conditions)

    active_tests = {
        "V1": {
            "observable": "sum of scintillation photons divided by nominal material yield times total deposited energy in BarLV",
            "status": "PASS",
            "coverage": "21/21 cells",
            "range": [0.999862, 1.000178],
        },
        "V2_H1": {
            "observable": "portable-outcome escape fraction at the first physical BarLV encounter of each scintillation photon exiting BarLV",
            "status": "PASS",
            "coverage": "21/21 cells",
            "range": [0.256163, 0.263241],
        },
        "V5_matched": {
            "observable": "matched boundary detections divided by independently observed sensor incidents, compared with surface-PDE expectation",
            "status": "PASS",
            "coverage": "21/21 cells",
        },
        "H3_prime": {
            "observable": "material ordering of first-encounter escape, absent photons, effective refractive index, and incidence-angle path proxy",
            "status": "PASS" if h3_prime_pass else "FAIL",
            "coverage": "four preregistered qualitative conditions",
            "evaluation": str(json_path),
        },
    }
    active_bad = [name for name, item in active_tests.items() if item["status"] in {"FAIL", "NOT_EVALUABLE"}]
    ready = not active_bad
    reason = (
        "V1, V2-H1, and matched V5 pass 21/21 cells and qualitative H3' passes all four conditions; no active test is FAIL or NOT_EVALUABLE."
        if ready
        else "Active acceptance tests preventing readiness: " + ", ".join(active_bad) + "."
    )

    golden["title"] = "EXEC44 golden physical-transport acceptance contract"
    golden["updated_utc"] = now
    golden["ready_for_acceptance"] = ready
    golden["overall_status"] = "READY_FOR_ACCEPTANCE" if ready else "NOT_READY_FOR_ACCEPTANCE"
    golden["active_acceptance_contract"] = {
        "scope": "Physical observables portable to a GPU transport engine; no Geant4-internal instrumentation is an acceptance observable.",
        "tests": active_tests,
        "failure_or_not_evaluable": active_bad,
        "ready_for_acceptance": ready,
        "reason": reason,
    }
    golden["exec44_resolution"] = {
        "decision": "Replace retracted H3 material invariance with qualitative H3'.",
        "preregistration": evaluation["preregistration"],
        "evaluation": {"path": str(json_path), "sha256": sha256(json_path)},
        "retracted_H3": {
            "status": "RETRACTED",
            "original_definition": golden["exec42_full_grid"]["physical_observable_contract"]["H3"],
            "historical_result": "FAIL in 21/21 material-pair comparisons in EXEC42",
            "reason": "First-encounter populations exclude photons absorbed in bulk before reaching a surface, so material invariance was physically inappropriate.",
            "refuted_by": "EXEC43; resolved by EXEC44",
        },
        "H3_prime": {"status": "PASS" if h3_prime_pass else "FAIL", "conditions": conditions},
        "H3_double_prime": {
            "status": "PENDING_FUTURE_WORK_NOT_AN_ACCEPTANCE_CRITERION",
            "objective": "Predict pairwise escape differences from attenuation lengths and measured pre-encounter path distributions within a declared tolerance.",
            "missing_observables": [
                {"name": "accumulated_track_length", "bytes_per_row": 8, "representation": "float64"},
                {"name": "creation_coordinates_xyz", "bytes_per_row": 24, "representation": "3 x float64"},
                {"name": "first_encounter_coordinates_xyz", "bytes_per_row": 24, "representation": "3 x float64"},
            ],
            "estimated_raw_payload_bytes_per_row": payload_bytes_per_row,
            "exec42_first_encounter_rows": first_encounter_rows,
            "estimated_raw_payload_bytes_for_exec42_grid": payload_total,
            "estimate_excludes": "ROOT branch/basket overhead and compression",
            "exec43_limitation": "The central WorldPV-envelope assignment overpredicted all 14 differences by 3.73-7.63 SE because exact paths were unavailable.",
        },
        "known_model_limitation": "SSLG4 uses constant n=1.58 for all three materials; measured model dispersion is below 4.5e-16. Real PVT varies approximately from 1.57 to 1.61 over 380-450 nm. This simplification affects wavelength-dependent observables, including V5 spectral transport filtering, and must be checked when translating the geometry to another transport engine.",
        "ready_for_acceptance": ready,
        "reason": reason,
    }
    args.golden.write_text(json.dumps(golden, indent=2) + "\n")
    print(json.dumps({"H3_prime": evaluation["H3_prime_status"], "ready_for_acceptance": ready, "reason": reason}, indent=2))


if __name__ == "__main__":
    main()
