#!/usr/bin/env python3
"""
validate_group_velocity_ej230.py
Port of mirror's validate_group_velocity.py for EJ-230 End-only+Mylar.

Validates that no GROUPVEL property is set in Materials.cc (confirming
Geant4 uses v_g = c/n from RINDEX). Computes and documents the
theoretical group velocity for EJ-230.

Outputs: analysis_ej230/csv/group_velocity_validation.csv
"""

import csv
import math
import pathlib
import re
import sys

HERE     = pathlib.Path(__file__).parent
REPO     = HERE.parent.parent
MAT_CC   = REPO / "src" / "Materials.cc"
OUT_CSV  = HERE.parent / "csv" / "group_velocity_validation.csv"

# EJ-230 constants (Materials.cc)
N_REFR      = 1.58          # flat RINDEX (Materials.cc:41)
C_MM_NS     = 299.792458    # speed of light mm/ns
V_G_THEORY  = C_MM_NS / N_REFR   # = 189.74 mm/ns

def check_groupvel_absent(mat_cc: pathlib.Path) -> tuple[bool, list[str]]:
    """Return (absent, matching_lines). absent=True means no GROUPVEL property."""
    if not mat_cc.exists():
        return (False, [f"Materials.cc not found at {mat_cc}"])
    text = mat_cc.read_text()
    pattern = re.compile(r"GROUPVEL", re.IGNORECASE)
    matches = []
    for i, line in enumerate(text.splitlines(), 1):
        if pattern.search(line):
            matches.append(f"  line {i}: {line.rstrip()}")
    return (len(matches) == 0, matches)

def main():
    absent, matches = check_groupvel_absent(MAT_CC)

    print("=== Group Velocity Validation: EJ-230 ===")
    print(f"Materials.cc: {MAT_CC}")
    print(f"GROUPVEL absent: {absent}")
    if not absent:
        print("  GROUPVEL found at:")
        for m in matches:
            print(m)
    else:
        print("  No GROUPVEL property set — Geant4 computes v_g = c/n.")

    print()
    print(f"n (RINDEX, flat)       = {N_REFR}")
    print(f"c                      = {C_MM_NS:.6f} mm/ns")
    print(f"v_g = c/n (theory)     = {V_G_THEORY:.4f} mm/ns")
    print()

    # Write CSV
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"quantity": "RINDEX_flat",
         "value": str(N_REFR),
         "unit": "dimensionless",
         "source": "Materials.cc:41",
         "note": "Flat refractive index used for all wavelengths"},
        {"quantity": "c_vacuum",
         "value": f"{C_MM_NS:.6f}",
         "unit": "mm/ns",
         "source": "PDG 2022",
         "note": "Speed of light in vacuum"},
        {"quantity": "v_g_theory",
         "value": f"{V_G_THEORY:.4f}",
         "unit": "mm/ns",
         "source": "c/n (derived)",
         "note": "Theoretical group velocity; no GROUPVEL property set in Geant4"},
        {"quantity": "GROUPVEL_set",
         "value": str(not absent),
         "unit": "bool",
         "source": str(MAT_CC),
         "note": "Whether GROUPVEL property was found in Materials.cc"},
        {"quantity": "validation_result",
         "value": "PASS" if absent else "FAIL",
         "unit": "status",
         "source": "validate_group_velocity_ej230.py",
         "note": "PASS = Geant4 uses v_g=c/n as expected; FAIL = explicit GROUPVEL found"},
    ]
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["quantity", "value", "unit", "source", "note"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written: {OUT_CSV}")
    if not absent:
        print("WARNING: GROUPVEL found — may override v_g = c/n!", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
