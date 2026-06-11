#!/usr/bin/env python3
"""EXEC_08 CTest guardrail — validate generated pair-scan macros.

Checks each .mac in the pairscan macro directory:
1. Contains /det/readout Top.
2. /muon/gunX value is in the expected range [-462, -422] mm.
3. /run/beamOn 3000.
4. Expected count of macros (41 positions).

Usage:
    python3 tests/check_pairscan_macros.py macros/pairscan
"""

from __future__ import annotations

import pathlib
import re
import sys

SCAN_LO   = -462.0
SCAN_HI   = -422.0
N_EVENTS  = 3000
N_MACROS_EXPECTED = 41
TOL       = 0.1


def fail(msg: str) -> None:
    print(f"check_pairscan_macros FAILED: {msg}", file=sys.stderr)
    sys.exit(1)


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else argv
    if not args:
        fail("Usage: check_pairscan_macros.py <macro_dir>")

    macro_dir = pathlib.Path(args[0])
    if not macro_dir.exists():
        fail(f"Macro directory not found: {macro_dir}")

    macros = sorted(macro_dir.glob("pairscan_x*.mac"))
    errors: list[str] = []

    if len(macros) != N_MACROS_EXPECTED:
        errors.append(f"Expected {N_MACROS_EXPECTED} macros, found {len(macros)}")

    for mac in macros:
        text = mac.read_text()

        # Check /det/readout Top
        if "/det/readout Top" not in text:
            errors.append(f"{mac.name}: missing '/det/readout Top'")

        # Check /run/beamOn N_EVENTS
        if f"/run/beamOn {N_EVENTS}" not in text:
            errors.append(f"{mac.name}: missing '/run/beamOn {N_EVENTS}'")

        # Check /muon/gunX value in range
        m = re.search(r"/muon/gunX\s+([+-]?\d+\.?\d*)\s+mm", text)
        if not m:
            errors.append(f"{mac.name}: missing /muon/gunX command")
        else:
            x = float(m.group(1))
            if not (SCAN_LO - TOL <= x <= SCAN_HI + TOL):
                errors.append(f"{mac.name}: gunX={x:.1f} mm outside "
                               f"[{SCAN_LO:.1f}, {SCAN_HI:.1f}] mm")

        # Check /gun/particle mu-
        if "/gun/particle mu-" not in text:
            errors.append(f"{mac.name}: missing '/gun/particle mu-'")

    if errors:
        for e in errors:
            print(f"  ERROR: {e}", file=sys.stderr)
        fail(f"{len(errors)} check(s) failed in {len(macros)} macros")

    print(f"check_pairscan_macros PASSED: {len(macros)} macros validated in {macro_dir}")


if __name__ == "__main__":
    main()
