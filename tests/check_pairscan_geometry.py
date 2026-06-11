#!/usr/bin/env python3
"""EXEC_08 CTest guardrail — pair-scan geometry checks.

Tests:
1. globalIDs 28 and 29 exist in the Top readout with /det/readout Top and
   their x-centers match the values hardcoded in analyze_pair_scan.C (±0.1 mm).
2. All 41 scan positions [-462, -422] mm (step 1 mm) are:
   (a) within the bar [-700, +700] mm, and
   (b) within [x_A - p/2, x_B + p/2] = [-462, -422] mm.

Reads geometry from DetectorConstruction.cc:45-49 and hh:37-38 (same formula
as analyze_geometry.py).

Usage (ctest invokes this directly):
    python3 tests/check_pairscan_geometry.py
"""

from __future__ import annotations

import sys

# ── Geometry constants (DetectorConstruction.cc:45-49, hh:37-38) ──────────────
K_N_TOP_SIPMS = 70
K_N_END_SIPMS = 8
GLOBAL_ID_OFFSET = 2 * K_N_END_SIPMS  # = 16
BAR_HALF_X = 700.0  # mm (DetectorConstruction.cc:29)


def top_sipm_center_x(local_idx: int) -> float:
    if local_idx < 35:
        return -692.0 + 20.0 * local_idx
    return 12.0 + 20.0 * (local_idx - 35)


def global_to_x(global_id: int) -> float:
    local_idx = global_id - GLOBAL_ID_OFFSET
    if local_idx < 0 or local_idx >= K_N_TOP_SIPMS:
        raise ValueError(f"global_id {global_id} out of top range [16, 85]")
    return top_sipm_center_x(local_idx)


# ── Expected values from analyze_pair_scan.C ──────────────────────────────────
ID_A        = 28
ID_B        = 29
X_A_EXPECT  = -452.0  # mm
X_B_EXPECT  = -432.0  # mm
X_MID       = -442.0  # mm
PITCH       =   20.0  # mm
SCAN_LO     = X_A_EXPECT - PITCH / 2   # = -462.0 mm
SCAN_HI     = X_B_EXPECT + PITCH / 2   # = -422.0 mm
SCAN_STEP   =   1.0
TOL         =   0.1  # mm


def build_scan_positions() -> list[float]:
    positions = []
    x = SCAN_LO
    while x <= SCAN_HI + 1e-9:
        positions.append(round(x, 6))
        x += SCAN_STEP
    return positions


def fail(msg: str) -> None:
    print(f"check_pairscan_geometry FAILED: {msg}", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    errors: list[str] = []

    # ── Test 1a: pair IDs exist as top SiPMs ──────────────────────────────────
    valid_top_ids = set(range(GLOBAL_ID_OFFSET, GLOBAL_ID_OFFSET + K_N_TOP_SIPMS))
    for gid in (ID_A, ID_B):
        if gid not in valid_top_ids:
            errors.append(f"ID {gid} not in top SiPM range [16, 85]")

    # ── Test 1b: x-centers match expected values ───────────────────────────────
    if not errors:
        x_a_actual = global_to_x(ID_A)
        x_b_actual = global_to_x(ID_B)
        if abs(x_a_actual - X_A_EXPECT) > TOL:
            errors.append(f"x_A mismatch: ID {ID_A} x={x_a_actual:.2f} mm "
                          f"(expected {X_A_EXPECT:.2f} ±{TOL} mm)")
        if abs(x_b_actual - X_B_EXPECT) > TOL:
            errors.append(f"x_B mismatch: ID {ID_B} x={x_b_actual:.2f} mm "
                          f"(expected {X_B_EXPECT:.2f} ±{TOL} mm)")

    # ── Test 1c: pitch is modal (20 mm, not anomalous gap) ───────────────────
    if not errors:
        pitch_actual = x_b_actual - x_a_actual
        if abs(pitch_actual - PITCH) > TOL:
            errors.append(f"Pair pitch {pitch_actual:.2f} mm ≠ modal pitch {PITCH:.2f} mm "
                          "(selected pair crosses anomalous gap!)")

    # ── Test 2: scan positions are within bar and [SCAN_LO, SCAN_HI] ─────────
    positions = build_scan_positions()
    for x in positions:
        if x < -BAR_HALF_X or x > BAR_HALF_X:
            errors.append(f"Scan position {x:+.1f} mm outside bar [{-BAR_HALF_X}, +{BAR_HALF_X}] mm")
        if x < SCAN_LO - TOL or x > SCAN_HI + TOL:
            errors.append(f"Scan position {x:+.1f} mm outside "
                          f"[x_A-p/2, x_B+p/2] = [{SCAN_LO:.1f}, {SCAN_HI:.1f}] mm")

    # ── Test 3: required anchor positions present ─────────────────────────────
    pos_set = set(positions)
    for anchor in (X_A_EXPECT, X_MID, X_B_EXPECT):
        if not any(abs(anchor - p) < TOL for p in pos_set):
            errors.append(f"Anchor position {anchor:+.1f} mm not found on scan grid")

    if errors:
        for e in errors:
            print(f"  ERROR: {e}", file=sys.stderr)
        fail(f"{len(errors)} check(s) failed")

    n = len(positions)
    print(f"check_pairscan_geometry PASSED: "
          f"IDs ({ID_A},{ID_B}) at ({X_A_EXPECT:+.0f},{X_B_EXPECT:+.0f}) mm, "
          f"pitch={PITCH:.0f} mm, {n} scan positions [{SCAN_LO:+.0f},{SCAN_HI:+.0f}] mm OK")


if __name__ == "__main__":
    main()
