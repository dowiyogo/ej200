#!/usr/bin/env python3
"""EXEC_08 Step 2 — extract top SiPM layout and select scan pair.

Reads the geometry from DetectorConstruction.cc:45-49 (hardcoded formula,
verified against kNTopSiPMs=70 in DetectorConstruction.hh:38).

Outputs pair_scan_config.json consumed by gen_pair_scan_macros.py.

Usage:
    python3 scripts/analyze_geometry.py [--output pair_scan_config.json]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

# ── Geometry constants (from DetectorConstruction.cc:45-49 + hh:37-38) ────────
K_N_TOP_SIPMS = 70
K_N_END_SIPMS = 8
GLOBAL_ID_TOP_OFFSET = 2 * K_N_END_SIPMS  # = 16


def top_sipm_center_x(local_idx: int) -> float:
    """Return x_center [mm] for top SiPM local index.

    DetectorConstruction.cc:45-49:
        if (idx < 35) return (-692.0 + 20.0 * idx) * mm;
        return (12.0 + 20.0 * (idx - 35)) * mm;
    """
    if local_idx < 0 or local_idx >= K_N_TOP_SIPMS:
        raise ValueError(f"local_idx {local_idx} out of [0, {K_N_TOP_SIPMS})")
    if local_idx < 35:
        return -692.0 + 20.0 * local_idx
    return 12.0 + 20.0 * (local_idx - 35)


def build_layout() -> list[dict]:
    """Return list of dicts: {global_id, local_idx, x_mm, pitch_to_next_mm}."""
    rows = []
    for i in range(K_N_TOP_SIPMS):
        gid = GLOBAL_ID_TOP_OFFSET + i
        x = top_sipm_center_x(i)
        rows.append({"global_id": gid, "local_idx": i, "x_mm": x, "pitch_to_next_mm": None})
    for k in range(len(rows) - 1):
        rows[k]["pitch_to_next_mm"] = rows[k + 1]["x_mm"] - rows[k]["x_mm"]
    return rows


def select_pair(layout: list[dict]) -> dict:
    """Apply pair-selection criteria (EXEC_08 Step 2).

    Criteria (in order):
    1. Pitch == p (modal; excludes central gap pair).
    2. Neither element is the first or last SiPM of its group.
    3. Located toward the middle third of its group.

    Returns {id_a, id_b, x_a_mm, x_b_mm, x_mid_mm, pitch_mm}.
    """
    pitches = [r["pitch_to_next_mm"] for r in layout if r["pitch_to_next_mm"] is not None]
    modal_pitch = max(set(pitches), key=pitches.count)

    # Left group: local indices 0-34 → global IDs 16-50
    left = [r for r in layout if r["local_idx"] <= 34]
    # Exclude first and last of the group
    valid_left = left[1:-1]  # exclude local_idx 0 (ID 16) and 34 (ID 50)

    # Middle third of left group (35 members; valid: 33; middle third ~ indices 11-21 of valid)
    n = len(valid_left)
    lo, hi = n // 3, 2 * n // 3
    mid_candidates = valid_left[lo:hi]

    # Select pair with pitch == modal_pitch from middle candidates
    for r in mid_candidates:
        if abs((r["pitch_to_next_mm"] or 0.0) - modal_pitch) < 0.1:
            next_gid = r["global_id"] + 1
            next_r = next((q for q in layout if q["global_id"] == next_gid), None)
            if next_r and next_r["local_idx"] <= 34:
                return {
                    "id_a": r["global_id"],
                    "id_b": next_r["global_id"],
                    "x_a_mm": r["x_mm"],
                    "x_b_mm": next_r["x_mm"],
                    "x_mid_mm": 0.5 * (r["x_mm"] + next_r["x_mm"]),
                    "pitch_mm": modal_pitch,
                }
    raise RuntimeError("No valid pair found — re-check selection criteria.")


def build_scan_positions(pair: dict) -> list[float]:
    """Step [x_A - p/2, x_B + p/2], step = max(p/20, 1 mm), ensure x_A/x_mid/x_B on grid."""
    p = pair["pitch_mm"]
    x_lo = pair["x_a_mm"] - p / 2.0
    x_hi = pair["x_b_mm"] + p / 2.0
    step = max(p / 20.0, 1.0)

    positions: list[float] = []
    x = x_lo
    while x <= x_hi + 1e-9:
        positions.append(round(x, 6))
        x += step

    required = {pair["x_a_mm"], pair["x_mid_mm"], pair["x_b_mm"]}
    missing = required - set(positions)
    if missing:
        raise RuntimeError(f"Required positions not on grid: {missing}")
    return positions


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="pair_scan_config.json",
                        type=pathlib.Path)
    parser.add_argument("--print-table", action="store_true",
                        help="Print full geometry table to stdout.")
    args = parser.parse_args(argv)

    layout = build_layout()

    if args.print_table:
        print(f"{'globalID':>8}  {'localIdx':>8}  {'x_mm':>8}  {'pitch_to_next_mm':>16}")
        for r in layout:
            pt = r["pitch_to_next_mm"]
            marker = "  <-- ANOMALOUS" if pt is not None and abs(pt - 20.0) > 0.5 else ""
            pt_str = f"{pt:>16.1f}" if pt is not None else f"{'—':>16}"
            print(f"{r['global_id']:>8}  {r['local_idx']:>8}  {r['x_mm']:>8.1f}  {pt_str}{marker}")

    pair = select_pair(layout)
    positions = build_scan_positions(pair)

    pitches = [r["pitch_to_next_mm"] for r in layout if r["pitch_to_next_mm"] is not None]
    modal_pitch = max(set(pitches), key=pitches.count)
    anomalous = [(layout[k]["global_id"], layout[k]["x_mm"],
                  layout[k + 1]["global_id"], layout[k + 1]["x_mm"],
                  layout[k]["pitch_to_next_mm"])
                 for k in range(len(layout) - 1)
                 if layout[k]["pitch_to_next_mm"] is not None
                 and abs(layout[k]["pitch_to_next_mm"] - modal_pitch) > 0.5]

    config = {
        "date": "2026-06-11",
        "k_n_top_sipms": K_N_TOP_SIPMS,
        "k_n_end_sipms": K_N_END_SIPMS,
        "modal_pitch_mm": modal_pitch,
        "anomalous_gaps": [
            {"id_lo": a[0], "x_lo_mm": a[1], "id_hi": a[2], "x_hi_mm": a[3], "pitch_mm": a[4]}
            for a in anomalous
        ],
        "pair": pair,
        "n_pe_threshold": 4,
        "events_per_position": 3000,
        "scan_positions_mm": positions,
        "n_positions": len(positions),
        "jitter_label": "INTRINSIC",
        "jitter_note": "Only /sipm/jitterSigma (electronics) exists; no SPTR hook. sigma_Dt does NOT include SPTR or electronic jitter.",
    }

    out = pathlib.Path(args.output)
    out.write_text(json.dumps(config, indent=2))
    print(f"Wrote {out}")
    print(f"Pair: IDs ({pair['id_a']}, {pair['id_b']}), "
          f"x_A={pair['x_a_mm']:+.1f} mm, x_B={pair['x_b_mm']:+.1f} mm, "
          f"x_mid={pair['x_mid_mm']:+.1f} mm, pitch={pair['pitch_mm']:.0f} mm")
    print(f"Scan: {len(positions)} positions from {positions[0]:+.1f} to {positions[-1]:+.1f} mm, "
          f"step={positions[1]-positions[0]:.1f} mm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
