"""EXEC_13 (EJ-230) geometry and analysis constants.

Geometry and channel layout are BYTE-IDENTICAL to exec07/common.py.
Only the scintillator constants differ (EJ-230 vs EJ-204).
"""

from __future__ import annotations

import os
import pathlib

# ── Geometry (identical to exec07) ──────────────────────────────────────────
N_EVENTS = 2000
N_CHANNELS = 86
EXPECTED_POSITIONS_MM = (
    -690, -670, *range(-650, 651, 50), 670, 690,
)
END_LEFT_IDS = tuple(range(0, 8))
END_RIGHT_IDS = tuple(range(8, 16))
TOP_IDS = tuple(range(16, 86))
TOP_LEFT_IDS = tuple(range(16, 51))
TOP_RIGHT_IDS = tuple(range(51, 86))
END_CLUSTERS = {
    "end_left_A_SUM4": (0, 1, 2, 3),
    "end_left_B_SUM4": (4, 5, 6, 7),
    "end_right_A_SUM4": (8, 9, 10, 11),
    "end_right_B_SUM4": (12, 13, 14, 15),
}
TOP_POSITIONS_MM = tuple(
    [-692.0 + 20.0 * i for i in range(35)]
    + [12.0 + 20.0 * i for i in range(35)]
)
BAR_HALF_LENGTH_MM = 700.0

# ── EJ-230 scintillator constants (Eljen datasheet Rev. Aug 2023) ────────────
TAU_R_NS = 0.5   # rise time
TAU_D_NS = 1.5   # decay time
SCINT_YIELD_PER_MEV = 9700
ABSLENGTH_CM = 120.0
EMISSION_PEAK_NM = 391.0
SIGMA_NUMERATOR_NS = (TAU_R_NS * TAU_D_NS) ** 0.5

SPR_RISE_NS = 0.5
SPR_FALL_NS = 5.0
LEADING_EDGE_THRESHOLD_PE = 4.0

KEY_POSITIONS = (-690, -650, -400, 0, 400, 650, 690)
N_THRESHOLDS = (4, 20)


def nearest_top_ids(x_mm: float, count: int = 3) -> tuple[int, ...]:
    ranked = sorted(
        range(len(TOP_POSITIONS_MM)),
        key=lambda index: (abs(TOP_POSITIONS_MM[index] - x_mm), index),
    )
    return tuple(16 + index for index in ranked[:count])


def nearest_end_name(x_mm: float) -> str:
    return "end_left_all" if x_mm <= 0 else "end_right_all"


def results_dir() -> pathlib.Path:
    """Return the canonical results directory (respects RESULTS_DIR env var)."""
    base = os.environ.get("RESULTS_DIR",
                          str(pathlib.Path.home() / "ej230" / "results_ej230"))
    return pathlib.Path(base)


def expected_file(data_dir: pathlib.Path, x_mm: int) -> pathlib.Path:
    return data_dir / f"photon_hits_x{x_mm}mm.root"
