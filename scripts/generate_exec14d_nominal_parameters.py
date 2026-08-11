#!/usr/bin/env python3
"""Export documented EJ-230 nominal timing parameters to a provenance CSV."""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys


REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from analysis.exec13.common13 import SIGMA_NUMERATOR_NS, TAU_D_NS, TAU_R_NS  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=REPO / "results_ej230_analysis/csv/exec14d_nominal_parameters.csv",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("tau_r_ns", "tau_d_ns", "sigma_numerator_ns", "reference"),
        )
        writer.writeheader()
        writer.writerow({
            "tau_r_ns": TAU_R_NS,
            "tau_d_ns": TAU_D_NS,
            "sigma_numerator_ns": SIGMA_NUMERATOR_NS,
            "reference": "Knoll, ch. 8",
        })
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
