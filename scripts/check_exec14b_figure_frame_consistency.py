#!/usr/bin/env python3
"""Check that position-specific EJ-230 figures are attached to the correct frames."""

from __future__ import annotations

import ast
import csv
import hashlib
import re
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
REPORT = REPO / "results_ej230_analysis/report/exec13_ej230_report_full.tex"
FIGS = REPO / "results_ej230_analysis/figs"
POSITIONS_CSV = REPO / "results_ej230_analysis/per_position_exec07.csv"
TN_CSV = REPO / "results_ej230_analysis/csv/exec13_tN_summary.csv"
KEY_POSITIONS = (-690, -650, -400, 0, 400, 650, 690)
DISPERSION_FRAMES = dict(zip((8, 12, 16, 20, 24, 28, 32), KEY_POSITIONS))
INCLUDE_RE = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*%?\s*\{([^{}]+)\}", re.DOTALL)


def frames(text: str) -> list[str]:
    output = []
    for match in re.finditer(r"\\begin\{frame\}", text):
        end = text.find(r"\end{frame}", match.end())
        if end < 0:
            raise RuntimeError("unterminated frame")
        output.append(text[match.end():end])
    return output


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    text = REPORT.read_text(encoding="utf-8")
    bodies = frames(text)
    errors: list[str] = []
    dispersion_paths: list[str] = []

    if len(bodies) < 119:
        errors.append(f"expected at least 119 frames, found {len(bodies)}")

    for frame_number, x_mm in DISPERSION_FRAMES.items():
        body = bodies[frame_number - 1]
        actual = INCLUDE_RE.findall(body)
        expected = [
            f"figs/exec13_tn_{x_mm}mm_top.png",
            f"figs/exec13_tn_{x_mm}mm_endL.png",
            f"figs/exec13_tn_{x_mm}mm_endR.png",
        ]
        if actual != expected:
            errors.append(
                f"frame {frame_number} x={x_mm}: ordered dispersion figures {actual}, "
                f"expected {expected}"
            )
        dispersion_paths.extend(actual)

    if len(dispersion_paths) != 21 or len(set(dispersion_paths)) != 21:
        errors.append("the seven dispersion frames do not use 21 distinct paths")
    existing_dispersion = [FIGS / Path(path).name for path in dispersion_paths]
    if all(path.is_file() for path in existing_dispersion):
        hashes = [digest(path) for path in existing_dispersion]
        if len(set(hashes)) != 21:
            errors.append("the 21 dispersion paths do not contain 21 distinct image hashes")

    with POSITIONS_CSV.open(newline="", encoding="utf-8") as stream:
        nearest = {int(row["x_beam_mm"]): int(row["nearest_top_id"]) for row in csv.DictReader(stream)}
    with TN_CSV.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            x_mm = int(row["x_mm"])
            if x_mm not in KEY_POSITIONS or row["group"] != "top_nearest" or int(row["N"]) != 4:
                continue
            channels = ast.literal_eval(row["channels"])
            if tuple(channels) != (nearest[x_mm],):
                errors.append(
                    f"x={x_mm}: top_nearest channels={channels}, expected {(nearest[x_mm],)}"
                )

    position_re = re.compile(r"\\begin\{frame\}\{Position x=(-?\d+) mm\}")
    for match in position_re.finditer(text):
        x_mm = int(match.group(1))
        end = text.find(r"\end{frame}", match.end())
        body = text[match.end():end]
        figures = INCLUDE_RE.findall(body)
        expected_prefix = f"figs/muon_{x_mm}mm_"
        wrong = [path for path in figures if not path.startswith(expected_prefix)]
        if wrong:
            errors.append(f"position frame x={x_mm} contains mismatched figures: {wrong}")
        table = (
            f"key_position_{x_mm}.tex" if x_mm in KEY_POSITIONS else f"position_{x_mm}.tex"
        )
        if rf"\input{{../tables/{table}}}" not in body:
            errors.append(f"position frame x={x_mm} does not import {table}")

    arrival_re = re.compile(
        r"\\begin\{frame\}\{Run \d+ \$\|\$ x = (?P<sign>\+?)(?P<x>-?\d+) mm "
        r"\$\|\$ Photon arrival .*?\}"
    )
    for match in arrival_re.finditer(text):
        x_mm = int(match.group("x"))
        end = text.find(r"\end{frame}", match.end())
        figures = INCLUDE_RE.findall(text[match.end():end])
        expected = f"figs/exec11_arrival_{x_mm}mm.png"
        if figures != [expected]:
            errors.append(f"arrival frame x={x_mm} contains {figures}, expected [{expected}]")

    tn_re = re.compile(
        r"\\begin\{frame\}\{Run \d+ \$\|\$ x = (?P<sign>\+?)(?P<x>-?\d+) mm "
        r"\$\|\$ Time-to-threshold \$t_N\$\}"
    )
    for match in tn_re.finditer(text):
        x_mm = int(match.group("x"))
        end = text.find(r"\end{frame}", match.end())
        figures = INCLUDE_RE.findall(text[match.end():end])
        expected = f"figs/exec13_tN_{x_mm}mm.png"
        if figures != [expected]:
            errors.append(f"tN frame x={x_mm} contains {figures}, expected [{expected}]")

    adaptive_expected = {
        "figs/exec14d_adaptive_tN_summary.png",
        *(f"figs/exec14d_adaptive_tN_{x_mm}mm.png" for x_mm in KEY_POSITIONS),
    }
    adaptive_actual = {
        path for body in bodies[119:] for path in INCLUDE_RE.findall(body)
        if path.startswith("figs/exec14d_adaptive_tN_")
    }
    if adaptive_actual != adaptive_expected:
        errors.append(
            f"appended adaptive figures={sorted(adaptive_actual)}, "
            f"expected={sorted(adaptive_expected)}"
        )

    print(
        f"frames={len(bodies)} dispersion_paths={len(dispersion_paths)} "
        f"position_frames={len(position_re.findall(text))} "
        f"arrival_frames={len(arrival_re.findall(text))} "
        f"tN_frames={len(tn_re.findall(text))} errors={len(errors)}"
    )
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
