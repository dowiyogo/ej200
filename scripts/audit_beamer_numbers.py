#!/usr/bin/env python3
"""Audit CSV provenance, generated TeX consistency, and Beamer numeric inputs."""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
import tempfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
REPORT = REPO / "results_ej230_analysis/report/exec13_ej230_report_full.tex"
TABLES = REPO / "results_ej230_analysis/tables"
SOURCE_RE = re.compile(r"^% source-sha256: (\S+) ([0-9a-f]{64})$", re.MULTILINE)
INPUT_RE = re.compile(r"\\input\{([^{}]+)\}")
NEWCOMMAND_RE = re.compile(r"^\\newcommand\{\\([^{}]+)\}", re.MULTILINE)

REQUIRED = {
    "provenance.tex",
    "key_position_-690.tex", "key_position_-650.tex", "key_position_-400.tex",
    "key_position_0.tex", "key_position_400.tex", "key_position_650.tex",
    "key_position_690.tex", "attenuation_fit.tex", "fano_fit.tex",
    "landau_key_positions.tex", "velocity_estimators.tex", "sum4_timing.tex",
    "tN_summary.tex", "numerical_conclusions.tex", "window_dip_test.tex",
    "endtop_endonly_tails.tex", "endtop_endonly_ratio.tex", "exec14_macros.tex",
    "top_localization_summary.tex", "endtop_physical_diagnosis.tex",
    "adaptive_tN.tex", "top_timing_definition.tex",
}
RESULT_FRAMES = {
    1, 33, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 51, 114, 116, 118,
    *range(5, 30, 4), *range(52, 113, 2),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frame_bodies(text: str) -> list[str]:
    bodies = []
    for match in re.finditer(r"\\begin\{frame\}", text):
        start = text.find("\n", match.end()) + 1
        end = text.find(r"\end{frame}", start)
        bodies.append(text[start:end])
    return bodies


def recursive_input_names(text: str, tables_dir: Path) -> set[str]:
    discovered: set[str] = set()
    pending = [Path(value).name for value in INPUT_RE.findall(text)]
    while pending:
        name = pending.pop()
        if name in discovered:
            continue
        discovered.add(name)
        path = tables_dir / name
        if path.is_file():
            pending.extend(Path(value).name for value in INPUT_RE.findall(path.read_text(encoding="utf-8")))
    return discovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--tables-dir", type=Path, default=TABLES)
    args = parser.parse_args()
    errors: list[str] = []
    report = args.report.read_text(encoding="utf-8")
    bodies = frame_bodies(report)

    if len(bodies) != 119:
        errors.append(f"expected 119 frames, found {len(bodies)}")
    macro_input = r"\input{../tables/exec14_macros.tex}"
    if report.count(macro_input) != 1:
        errors.append(f"report must import {macro_input} exactly once")
    for forbidden in (
        "material comparison",
        "mirror of EXEC\\_12",
        "HOOK_",
        "HOOK\\_",
        "N/A",
        "88 ps",
        "EJ-204",
        "/home/",
        "/mnt/",
    ):
        if forbidden in report:
            errors.append(f"visible standalone report still contains forbidden phrase: {forbidden}")

    report_inputs = recursive_input_names(report, args.tables_dir)
    missing_inputs = sorted(REQUIRED - report_inputs)
    if missing_inputs:
        errors.append(f"required generated tables are not imported: {missing_inputs}")
    absent_inputs = sorted(name for name in report_inputs if not (args.tables_dir / name).is_file())
    if absent_inputs:
        errors.append(f"imported generated tables are absent: {absent_inputs}")

    for number in sorted(RESULT_FRAMES):
        body = bodies[number - 1] if number <= len(bodies) else ""
        if r"\input{" not in body and number != 38:
            errors.append(f"result frame {number} has no generated input")
    for number, body in enumerate(bodies, 1):
        if re.search(r"\d+(?:\.\d+)?\\pm\d", body):
            errors.append(f"frame {number} contains a hardcoded value with uncertainty")

    macro_owners: dict[str, list[str]] = {}
    for path in sorted(args.tables_dir.glob("*.tex")):
        text = path.read_text(encoding="utf-8")
        if re.search(r"\\pm0(?:\.0+)?(?=[^0-9.]|$)", text):
            errors.append(f"{path.name}: contains an uncertainty rounded to zero")
        for macro in NEWCOMMAND_RE.findall(text):
            macro_owners.setdefault(macro, []).append(path.name)
            if not macro.isalpha():
                errors.append(f"{path.name}: invalid non-alphabetic newcommand name: {macro}")
        sources = SOURCE_RE.findall(text)
        if not sources:
            errors.append(f"{path.name}: missing source-sha256 metadata")
        for relative, expected in sources:
            source = REPO / relative
            if not source.is_file():
                errors.append(f"{path.name}: missing source {relative}")
            elif sha256(source) != expected:
                errors.append(f"{path.name}: stale source hash for {relative}")
    duplicates = {name: owners for name, owners in macro_owners.items() if len(owners) > 1}
    if duplicates:
        errors.append(f"duplicate newcommand definitions: {duplicates}")

    with tempfile.TemporaryDirectory(prefix="exec14b_tables_") as temp:
        result = subprocess.run(
            [
                sys.executable,
                str(REPO / "scripts/generate_exec14b_tables.py"),
                "--tables-dir",
                temp,
            ],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            errors.append(f"table regeneration failed: {result.stderr.strip() or result.stdout.strip()}")
        else:
            generated = Path(temp)
            for path in generated.glob("*.tex"):
                if path.name == "provenance.tex":
                    continue
                committed = args.tables_dir / path.name
                if not committed.is_file():
                    errors.append(f"generated table absent from results: {path.name}")
                elif path.read_bytes() != committed.read_bytes():
                    errors.append(f"generated table differs from CSV regeneration: {path.name}")

    print(
        f"frames={len(bodies)} report_inputs={len(report_inputs)} "
        f"tables={len(list(args.tables_dir.glob('*.tex')))} errors={len(errors)}"
    )
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
