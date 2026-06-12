#!/usr/bin/env python3
"""Audit every image reference in the standalone EJ-230 Beamer report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

try:
    from PIL import Image
except ImportError:
    Image = None


REPO = Path(__file__).resolve().parents[1]
DEFAULT_TEX = REPO / "results_ej230_analysis/report/exec13_ej230_report_full.tex"
DEFAULT_CSV = REPO / "results_ej230_analysis/report/figure_manifest.csv"
DEFAULT_MD = REPO / "results_ej230_analysis/report/figure_manifest.md"
FORBIDDEN_DIR = Path("/home/reriosto/SHiP/ej200/analysis/exec07/figs")

INCLUDE_RE = re.compile(
    r"\\includegraphics\s*(?:\[[^\]]*\])?\s*%?\s*\{([^{}]+)\}", re.DOTALL
)
FRAME_RE = re.compile(r"\\begin\{frame\}(?:\[[^\]]*\])?\{([^{}]*)\}")
BACKREF_RE = re.compile(r"\\[1-9]")
GRAPHICSPATH_RE = re.compile(r"\\graphicspath\s*\{((?:\{[^{}]*\}\s*)+)\}")

GENERATORS = {
    "muon_": "analysis/exec07_photon_budget.py",
    "P1_": "analysis/exec07_photon_budget.py",
    "P2_": "analysis/exec07_photon_budget.py",
    "P3_": "analysis/exec07_photon_budget.py",
    "P4_": "analysis/exec07_photon_budget.py",
    "P5_": "analysis/exec07/exec10_time_figures.py",
    "P6_": "analysis/exec07_photon_budget.py",
    "P7_": "analysis/exec07_photon_budget.py",
    "fits_attenuation_velocity": "analysis/exec07_photon_budget.py",
    "exec08b_": "analysis/exec07/exec08b_window_dip.py",
    "exec09_": "analysis/exec07/exec09_timing_mechanism.py",
    "exec10_fano": "analysis/exec07/exec10_landau_analysis.py",
    "exec10_landau": "analysis/exec07/exec10_landau_analysis.py",
    "exec10_velocity": "analysis/exec07/exec10_veff_diagnosis.py",
    "exec11_": "analysis/exec07/exec11_time_arrival.py",
    "exec13_tN": "analysis/exec13/exec13_tN_analysis.py",
    "exec13_tn": "analysis/exec07/exec12b_tn_dispersion.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def image_valid(path: Path) -> tuple[bool, str]:
    if not path.is_file():
        return False, "missing"
    if Image is not None:
        try:
            with Image.open(path) as image:
                image.verify()
                return True, f"Pillow:{image.format}"
        except Exception as exc:
            return False, f"Pillow:{type(exc).__name__}:{exc}"
    result = subprocess.run(
        ["file", "--brief", str(path)], text=True, capture_output=True, check=False
    )
    valid = result.returncode == 0 and "image" in result.stdout.lower()
    return valid, f"file:{result.stdout.strip()}"


def forbidden_hashes() -> dict[str, list[str]]:
    hashes: dict[str, list[str]] = defaultdict(list)
    if FORBIDDEN_DIR.is_dir():
        for path in FORBIDDEN_DIR.rglob("*"):
            if path.is_file():
                try:
                    hashes[sha256(path)].append(str(path))
                except OSError:
                    pass
    return hashes


def graphicspaths(text: str, report_dir: Path) -> list[Path]:
    match = GRAPHICSPATH_RE.search(text)
    if not match:
        return []
    return [(report_dir / item).resolve() for item in re.findall(r"\{([^{}]*)\}", match.group(1))]


def resolve_reference(raw: str, report_dir: Path, paths: list[Path]) -> tuple[Path, list[Path]]:
    candidates = [(report_dir / raw).resolve()]
    candidates.extend((base / raw).resolve() for base in paths)
    for candidate in candidates:
        if candidate.is_file():
            return candidate, candidates
    # Canonical convention for this report: report/../<raw>.
    return (report_dir.parent / raw).resolve(), candidates


def classify(name: str, digest: str, forbidden: dict[str, list[str]]) -> str:
    if digest and digest in forbidden and not name.endswith("_geometry.png"):
        return "FORBIDDEN_EJ204_COPY"
    base = Path(name).name
    if base.startswith(("exec08b_", "exec09_")):
        return "EJ230_SPECIAL_SIM"
    if base.startswith("muon_") and base.endswith("_geometry.png"):
        return "STATIC_GEOMETRY"
    if base.endswith((".png", ".pdf", ".jpg", ".jpeg")):
        return "EJ230_ANALYSIS"
    return "UNKNOWN"


def generator_for(name: str) -> str:
    base = Path(name).name
    for prefix, generator in GENERATORS.items():
        if base.startswith(prefix):
            return generator
    return ""


def frame_at(text: str, position: int) -> tuple[int, str]:
    frames = list(FRAME_RE.finditer(text, 0, position))
    if not frames:
        return 0, "<preamble>"
    return len(frames), frames[-1].group(1)


def write_manifest(rows: list[dict[str, str]], csv_path: Path, md_path: Path) -> None:
    fields = list(rows[0]) if rows else [
        "frame", "frametitle", "include_exact", "relative_path", "expected_absolute",
        "exists", "size_bytes", "extension", "valid", "validation", "sha256",
        "possible_generator", "classification", "search_candidates",
    ]
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with md_path.open("w", encoding="utf-8") as stream:
        stream.write("# EJ-230 Beamer figure manifest\n\n")
        stream.write("| Frame | Frametime | Figure | Exists | Valid | Classification | Generator |\n")
        stream.write("|---:|---|---|---|---|---|---|\n")
        for row in rows:
            title = row["frametitle"].replace("|", r"\|")
            stream.write(
                f"| {row['frame']} | {title} | `{row['relative_path']}` | "
                f"{row['exists']} | {row['valid']} | {row['classification']} | "
                f"`{row['possible_generator']}` |\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tex", type=Path, default=DEFAULT_TEX)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    args = parser.parse_args()

    text = args.tex.read_text(encoding="utf-8")
    report_dir = args.tex.resolve().parent
    paths = graphicspaths(text, report_dir)
    forbidden = forbidden_hashes()
    rows: list[dict[str, str]] = []
    errors: list[str] = []
    resolved_by_path: dict[Path, set[str]] = defaultdict(set)

    backrefs = sorted(set(BACKREF_RE.findall(text)))
    if backrefs:
        errors.append(f"residual regex backreferences in TeX: {', '.join(backrefs)}")

    raw_graphicspath = GRAPHICSPATH_RE.search(text)
    if raw_graphicspath and any(
        match.group(1).strip().startswith("figs/") for match in INCLUDE_RE.finditer(text)
    ):
        entries = re.findall(r"\{([^{}]*)\}", raw_graphicspath.group(1))
        if any(entry.rstrip("/").endswith("figs") for entry in entries):
            errors.append("graphicspath points at figs/ while includegraphics paths already contain figs/")

    for match in INCLUDE_RE.finditer(text):
        raw = match.group(1).strip()
        frame, title = frame_at(text, match.start())
        expected, candidates = resolve_reference(raw, report_dir, paths)
        exists = expected.is_file()
        valid, validation = image_valid(expected)
        digest = sha256(expected) if exists else ""
        category = classify(raw, digest, forbidden)
        if not exists:
            errors.append(f"frame {frame}: missing {raw}")
        elif not valid:
            errors.append(f"frame {frame}: invalid image {raw}: {validation}")
        if category == "FORBIDDEN_EJ204_COPY":
            errors.append(f"frame {frame}: forbidden EJ-204 copy {raw}")
        if exists:
            resolved_by_path[expected].add(raw)
        rows.append({
            "frame": str(frame),
            "frametitle": title,
            "include_exact": match.group(0).replace("\n", "\\n"),
            "relative_path": raw,
            "expected_absolute": str(expected),
            "exists": str(exists),
            "size_bytes": str(expected.stat().st_size if exists else 0),
            "extension": expected.suffix.lower(),
            "valid": str(valid),
            "validation": validation,
            "sha256": digest,
            "possible_generator": generator_for(raw),
            "classification": category,
            "search_candidates": ";".join(str(path) for path in candidates),
        })

    for path, aliases in resolved_by_path.items():
        if len(aliases) > 1:
            errors.append(f"distinct references resolve to one file: {sorted(aliases)} -> {path}")

    write_manifest(rows, args.csv, args.md)
    found = sum(row["exists"] == "True" for row in rows)
    print(f"references={len(rows)} found={found} missing={len(rows) - found} errors={len(errors)}")
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
