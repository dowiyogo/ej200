#!/usr/bin/env python3
"""OCR referenced raster figures and reject inherited visible EJ-204 labels."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "results_ej230_analysis/report/figure_manifest.csv"
FORBIDDEN = re.compile(
    r"EJ[\s_-]*204|OPSC[\s_-]*101|EXEC[\s_\\-]*12|material comparison|mirror of",
    re.IGNORECASE,
)


def ocr(path: Path) -> tuple[Path, str, str]:
    result = subprocess.run(
        ["tesseract", str(path), "stdout", "--psm", "11"],
        text=True,
        capture_output=True,
        check=False,
        timeout=180,
    )
    return path, result.stdout, result.stderr.strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as stream:
        paths = {
            Path(row["expected_absolute"])
            for row in csv.DictReader(stream)
            if row["exists"] == "True" and Path(row["expected_absolute"]).suffix.lower() == ".png"
        }

    errors: list[str] = []
    matches: list[str] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(ocr, path): path for path in sorted(paths)}
        for future in as_completed(futures):
            path = futures[future]
            try:
                _, text, stderr = future.result()
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{path}: OCR failed: {type(exc).__name__}: {exc}")
                continue
            if stderr and "Error" in stderr:
                errors.append(f"{path}: tesseract: {stderr}")
            match = FORBIDDEN.search(text)
            if match:
                matches.append(f"{path}: {match.group(0)!r}")

    print(f"unique_png={len(paths)} forbidden_matches={len(matches)} ocr_errors={len(errors)}")
    for item in matches:
        print(f"ERROR: forbidden raster text: {item}", file=sys.stderr)
    for item in errors:
        print(f"ERROR: {item}", file=sys.stderr)
    return 1 if matches or errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
