#!/usr/bin/env python3
"""Machine-readable PDF preflight for the standalone EJ-230 report."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import fitz


REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results_ej230_analysis"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdf", type=Path, default=RESULTS / "report/exec13_ej230_report_full.pdf"
    )
    parser.add_argument(
        "--manifest", type=Path, default=RESULTS / "report/figure_manifest.csv"
    )
    parser.add_argument(
        "--output", type=Path, default=RESULTS / "report/pdf_preflight.json"
    )
    args = parser.parse_args()

    expected_image_pages: set[int] = set()
    with args.manifest.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            expected_image_pages.add(int(row["frame"]))

    document = fitz.open(args.pdf)
    literal_figs = []
    filename_text = []
    pages_with_images = []
    image_counts = {}
    suspicious = []
    for number, page in enumerate(document, 1):
        text = page.get_text()
        images = page.get_images(full=True)
        image_counts[str(number)] = len(images)
        if images:
            pages_with_images.append(number)
        if "figs/" in text:
            literal_figs.append(number)
        if re.search(r"(?:exec|muon|P[1-9]_)[^\s]*\.(?:png|pdf|jpg)", text):
            filename_text.append(number)
        if number in expected_image_pages and not images:
            suspicious.append(number)

    result = {
        "pdf": str(args.pdf.resolve()),
        "pages_total": len(document),
        "pages_with_raster_images": len(pages_with_images),
        "total_page_image_references": sum(image_counts.values()),
        "pages_with_literal_figs": literal_figs,
        "pages_with_filename_text": filename_text,
        "pages_with_suspicious_empty_image_boxes": suspicious,
        "image_counts_by_page": image_counts,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items() if key != "image_counts_by_page"}, indent=2))
    passed = (
        result["pages_total"] >= 119
        and not literal_figs
        and not filename_text
        and not suspicious
        and bool(pages_with_images)
    )
    print(f"pdf_preflight={'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
