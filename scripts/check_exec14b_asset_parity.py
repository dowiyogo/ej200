#!/usr/bin/env python3
"""Compare the ordered EXEC_12/EJ-230 image-reference manifests."""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
REFERENCE = REPO / "analysis/exec07/exec12_report_full.tex"
REPORT = REPO / "results_ej230_analysis/report/exec13_ej230_report_full.tex"
INCLUDE_RE = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*%?\s*\{([^{}]+)\}", re.DOTALL)


def references(path: Path, frame_limit: int | None = None) -> list[str]:
    text = path.read_text(encoding="utf-8")
    if frame_limit is not None:
        frames = []
        for match in re.finditer(r"\\begin\{frame\}", text):
            end = text.find(r"\end{frame}", match.end())
            frames.append(text[match.start():end + len(r"\end{frame}")])
        text = "\n".join(frames[:frame_limit])
    return [item.strip() for item in INCLUDE_RE.findall(text)]


def normalize(path: str) -> str:
    path = path.replace("exec12_tN", "execXX_tN").replace("exec13_tN", "execXX_tN")
    path = re.sub(r"figs/(?:exec13_)?tn_", "figs/execXX_tn_", path)
    return path


def main() -> int:
    expected = references(REFERENCE)
    observed = references(REPORT, frame_limit=119)
    expected_normalized = [normalize(path) for path in expected]
    observed_normalized = [normalize(path) for path in observed]
    errors: list[str] = []

    if len(expected) != len(observed):
        errors.append(f"reference count={len(expected)}, report count={len(observed)}")
    if Counter(expected_normalized) != Counter(observed_normalized):
        missing = Counter(expected_normalized) - Counter(observed_normalized)
        extra = Counter(observed_normalized) - Counter(expected_normalized)
        errors.append(f"missing normalized references: {dict(missing)}")
        errors.append(f"extra normalized references: {dict(extra)}")

    print(
        f"asset_parity={'PASS' if not errors else 'FAIL'} "
        f"reference_count={len(expected)} report_count={len(observed)} "
        f"reference_unique={len(set(expected_normalized))} "
        f"report_unique={len(set(observed_normalized))}"
    )
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
