#!/usr/bin/env python3
"""Check ordered Beamer frametitle parity against the EXEC_12 reference."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ALLOW_RE = re.compile(r"<!--\s*parity-allow:\s*(\d+)\s*-->")


def balanced_argument(text: str, start: int) -> tuple[str, int]:
    if text[start] != "{":
        raise ValueError("argument does not start with {")
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "{" and (index == 0 or text[index - 1] != "\\"):
            depth += 1
        elif text[index] == "}" and (index == 0 or text[index - 1] != "\\"):
            depth -= 1
            if depth == 0:
                return text[start + 1:index], index + 1
    raise ValueError("unterminated balanced argument")


def titles(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    output = []
    for match in re.finditer(r"\\begin\{frame\}(?:\[[^\]]*\])?", text):
        cursor = match.end()
        while cursor < len(text) and text[cursor].isspace():
            cursor += 1
        title, _ = balanced_argument(text, cursor)
        output.append(title)
    return output


def normalize(title: str) -> str:
    return (
        title.replace("EJ-204", "EJ-MATERIAL")
        .replace("EJ-230", "EJ-MATERIAL")
        .replace("exec12", "execXX")
        .replace("exec13", "execXX")
        .replace("EXEC\\_12", "EXEC\\_XX")
        .replace("EXEC\\_13", "EXEC\\_XX")
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reference", type=Path, default=REPO / "analysis/exec07/exec12_report_full.tex"
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=REPO / "results_ej230_analysis/report/exec13_ej230_report_full.tex",
    )
    parser.add_argument(
        "--allowlist",
        type=Path,
        default=REPO / "docs/exec14d_title_exceptions.md",
    )
    args = parser.parse_args()
    expected, observed = titles(args.reference), titles(args.report)
    allowlisted = {
        int(value)
        for value in ALLOW_RE.findall(args.allowlist.read_text(encoding="utf-8"))
    }
    mismatches = []
    for index, (left, right) in enumerate(zip(expected, observed), 1):
        if normalize(left) != normalize(right):
            mismatches.append((index, left, right))
    if len(observed) < len(expected):
        print(f"ERROR: frame count reference={len(expected)} report={len(observed)}")
        return 1
    appended = observed[len(expected):]
    invalid_appended = [title for title in appended if not title.startswith("Backup: far-End adaptive")]
    if invalid_appended:
        print(f"ERROR: non-adaptive frames appended after primary deck: {invalid_appended}")
        return 1
    unexpected = [item for item in mismatches if item[0] not in allowlisted]
    unused = sorted(allowlisted - {item[0] for item in mismatches})
    if unexpected or unused:
        for index, left, right in unexpected:
            print(f"ERROR frame {index}: reference={left!r} report={right!r}")
        for index in unused:
            print(f"ERROR frame {index}: allowlisted title exception is not present")
        return 1
    exact = len(expected) - len(mismatches)
    print(
        f"frame_parity=PASS primary_frames={len(expected)} appended_frames={len(appended)} "
        f"exact_titles={exact} allowlisted_exceptions={len(mismatches)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
