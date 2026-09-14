#!/usr/bin/env python3
"""Preserve non-ROOT EXEC40--43 analysis products in the repository."""

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
ARCHIVE = Path(__file__).resolve().parent
CAMPAIGNS = {
    "exec40": (Path("/home/rrios/exec40_20260913"), [Path("/home/rrios/REPORT_EXEC40_20260913.md")]),
    "exec41": (Path("/home/rrios/exec41_20260913"), [Path("/home/rrios/REPORT_EXEC41_20260913.md")]),
    "exec42": (
        Path("/home/rrios/exec42_20260913"),
        [Path("/home/rrios/REPORT_EXEC42A_20260913.md"), Path("/home/rrios/REPORT_EXEC42_20260914.md")],
    ),
    "exec43": (Path("/home/rrios/exec43_20260914"), [Path("/home/rrios/REPORT_EXEC43_20260914.md")]),
}
ALLOWED_SUFFIXES = {".csv", ".json", ".png", ".pdf", ".svg"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    summary = {
        "schema": "ej200.exec40_43.nonroot_archive.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "measurement_command": "du -sh --exclude='*.root' /home/rrios/exec4[0-3]_2026*/",
        "selection": "All CSV, JSON/metadata, PNG, PDF, and SVG files under each campaign root, plus published Markdown reports.",
        "excluded": {
            "ROOT": "Excluded from git; the transport and analysis ROOT corpus is approximately 147 GB.",
            "other_nonroot": "Build products, source mirrors, macros, and logs are outside the requested reports/analysis-output set.",
        },
        "campaigns": {},
    }
    for name, (source, reports) in CAMPAIGNS.items():
        destination = ARCHIVE / name
        if destination.exists():
            raise SystemExit(f"refusing to overwrite existing archive: {destination}")
        all_nonroot = [p for p in source.rglob("*") if p.is_file() and p.suffix.lower() != ".root"]
        selected = [p for p in all_nonroot if p.suffix.lower() in ALLOWED_SUFFIXES]
        records = []
        for path in sorted(selected):
            relative = path.relative_to(source)
            target = destination / "campaign" / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            records.append({
                "source": str(path),
                "archived": str(target.relative_to(REPO)),
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            })
        for report in reports:
            target = destination / report.name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(report, target)
            records.append({
                "source": str(report),
                "archived": str(target.relative_to(REPO)),
                "bytes": report.stat().st_size,
                "sha256": digest(report),
            })
        entry = {
            "source_root": str(source),
            "nonroot_disk_file_count": len(all_nonroot),
            "nonroot_disk_bytes": sum(p.stat().st_size for p in all_nonroot),
            "archived_file_count": len(records),
            "archived_bytes": sum(row["bytes"] for row in records),
            "files": records,
        }
        (destination / "INDEX.json").write_text(json.dumps(entry, indent=2) + "\n")
        summary["campaigns"][name] = {k: v for k, v in entry.items() if k != "files"}
    summary["nonroot_disk_total_bytes"] = sum(v["nonroot_disk_bytes"] for v in summary["campaigns"].values())
    summary["archived_total_bytes"] = sum(v["archived_bytes"] for v in summary["campaigns"].values())
    (ARCHIVE / "EXEC40_43_ARCHIVE.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
