#!/usr/bin/env python3
"""
Generate a Beamer presentation from the manifest written by
fpt_vs_n_profile_batch_slides.C.

The manifest carries scan metadata explicitly, so this script does not infer
captions or ordering from PNG filenames.

Usage:
    python3 fpt_manifest_to_beamer.py fpt_png/fpt_beamer_manifest.tsv \
        presentacion_timing_detector_scan.tex
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path


REQUIRED_COLUMNS = {
    "slide_index",
    "run",
    "gun_x_mm",
    "face",
    "image_path",
    "caption",
}


def latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def latex_path(path: Path) -> str:
    return str(path).replace(os.sep, "/")


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Manifest vacio: {path}")

        missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValueError(f"Manifest sin columnas requeridas: {missing_list}")

        rows = [row for row in reader if row.get("image_path")]

    rows.sort(key=lambda row: int(row["slide_index"]))
    return rows


def resolve_image_path(raw_path: str, manifest_path: Path) -> Path:
    image_path = Path(raw_path)
    if image_path.is_absolute():
        return image_path

    candidates = [
        Path.cwd() / image_path,
        manifest_path.parent / image_path,
        manifest_path.parent / image_path.name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    return candidates[0]


def image_path_for_tex(raw_path: str, manifest_path: Path, output_tex: Path) -> str:
    image_path = resolve_image_path(raw_path, manifest_path)

    rel_path = os.path.relpath(image_path, output_tex.parent)
    return latex_path(Path(rel_path))


def frame_title(row: dict[str, str]) -> str:
    return (
        f"Run {int(row['run']):03d} | "
        f"x = {int(round(float(row['gun_x_mm']))):+d} mm | "
        f"{row['face']}"
    )


def generate_beamer(
    rows: list[dict[str, str]],
    manifest_path: Path,
    output_tex: Path,
    title: str,
    author: str,
) -> None:
    if not rows:
        raise ValueError("Manifest sin slides")

    lines: list[str] = []
    lines.append(
        r"""\documentclass{beamer}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[spanish]{babel}
\usepackage{graphicx}

\usetheme{Madrid}
\usecolortheme{default}

"""
    )
    lines.append(f"\\title{{{latex_escape(title)}}}\n")
    lines.append(f"\\author{{{latex_escape(author)}}}\n")
    lines.append("\\date{\\today}\n\n")
    lines.append(
        r"""\begin{document}

\begin{frame}
  \titlepage
\end{frame}

\begin{frame}
  \frametitle{Contenido}
  \tableofcontents
\end{frame}

\section{FPT vs n}

"""
    )

    for row in rows:
        img = image_path_for_tex(row["image_path"], manifest_path, output_tex)
        title_text = latex_escape(frame_title(row))
        caption = latex_escape(row["caption"])
        lines.append(
            f"""\\begin{{frame}}[fragile]
  \\frametitle{{{title_text}}}
  \\begin{{center}}
    \\includegraphics[width=0.90\\textwidth,height=0.80\\textheight,keepaspectratio]{{{img}}}

    \\vspace{{0.1cm}}
    \\small {caption}
  \\end{{center}}
\\end{{frame}}

"""
        )

    lines.append("\\end{document}\n")
    output_tex.parent.mkdir(parents=True, exist_ok=True)
    output_tex.write_text("".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a Beamer .tex file from an FPT scan manifest."
    )
    parser.add_argument("manifest", type=Path, help="TSV manifest from ROOT batch macro")
    parser.add_argument("output_tex", type=Path, help="Output Beamer .tex file")
    parser.add_argument(
        "--title",
        default="Timing Detector Analysis",
        help="Presentation title",
    )
    parser.add_argument(
        "--author",
        default="SHiP Experiment",
        help="Presentation author",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        rows = read_manifest(args.manifest)
        generate_beamer(rows, args.manifest, args.output_tex, args.title, args.author)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Presentacion generada: {args.output_tex}")
    print(f"Slides creadas: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
