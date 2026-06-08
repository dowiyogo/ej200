#!/usr/bin/env python3
"""Extract diagnostic photon budgets and the cheap 31-position efficiencies."""

import argparse
import csv
import re
from pathlib import Path


RUN_SUMMARY = re.compile(
    r"Run ID\s*: (?P<run>\d+).*?"
    r"Events run\s*: (?P<events>\d+).*?"
    r"End-left  photons\s*: (?P<left>\d+).*?"
    r"End-right photons\s*: (?P<right>\d+).*?"
    r"Top SiPM  photons\s*: (?P<top>\d+).*?"
    r"Scint photons generated: (?P<scint>\d+)",
    re.S,
)

BUDGET = re.compile(
    r"=== Photon Budget Diagnostic ==="
    r".*?Generated optical photons\s*: (?P<generated>\d+)"
    r".*?Generated scintillation\s*: (?P<scint>\d+)"
    r".*?Generated Cherenkov\s*: (?P<cherenkov>\d+)"
    r".*?First reach of End face\s*: (?P<end_face>\d+)"
    r".*?Entered End SiPM volume\s*: (?P<entered_end>\d+)"
    r".*?Entered Top SiPM volume\s*: (?P<entered_top>\d+)"
    r".*?Detected End PE\s*: (?P<detected_end>\d+)"
    r".*?Detected Top PE\s*: (?P<detected_top>\d+)"
    r".*?PDE rejected End/Top\s*: (?P<rejected_end>\d+) / (?P<rejected_top>\d+)"
    r".*?Bulk OpAbsorption\s*: (?P<bulk>\d+)"
    r".*?Surface absorption\s*: (?P<surface>\d+)"
    r".*?Direct Bar -> World escape: (?P<escaped>\d+)"
    r".*?Wavelength-filter killed\s*: (?P<wavelength>\d+)",
    re.S,
)


def percent(value, total):
    return 100.0 * value / total if total else 0.0


def master_summaries(text, expected_events):
    return [
        {key: int(value) for key, value in match.groupdict().items()}
        for match in RUN_SUMMARY.finditer(text)
        if int(match.group("events")) == expected_events
    ]


def write_diagnostic(log_path, output_path, positions):
    text = log_path.read_text(errors="replace")
    summaries = master_summaries(text, 1000)
    budgets = [
        {key: int(value) for key, value in match.groupdict().items()}
        for match in BUDGET.finditer(text)
    ]
    if len(summaries) != len(positions) or len(budgets) != len(positions):
        raise RuntimeError(
            f"Expected {len(positions)} summaries/budgets, got "
            f"{len(summaries)}/{len(budgets)}"
        )

    fields = [
        "x_mm", "generated_optical", "generated_scintillation", "generated_cherenkov",
        "first_end_face", "first_end_face_percent", "entered_end", "entered_end_percent",
        "entered_top", "entered_top_percent", "detected_end", "detected_end_percent",
        "detected_top", "detected_top_percent", "rejected_pde_end", "rejected_pde_top",
        "bulk_absorption", "bulk_absorption_percent", "surface_absorption",
        "surface_absorption_percent", "escaped_world", "escaped_world_percent",
        "wavelength_killed", "terminal_accounted_percent", "mean_end_pe_per_event",
        "mean_top_pe_per_event",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for x, summary, budget in zip(positions, summaries, budgets):
            generated = budget["generated"]
            terminal = (
                budget["bulk"] + budget["surface"] + budget["escaped"]
                + budget["wavelength"] + budget["detected_end"] + budget["detected_top"]
                + budget["rejected_end"] + budget["rejected_top"]
            )
            writer.writerow({
                "x_mm": x,
                "generated_optical": generated,
                "generated_scintillation": budget["scint"],
                "generated_cherenkov": budget["cherenkov"],
                "first_end_face": budget["end_face"],
                "first_end_face_percent": f"{percent(budget['end_face'], generated):.6f}",
                "entered_end": budget["entered_end"],
                "entered_end_percent": f"{percent(budget['entered_end'], generated):.6f}",
                "entered_top": budget["entered_top"],
                "entered_top_percent": f"{percent(budget['entered_top'], generated):.6f}",
                "detected_end": budget["detected_end"],
                "detected_end_percent": f"{percent(budget['detected_end'], generated):.6f}",
                "detected_top": budget["detected_top"],
                "detected_top_percent": f"{percent(budget['detected_top'], generated):.6f}",
                "rejected_pde_end": budget["rejected_end"],
                "rejected_pde_top": budget["rejected_top"],
                "bulk_absorption": budget["bulk"],
                "bulk_absorption_percent": f"{percent(budget['bulk'], generated):.6f}",
                "surface_absorption": budget["surface"],
                "surface_absorption_percent": f"{percent(budget['surface'], generated):.6f}",
                "escaped_world": budget["escaped"],
                "escaped_world_percent": f"{percent(budget['escaped'], generated):.6f}",
                "wavelength_killed": budget["wavelength"],
                "terminal_accounted_percent": f"{percent(terminal, generated):.6f}",
                "mean_end_pe_per_event": f"{(summary['left'] + summary['right']) / summary['events']:.6f}",
                "mean_top_pe_per_event": f"{summary['top'] / summary['events']:.6f}",
            })


def write_cheap_efficiency(log_path, output_path, positions):
    text = log_path.read_text(errors="replace")
    summaries = master_summaries(text, 10000)
    if len(summaries) != len(positions):
        raise RuntimeError(f"Expected {len(positions)} master summaries, got {len(summaries)}")
    fields = [
        "x_mm", "generated_scintillation", "end_left_pe_per_event",
        "end_right_pe_per_event", "end_total_pe_per_event", "top_pe_per_event",
        "end_efficiency_percent", "top_efficiency_percent", "total_efficiency_percent",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for x, row in zip(positions, summaries):
            end = row["left"] + row["right"]
            total = end + row["top"]
            events = row["events"]
            writer.writerow({
                "x_mm": x,
                "generated_scintillation": row["scint"],
                "end_left_pe_per_event": f"{row['left'] / events:.6f}",
                "end_right_pe_per_event": f"{row['right'] / events:.6f}",
                "end_total_pe_per_event": f"{end / events:.6f}",
                "top_pe_per_event": f"{row['top'] / events:.6f}",
                "end_efficiency_percent": f"{percent(end, row['scint']):.6f}",
                "top_efficiency_percent": f"{percent(row['top'], row['scint']):.6f}",
                "total_efficiency_percent": f"{percent(total, row['scint']):.6f}",
            })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--diag-log", type=Path, required=True)
    parser.add_argument("--scan-log", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    diag_positions = [0, 350, 650]
    scan_positions = list(range(-650, 651, 50)) + [-690, -670, 670, 690]
    scan_positions.sort()
    write_diagnostic(args.diag_log, args.output_dir / "photon_budget.csv", diag_positions)
    write_cheap_efficiency(
        args.scan_log, args.output_dir / "cheap_efficiency_31_positions.csv", scan_positions
    )


if __name__ == "__main__":
    main()
