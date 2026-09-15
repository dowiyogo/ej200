#!/usr/bin/env python3
"""Validate EXEC_46 per-photon instrumentation and write auditable sidecars."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


TREE_NAME = "sipm_hits"
STEP_SIZE = "64 MB"
PATH_TOLERANCE_MM = 1.0e-6
TIME_TOLERANCE_NS = 1.0e-12
SCINTILLATION_SOURCE = 1
LEFT_FACE = 0
BOUNDARY_EDGES = np.arange(-0.5, 1600.5 + 5.0, 5.0)
FIRST_BOUNDARY_EDGES = np.arange(-0.5, 50.5 + 1.0, 1.0)
PATH_EDGES_MM = np.arange(0.0, 22000.0 + 50.0, 50.0)
WAVELENGTH_EDGES_NM = np.arange(379.5, 500.5 + 1.0, 1.0)
REQUIRED_BRANCHES = {
    "event_id", "track_id", "face_type", "t_creation_ns",
    "t_detection_ns", "x_creation_mm", "y_creation_mm", "z_creation_mm",
    "x_mm", "y_mm", "z_mm", "wl_nm", "wl_nm_created",
    "path_length_mm", "n_boundary_encounters", "exit_angle_deg",
    "source_type",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def distribution(values):
    q16, q50, q84 = np.quantile(values, [0.16, 0.50, 0.84])
    return {
        "n": int(values.size), "mean": float(np.mean(values)),
        "q16": float(q16), "median": float(q50), "q84": float(q84),
        "minimum": float(np.min(values)), "maximum": float(np.max(values)),
    }


def analyze_cell(label, root_path):
    keys = []
    scint_wavelengths = []
    face_boundaries = {face: [] for face in (0, 1, 2)}
    face_paths = {face: [] for face in (0, 1, 2)}
    left_paths = []
    left_boundaries = []
    first_left = {}
    rows = 0
    minimum_path_margin = np.inf
    minimum_creation_time = np.inf
    minimum_propagation_time = np.inf
    minimum_boundaries = np.inf
    maximum_wavelength_change = 0.0
    violations = {name: 0 for name in
                  ("duplicate_track_id", "path", "creation_time", "time_order", "boundary")}

    with uproot.open(root_path) as root_file:
        tree = root_file[TREE_NAME]
        missing = REQUIRED_BRANCHES - set(tree.keys())
        if missing:
            raise RuntimeError(f"{label}: missing branches {sorted(missing)}")
        for arrays in tree.iterate(sorted(REQUIRED_BRANCHES), step_size=STEP_SIZE,
                                   library="np"):
            count = len(arrays["event_id"])
            rows += count
            key = ((arrays["event_id"].astype(np.uint64) << np.uint64(32)) |
                   arrays["track_id"].astype(np.uint32).astype(np.uint64))
            keys.append(key)

            chord = np.sqrt(
                (arrays["x_mm"] - arrays["x_creation_mm"]) ** 2 +
                (arrays["y_mm"] - arrays["y_creation_mm"]) ** 2 +
                (arrays["z_mm"] - arrays["z_creation_mm"]) ** 2)
            margin = arrays["path_length_mm"] - chord
            propagation = arrays["t_detection_ns"] - arrays["t_creation_ns"]
            boundaries = arrays["n_boundary_encounters"]
            minimum_path_margin = min(minimum_path_margin, float(np.min(margin)))
            minimum_creation_time = min(minimum_creation_time,
                                        float(np.min(arrays["t_creation_ns"])))
            minimum_propagation_time = min(minimum_propagation_time,
                                           float(np.min(propagation)))
            minimum_boundaries = min(minimum_boundaries, float(np.min(boundaries)))
            maximum_wavelength_change = max(
                maximum_wavelength_change,
                float(np.max(np.abs(arrays["wl_nm"] - arrays["wl_nm_created"]))))
            violations["path"] += int(np.count_nonzero(
                ~np.isfinite(margin) | (margin < -PATH_TOLERANCE_MM)))
            violations["creation_time"] += int(np.count_nonzero(
                ~np.isfinite(arrays["t_creation_ns"]) |
                (arrays["t_creation_ns"] < -TIME_TOLERANCE_NS)))
            violations["time_order"] += int(np.count_nonzero(
                ~np.isfinite(propagation) | (propagation < -TIME_TOLERANCE_NS)))
            violations["boundary"] += int(np.count_nonzero(
                ~np.isfinite(boundaries.astype(float)) | (boundaries < 0)))

            scint = arrays["source_type"] == SCINTILLATION_SOURCE
            scint_wavelengths.append(arrays["wl_nm_created"][scint])
            for face in (0, 1, 2):
                selected = arrays["face_type"] == face
                face_boundaries[face].append(boundaries[selected])
                face_paths[face].append(arrays["path_length_mm"][selected])

            selected = arrays["face_type"] == LEFT_FACE
            left_paths.append(arrays["path_length_mm"][selected])
            left_boundaries.append(boundaries[selected])
            for event, arrival, boundary, path in zip(
                    arrays["event_id"][selected], arrays["t_detection_ns"][selected],
                    boundaries[selected], arrays["path_length_mm"][selected]):
                old = first_left.get(int(event))
                if old is None or arrival < old[0]:
                    first_left[int(event)] = (float(arrival), int(boundary), float(path))

        if rows != tree.num_entries:
            raise RuntimeError(f"{label}: incomplete read {rows}/{tree.num_entries}")

    all_keys = np.concatenate(keys)
    unique_count = np.unique(all_keys).size
    violations["duplicate_track_id"] = int(rows - unique_count)
    if any(violations.values()):
        raise RuntimeError(f"{label}: hard guardrail failure {violations}")

    scint_wavelengths = np.concatenate(scint_wavelengths)
    left_boundaries = np.concatenate(left_boundaries)
    left_paths = np.concatenate(left_paths)
    first_values = np.asarray(list(first_left.values()), dtype=float)
    face_summary = {}
    for face in (0, 1, 2):
        face_summary[str(face)] = {
            "boundary": distribution(np.concatenate(face_boundaries[face])),
            "path_length_mm": distribution(np.concatenate(face_paths[face])),
        }
    return {
        "label": label, "root": str(Path(root_path).resolve()),
        "root_sha256": sha256(root_path), "rows": rows,
        "violations": violations, "minimum_path_margin_mm": minimum_path_margin,
        "minimum_creation_time_ns": minimum_creation_time,
        "minimum_propagation_time_ns": minimum_propagation_time,
        "minimum_boundary_encounters": minimum_boundaries,
        "maximum_abs_wavelength_change_nm": maximum_wavelength_change,
        "face_summary": face_summary,
        "first_left_boundary": distribution(first_values[:, 1]),
        "first_left_path_mm": distribution(first_values[:, 2]),
        "scintillation_wavelength": distribution(scint_wavelengths),
        "left_boundaries": left_boundaries, "left_paths": left_paths,
        "first_left_boundaries": first_values[:, 1],
        "scint_wavelengths": scint_wavelengths,
    }


def write_outputs(results, spectrum_path, output_dir, metadata):
    output_dir.mkdir(parents=True, exist_ok=True)
    spectrum = np.loadtxt(spectrum_path)
    spectrum_wavelength = spectrum[:, 0]
    spectrum_intensity = spectrum[:, 1]
    configured_peak = float(spectrum_wavelength[np.argmax(spectrum_intensity)])
    configured_min = float(spectrum_wavelength[spectrum_intensity > 0].min())
    configured_max = float(spectrum_wavelength[spectrum_intensity > 0].max())

    for result in results:
        wavelengths = result["scint_wavelengths"]
        counts, edges = np.histogram(wavelengths, bins=WAVELENGTH_EDGES_NM)
        peak = float(0.5 * (edges[np.argmax(counts)] + edges[np.argmax(counts) + 1]))
        result["emission_sanity"] = {
            "configured_peak_nm": configured_peak,
            "detected_scintillation_modal_bin_nm": peak,
            "peak_difference_nm": peak - configured_peak,
            "configured_nonzero_support_nm": [configured_min, configured_max],
            "fraction_inside_configured_support": float(np.mean(
                (wavelengths >= configured_min) & (wavelengths <= configured_max))),
            "pass_support_and_peak": bool(
                np.all((wavelengths >= configured_min) & (wavelengths <= configured_max)) and
                abs(peak - configured_peak) <= 2.0),
            "scope": "Detected scintillation subset; PDE and transport condition the shape.",
        }

    panels = {
        "boundary_left_all": (BOUNDARY_EDGES,
            [np.histogram(r["left_boundaries"], BOUNDARY_EDGES)[0] for r in results]),
        "path_left_all_mm": (PATH_EDGES_MM,
            [np.histogram(r["left_paths"], PATH_EDGES_MM)[0] for r in results]),
        "boundary_left_first": (FIRST_BOUNDARY_EDGES,
            [np.histogram(r["first_left_boundaries"], FIRST_BOUNDARY_EDGES)[0]
             for r in results]),
        "wavelength_scint_nm": (WAVELENGTH_EDGES_NM,
            [np.histogram(r["scint_wavelengths"], WAVELENGTH_EDGES_NM)[0]
             for r in results]),
    }

    csv_path = output_dir / "validation_diagnostics.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["panel", "cell", "bin_low", "bin_high", "count", "density"])
        for panel, (edges, histograms) in panels.items():
            widths = np.diff(edges)
            for result, counts in zip(results, histograms):
                density = counts / (counts.sum() * widths)
                for low, high, count, value in zip(edges[:-1], edges[1:], counts, density):
                    writer.writerow([panel, result["label"], low, high, int(count), value])
        for wavelength, intensity in spectrum:
            writer.writerow(["configured_emission", "OPSC-100", wavelength, wavelength,
                             "", intensity / spectrum_intensity.max()])

    root_path = output_dir / "validation_diagnostics.root"
    with uproot.recreate(root_path) as root_file:
        for panel, (edges, histograms) in panels.items():
            for result, counts in zip(results, histograms):
                root_file[f"{panel}_{result['label']}"] = (counts, edges)
        root_file["configured_emission"] = {
            "wavelength_nm": spectrum_wavelength,
            "relative_intensity": spectrum_intensity,
        }

    colors = ["#255F85", "#D1495B"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)
    for result_index, (result, color) in enumerate(zip(results, colors)):
        label = result["label"].replace("xm", "x=-").replace("x0", "x=0") + " mm"
        for axis, panel in zip(axes.flat[:2], ("boundary_left_all", "path_left_all_mm")):
            edges, histograms = panels[panel]
            counts = histograms[result_index]
            axis.stairs(counts / counts.sum(), edges, label=label, color=color)
        edges, histograms = panels["boundary_left_first"]
        counts = histograms[result_index]
        axes[1, 0].stairs(counts / counts.sum(), edges, label=label, color=color)
        edges, histograms = panels["wavelength_scint_nm"]
        counts = histograms[result_index]
        axes[1, 1].stairs(counts / counts.max(), edges, label=label, color=color)
    axes[0, 0].set(xlabel="Boundary encounters before left detection",
                   ylabel="Fraction / bin", yscale="log")
    axes[0, 0].set_xlim(-0.5, 250.0)
    axes[0, 1].set(xlabel="Path length to left detection [mm]", ylabel="Fraction / bin",
                   yscale="log")
    axes[0, 1].set_xlim(0.0, 2500.0)
    axes[1, 0].set(xlabel="Boundary encounters of first left hit/event",
                   ylabel="Fraction / bin", yscale="log")
    axes[1, 1].plot(spectrum_wavelength, spectrum_intensity / spectrum_intensity.max(),
                    color="black", linestyle="--", label="configured emission")
    axes[1, 1].set(xlabel="Created wavelength [nm]", ylabel="Peak-normalized shape")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    fig.suptitle("EXEC_46 validation: detected-photon trajectory observables")
    fig.savefig(output_dir / "validation_diagnostics.pdf")
    fig.savefig(output_dir / "validation_diagnostics.png", dpi=180)
    plt.close(fig)

    serializable = []
    for result in results:
        serializable.append({key: value for key, value in result.items()
                             if not isinstance(value, np.ndarray)})
    metrics = dict(metadata, spectrum_file=str(Path(spectrum_path).resolve()),
                   spectrum_sha256=sha256(spectrum_path), cells=serializable)
    (output_dir / "validation_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    figure_meta = dict(metadata,
        figure="validation_diagnostics",
        question="Are the EXEC_46 branches populated and physically ordered at center and near end?",
        method=("Full sipm_hits scan; hard guardrails; full-range histograms with display zoom; "
                "detected scintillation selected by source_type=1."),
        inputs=[r["root"] for r in results], csv=str(csv_path), root=str(root_path),
        note=("The configured curve is the source property. Detected photons are conditioned by "
              "transport and PDE, so only support and modal-peak sanity are treated as pass/fail."))
    (output_dir / "validation_diagnostics.meta.json").write_text(
        json.dumps(figure_meta, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--x0", type=Path, required=True)
    parser.add_argument("--xm650", type=Path, required=True)
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    results = [analyze_cell("x0", args.x0), analyze_cell("xm650", args.xm650)]
    metadata = {
        "campaign": "EXEC_46", "created_utc": datetime.now(timezone.utc).isoformat(),
        "simulation_commit": args.commit, "events_per_cell": 500,
        "workers": 4, "eventModulo": 1, "seeds": [26092601, 8349041],
        "configuration": "EndTop", "material": "EJ-200", "opsc": "OPSC-100",
        "jitter_ns": 0.0,
        "analysis_command": " ".join(__import__("sys").argv),
    }
    write_outputs(results, args.spectrum, args.output, metadata)
    print(json.dumps({"status": "PASS", "cells": [
        {"label": r["label"], "rows": r["rows"],
         "violations": r["violations"],
         "first_left_boundary_median": r["first_left_boundary"]["median"]}
        for r in results]}, indent=2))


if __name__ == "__main__":
    main()
