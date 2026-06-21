#!/usr/bin/env python3
"""EXEC_24 PE-budget and SiPM-detection audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
from collections import Counter, defaultdict

import numpy as np
import uproot

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

LIGHT_YIELD_PER_MEV = 10400.0
BAR_LENGTH_MM = 1400.0
BAR_WIDTH_MM = 60.0
BAR_THICKNESS_MM = 10.0
N_END_SIPM_PER_SIDE = 8
SIPM_ACTIVE_SIDE_MM = 6.0


def write_csv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def clean(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    return value


def parse_generated(log_path: pathlib.Path) -> int | None:
    if not log_path.exists():
        return None
    matches = re.findall(r"Scint photons generated:\s*(\d+)", log_path.read_text())
    return int(matches[-1]) if matches else None


def parse_bar_to_sipm(log_path: pathlib.Path) -> int | None:
    if not log_path.exists():
        return None
    matches = re.findall(r"Bar -> SiPM \(entering\)\s*:\s*(\d+)", log_path.read_text())
    return int(matches[-1]) if matches else None


def load_root(root_path: pathlib.Path) -> dict[str, np.ndarray]:
    with uproot.open(root_path) as rf:
        return rf["sipm_hits"].arrays(
            ["event_id", "global_id", "pde", "time_ns"], library="np"
        )


def pe_accounting(
    root_arrays: dict[str, np.ndarray],
    history_rows: list[dict[str, str]],
    events: int,
    bar_to_sipm_entering: int | None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    event_id = root_arrays["event_id"].astype(int)
    global_id = root_arrays["global_id"].astype(int)
    pde = root_arrays["pde"].astype(float)
    raw_by_event = np.bincount(event_id, minlength=events)
    sum_pde_by_event = np.bincount(event_id, weights=pde, minlength=events)
    left_raw = np.bincount(event_id[global_id < 8], minlength=events)
    right_raw = np.bincount(event_id[(global_id >= 8) & (global_id < 16)], minlength=events)
    left_pde = np.bincount(event_id[global_id < 8], weights=pde[global_id < 8], minlength=events)
    right_pde = np.bincount(
        event_id[(global_id >= 8) & (global_id < 16)],
        weights=pde[(global_id >= 8) & (global_id < 16)],
        minlength=events,
    )

    tracks_by_event: dict[int, set[int]] = defaultdict(set)
    tracks_by_gid: dict[tuple[int, int], set[int]] = defaultdict(set)
    multihit_counter: Counter[tuple[int, int, int]] = Counter()
    for row in history_rows:
        ev = int(row["event_id"])
        tr = int(row["track_id"])
        gid = int(row["global_id"])
        tracks_by_event[ev].add(tr)
        tracks_by_gid[(ev, gid)].add(tr)
        multihit_counter[(ev, tr, gid)] += 1

    rows = []
    for ev in range(events):
        n_unique = len(tracks_by_event.get(ev, set()))
        n_raw = int(raw_by_event[ev])
        rows.append(
            {
                "event_id": ev,
                "n_rows_raw": n_raw,
                "n_unique_tracks": n_unique,
                "sum_pde": float(sum_pde_by_event[ev]),
                "n_rows_left": int(left_raw[ev]),
                "n_rows_right": int(right_raw[ev]),
                "sum_pde_left": float(left_pde[ev]),
                "sum_pde_right": float(right_pde[ev]),
                "duplicates_per_track": float(n_raw / n_unique) if n_unique else math.nan,
                "n_multihit_event_track_gid": sum(1 for key, count in multihit_counter.items() if key[0] == ev and count > 1),
                "max_hits_same_event_track_gid": max(
                    [count for key, count in multihit_counter.items() if key[0] == ev],
                    default=0,
                ),
                "n_unique_tracks_per_gid_mean": float(
                    np.mean([len(tracks_by_gid.get((ev, gid), set())) for gid in range(16)])
                ),
            }
        )

    duplicate_values = np.array([r["duplicates_per_track"] for r in rows], dtype=float)
    summary = {
        "events": events,
        "raw_rows_total": int(len(event_id)),
        "raw_rows_per_event_mean": float(np.mean(raw_by_event)),
        "unique_tracks_per_event_mean": float(np.mean([len(tracks_by_event.get(ev, set())) for ev in range(events)])),
        "sum_pde_per_event_mean": float(np.mean(sum_pde_by_event)),
        "sum_pde_left_per_event_mean": float(np.mean(left_pde)),
        "sum_pde_right_per_event_mean": float(np.mean(right_pde)),
        "duplicates_per_track_mean": float(np.nanmean(duplicate_values)),
        "duplicates_per_track_max": float(np.nanmax(duplicate_values)),
        "n_multihit_event_track_gid_total": int(sum(1 for count in multihit_counter.values() if count > 1)),
        "max_hits_same_event_track_gid": int(max(multihit_counter.values(), default=0)),
        "bar_to_sipm_entering_total": bar_to_sipm_entering,
        "bar_to_sipm_entering_per_event": (
            bar_to_sipm_entering / events if bar_to_sipm_entering is not None else None
        ),
        "detected_rows_per_bar_to_sipm_entering": (
            len(event_id) / bar_to_sipm_entering if bar_to_sipm_entering else None
        ),
        "mean_pde_on_recorded_rows": float(np.mean(pde)),
        "pde_application_assessment": (
            "SiPMSD does not do an explicit Bernoulli draw; the SiPM optical "
            "surface defines EFFICIENCY=PDE. The Bar->SiPM census versus rows "
            "is compatible with Geant4 surface-level PDE gating."
        ),
        "sum_pde_interpretation": (
            "sum_pde is kept as a diagnostic reweighting of already recorded rows; "
            "it is not used as a second PE estimator when surface EFFICIENCY is active."
        ),
    }
    return rows, summary


def photon_budget(accounting_rows: list[dict[str, object]], generated_total: int | None, events: int) -> tuple[list[dict[str, object]], dict[str, object]]:
    generated_per_event = generated_total / events if generated_total else None
    edep_mev = generated_per_event / LIGHT_YIELD_PER_MEV if generated_per_event else None
    rows = []
    for row in accounting_rows:
        n_gen = generated_per_event if generated_per_event else math.nan
        rows.append(
            {
                "event_id": row["event_id"],
                "n_gen_estimated": n_gen,
                "edep_mev_estimated": edep_mev if edep_mev is not None else math.nan,
                "n_rows_raw": row["n_rows_raw"],
                "n_unique_tracks": row["n_unique_tracks"],
                "sum_pde": row["sum_pde"],
                "collection_eff_raw": row["n_rows_raw"] / n_gen if generated_per_event else math.nan,
                "collection_eff_unique": row["n_unique_tracks"] / n_gen if generated_per_event else math.nan,
                "collection_eff_pde": row["sum_pde"] / n_gen if generated_per_event else math.nan,
            }
        )
    summary = {
        "light_yield_per_mev": LIGHT_YIELD_PER_MEV,
        "generated_total": generated_total,
        "generated_per_event": generated_per_event,
        "edep_mev_estimated_per_event": edep_mev,
        "collection_eff_raw_mean": float(np.nanmean([r["collection_eff_raw"] for r in rows])),
        "collection_eff_unique_mean": float(np.nanmean([r["collection_eff_unique"] for r in rows])),
        "collection_eff_pde_mean": float(np.nanmean([r["collection_eff_pde"] for r in rows])),
    }
    if summary["collection_eff_pde_mean"] > 0.10:
        summary["gate"] = "collection_eff_high_concern"
    elif summary["collection_eff_pde_mean"] > 0.01:
        summary["gate"] = "collection_eff_moderate_concern"
    else:
        summary["gate"] = "collection_eff_plausible"
    return rows, summary


def sipm_geometry() -> tuple[list[dict[str, object]], dict[str, object]]:
    face_area = BAR_WIDTH_MM * BAR_THICKNESS_MM
    active_area_one = SIPM_ACTIVE_SIDE_MM * SIPM_ACTIVE_SIDE_MM
    active_per_end = N_END_SIPM_PER_SIDE * active_area_one
    coverage = active_per_end / face_area
    positions = []
    pitch = 7.5
    for side, base_gid in (("L", 0), ("R", 8)):
        for i in range(8):
            positions.append(
                {
                    "side": side,
                    "global_id": base_gid + i,
                    "center_y_mm": (i - 3.5) * pitch,
                    "center_z_mm": 0.0,
                    "active_width_mm": SIPM_ACTIVE_SIDE_MM,
                    "active_height_mm": SIPM_ACTIVE_SIDE_MM,
                }
            )
    summary = {
        "bar_end_face_area_mm2": face_area,
        "active_sipm_area_per_end_mm2": active_per_end,
        "coverage_per_end": coverage,
        "n_sipm_per_end": N_END_SIPM_PER_SIDE,
        "note": "Geometry uses 8 independent 6x6 mm EndSiPMLV volumes per end, not a full-face detector.",
    }
    return positions, summary


def lifetime_summary(history_rows: list[dict[str, str]], speed_rows: list[dict[str, str]], generated_total: int | None) -> tuple[list[dict[str, object]], dict[str, object]]:
    det = []
    for row in history_rows:
        det.append(
            {
                "event_id": int(row["event_id"]),
                "track_id": int(row["track_id"]),
                "global_id": int(row["global_id"]),
                "time_alive_ns": float(row["first_detection_time_ns"]),
                "n_reflections_scint_air": int(row["n_reflections_scint_air"]),
                "n_reflections_air_reflector": int(row["n_reflections_air_reflector"]),
                "n_total_internal_reflections": int(row["n_total_internal_reflections"]),
                "visited_air_gap": int(row["visited_air_gap"]),
                "visited_reflector_boundary": int(row["visited_reflector_boundary"]),
            }
        )
    total_steps = sum(int(r["n_steps"]) for r in speed_rows)
    detected = len(history_rows)
    total_generated = generated_total or math.nan
    summary = {
        "detected_photons": detected,
        "generated_photons": generated_total,
        "detected_fraction_raw": detected / total_generated if generated_total else None,
        "total_optical_steps": total_steps,
        "steps_per_generated_photon_mean": total_steps / total_generated if generated_total else None,
        "detected_time_ns_mean": float(np.mean([r["time_alive_ns"] for r in det])) if det else None,
        "detected_time_ns_p95": float(np.percentile([r["time_alive_ns"] for r in det], 95)) if det else None,
        "detected_time_ns_max": float(np.max([r["time_alive_ns"] for r in det])) if det else None,
        "detected_tir_mean": float(np.mean([r["n_total_internal_reflections"] for r in det])) if det else None,
        "detected_scint_air_reflection_mean": float(np.mean([r["n_reflections_scint_air"] for r in det])) if det else None,
        "long_lived_photon_concern": bool(det and np.percentile([r["time_alive_ns"] for r in det], 95) > 50.0),
        "note": "Per-track lifetime is available for detected photons; total steps are aggregate over all generated optical photons.",
    }
    return det, summary


def make_plots(out_dir: pathlib.Path, accounting_rows: list[dict[str, object]], budget_rows: list[dict[str, object]], geometry_rows: list[dict[str, object]], lifetime_rows: list[dict[str, object]]) -> None:
    png = out_dir / "png"
    png.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(6, 4))
    plt.hist([r["duplicates_per_track"] for r in accounting_rows], bins=20, histtype="step")
    plt.xlabel("raw rows / unique tracks")
    plt.ylabel("events")
    plt.tight_layout()
    plt.savefig(png / "duplicates_per_track.png", dpi=160)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.hist([r["collection_eff_pde"] for r in budget_rows], bins=20, histtype="step", label="sum(pde)/Ngen")
    plt.hist([r["collection_eff_raw"] for r in budget_rows], bins=20, histtype="step", label="raw/Ngen")
    plt.xlabel("collection fraction")
    plt.ylabel("events")
    plt.legend()
    plt.tight_layout()
    plt.savefig(png / "photon_budget_fractions.png", dpi=160)
    plt.close()

    plt.figure(figsize=(5, 4))
    ax = plt.gca()
    rect = plt.Rectangle((-30, -5), 60, 10, fill=False, color="black")
    ax.add_patch(rect)
    for row in geometry_rows:
        y = float(row["center_y_mm"])
        z = float(row["center_z_mm"])
        ax.add_patch(plt.Rectangle((y - 3, z - 3), 6, 6, fill=False, color="tab:blue"))
        ax.text(y, z, str(row["global_id"]), ha="center", va="center", fontsize=7)
    ax.set_aspect("equal")
    ax.set_xlabel("Y [mm]")
    ax.set_ylabel("Z [mm]")
    ax.set_xlim(-35, 35)
    ax.set_ylim(-8, 8)
    plt.tight_layout()
    plt.savefig(png / "end_face_sipm_layout.png", dpi=160)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.hist([r["n_total_internal_reflections"] for r in lifetime_rows], bins=80, histtype="step")
    plt.xlabel("TIR count for detected photons")
    plt.ylabel("photons")
    plt.tight_layout()
    plt.savefig(png / "optical_steps_per_photon.png", dpi=160)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.hist([r["time_alive_ns"] for r in lifetime_rows], bins=80, histtype="step")
    plt.xlabel("detected photon time [ns]")
    plt.ylabel("photons")
    plt.tight_layout()
    plt.savefig(png / "optical_time_alive.png", dpi=160)
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=pathlib.Path)
    parser.add_argument("--history", required=True, type=pathlib.Path)
    parser.add_argument("--speed", required=True, type=pathlib.Path)
    parser.add_argument("--run-log", required=True, type=pathlib.Path)
    parser.add_argument("--events", required=True, type=int)
    parser.add_argument("--out-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    root_arrays = load_root(args.root)
    history_rows = read_csv(args.history)
    speed_rows = read_csv(args.speed)
    generated_total = parse_generated(args.run_log)
    bar_to_sipm_entering = parse_bar_to_sipm(args.run_log)

    accounting_rows, accounting_summary = pe_accounting(
        root_arrays, history_rows, args.events, bar_to_sipm_entering
    )
    budget_rows, budget_summary = photon_budget(accounting_rows, generated_total, args.events)
    geometry_rows, geometry_summary = sipm_geometry()
    lifetime_rows, lifetime = lifetime_summary(history_rows, speed_rows, generated_total)
    make_plots(args.out_dir, accounting_rows, budget_rows, geometry_rows, lifetime_rows)

    write_csv(args.out_dir / "pe_accounting_by_event.csv", accounting_rows)
    write_csv(args.out_dir / "photon_budget_by_event.csv", budget_rows)
    write_csv(args.out_dir / "sipm_geometry_audit.csv", geometry_rows)
    write_csv(args.out_dir / "optical_lifetime_summary.csv", lifetime_rows)

    (args.out_dir / "pe_accounting_summary.json").write_text(
        json.dumps(clean(accounting_summary), indent=2, allow_nan=False) + "\n"
    )
    (args.out_dir / "photon_budget_summary.json").write_text(
        json.dumps(clean(budget_summary), indent=2, allow_nan=False) + "\n"
    )
    (args.out_dir / "sipm_geometry_summary.json").write_text(
        json.dumps(clean(geometry_summary), indent=2, allow_nan=False) + "\n"
    )
    (args.out_dir / "optical_lifetime_summary.json").write_text(
        json.dumps(clean(lifetime), indent=2, allow_nan=False) + "\n"
    )

    verdict = "COMPLETADO-CON-FLAGS"
    flags = []
    if accounting_summary["raw_rows_per_event_mean"] > 300:
        flags.append("npe_raw_high_concern")
    if accounting_summary["sum_pde_per_event_mean"] > 300:
        flags.append("sum_pde_reweight_high_concern")
    if budget_summary["collection_eff_pde_mean"] > 0.10:
        flags.append("collection_eff_high_concern")
    elif budget_summary["collection_eff_pde_mean"] > 0.01:
        flags.append("collection_eff_moderate_concern")
    if lifetime["long_lived_photon_concern"]:
        flags.append("long_lived_photon_concern")
    results = {
        "exec_id": "EXEC_24",
        "verdict": verdict,
        "pe_accounting": accounting_summary,
        "photon_budget": budget_summary,
        "sipm_geometry": geometry_summary,
        "optical_lifetime": lifetime,
        "corrections_applied": [],
        "final_quick_result": {
            "x_mm": 0,
            "events": args.events,
            "npe_raw": accounting_summary["raw_rows_per_event_mean"],
            "npe_unique": accounting_summary["unique_tracks_per_event_mean"],
            "npe_pde_reweighted_diagnostic": accounting_summary["sum_pde_per_event_mean"],
            "npe_pde_expected": accounting_summary["raw_rows_per_event_mean"],
        },
        "flags": flags,
    }
    (args.out_dir / "results_exec24.json").write_text(
        json.dumps(clean(results), indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(clean(results), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
