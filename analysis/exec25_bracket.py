#!/usr/bin/env python3
"""EXEC_25 PDE semantics and optical realism bracket runner."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import time
import zlib
from dataclasses import dataclass

import numpy as np
import uproot

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

LIGHT_YIELD_PER_MEV = 10400.0
DEFAULT_THREADS = 16
WINDOWS_NS = (5.0, 10.0, 20.0, 50.0)


@dataclass(frozen=True)
class Config:
    name: str
    stage: str
    events: int
    sipm_eff_mode: str = "nominal"
    reflectivity: float = 0.95
    finish: str = "polished"
    sigma_alpha_deg: float = 0.1
    air_gap_mm: float = 0.10
    specular_lobe: float = 1.0
    x_mm: float = 0.0


def write_csv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: pathlib.Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(payload), indent=2, allow_nan=False) + "\n")


def clean(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    return value


def parse_last_int(pattern: str, text: str) -> int | None:
    matches = re.findall(pattern, text)
    return int(matches[-1]) if matches else None


def macro_text(config: Config, seed_base: int) -> str:
    seed1 = seed_base + zlib.crc32(config.name.encode("utf-8")) % 100000
    seed2 = seed_base + 200000 + zlib.crc32(config.name[::-1].encode("utf-8")) % 100000
    return f"""/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/det/readout EndTop
/ship/geom/topSurface mylar
/ship/geom/sipmEfficiencyMode {config.sipm_eff_mode}
/ship/geom/mylar/reflectivity {config.reflectivity}
/ship/geom/mylar/specularLobe {config.specular_lobe}
/ship/geom/mylar/sigmaAlpha {config.sigma_alpha_deg} deg
/ship/geom/mylar/finish {config.finish}
/ship/geom/mylar/airGap {config.air_gap_mm} mm
/det/scintillator OPSC-101
/sipm/model AFBR-S4N66P024M
/run/initialize
/sipm/jitterSigma 0 ns
/gun/particle mu-
/gun/energy 1 GeV
/muon/angle 0
/random/setSeeds {seed1} {seed2}
/muon/gunX {config.x_mm} mm
/run/beamOn {config.events}
"""


def run_config(config: Config, repo: pathlib.Path, out_dir: pathlib.Path, threads: int) -> dict[str, object]:
    sim = repo / "build" / "ej200_bar_sim"
    raw_dir = out_dir / "raw" / config.name
    work = raw_dir / "work"
    audit = raw_dir / "audit"
    for path in (work, audit):
        path.mkdir(parents=True, exist_ok=True)
    sslg4_link = work / "sslg4"
    if not sslg4_link.exists():
        sslg4_link.symlink_to(repo / "build" / "sslg4")
    macro = work / "run.mac"
    macro.write_text(macro_text(config, seed_base=810000))
    log = raw_dir / "run.log"
    env = os.environ.copy()
    env["SHIP_OPTICAL_AUDIT_DIR"] = str(audit)
    env["SHIP_OPTICAL_AUDIT_TAG"] = config.name
    start = time.monotonic()
    with log.open("w") as handle:
        proc = subprocess.run(
            [str(sim), "-t", str(threads), "-m", str(macro.name)],
            cwd=work,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    runtime_s = time.monotonic() - start
    if proc.returncode != 0:
        raise RuntimeError(f"{config.name}: simulation failed, see {log}")
    root_src = work / "photon_hits_run000.root"
    root_final = raw_dir / "photon_hits.root"
    if root_src.exists():
        if root_final.exists():
            root_final.unlink()
        shutil.move(str(root_src), str(root_final))
    return analyze_config(config, root_final, log, audit, runtime_s)


def read_root(root_path: pathlib.Path) -> dict[str, np.ndarray]:
    if not root_path.exists() or root_path.stat().st_size == 0:
        return {
            "event_id": np.array([], dtype=np.int64),
            "global_id": np.array([], dtype=np.int64),
            "time_ns": np.array([], dtype=float),
            "pde": np.array([], dtype=float),
        }
    with uproot.open(root_path) as handle:
        tree = handle["sipm_hits"]
        return tree.arrays(["event_id", "global_id", "time_ns", "pde"], library="np")


def read_csv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sigma_end_ps(event_id: np.ndarray, global_id: np.ndarray, time_ns: np.ndarray, events: int) -> float | None:
    values = []
    for event in range(events):
        mask = event_id == event
        left = time_ns[mask & (global_id < 8)]
        right = time_ns[mask & (global_id >= 8) & (global_id < 16)]
        if left.size and right.size:
            values.append(0.5 * (float(np.min(left)) + float(np.min(right))))
    if len(values) < 5:
        return None
    arr = np.asarray(values, dtype=float)
    med = float(np.median(arr))
    mad = 1.4826 * float(np.median(np.abs(arr - med)))
    core = arr[np.abs(arr - med) <= 2.0 * max(mad, 1.0e-9)]
    if core.size < 2:
        return None
    return 1000.0 * float(np.std(core, ddof=1))


def analyze_config(
    config: Config,
    root_path: pathlib.Path,
    log_path: pathlib.Path,
    audit_dir: pathlib.Path,
    runtime_s: float,
) -> dict[str, object]:
    text = log_path.read_text(errors="replace") if log_path.exists() else ""
    arrays = read_root(root_path)
    event_id = arrays["event_id"].astype(np.int64)
    global_id = arrays["global_id"].astype(np.int64)
    time_ns = arrays["time_ns"].astype(float)
    pde = arrays["pde"].astype(float)
    end_mask = (global_id >= 0) & (global_id < 16)
    left_counts = np.bincount(event_id[end_mask & (global_id < 8)], minlength=config.events)
    right_counts = np.bincount(event_id[end_mask & (global_id >= 8)], minlength=config.events)
    generated = parse_last_int(r"Scint photons generated:\s*(\d+)", text)
    entering = parse_last_int(r"Bar -> SiPM \(entering\)\s*:\s*(\d+)", text)
    bar_to_reflector = parse_last_int(r"Bar -> reflector panel\s*:\s*(\d+)", text)
    killed_world = parse_last_int(r"Killed in WorldLV\s*:\s*(\d+)", text)
    history = read_csv(audit_dir / "detected_photon_history.csv")
    track_keys = {(row["event_id"], row["track_id"]) for row in history}
    multi = {}
    for row in history:
        key = (row["event_id"], row["track_id"], row["global_id"])
        multi[key] = multi.get(key, 0) + 1
    speed_rows = read_csv(audit_dir / "speed_audit_steps.csv")
    total_steps = sum(int(row["n_steps"]) for row in speed_rows)
    missing_rindex = sum(int(row["missing_rindex"]) for row in speed_rows)
    max_super = max((float(row["superluminal_fraction"]) for row in speed_rows), default=0.0)
    recorded = int(end_mask.sum())
    paired = np.count_nonzero((left_counts > 0) & (right_counts > 0))
    sigma_ps = sigma_end_ps(event_id[end_mask], global_id[end_mask], time_ns[end_mask], config.events)
    detected_times = time_ns[end_mask]
    return {
        "name": config.name,
        "stage": config.stage,
        "events": config.events,
        "x_mm": config.x_mm,
        "sipm_eff_mode": config.sipm_eff_mode,
        "reflectivity": config.reflectivity,
        "finish": config.finish,
        "sigma_alpha_deg": config.sigma_alpha_deg,
        "air_gap_mm": config.air_gap_mm,
        "specular_lobe": config.specular_lobe,
        "runtime_s": runtime_s,
        "n_gen_total": generated,
        "n_gen_per_event": generated / config.events if generated else None,
        "edep_mev_per_event": generated / config.events / LIGHT_YIELD_PER_MEV if generated else None,
        "bar_to_sipm_entering_total": entering,
        "bar_to_sipm_entering_per_event": entering / config.events if entering else None,
        "recorded_hits_total": recorded,
        "npe_end": recorded / config.events,
        "npe_l": float(np.mean(left_counts)),
        "npe_r": float(np.mean(right_counts)),
        "recorded_over_generated": recorded / generated if generated else None,
        "recorded_over_entering": recorded / entering if entering else None,
        "mean_pde_branch": float(np.mean(pde)) if pde.size else None,
        "unique_tracks_total": len(track_keys),
        "duplicates_per_track": recorded / len(track_keys) if track_keys else (0.0 if recorded == 0 else None),
        "multihit_event_track_gid": sum(1 for value in multi.values() if value > 1),
        "paired_fraction": paired / config.events,
        "sigma_end_ps": sigma_ps,
        "time_mean_ns": float(np.mean(detected_times)) if detected_times.size else None,
        "time_p95_ns": float(np.percentile(detected_times, 95)) if detected_times.size else None,
        "time_max_ns": float(np.max(detected_times)) if detected_times.size else None,
        "steps_per_generated_photon": total_steps / generated if generated else None,
        "total_optical_steps": total_steps,
        "missing_rindex": missing_rindex,
        "max_superluminal_fraction": max_super,
        "speed_gate": "PASS" if missing_rindex == 0 and max_super <= 0.0 else "FAIL",
        "bar_to_reflector_total": bar_to_reflector,
        "killed_world_total": killed_world,
        "root_path": str(root_path),
        "log_path": str(log_path),
        "audit_dir": str(audit_dir),
    }


def pde_configs() -> list[Config]:
    return [
        Config("pde_nominal", "pde_semantics", 10, "nominal"),
        Config("pde_zero", "pde_semantics", 10, "zero"),
        Config("pde_one", "pde_semantics", 10, "one"),
    ]


def bracket_configs() -> list[Config]:
    configs = []
    for r in (0.00, 0.50, 0.70, 0.85, 0.90, 0.95, 0.98):
        configs.append(Config(f"R{int(r * 100):03d}_polished", "reflectivity", 20, reflectivity=r))
    for name, finish, sigma in (
        ("rough_polished_s0p0", "polished", 0.0),
        ("rough_polished_s0p1", "polished", 0.1),
        ("rough_polished_s0p3", "polished", 0.3),
        ("rough_ground_s0p1", "ground", 0.1),
    ):
        configs.append(Config(name, "roughness", 20, reflectivity=0.95, finish=finish, sigma_alpha_deg=sigma))
    for gap in (0.02, 0.05, 0.10, 0.20):
        configs.append(Config(f"gap_{str(gap).replace('.', 'p')}mm", "gap", 20, air_gap_mm=gap))
    return configs


def pde_pass(rows: list[dict[str, object]]) -> tuple[bool, dict[str, object]]:
    by_name = {row["name"]: row for row in rows}
    zero = by_name["pde_zero"]
    one = by_name["pde_one"]
    nominal = by_name["pde_nominal"]
    zero_ratio = zero["recorded_over_entering"] or 0.0
    one_ratio = one["recorded_over_entering"] or 0.0
    nominal_ratio = nominal["recorded_over_entering"] or 0.0
    nominal_pde = nominal["mean_pde_branch"] or 0.0
    passed = (
        zero["recorded_hits_total"] <= 1
        and one_ratio >= 0.95
        and abs(nominal_ratio - nominal_pde) <= 0.10
    )
    return passed, {
        "pass": passed,
        "zero_recorded_over_entering": zero_ratio,
        "one_recorded_over_entering": one_ratio,
        "nominal_recorded_over_entering": nominal_ratio,
        "nominal_mean_pde_branch": nominal_pde,
        "zero_recorded_hits": zero["recorded_hits_total"],
        "one_recorded_hits": one["recorded_hits_total"],
        "nominal_recorded_hits": nominal["recorded_hits_total"],
    }


def candidate_from(rows: list[dict[str, object]]) -> dict[str, object]:
    eligible = [
        row for row in rows
        if row["stage"] in {"reflectivity", "roughness", "gap"}
        and row["speed_gate"] == "PASS"
        and row["missing_rindex"] == 0
        and row["duplicates_per_track"] in (0.0, 1.0)
        and 20.0 <= float(row["npe_end"]) <= 300.0
        and float(row["paired_fraction"]) >= 0.8
    ]
    if not eligible:
        return {
            "exists": False,
            "name": None,
            "reflectivity": None,
            "finish": None,
            "sigma_alpha": None,
            "air_gap_mm": None,
            "npe_end_center": None,
            "sigma_end_center_ps": None,
            "runtime_100ev_s": None,
        }
    eligible.sort(key=lambda row: (abs(float(row["npe_end"]) - 150.0), row["runtime_s"]))
    best = eligible[0]
    return {
        "exists": True,
        "name": best["name"],
        "reflectivity": best["reflectivity"],
        "finish": best["finish"],
        "sigma_alpha": best["sigma_alpha_deg"],
        "air_gap_mm": best["air_gap_mm"],
        "npe_end_center": best["npe_end"],
        "sigma_end_center_ps": best["sigma_end_ps"],
        "runtime_100ev_s": None,
    }


def time_window_rows(rows: list[dict[str, object]], out_dir: pathlib.Path) -> list[dict[str, object]]:
    candidates = [row for row in rows if row["stage"] in {"reflectivity", "roughness"}]
    candidates.sort(key=lambda row: float(row["npe_end"]))
    selected_names = ["R095_polished"]
    selected_names.extend(row["name"] for row in candidates[:2] if row["name"] not in selected_names)
    selected = [row for row in rows if row["name"] in selected_names]
    output = []
    for row in selected:
        arrays = read_root(pathlib.Path(str(row["root_path"])))
        event_id = arrays["event_id"].astype(np.int64)
        global_id = arrays["global_id"].astype(np.int64)
        time_ns = arrays["time_ns"].astype(float)
        end = (global_id >= 0) & (global_id < 16)
        event_id = event_id[end]
        time_ns = time_ns[end]
        counts = {window: [] for window in WINDOWS_NS}
        total = []
        for event in range(int(row["events"])):
            times = time_ns[event_id == event]
            total.append(int(times.size))
            if times.size == 0:
                for window in WINDOWS_NS:
                    counts[window].append(0)
                continue
            t0 = float(np.min(times))
            for window in WINDOWS_NS:
                counts[window].append(int(np.count_nonzero(times - t0 <= window)))
        output.append({
            "config": row["name"],
            "total_npe": float(np.mean(total)),
            "npe_5ns": float(np.mean(counts[5.0])),
            "npe_10ns": float(np.mean(counts[10.0])),
            "npe_20ns": float(np.mean(counts[20.0])),
            "npe_50ns": float(np.mean(counts[50.0])),
            "t0_reference": "first END PE per event",
        })
    if output:
        png = out_dir / "png"
        png.mkdir(parents=True, exist_ok=True)
        xs = list(WINDOWS_NS)
        plt.figure(figsize=(6, 4))
        for row in output:
            ys = [row[f"npe_{int(x)}ns"] for x in xs]
            plt.plot(xs, ys, marker="o", label=str(row["config"]))
        plt.xlabel("time window from first END PE [ns]")
        plt.ylabel("mean Npe")
        plt.legend()
        plt.tight_layout()
        plt.savefig(png / "npe_vs_time_window.png", dpi=160)
        plt.close()
    return output


def make_plots(rows: list[dict[str, object]], out_dir: pathlib.Path) -> None:
    png = out_dir / "png"
    png.mkdir(parents=True, exist_ok=True)
    refl = [row for row in rows if row["stage"] == "reflectivity"]
    if refl:
        plt.figure(figsize=(6, 4))
        plt.plot([row["reflectivity"] for row in refl], [row["npe_end"] for row in refl], marker="o")
        plt.axhspan(20, 300, color="tab:green", alpha=0.15)
        plt.xlabel("reflectivity")
        plt.ylabel("mean END Npe")
        plt.tight_layout()
        plt.savefig(png / "reflectivity_bracket.png", dpi=160)
        plt.close()
    gap = [row for row in rows if row["stage"] == "gap"]
    if gap:
        plt.figure(figsize=(6, 4))
        plt.plot([row["air_gap_mm"] for row in gap], [row["npe_end"] for row in gap], marker="o")
        plt.xlabel("air gap [mm]")
        plt.ylabel("mean END Npe")
        plt.tight_layout()
        plt.savefig(png / "gap_bracket.png", dpi=160)
        plt.close()
    rough = [row for row in rows if row["stage"] == "roughness"]
    if rough:
        labels = [f"{row['finish']} s={row['sigma_alpha_deg']}" for row in rough]
        plt.figure(figsize=(7, 4))
        plt.bar(labels, [row["npe_end"] for row in rough])
        plt.xticks(rotation=25, ha="right")
        plt.ylabel("mean END Npe")
        plt.tight_layout()
        plt.savefig(png / "roughness_bracket.png", dpi=160)
        plt.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("out/EXEC_25"))
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--skip-runs", action="store_true")
    args = parser.parse_args()

    repo = args.repo.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    if not args.skip_runs:
        for config in pde_configs():
            print(f"[EXEC_25] running {config.name}", flush=True)
            rows.append(run_config(config, repo, out_dir, args.threads))
        write_csv(out_dir / "pde_semantics_test.csv", rows)
        passed, semantics = pde_pass(rows)
        write_json(out_dir / "pde_semantics_summary.json", semantics)
        if not passed:
            results = {
                "exec_id": "EXEC_25",
                "verdict": "ABORTADO-EN-T1",
                "pde_semantics": semantics,
                "exec23_reference": {
                    "npe_end_center": 2035.26,
                    "collection_eff_raw": 0.097205,
                },
                "best_candidate": candidate_from([]),
                "next_step": "Fix or re-audit SiPM EFFICIENCY semantics before optical bracketing.",
            }
            write_json(out_dir / "results_exec25.json", results)
            print(json.dumps(clean(results), indent=2, allow_nan=False))
            return 2
        bracket_rows = []
        for config in bracket_configs():
            print(f"[EXEC_25] running {config.name}", flush=True)
            bracket_rows.append(run_config(config, repo, out_dir, args.threads))
        rows.extend(bracket_rows)
    else:
        raise RuntimeError("--skip-runs is not implemented for this one-shot runner")

    pde_rows = [row for row in rows if row["stage"] == "pde_semantics"]
    bracket_rows = [row for row in rows if row["stage"] != "pde_semantics"]
    passed, semantics = pde_pass(pde_rows)
    write_csv(out_dir / "pde_semantics_test.csv", pde_rows)
    write_json(out_dir / "pde_semantics_summary.json", semantics)
    write_csv(out_dir / "bracket_results.csv", bracket_rows)
    write_json(out_dir / "bracket_results.json", bracket_rows)
    windows = time_window_rows(rows, out_dir)
    write_csv(out_dir / "time_window_diagnostic.csv", windows)
    make_plots(rows, out_dir)

    best = candidate_from(bracket_rows)
    flags = []
    if not best["exists"]:
        flags.append("no_candidate_in_20_to_300_pe_range")
    if any(row["speed_gate"] != "PASS" for row in rows):
        flags.append("speed_gate_failure")
    if any(row["missing_rindex"] != 0 for row in rows):
        flags.append("missing_rindex")
    if any(float(row["npe_end"]) > 300.0 for row in bracket_rows):
        flags.append("high_collection_configs_remain")
    if any(row["name"] == "R000_polished" and float(row["npe_end"]) > 300.0 for row in bracket_rows):
        flags.append("reflectivity_not_primary_driver")
    if all(row["stage"] != "gap" or float(row["npe_end"]) > 300.0 for row in bracket_rows):
        flags.append("gap_thickness_not_primary_driver")
    if windows and min(float(row["npe_5ns"]) for row in windows) > 300.0:
        flags.append("prompt_collection_excess")
    results = {
        "exec_id": "EXEC_25",
        "verdict": "COMPLETADO-CON-FLAGS" if flags else "COMPLETADO",
        "pde_semantics": semantics,
        "exec23_reference": {
            "npe_end_center": 2035.26,
            "collection_eff_raw": 0.097205,
        },
        "best_candidate": best,
        "n_configs_run": len(rows),
        "flags": flags,
        "next_step": (
            "EXEC_26 scan END-only with selected candidate"
            if best["exists"]
            else "Implement partial-contact or non-ideal wrapping model before scan."
        ),
    }
    write_json(out_dir / "results_exec25.json", results)
    print(json.dumps(clean(results), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
