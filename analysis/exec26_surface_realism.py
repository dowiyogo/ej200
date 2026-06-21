#!/usr/bin/env python3
"""EXEC_26 scintillator-air surface realism bracket."""

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


@dataclass(frozen=True)
class Config:
    name: str
    stage: str
    events: int
    x_mm: float = 0.0
    scint_air_finish: str = "polished"
    scint_air_sigma_rad: float = 0.0
    reflector_r: float = 0.0
    surface_loss: float = 0.0
    mylar_finish: str = "polished"
    mylar_sigma_deg: float = 0.1
    air_gap_mm: float = 0.10


def clean(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    return value


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


def read_csv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def parse_last_int(pattern: str, text: str) -> int | None:
    matches = re.findall(pattern, text)
    return int(matches[-1]) if matches else None


def macro_text(config: Config) -> str:
    seed1 = 910000 + zlib.crc32(config.name.encode("utf-8")) % 100000
    seed2 = 1110000 + zlib.crc32(config.name[::-1].encode("utf-8")) % 100000
    return f"""/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/det/readout EndTop
/ship/geom/topSurface mylar
/ship/geom/sipmEfficiencyMode nominal
/ship/geom/scintAir/finish {config.scint_air_finish}
/ship/geom/scintAir/sigmaAlpha {config.scint_air_sigma_rad} rad
/ship/geom/mylar/reflectivity {config.reflector_r}
/ship/geom/mylar/specularLobe 1.0
/ship/geom/mylar/sigmaAlpha {config.mylar_sigma_deg} deg
/ship/geom/mylar/finish {config.mylar_finish}
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
    work.mkdir(parents=True, exist_ok=True)
    audit.mkdir(parents=True, exist_ok=True)
    sslg4_link = work / "sslg4"
    if not sslg4_link.exists():
        sslg4_link.symlink_to(repo / "build" / "sslg4")
    (work / "run.mac").write_text(macro_text(config))
    env = os.environ.copy()
    env["SHIP_OPTICAL_AUDIT_DIR"] = str(audit)
    env["SHIP_OPTICAL_AUDIT_TAG"] = config.name
    env["EXEC26_SCINT_AIR_SURFACE_LOSS"] = str(config.surface_loss)
    log = raw_dir / "run.log"
    start = time.monotonic()
    with log.open("w") as handle:
        proc = subprocess.run(
            [str(sim), "-t", str(threads), "-m", "run.mac"],
            cwd=work,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    runtime_s = time.monotonic() - start
    if proc.returncode != 0:
        raise RuntimeError(f"{config.name}: simulation failed; see {log}")
    root_src = work / "photon_hits_run000.root"
    root_final = raw_dir / "photon_hits.root"
    if root_src.exists():
        if root_final.exists():
            root_final.unlink()
        shutil.move(str(root_src), str(root_final))
    return analyze_config(config, root_final, log, audit, runtime_s)


def root_arrays(path: pathlib.Path) -> dict[str, np.ndarray]:
    if not path.exists() or path.stat().st_size == 0:
        return {
            "event_id": np.array([], dtype=np.int64),
            "global_id": np.array([], dtype=np.int64),
            "time_ns": np.array([], dtype=float),
        }
    with uproot.open(path) as handle:
        return handle["sipm_hits"].arrays(["event_id", "global_id", "time_ns"], library="np")


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
    arr = np.asarray(values)
    median = float(np.median(arr))
    mad = 1.4826 * float(np.median(np.abs(arr - median)))
    core = arr[np.abs(arr - median) <= 2.0 * max(mad, 1.0e-9)]
    if core.size < 2:
        return None
    return 1000.0 * float(np.std(core, ddof=1))


def boundary_statuses(audit_dir: pathlib.Path, run_tag: str) -> dict[str, int]:
    counts = {
        "TotalInternalReflection": 0,
        "FresnelReflection": 0,
        "FresnelRefraction": 0,
        "Absorption": 0,
        "NoRINDEX": 0,
        "SurfaceLossKilled": 0,
    }
    for row in read_csv(audit_dir / "optical_boundary_counts.csv"):
        if row.get("run_tag") != run_tag:
            continue
        boundary = row.get("boundary", "")
        status = row.get("status", "")
        count = int(row.get("count", "0"))
        if status in counts:
            counts[status] += count
        if boundary == "BarLV->ScintAirSurfaceLoss" and status == "Killed":
            counts["SurfaceLossKilled"] += count
    return counts


def analyze_config(
    config: Config,
    root_path: pathlib.Path,
    log_path: pathlib.Path,
    audit_dir: pathlib.Path,
    runtime_s: float,
) -> dict[str, object]:
    text = log_path.read_text(errors="replace")
    arrays = root_arrays(root_path)
    event_id = arrays["event_id"].astype(np.int64)
    global_id = arrays["global_id"].astype(np.int64)
    time_ns = arrays["time_ns"].astype(float)
    end_mask = (global_id >= 0) & (global_id < 16)
    left_counts = np.bincount(event_id[end_mask & (global_id < 8)], minlength=config.events)
    right_counts = np.bincount(event_id[end_mask & (global_id >= 8)], minlength=config.events)
    generated = parse_last_int(r"Scint photons generated:\s*(\d+)", text)
    entering = parse_last_int(r"Bar -> SiPM \(entering\)\s*:\s*(\d+)", text)
    speed_rows = read_csv(audit_dir / "speed_audit_steps.csv")
    total_steps = sum(int(row["n_steps"]) for row in speed_rows)
    missing_rindex = sum(int(row["missing_rindex"]) for row in speed_rows)
    max_super = max((float(row["superluminal_fraction"]) for row in speed_rows), default=0.0)
    recorded = int(end_mask.sum())
    paired = int(np.count_nonzero((left_counts > 0) & (right_counts > 0)))
    times = time_ns[end_mask]
    statuses = boundary_statuses(audit_dir, config.name)
    row = {
        "name": config.name,
        "stage": config.stage,
        "events": config.events,
        "x_mm": config.x_mm,
        "scint_air_finish": config.scint_air_finish,
        "scint_air_sigma_rad": config.scint_air_sigma_rad,
        "reflector_R": config.reflector_r,
        "surface_loss": config.surface_loss,
        "runtime_s": runtime_s,
        "n_gen_total": generated,
        "n_gen_per_event": generated / config.events if generated else None,
        "edep_mev_per_event": generated / config.events / LIGHT_YIELD_PER_MEV if generated else None,
        "bar_to_sipm_entering_total": entering,
        "bar_to_sipm_entering_per_event": entering / config.events if entering else None,
        "npe_end": recorded / config.events,
        "npe_l": float(np.mean(left_counts)),
        "npe_r": float(np.mean(right_counts)),
        "recorded_over_generated": recorded / generated if generated else None,
        "recorded_over_entering": recorded / entering if entering else None,
        "paired_fraction": paired / config.events,
        "sigma_end_ps": sigma_end_ps(event_id[end_mask], global_id[end_mask], times, config.events),
        "detected_time_p95_ns": float(np.percentile(times, 95)) if times.size else None,
        "detected_time_max_ns": float(np.max(times)) if times.size else None,
        "steps_per_generated_photon": total_steps / generated if generated else None,
        "speed_gate": "PASS" if missing_rindex == 0 and max_super < 0.001 else "FAIL",
        "missing_rindex": missing_rindex,
        "max_superluminal_fraction": max_super,
        "root_path": str(root_path),
        "log_path": str(log_path),
        "audit_dir": str(audit_dir),
    }
    row.update({f"boundary_{key}": value for key, value in statuses.items()})
    return row


def roughness_configs() -> list[Config]:
    return [
        Config("rough_R0_polished_s0p00", "roughness_R0", 20,
               scint_air_finish="polished", scint_air_sigma_rad=0.00, reflector_r=0.0),
        Config("rough_R0_ground_s0p05", "roughness_R0", 20,
               scint_air_finish="ground", scint_air_sigma_rad=0.05, reflector_r=0.0),
        Config("rough_R0_ground_s0p10", "roughness_R0", 20,
               scint_air_finish="ground", scint_air_sigma_rad=0.10, reflector_r=0.0),
        Config("rough_R0_ground_s0p20", "roughness_R0", 20,
               scint_air_finish="ground", scint_air_sigma_rad=0.20, reflector_r=0.0),
        Config("rough_R0_ground_s0p30", "roughness_R0", 20,
               scint_air_finish="ground", scint_air_sigma_rad=0.30, reflector_r=0.0),
        Config("rough_R0_ground_s0p50", "roughness_R0", 20,
               scint_air_finish="ground", scint_air_sigma_rad=0.50, reflector_r=0.0),
    ]


def surface_loss_configs(base: dict[str, object], reflector_r: float) -> list[Config]:
    finish = str(base["scint_air_finish"])
    sigma = float(base["scint_air_sigma_rad"])
    return [
        Config(f"loss_R{int(reflector_r * 100):03d}_{str(loss).replace('.', 'p')}",
               "surface_loss", 20,
               scint_air_finish=finish, scint_air_sigma_rad=sigma,
               reflector_r=reflector_r, surface_loss=loss)
        for loss in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9)
    ]


def in_candidate_range(row: dict[str, object]) -> bool:
    return (
        row["speed_gate"] == "PASS"
        and int(row["missing_rindex"]) == 0
        and 20.0 <= float(row["npe_end"]) <= 300.0
        and float(row["paired_fraction"]) >= 0.8
        and (row["sigma_end_ps"] is None or float(row["sigma_end_ps"]) < 250.0)
    )


def choose_candidate(rows: list[dict[str, object]]) -> dict[str, object] | None:
    eligible = [row for row in rows if in_candidate_range(row)]
    if not eligible:
        return None
    eligible.sort(key=lambda row: (float(row["surface_loss"]), abs(float(row["npe_end"]) - 150.0)))
    return eligible[0]


def make_plots(out_dir: pathlib.Path, rough: list[dict[str, object]],
               loss: list[dict[str, object]], refl: list[dict[str, object]]) -> None:
    png = out_dir / "png"
    png.mkdir(parents=True, exist_ok=True)
    if rough:
        plt.figure(figsize=(6, 4))
        xs = [row["scint_air_sigma_rad"] for row in rough]
        ys = [row["npe_end"] for row in rough]
        labels = [row["scint_air_finish"] for row in rough]
        plt.plot(xs, ys, marker="o")
        for x, y, label in zip(xs, ys, labels):
            plt.text(x, y, label, fontsize=7)
        plt.axhspan(20, 300, color="tab:green", alpha=0.15)
        plt.xlabel("scint-air sigma alpha [rad]")
        plt.ylabel("mean END Npe")
        plt.tight_layout()
        plt.savefig(png / "scint_air_roughness_R0.png", dpi=160)
        plt.close()
    if loss:
        plt.figure(figsize=(6, 4))
        plt.plot([row["surface_loss"] for row in loss],
                 [row["npe_end"] for row in loss], marker="o")
        plt.axhspan(20, 300, color="tab:green", alpha=0.15)
        plt.xlabel("surface loss probability")
        plt.ylabel("mean END Npe")
        plt.tight_layout()
        plt.savefig(png / "surface_loss_bracket.png", dpi=160)
        plt.close()
    if refl:
        plt.figure(figsize=(6, 4))
        plt.plot([row["reflector_R"] for row in refl],
                 [row["npe_end"] for row in refl], marker="o")
        plt.axhspan(20, 300, color="tab:green", alpha=0.15)
        plt.xlabel("reflector R")
        plt.ylabel("mean END Npe")
        plt.tight_layout()
        plt.savefig(png / "scint_air_plus_reflector.png", dpi=160)
        plt.close()


def make_candidate_time_pdf(out_dir: pathlib.Path, validation_rows: list[dict[str, object]]) -> None:
    x0_rows = [row for row in validation_rows if float(row["x_mm"]) == 0.0]
    if not x0_rows:
        return
    arrays = root_arrays(pathlib.Path(str(x0_rows[0]["root_path"])))
    global_id = arrays["global_id"].astype(np.int64)
    time_ns = arrays["time_ns"].astype(float)
    times = time_ns[(global_id >= 0) & (global_id < 16)]
    if times.size == 0:
        return
    png = out_dir / "png"
    png.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6, 4))
    plt.hist(times, bins=80, histtype="step", density=True)
    plt.xlabel("END photon detection time [ns]")
    plt.ylabel("density")
    plt.tight_layout()
    plt.savefig(png / "candidate_time_pdf_x0.png", dpi=160)
    plt.close()


def candidate_json(row: dict[str, object] | None) -> dict[str, object]:
    if row is None:
        return {
            "exists": False,
            "name": None,
            "finish": None,
            "sigma_alpha": None,
            "reflector_R": None,
            "surface_loss": None,
            "npe_center": None,
            "sigma_center_ps": None,
            "paired_fraction": None,
            "speed_gate_pass": None,
        }
    return {
        "exists": True,
        "name": row["name"],
        "finish": row["scint_air_finish"],
        "sigma_alpha": row["scint_air_sigma_rad"],
        "reflector_R": row["reflector_R"],
        "surface_loss": row["surface_loss"],
        "npe_center": row["npe_end"],
        "sigma_center_ps": row["sigma_end_ps"],
        "paired_fraction": row["paired_fraction"],
        "speed_gate_pass": row["speed_gate"] == "PASS",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("out/EXEC_26"))
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()

    repo = args.repo.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rough_rows = []
    for config in roughness_configs():
        print(f"[EXEC_26] running {config.name}", flush=True)
        rough_rows.append(run_config(config, repo, out_dir, args.threads))
    write_csv(out_dir / "scint_air_roughness_R0.csv", rough_rows)
    write_json(out_dir / "scint_air_roughness_R0.json", rough_rows)

    rough_candidate = choose_candidate(rough_rows)
    candidate = rough_candidate
    reflector_rows: list[dict[str, object]] = []
    loss_rows: list[dict[str, object]] = []
    validation_rows: list[dict[str, object]] = []

    if candidate is not None:
        for r in (0.85, 0.95):
            config = Config(
                f"{candidate['name']}_R{int(r * 100):03d}",
                "plus_reflector",
                20,
                scint_air_finish=str(candidate["scint_air_finish"]),
                scint_air_sigma_rad=float(candidate["scint_air_sigma_rad"]),
                reflector_r=r,
                surface_loss=float(candidate["surface_loss"]),
            )
            print(f"[EXEC_26] running {config.name}", flush=True)
            reflector_rows.append(run_config(config, repo, out_dir, args.threads))
        candidate = choose_candidate(reflector_rows)

    if candidate is None:
        loss_base = rough_candidate or min(rough_rows, key=lambda row: float(row["npe_end"]))
        for config in surface_loss_configs(loss_base, reflector_r=0.95):
            print(f"[EXEC_26] running {config.name}", flush=True)
            loss_rows.append(run_config(config, repo, out_dir, args.threads))
        candidate = choose_candidate(loss_rows)
        if candidate is not None:
            for r in (0.85, 0.95):
                config = Config(
                    f"{candidate['name']}_R{int(r * 100):03d}",
                    "plus_reflector_loss",
                    20,
                    scint_air_finish=str(candidate["scint_air_finish"]),
                    scint_air_sigma_rad=float(candidate["scint_air_sigma_rad"]),
                    reflector_r=r,
                    surface_loss=float(candidate["surface_loss"]),
                )
                print(f"[EXEC_26] running {config.name}", flush=True)
                reflector_rows.append(run_config(config, repo, out_dir, args.threads))
            candidate = choose_candidate(reflector_rows) or candidate

    if reflector_rows:
        write_csv(out_dir / "scint_air_plus_reflector.csv", reflector_rows)
        write_json(out_dir / "scint_air_plus_reflector.json", reflector_rows)
    if loss_rows:
        write_csv(out_dir / "surface_loss_bracket.csv", loss_rows)
        write_json(out_dir / "surface_loss_bracket.json", loss_rows)

    if candidate is not None:
        write_csv(out_dir / "candidate_speed_audit.csv", [candidate])
        write_json(out_dir / "candidate_speed_audit.json", candidate)
        validation_config = Config(
            f"{candidate['name']}_val100_x0",
            "candidate_validation",
            100,
            scint_air_finish=str(candidate["scint_air_finish"]),
            scint_air_sigma_rad=float(candidate["scint_air_sigma_rad"]),
            reflector_r=float(candidate["reflector_R"]),
            surface_loss=float(candidate["surface_loss"]),
        )
        print(f"[EXEC_26] running {validation_config.name}", flush=True)
        val0 = run_config(validation_config, repo, out_dir, args.threads)
        validation_rows.append(val0)
        if in_candidate_range(val0):
            for x in (-690.0, 690.0):
                edge_config = Config(
                    f"{candidate['name']}_val50_x{int(x)}",
                    "candidate_validation",
                    50,
                    x_mm=x,
                    scint_air_finish=str(candidate["scint_air_finish"]),
                    scint_air_sigma_rad=float(candidate["scint_air_sigma_rad"]),
                    reflector_r=float(candidate["reflector_R"]),
                    surface_loss=float(candidate["surface_loss"]),
                )
                print(f"[EXEC_26] running {edge_config.name}", flush=True)
                validation_rows.append(run_config(edge_config, repo, out_dir, args.threads))
        write_csv(out_dir / "candidate_validation.csv", validation_rows)
        write_json(out_dir / "candidate_validation.json", validation_rows)
        make_candidate_time_pdf(out_dir, validation_rows)

    make_plots(out_dir, rough_rows, loss_rows, reflector_rows)

    all_rows = rough_rows + loss_rows + reflector_rows + validation_rows
    flags = []
    if candidate is None:
        flags.append("no_candidate")
    if any(row["speed_gate"] != "PASS" for row in all_rows):
        flags.append("speed_gate_failure")
    if any(int(row["missing_rindex"]) != 0 for row in all_rows):
        flags.append("missing_rindex")
    if loss_rows:
        flags.append("surface_loss_model_used")
    if min(float(row["npe_end"]) for row in rough_rows) > 300.0:
        flags.append("roughness_alone_high_collection")
    if validation_rows and any(
        float(row["npe_end"]) > 300.0
        for row in validation_rows
        if float(row["x_mm"]) != 0.0
    ):
        flags.append("edge_validation_high_collection")

    results = {
        "exec_id": "EXEC_26",
        "verdict": "COMPLETADO-CON-FLAGS" if flags else "COMPLETADO",
        "exec25_reference": {
            "npe_exec23": 2035.26,
            "lowest_exec25_npe": 1781.95,
            "dominant_driver": "scintillator-air TIR / effective prompt acceptance",
        },
        "roughness_bracket": {
            "ran": True,
            "lowest_npe": min(float(row["npe_end"]) for row in rough_rows),
            "best_sigma_alpha": min(rough_rows, key=lambda row: float(row["npe_end"]))["scint_air_sigma_rad"],
        },
        "surface_loss_model": {
            "used": bool(loss_rows),
            "loss_fraction": candidate["surface_loss"] if candidate is not None else None,
            "npe": candidate["npe_end"] if candidate is not None else None,
        },
        "candidate": candidate_json(candidate),
        "flags": flags,
        "next_step": (
            "EXEC_27 scan END-only 11 positions"
            if candidate is not None and "edge_validation_high_collection" not in flags
            else "Review edge response before scan or add position-dependent coupling/contact model"
        ),
    }
    write_json(out_dir / "results_exec26.json", results)
    print(json.dumps(clean(results), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
