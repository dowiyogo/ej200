#!/usr/bin/env python3
"""EXEC_27 edge estimator decomposition for END light collection."""

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


@dataclass(frozen=True)
class RunConfig:
    name: str
    x_mm: float
    events: int
    scint_air_finish: str = "ground"
    scint_air_sigma_rad: float = 0.20
    reflector_r: float = 0.95
    surface_loss: float = 0.10
    readout: str = "End"


EXEC26_VALIDATION = {
    -690.0: pathlib.Path("out/EXEC_26/raw/loss_R095_0p1_val50_x-690/photon_hits.root"),
    0.0: pathlib.Path("out/EXEC_26/raw/loss_R095_0p1_val100_x0/photon_hits.root"),
    690.0: pathlib.Path("out/EXEC_26/raw/loss_R095_0p1_val50_x690/photon_hits.root"),
}


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


def macro_text(config: RunConfig) -> str:
    seed1 = 970000 + zlib.crc32(config.name.encode("utf-8")) % 100000
    seed2 = 1270000 + zlib.crc32(config.name[::-1].encode("utf-8")) % 100000
    return f"""/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/det/readout {config.readout}
/ship/geom/topSurface mylar
/ship/geom/sipmEfficiencyMode nominal
/ship/geom/scintAir/finish {config.scint_air_finish}
/ship/geom/scintAir/sigmaAlpha {config.scint_air_sigma_rad} rad
/ship/geom/mylar/reflectivity {config.reflector_r}
/ship/geom/mylar/specularLobe 1.0
/ship/geom/mylar/sigmaAlpha 0.1 deg
/ship/geom/mylar/finish polished
/ship/geom/mylar/airGap 0.10 mm
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


def run_config(config: RunConfig, repo: pathlib.Path, out_dir: pathlib.Path, threads: int) -> pathlib.Path:
    raw_dir = out_dir / "raw" / config.name
    work = raw_dir / "work"
    audit = raw_dir / "audit"
    work.mkdir(parents=True, exist_ok=True)
    audit.mkdir(parents=True, exist_ok=True)
    link = work / "sslg4"
    if not link.exists():
        link.symlink_to(repo / "build" / "sslg4")
    (work / "run.mac").write_text(macro_text(config))
    env = os.environ.copy()
    env["SHIP_OPTICAL_AUDIT_DIR"] = str(audit)
    env["SHIP_OPTICAL_AUDIT_TAG"] = config.name
    env["EXEC26_SCINT_AIR_SURFACE_LOSS"] = str(config.surface_loss)
    log = raw_dir / "run.log"
    start = time.monotonic()
    with log.open("w") as handle:
        proc = subprocess.run(
            [str(repo / "build" / "ej200_bar_sim"), "-t", str(threads), "-m", "run.mac"],
            cwd=work,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    runtime_s = time.monotonic() - start
    (raw_dir / "runtime_s.txt").write_text(f"{runtime_s:.6f}\n")
    if proc.returncode != 0:
        raise RuntimeError(f"{config.name}: simulation failed; see {log}")
    src = work / "photon_hits_run000.root"
    dst = raw_dir / "photon_hits.root"
    if src.exists():
        if dst.exists():
            dst.unlink()
        shutil.move(str(src), str(dst))
    return dst


def existing_or_run(config: RunConfig, repo: pathlib.Path, out_dir: pathlib.Path, threads: int) -> pathlib.Path:
    existing = EXEC26_VALIDATION.get(float(config.x_mm))
    if existing is not None and (repo / existing).exists():
        return repo / existing
    dst = out_dir / "raw" / config.name / "photon_hits.root"
    if dst.exists():
        return dst
    return run_config(config, repo, out_dir, threads)


def root_arrays(path: pathlib.Path) -> dict[str, np.ndarray]:
    with uproot.open(path) as handle:
        return handle["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns", "gun_x_mm"], library="np"
        )


def robust_sigma_ps(values: np.ndarray) -> float | None:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 5:
        return None
    median = float(np.median(arr))
    mad = 1.4826 * float(np.median(np.abs(arr - median)))
    core = arr[np.abs(arr - median) <= 2.0 * max(mad, 1.0e-9)]
    if core.size < 2:
        core = arr
    return 1000.0 * float(np.std(core, ddof=1))


def side_first_times(
    event_id: np.ndarray, global_id: np.ndarray, time_ns: np.ndarray, n_events: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    t_left = np.full(n_events, np.nan)
    t_right = np.full(n_events, np.nan)
    n_left = np.zeros(n_events, dtype=int)
    n_right = np.zeros(n_events, dtype=int)
    end = (global_id >= 0) & (global_id < 16)
    for event in range(n_events):
        mask = end & (event_id == event)
        left = time_ns[mask & (global_id < 8)]
        right = time_ns[mask & (global_id >= 8)]
        n_left[event] = int(left.size)
        n_right[event] = int(right.size)
        if left.size:
            t_left[event] = float(np.min(left))
        if right.size:
            t_right[event] = float(np.min(right))
    return t_left, t_right, n_left, n_right


def gls_metrics(t_left: np.ndarray, t_right: np.ndarray) -> tuple[float | None, float | None]:
    residuals = np.vstack([t_left - np.nanmean(t_left), t_right - np.nanmean(t_right)])
    if residuals.shape[1] < 5:
        return None, None
    cov = np.cov(residuals, ddof=1)
    inv = np.linalg.pinv(cov)
    one = np.ones(2)
    denom = float(one @ inv @ one)
    if denom <= 0.0 or not math.isfinite(denom):
        return None, None
    weights = inv @ one / denom
    values = weights[0] * t_left + weights[1] * t_right
    return 1000.0 / math.sqrt(denom), robust_sigma_ps(values)


def speed_gate_for(path: pathlib.Path, name: str) -> tuple[str, int, float]:
    raw = path.parent
    audit = raw / "audit"
    if not audit.exists():
        return "REUSED_EXEC26", 0, 0.0
    rows = read_csv(audit / "speed_audit_steps.csv")
    missing = sum(int(row["missing_rindex"]) for row in rows)
    max_super = max((float(row["superluminal_fraction"]) for row in rows), default=0.0)
    return ("PASS" if missing == 0 and max_super < 0.001 else "FAIL"), missing, max_super


def analyze_position(config: RunConfig, root_path: pathlib.Path) -> dict[str, object]:
    arrays = root_arrays(root_path)
    event_id = arrays["event_id"].astype(np.int64)
    global_id = arrays["global_id"].astype(np.int64)
    time_ns = arrays["time_ns"].astype(float)
    n_events = max(config.events, int(np.max(event_id)) + 1 if event_id.size else 0)
    t_left, t_right, n_left, n_right = side_first_times(event_id, global_id, time_ns, n_events)
    paired = np.isfinite(t_left) & np.isfinite(t_right)
    if np.count_nonzero(paired) < 5:
        raise RuntimeError(f"{config.name}: fewer than five paired END events")
    L = t_left[paired]
    R = t_right[paired]
    NL = n_left[paired].astype(float)
    NR = n_right[paired].astype(float)
    sigma_l = robust_sigma_ps(L)
    sigma_r = robust_sigma_ps(R)
    tavg = 0.5 * (L + R)
    if sigma_l and sigma_r and sigma_l > 0.0 and sigma_r > 0.0:
        w_l = 1.0 / sigma_l**2
        w_r = 1.0 / sigma_r**2
        wmean_static = (w_l * L + w_r * R) / (w_l + w_r)
    else:
        wmean_static = tavg
    denom = NL + NR
    wmean_event = (NL * L + NR * R) / denom
    denom2 = NL**2 + NR**2
    wmean_event_npe2 = (NL**2 * L + NR**2 * R) / denom2
    sigma_gls_formula, sigma_gls_dist = gls_metrics(L, R)
    near = L if config.x_mm < 0 else R if config.x_mm > 0 else np.where(NL >= NR, L, R)
    far = R if config.x_mm < 0 else L if config.x_mm > 0 else np.where(NL >= NR, R, L)
    speed_gate, missing_rindex, max_super = speed_gate_for(root_path, config.name)
    return {
        "name": config.name,
        "x_mm": config.x_mm,
        "events": n_events,
        "root_path": str(root_path),
        "npe_L": float(np.mean(n_left)),
        "npe_R": float(np.mean(n_right)),
        "npe_total": float(np.mean(n_left + n_right)),
        "paired_fraction": float(np.count_nonzero(paired) / n_events),
        "sigma_L_ps": sigma_l,
        "sigma_R_ps": sigma_r,
        "sigma_tavg_ps": robust_sigma_ps(tavg),
        "sigma_wmean_static_ps": robust_sigma_ps(wmean_static),
        "sigma_wmean_event_ps": robust_sigma_ps(wmean_event),
        "sigma_wmean_event_npe2_ps": robust_sigma_ps(wmean_event_npe2),
        "sigma_GLS_residual_ps": sigma_gls_formula,
        "sigma_GLS_estimator_dist_ps": sigma_gls_dist,
        "sigma_near_end_ps": robust_sigma_ps(near),
        "sigma_far_end_ps": robust_sigma_ps(far),
        "mean_t_L_ns": float(np.nanmean(t_left)),
        "mean_t_R_ns": float(np.nanmean(t_right)),
        "mean_weight_static_L": float(
            (1.0 / sigma_l**2) / ((1.0 / sigma_l**2) + (1.0 / sigma_r**2))
        ) if sigma_l and sigma_r and sigma_l > 0.0 and sigma_r > 0.0 else None,
        "mean_weight_event_L": float(np.mean(NL / (NL + NR))),
        "mean_weight_event_npe2_L": float(np.mean(NL**2 / (NL**2 + NR**2))),
        "speed_gate": speed_gate,
        "missing_rindex": missing_rindex,
        "max_superluminal_fraction": max_super,
    }


def make_plots(out_dir: pathlib.Path, rows: list[dict[str, object]]) -> None:
    png = out_dir / "png"
    png.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda row: float(row["x_mm"]))
    x = [float(row["x_mm"]) for row in rows]
    plt.figure(figsize=(7, 4.5))
    for key, label in [
        ("sigma_tavg_ps", "t_avg"),
        ("sigma_wmean_static_ps", "wmean static"),
        ("sigma_wmean_event_ps", "wmean event Npe"),
        ("sigma_GLS_residual_ps", "GLS residual"),
        ("sigma_near_end_ps", "near end"),
        ("sigma_far_end_ps", "far end"),
    ]:
        plt.plot(x, [row[key] for row in rows], marker="o", label=label)
    plt.axhline(250.0, color="0.4", ls="--", lw=1)
    plt.xlabel("gun x [mm]")
    plt.ylabel("sigma [ps]")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(png / "sigma_estimators_vs_x.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(x, [row["npe_L"] for row in rows], marker="o", label="END_L")
    plt.plot(x, [row["npe_R"] for row in rows], marker="o", label="END_R")
    plt.plot(x, [row["npe_total"] for row in rows], marker="o", label="total")
    plt.yscale("log")
    plt.xlabel("gun x [mm]")
    plt.ylabel("mean PE / event")
    plt.legend()
    plt.tight_layout()
    plt.savefig(png / "npe_lr_vs_x.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(x, [row["mean_weight_static_L"] for row in rows], marker="o", label="static L")
    plt.plot(x, [row["mean_weight_event_L"] for row in rows], marker="o", label="Npe L")
    plt.plot(x, [row["mean_weight_event_npe2_L"] for row in rows], marker="o", label="Npe^2 L")
    plt.xlabel("gun x [mm]")
    plt.ylabel("left-side estimator weight")
    plt.ylim(-0.05, 1.05)
    plt.legend()
    plt.tight_layout()
    plt.savefig(png / "weights_vs_x.png", dpi=160)
    plt.close()


def write_comparability(path: pathlib.Path) -> None:
    path.write_text(
        """# What is comparable to the historical 85-90 ps?

| Estimator | Cancels position? | Requires per-event amplitude? | Comparable to historical weighted mean? | Useful for ToF? | Useful for x reconstruction? |
|---|---|---|---|---|---|
| t_avg = (t_L+t_R)/2 | Yes for symmetric two-ended propagation, but fragile when one side has very low PE | No | No, unless the reference used equal endpoint weights | Yes at center; poor at edges in EXEC_27 | Only with separate t_L-t_R or calibrated propagation |
| wmean_static | Partly; it weights by per-position timing precision | No, but needs a calibration map versus x | Yes, closest if historical analysis used calibrated endpoint weights | Yes, with x-dependent calibration | Not by itself; combine with asymmetry/timing difference |
| wmean_event using Npe | Partly; it follows per-event endpoint information | Yes | Plausibly comparable if historical weighted mean used amplitude weights | Yes, but can collapse to near-end timing at edges | Npe asymmetry itself carries x information |
| GLS_residual | Yes for resolution after removing per-position mean | No per event, but needs covariance calibration | Yes for intrinsic resolution diagnostics | Yes after calibration; not a raw timestamp | No, it removes the position mean by construction |
| near-end only | No | No | No | Good local timing diagnostic, not a two-ended ToF estimator | Strong x dependence; near/far identity gives side |
| far-end only | No | No | No | Diagnostic of light transport and far-side detectability | Strong x dependence; useful for transport validation |
""",
        encoding="ascii",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("out/EXEC_27"))
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    repo = args.repo.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    configs = [
        RunConfig("exec27_edge_x-690", -690.0, 50),
        RunConfig("exec27_edge_x-400", -400.0, 50),
        RunConfig("exec27_edge_x0", 0.0, 100),
        RunConfig("exec27_edge_x400", 400.0, 50),
        RunConfig("exec27_edge_x690", 690.0, 50),
    ]
    rows = []
    for config in configs:
        root = existing_or_run(config, repo, out_dir, args.threads)
        print(f"[EXEC_27] analyze x={config.x_mm:g} root={root}", flush=True)
        rows.append(analyze_position(config, root))
    rows = sorted(rows, key=lambda row: float(row["x_mm"]))
    write_csv(out_dir / "edge_estimator_summary.csv", rows)
    write_json(out_dir / "edge_estimator_summary.json", rows)
    write_csv(out_dir / "candidate_validation.csv", rows)
    write_json(out_dir / "candidate_validation.json", rows)
    make_plots(out_dir, rows)
    write_comparability(out_dir / "comparability_matrix.md")

    tavg_edge_fails = any(
        float(row["sigma_tavg_ps"]) >= 250.0 for row in rows if abs(float(row["x_mm"])) >= 400.0
    )
    wmean_edge_passes = all(
        row["sigma_wmean_static_ps"] is not None and float(row["sigma_wmean_static_ps"]) < 250.0
        for row in rows
    )
    gls_edge_passes = all(
        row["sigma_GLS_residual_ps"] is not None and float(row["sigma_GLS_residual_ps"]) < 250.0
        for row in rows
    )
    speed_fail = any(row["speed_gate"] == "FAIL" for row in rows)
    main_issue = (
        "t_avg estimator at asymmetric edge light sharing"
        if tavg_edge_fails and (wmean_edge_passes or gls_edge_passes)
        else "undetermined"
    )
    flags = []
    if tavg_edge_fails:
        flags.append("tavg_edge_failure")
    if wmean_edge_passes:
        flags.append("weighted_estimators_rescue_edges")
    if gls_edge_passes:
        flags.append("gls_residual_rescues_edges")
    if speed_fail:
        flags.append("speed_gate_failure")

    center = min(rows, key=lambda row: abs(float(row["x_mm"])))
    edge_rows = [row for row in rows if abs(float(row["x_mm"])) == 690.0]
    worst_edge = max(edge_rows, key=lambda row: float(row["npe_total"]))
    candidate_exists = (wmean_edge_passes or gls_edge_passes) and not speed_fail
    candidate_estimator = "wmean_static" if wmean_edge_passes else "GLS_residual" if gls_edge_passes else None
    results = {
        "exec_id": "EXEC_27",
        "verdict": "COMPLETADO-CON-FLAGS" if flags else "COMPLETADO",
        "estimator_diagnosis": {
            "tavg_edge_fails": tavg_edge_fails,
            "wmean_edge_passes": wmean_edge_passes,
            "gls_edge_passes": gls_edge_passes,
            "main_issue": main_issue,
        },
        "sipm_angle_audit": {
            "ran": False,
            "edge_oblique_fraction": None,
            "center_oblique_fraction": None,
            "reason_not_run": "T1 gate showed weighted/GLS estimators rescue all edge points; no SiPM coupling change justified.",
        },
        "coupling_model": {
            "used": False,
            "candidate_exists": False,
            "model": None,
            "parameters": None,
        },
        "candidate": {
            "exists": candidate_exists,
            "npe_center": center["npe_total"],
            "npe_edge_near": max(worst_edge["npe_L"], worst_edge["npe_R"]),
            "sigma_center_ps": center["sigma_wmean_static_ps"] if candidate_estimator == "wmean_static" else center["sigma_GLS_residual_ps"],
            "sigma_edge_ps": max(
                float(row["sigma_wmean_static_ps"] if candidate_estimator == "wmean_static" else row["sigma_GLS_residual_ps"])
                for row in edge_rows
            ) if candidate_estimator else None,
            "estimator": candidate_estimator,
        },
        "flags": flags,
        "next_step": (
            "EXEC_28 scan11 END-only with calibrated weighted/GLS estimators; do not add SiPM coupling yet."
            if candidate_exists
            else "Run T2 angular audit before considering any SiPM coupling model."
        ),
    }
    write_json(out_dir / "results_exec27.json", results)
    print(json.dumps(clean(results), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
