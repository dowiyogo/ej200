#!/usr/bin/env python3
"""EXEC_28 END-only scan11 with weighted/GLS estimators."""

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


CONFIG = {
    "scint_air_finish": "ground",
    "sigma_alpha_rad": 0.20,
    "reflector_R": 0.95,
    "surface_loss": 0.10,
    "air_gap_mm": 0.10,
    "pde_mode": "nominal",
    "scintillator": "OPSC-101",
    "readout": "End",
    "top_surface": "mylar",
}

POSITIONS = [-690.0, -552.0, -414.0, -276.0, -138.0, 0.0, 138.0, 276.0, 414.0, 552.0, 690.0]
RUN_ORDER = [0.0, -690.0, 690.0, -414.0, 414.0, -552.0, 552.0, -276.0, 276.0, -138.0, 138.0]
REAL_SINGLE_FILT_PS = 250.0
REAL_BOTHENDS_WMEAN_PS = (85.0, 90.0)
REAL_SPATIAL_CM = 1.3


@dataclass(frozen=True)
class RunConfig:
    name: str
    x_mm: float
    events: int


def clean(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, np.floating):
        return clean(float(value))
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    return value


def write_csv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(clean(rows))


def write_json(path: pathlib.Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(payload), indent=2, allow_nan=False) + "\n", encoding="ascii")


def read_csv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def parse_last_int(pattern: str, text: str) -> int | None:
    matches = re.findall(pattern, text)
    return int(matches[-1]) if matches else None


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


def rms(values: np.ndarray) -> float | None:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 2:
        return None
    return float(np.std(arr, ddof=1))


def macro_text(config: RunConfig) -> str:
    seed1 = 980000 + zlib.crc32(config.name.encode("utf-8")) % 100000
    seed2 = 1380000 + zlib.crc32(config.name[::-1].encode("utf-8")) % 100000
    return f"""/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/det/readout {CONFIG["readout"]}
/ship/geom/topSurface {CONFIG["top_surface"]}
/ship/geom/sipmEfficiencyMode {CONFIG["pde_mode"]}
/ship/geom/scintAir/finish {CONFIG["scint_air_finish"]}
/ship/geom/scintAir/sigmaAlpha {CONFIG["sigma_alpha_rad"]} rad
/ship/geom/mylar/reflectivity {CONFIG["reflector_R"]}
/ship/geom/mylar/specularLobe 1.0
/ship/geom/mylar/sigmaAlpha 0.1 deg
/ship/geom/mylar/finish polished
/ship/geom/mylar/airGap {CONFIG["air_gap_mm"]} mm
/det/scintillator {CONFIG["scintillator"]}
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


def run_config(
    config: RunConfig,
    repo: pathlib.Path,
    out_dir: pathlib.Path,
    threads: int,
    timeout_s: float | None,
) -> dict[str, object]:
    raw_dir = out_dir / "raw" / config.name
    work = raw_dir / "work"
    audit = raw_dir / "audit"
    work.mkdir(parents=True, exist_ok=True)
    audit.mkdir(parents=True, exist_ok=True)
    link = work / "sslg4"
    if not link.exists():
        link.symlink_to(repo / "build" / "sslg4")
    (work / "run.mac").write_text(macro_text(config), encoding="ascii")
    env = os.environ.copy()
    env["SHIP_OPTICAL_AUDIT_DIR"] = str(audit)
    env["SHIP_OPTICAL_AUDIT_TAG"] = config.name
    env["EXEC26_SCINT_AIR_SURFACE_LOSS"] = str(CONFIG["surface_loss"])
    log = raw_dir / "run.log"
    start = time.monotonic()
    timed_out = False
    returncode: int | None = None
    with log.open("w") as handle:
        try:
            proc = subprocess.run(
                [str(repo / "build" / "ej200_bar_sim"), "-t", str(threads), "-m", "run.mac"],
                cwd=work,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                timeout=timeout_s,
                check=False,
            )
            returncode = proc.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
            returncode = 124
    runtime_s = time.monotonic() - start
    (raw_dir / "runtime_s.txt").write_text(f"{runtime_s:.6f}\n", encoding="ascii")
    root_src = work / "photon_hits_run000.root"
    root_dst = raw_dir / "photon_hits.root"
    if root_src.exists():
        if root_dst.exists():
            root_dst.unlink()
        shutil.move(str(root_src), str(root_dst))
    return {
        "name": config.name,
        "x_mm_nominal": config.x_mm,
        "n_events_requested": config.events,
        "runtime_s": runtime_s,
        "root_path": str(root_dst),
        "log_path": str(log),
        "audit_dir": str(audit),
        "returncode": returncode,
        "timed_out": timed_out,
        "completed": (returncode == 0 and root_dst.exists()),
    }


def root_arrays(path: pathlib.Path) -> dict[str, np.ndarray]:
    with uproot.open(path) as handle:
        return handle["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns", "gun_x_mm"], library="np"
        )


def event_series(root_path: pathlib.Path) -> dict[str, np.ndarray]:
    arrays = root_arrays(root_path)
    event_id = arrays["event_id"].astype(np.int64)
    global_id = arrays["global_id"].astype(np.int64)
    time_ns = arrays["time_ns"].astype(float)
    gun_x = arrays["gun_x_mm"].astype(float)
    if event_id.size == 0:
        raise RuntimeError(f"{root_path}: no hits")
    n_events = int(np.max(event_id)) + 1
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
    return {
        "t_left": t_left,
        "t_right": t_right,
        "n_left": n_left,
        "n_right": n_right,
        "gun_x_median": np.array([float(np.median(gun_x)) if gun_x.size else math.nan]),
        "n_events": np.array([n_events]),
    }


def gls(t_left: np.ndarray, t_right: np.ndarray) -> tuple[float | None, float | None, float | None]:
    if t_left.size < 5:
        return None, None, None
    residuals = np.vstack([t_left - np.mean(t_left), t_right - np.mean(t_right)])
    cov = np.cov(residuals, ddof=1)
    rho = None
    if cov[0, 0] > 0.0 and cov[1, 1] > 0.0:
        rho = float(cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1]))
    inv = np.linalg.pinv(cov)
    one = np.ones(2)
    denom = float(one @ inv @ one)
    if denom <= 0.0 or not math.isfinite(denom):
        return None, None, rho
    weights = inv @ one / denom
    estimator = weights[0] * t_left + weights[1] * t_right
    return 1000.0 / math.sqrt(denom), robust_sigma_ps(estimator), rho


def speed_summary(audit_dir: pathlib.Path) -> dict[str, object]:
    rows = read_csv(audit_dir / "speed_audit_steps.csv")
    if not rows:
        return {"speed_gate": "MISSING", "missing_rindex": None, "max_ratio": None, "superluminal_fraction": None}
    missing = sum(int(row["missing_rindex"]) for row in rows)
    max_super = max((float(row["superluminal_fraction"]) for row in rows), default=0.0)
    max_ratio = max((float(row.get("max_ratio", row.get("max_v_over_c_over_n", "0")) or 0.0) for row in rows), default=0.0)
    return {
        "speed_gate": "PASS" if missing == 0 and max_super < 0.001 and max_ratio <= 1.001 else "FAIL",
        "missing_rindex": missing,
        "max_ratio": max_ratio,
        "superluminal_fraction": max_super,
    }


def analyze_completed(inventory: list[dict[str, object]], out_dir: pathlib.Path) -> dict[str, object]:
    npe_rows: list[dict[str, object]] = []
    light_rows: list[dict[str, object]] = []
    endpoint_rows: list[dict[str, object]] = []
    combined_rows: list[dict[str, object]] = []
    spatial_rows: list[dict[str, object]] = []
    series_by_x: dict[float, dict[str, np.ndarray]] = {}

    for item in inventory:
        if not item.get("completed"):
            continue
        root_path = pathlib.Path(str(item["root_path"]))
        series = event_series(root_path)
        x = float(series["gun_x_median"][0])
        series_by_x[x] = series
        n_left = series["n_left"]
        n_right = series["n_right"]
        n_total = n_left + n_right
        paired = (n_left > 0) & (n_right > 0)
        far_empty = n_right == 0 if x < 0 else n_left == 0 if x > 0 else (n_left == 0) | (n_right == 0)
        near = n_left if x < 0 else n_right if x > 0 else np.maximum(n_left, n_right)
        far = n_right if x < 0 else n_left if x > 0 else np.minimum(n_left, n_right)
        light_share = np.divide(n_left, n_total, out=np.full_like(n_total, np.nan, dtype=float), where=n_total > 0)
        npe_row = {
            "x_mm": x,
            "Npe_L_mean": float(np.mean(n_left)),
            "Npe_L_median": float(np.median(n_left)),
            "Npe_L_rms": rms(n_left),
            "Npe_R_mean": float(np.mean(n_right)),
            "Npe_R_median": float(np.median(n_right)),
            "Npe_R_rms": rms(n_right),
            "Npe_END_mean": float(np.mean(n_total)),
            "Npe_near_mean": float(np.mean(near)),
            "Npe_far_mean": float(np.mean(far)),
            "light_share_mean": float(np.nanmean(light_share)),
            "paired_fraction": float(np.mean(paired)),
            "fraction_far_end_empty": float(np.mean(far_empty)),
        }
        npe_rows.append(npe_row)
        light_rows.append({
            "x_mm": x,
            "Npe_L": npe_row["Npe_L_mean"],
            "Npe_R": npe_row["Npe_R_mean"],
            "Npe_near": npe_row["Npe_near_mean"],
            "Npe_far": npe_row["Npe_far_mean"],
            "light_share": npe_row["light_share_mean"],
            "paired_fraction": npe_row["paired_fraction"],
            "fraction_far_end_empty": npe_row["fraction_far_end_empty"],
        })

        t_left_all = series["t_left"]
        t_right_all = series["t_right"]
        pair_mask = np.isfinite(t_left_all) & np.isfinite(t_right_all)
        if np.count_nonzero(pair_mask) < 5:
            raise RuntimeError(f"x={x}: cannot form t_L/t_R")
        t_left = t_left_all[pair_mask]
        t_right = t_right_all[pair_mask]
        nl = n_left[pair_mask].astype(float)
        nr = n_right[pair_mask].astype(float)
        sigma_l = robust_sigma_ps(t_left)
        sigma_r = robust_sigma_ps(t_right)
        near_t = t_left if x < 0 else t_right if x > 0 else np.where(nl >= nr, t_left, t_right)
        far_t = t_right if x < 0 else t_left if x > 0 else np.where(nl >= nr, t_right, t_left)
        endpoint_rows.append({
            "x_mm": x,
            "sigma_L_ps": sigma_l,
            "sigma_R_ps": sigma_r,
            "sigma_near_ps": robust_sigma_ps(near_t),
            "sigma_far_ps": robust_sigma_ps(far_t),
            "mean_tL_ns": float(np.mean(t_left)),
            "mean_tR_ns": float(np.mean(t_right)),
        })
        tavg = 0.5 * (t_left + t_right)
        if sigma_l and sigma_r and sigma_l > 0.0 and sigma_r > 0.0:
            wl = 1.0 / sigma_l**2
            wr = 1.0 / sigma_r**2
            wmean_static = (wl * t_left + wr * t_right) / (wl + wr)
            weight_l = wl / (wl + wr)
        else:
            wmean_static = tavg
            weight_l = None
        wmean_event = (nl * t_left + nr * t_right) / (nl + nr)
        wmean_event_npe2 = (nl**2 * t_left + nr**2 * t_right) / (nl**2 + nr**2)
        sigma_gls, sigma_gls_dist, rho = gls(t_left, t_right)
        delta_t = t_left - t_right
        combined_rows.append({
            "x_mm": x,
            "sigma_tavg_ps": robust_sigma_ps(tavg),
            "sigma_wmean_static_ps": robust_sigma_ps(wmean_static),
            "sigma_wmean_event_ps": robust_sigma_ps(wmean_event),
            "sigma_wmean_event_npe2_ps": robust_sigma_ps(wmean_event_npe2),
            "sigma_GLS_residual_ps": sigma_gls,
            "sigma_GLS_estimator_dist_ps": sigma_gls_dist,
            "sigma_near_ps": robust_sigma_ps(near_t),
            "sigma_far_ps": robust_sigma_ps(far_t),
            "weight_L_fraction": weight_l,
            "weight_R_fraction": None if weight_l is None else 1.0 - weight_l,
            "rho_LR": rho,
            "mean_Delta_t_ns": float(np.mean(delta_t)),
            "sigma_Delta_t_ns": robust_sigma_ps(delta_t) / 1000.0 if robust_sigma_ps(delta_t) else None,
        })

    npe_rows.sort(key=lambda row: float(row["x_mm"]))
    light_rows.sort(key=lambda row: float(row["x_mm"]))
    endpoint_rows.sort(key=lambda row: float(row["x_mm"]))
    combined_rows.sort(key=lambda row: float(row["x_mm"]))
    write_csv(out_dir / "scan11_npe_summary.csv", npe_rows)
    write_csv(out_dir / "scan11_light_sharing.csv", light_rows)
    write_csv(out_dir / "scan11_endpoint_timing.csv", endpoint_rows)
    write_csv(out_dir / "scan11_combined_estimators.csv", combined_rows)

    fit_rows = [
        row for row in combined_rows
        if row["sigma_Delta_t_ns"] is not None and
        next(n for n in npe_rows if float(n["x_mm"]) == float(row["x_mm"]))["Npe_far_mean"] >= 5.0
    ]
    fit_ok = len(fit_rows) >= 5
    fit = {"fit_ok": fit_ok, "n_fit_points": len(fit_rows)}
    if fit_ok:
        xs = np.array([float(row["x_mm"]) for row in fit_rows])
        ys = np.array([float(row["mean_Delta_t_ns"]) for row in fit_rows])
        b, a = np.polyfit(xs, ys, 1)
        pred = a + b * xs
        residual = ys - pred
        v_eff_cm_per_ns = abs(2.0 / b) / 10.0 if b != 0 else None
        fit.update({
            "intercept_ns": float(a),
            "slope_ns_per_mm": float(b),
            "v_eff_cm_per_ns": v_eff_cm_per_ns,
            "residual_rms_ns": rms(residual),
        })
        for row in combined_rows:
            sx = None
            if row["sigma_Delta_t_ns"] is not None and b != 0:
                sx = float(row["sigma_Delta_t_ns"]) / abs(float(b)) / 10.0
            spatial_rows.append({
                "x_mm": row["x_mm"],
                "mean_Delta_t_ns": row["mean_Delta_t_ns"],
                "sigma_Delta_t_ns": row["sigma_Delta_t_ns"],
                "sigma_x_cm": sx,
                "included_in_fit": any(float(fr["x_mm"]) == float(row["x_mm"]) for fr in fit_rows),
            })
    else:
        for row in combined_rows:
            spatial_rows.append({
                "x_mm": row["x_mm"],
                "mean_Delta_t_ns": row["mean_Delta_t_ns"],
                "sigma_Delta_t_ns": row["sigma_Delta_t_ns"],
                "sigma_x_cm": None,
                "included_in_fit": False,
            })
    write_csv(out_dir / "scan11_spatial.csv", spatial_rows)
    write_json(out_dir / "scan11_veff_fit.json", fit)

    make_plots(out_dir, npe_rows, endpoint_rows, combined_rows, spatial_rows)
    anchor = anchor_comparison(combined_rows, spatial_rows)
    write_csv(out_dir / "real_anchor_comparison.csv", [anchor])
    write_json(out_dir / "real_anchor_comparison.json", anchor)
    return {
        "npe": npe_rows,
        "light": light_rows,
        "endpoint": endpoint_rows,
        "combined": combined_rows,
        "spatial": spatial_rows,
        "fit": fit,
        "anchor": anchor,
    }


def anchor_comparison(combined_rows: list[dict[str, object]], spatial_rows: list[dict[str, object]]) -> dict[str, object]:
    center = min(combined_rows, key=lambda row: abs(float(row["x_mm"])))
    w = np.array([float(row["sigma_wmean_static_ps"]) for row in combined_rows if row["sigma_wmean_static_ps"] is not None])
    g = np.array([float(row["sigma_GLS_residual_ps"]) for row in combined_rows if row["sigma_GLS_residual_ps"] is not None])
    sx = np.array([float(row["sigma_x_cm"]) for row in spatial_rows if row["sigma_x_cm"] is not None])
    center_spatial = min(spatial_rows, key=lambda row: abs(float(row["x_mm"])))
    return {
        "REAL_SINGLE_FILT_PS": REAL_SINGLE_FILT_PS,
        "REAL_BOTHENDS_WMEAN_LOW_PS": REAL_BOTHENDS_WMEAN_PS[0],
        "REAL_BOTHENDS_WMEAN_HIGH_PS": REAL_BOTHENDS_WMEAN_PS[1],
        "REAL_SPATIAL_CM": REAL_SPATIAL_CM,
        "max_sigma_wmean_static_ps": float(np.max(w)) if w.size else None,
        "median_sigma_wmean_static_ps": float(np.median(w)) if w.size else None,
        "center_sigma_wmean_static_ps": center["sigma_wmean_static_ps"],
        "max_sigma_GLS_ps": float(np.max(g)) if g.size else None,
        "median_sigma_GLS_ps": float(np.median(g)) if g.size else None,
        "center_sigma_GLS_ps": center["sigma_GLS_residual_ps"],
        "center_sigma_tavg_ps": center["sigma_tavg_ps"],
        "spatial_sigma_x_center_cm": center_spatial["sigma_x_cm"],
        "spatial_sigma_x_median_cm": float(np.median(sx)) if sx.size else None,
    }


def make_plots(
    out_dir: pathlib.Path,
    npe_rows: list[dict[str, object]],
    endpoint_rows: list[dict[str, object]],
    combined_rows: list[dict[str, object]],
    spatial_rows: list[dict[str, object]],
) -> None:
    png = out_dir / "png"
    png.mkdir(parents=True, exist_ok=True)
    x = [float(row["x_mm"]) for row in npe_rows]
    plt.figure(figsize=(7, 4.5))
    plt.plot(x, [row["Npe_L_mean"] for row in npe_rows], marker="o", label="END_L")
    plt.plot(x, [row["Npe_R_mean"] for row in npe_rows], marker="o", label="END_R")
    plt.yscale("log")
    plt.xlabel("gun x [mm]")
    plt.ylabel("mean PE / event")
    plt.legend()
    plt.tight_layout()
    plt.savefig(png / "npe_lr_vs_x.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(x, [row["Npe_near_mean"] for row in npe_rows], marker="o", label="near")
    plt.plot(x, [row["Npe_far_mean"] for row in npe_rows], marker="o", label="far")
    plt.yscale("log")
    plt.xlabel("gun x [mm]")
    plt.ylabel("mean PE / event")
    plt.legend()
    plt.tight_layout()
    plt.savefig(png / "npe_near_far_vs_x.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(x, [row["light_share_mean"] for row in npe_rows], marker="o")
    plt.xlabel("gun x [mm]")
    plt.ylabel("END_L / END total")
    plt.ylim(-0.05, 1.05)
    plt.tight_layout()
    plt.savefig(png / "light_share_vs_x.png", dpi=160)
    plt.close()

    ex = [float(row["x_mm"]) for row in endpoint_rows]
    plt.figure(figsize=(7, 4.5))
    plt.plot(ex, [row["sigma_L_ps"] for row in endpoint_rows], marker="o", label="L")
    plt.plot(ex, [row["sigma_R_ps"] for row in endpoint_rows], marker="o", label="R")
    plt.axhline(REAL_SINGLE_FILT_PS, color="0.4", ls="--", lw=1)
    plt.xlabel("gun x [mm]")
    plt.ylabel("sigma [ps]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(png / "sigma_LR_vs_x.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(ex, [row["sigma_near_ps"] for row in endpoint_rows], marker="o", label="near")
    plt.plot(ex, [row["sigma_far_ps"] for row in endpoint_rows], marker="o", label="far")
    plt.axhline(REAL_SINGLE_FILT_PS, color="0.4", ls="--", lw=1)
    plt.xlabel("gun x [mm]")
    plt.ylabel("sigma [ps]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(png / "sigma_near_far_vs_x.png", dpi=160)
    plt.close()

    cx = [float(row["x_mm"]) for row in combined_rows]
    plt.figure(figsize=(7, 4.5))
    for key, label in [
        ("sigma_tavg_ps", "t_avg"),
        ("sigma_wmean_static_ps", "wmean static"),
        ("sigma_wmean_event_ps", "wmean event"),
        ("sigma_GLS_residual_ps", "GLS"),
    ]:
        plt.plot(cx, [row[key] for row in combined_rows], marker="o", label=label)
    plt.axhspan(REAL_BOTHENDS_WMEAN_PS[0], REAL_BOTHENDS_WMEAN_PS[1], color="tab:green", alpha=0.15)
    plt.xlabel("gun x [mm]")
    plt.ylabel("sigma [ps]")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(png / "sigma_estimators_vs_x.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(cx, [row["weight_L_fraction"] for row in combined_rows], marker="o", label="L")
    plt.plot(cx, [row["weight_R_fraction"] for row in combined_rows], marker="o", label="R")
    plt.xlabel("gun x [mm]")
    plt.ylabel("static weight fraction")
    plt.ylim(-0.05, 1.05)
    plt.legend()
    plt.tight_layout()
    plt.savefig(png / "weights_static_vs_x.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(cx, [row["rho_LR"] for row in combined_rows], marker="o")
    plt.xlabel("gun x [mm]")
    plt.ylabel("rho(L,R)")
    plt.ylim(-1.05, 1.05)
    plt.tight_layout()
    plt.savefig(png / "gls_correlation_vs_x.png", dpi=160)
    plt.close()

    sx = [float(row["x_mm"]) for row in spatial_rows]
    plt.figure(figsize=(7, 4.5))
    plt.plot(sx, [row["mean_Delta_t_ns"] for row in spatial_rows], marker="o")
    plt.xlabel("gun x [mm]")
    plt.ylabel("mean Delta_t [ns]")
    plt.tight_layout()
    plt.savefig(png / "delta_t_mean_vs_x.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(sx, [row["sigma_x_cm"] for row in spatial_rows], marker="o")
    plt.axhline(REAL_SPATIAL_CM, color="0.4", ls="--", lw=1)
    plt.xlabel("gun x [mm]")
    plt.ylabel("sigma_x [cm]")
    plt.tight_layout()
    plt.savefig(png / "sigma_x_vs_x.png", dpi=160)
    plt.close()


def mirror_flags(npe_rows: list[dict[str, object]]) -> list[str]:
    flags = []
    by_x = {float(row["x_mm"]): row for row in npe_rows}
    for x in sorted(abs(float(row["x_mm"])) for row in npe_rows if float(row["x_mm"]) > 0):
        pos = by_x.get(x)
        neg = by_x.get(-x)
        if not pos or not neg:
            continue
        pairs = [
            (float(pos["Npe_R_mean"]), float(neg["Npe_L_mean"])),
            (float(pos["Npe_L_mean"]), float(neg["Npe_R_mean"])),
        ]
        for a, b in pairs:
            denom = max((a + b) / 2.0, 1.0e-9)
            if abs(a - b) / denom > 0.15:
                flags.append("asymmetry_concern")
                return flags
    return flags


def write_report(
    path: pathlib.Path,
    tag: str,
    inventory: list[dict[str, object]],
    analysis: dict[str, object],
    speed_rows: list[dict[str, object]],
    flags: list[str],
    results: dict[str, object],
) -> None:
    completed = [row for row in inventory if row.get("completed")]
    lines = [
        "# EXEC_28 - END-only scan11 weighted/GLS baseline",
        "",
        "## Veredicto",
        "",
        str(results["verdict"]),
        "",
        "## Resumen ejecutivo",
        "",
        f"- Configuracion optica: ground sigma_alpha=0.20 rad, reflector R=0.95, surface_loss=0.10, PDE nominal.",
        f"- Posiciones completadas: {len(completed)} de 11; eventos target por posicion: {results['scan']['events_per_position_target']}.",
        f"- Centro: t_avg={results['timing']['center_sigma_tavg_ps']:.2f} ps, wmean_static={results['timing']['center_sigma_wmean_static_ps']:.2f} ps, GLS={results['timing']['center_sigma_GLS_ps']:.2f} ps.",
        f"- Mediana wmean_static={results['timing']['median_sigma_wmean_static_ps']:.2f} ps; max={results['timing']['max_sigma_wmean_static_ps']:.2f} ps.",
        f"- v_eff={results['spatial']['v_eff_cm_per_ns']:.2f} cm/ns; sigma_x centro={results['spatial']['center_sigma_x_cm']:.2f} cm; sigma_x mediana={results['spatial']['median_sigma_x_cm']:.2f} cm.",
        "- Los estimadores weighted pueden quedar dominados por el extremo cercano en los bordes; los valores de pocos ps en borde no son una resolucion global two-ended.",
        "",
        "## Configuracion",
        "",
        f"- scint-air finish: {CONFIG['scint_air_finish']}",
        f"- sigma_alpha: {CONFIG['sigma_alpha_rad']} rad",
        f"- reflector R: {CONFIG['reflector_R']}",
        f"- surface_loss: {CONFIG['surface_loss']}",
        f"- PDE mode: {CONFIG['pde_mode']}",
        f"- air gap: {CONFIG['air_gap_mm']} mm",
        "- Branch: `exec28-endonly-scan11-weighted`",
        "- Commit: pending at report-writing time",
        "",
        "## Scan inventory",
        "",
        "| x nominal | events completed | runtime s | root |",
        "|---:|---:|---:|---|",
    ]
    for row in inventory:
        lines.append(
            f"| {float(row['x_mm_nominal']):.0f} | {int(row.get('n_events_completed') or 0)} | {float(row['runtime_s']):.1f} | `{row['root_path']}` |"
        )
    lines += [
        "",
        "## Light sharing",
        "",
        "| x | Npe_L | Npe_R | Npe_near | Npe_far | light_share | paired_fraction |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in analysis["light"]:
        lines.append(
            f"| {float(row['x_mm']):.0f} | {float(row['Npe_L']):.2f} | {float(row['Npe_R']):.2f} | {float(row['Npe_near']):.2f} | {float(row['Npe_far']):.2f} | {float(row['light_share']):.3f} | {float(row['paired_fraction']):.3f} |"
        )
    lines += [
        "",
        "## Endpoint timing",
        "",
        "| x | sigma_L | sigma_R | sigma_near | sigma_far |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in analysis["endpoint"]:
        lines.append(
            f"| {float(row['x_mm']):.0f} | {float(row['sigma_L_ps']):.2f} | {float(row['sigma_R_ps']):.2f} | {float(row['sigma_near_ps']):.2f} | {float(row['sigma_far_ps']):.2f} |"
        )
    lines += [
        "",
        "## Combined timing estimators",
        "",
        "| x | sigma_tavg | sigma_wmean_static | sigma_wmean_event | sigma_GLS | weight_L_fraction | rho_LR |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in analysis["combined"]:
        lines.append(
            f"| {float(row['x_mm']):.0f} | {float(row['sigma_tavg_ps']):.2f} | {float(row['sigma_wmean_static_ps']):.2f} | {float(row['sigma_wmean_event_ps']):.2f} | {float(row['sigma_GLS_residual_ps']):.2f} | {float(row['weight_L_fraction']):.3f} | {float(row['rho_LR']):.3f} |"
        )
    lines += [
        "",
        "Mandatory caveat: weighted estimators can become near-end dominated at the edges. Do not interpret edge few-ps values as global two-ended resolution.",
        "",
        "## Spatial reconstruction",
        "",
        f"- v_eff: {results['spatial']['v_eff_cm_per_ns']:.3f} cm/ns",
        f"- fit_ok: {results['spatial']['fit_ok']}",
        f"- center sigma_x: {results['spatial']['center_sigma_x_cm']:.3f} cm",
        f"- median sigma_x: {results['spatial']['median_sigma_x_cm']:.3f} cm",
        "",
        "| x | mean Delta_t ns | sigma Delta_t ns | sigma_x cm | fit point |",
        "|---:|---:|---:|---:|---|",
    ]
    for row in analysis["spatial"]:
        sx = "n/a" if row["sigma_x_cm"] is None else f"{float(row['sigma_x_cm']):.3f}"
        lines.append(
            f"| {float(row['x_mm']):.0f} | {float(row['mean_Delta_t_ns']):.4f} | {float(row['sigma_Delta_t_ns']):.4f} | {sx} | {row['included_in_fit']} |"
        )
    lines += [
        "",
        "## Real anchor comparison",
        "",
        f"- REAL_SINGLE_FILT_PS={REAL_SINGLE_FILT_PS:.0f}",
        f"- REAL_BOTHENDS_WMEAN_PS={REAL_BOTHENDS_WMEAN_PS[0]:.0f}-{REAL_BOTHENDS_WMEAN_PS[1]:.0f}",
        f"- REAL_SPATIAL_CM={REAL_SPATIAL_CM:.1f}",
        f"- center wmean_static={results['timing']['center_sigma_wmean_static_ps']:.2f} ps",
        f"- center GLS={results['timing']['center_sigma_GLS_ps']:.2f} ps",
        f"- center t_avg={results['timing']['center_sigma_tavg_ps']:.2f} ps",
        f"- center sigma_x={results['spatial']['center_sigma_x_cm']:.2f} cm",
        "",
        "## Speed/RINDEX gate",
        "",
        "| x | gate | max_ratio | superluminal_fraction | missing_rindex |",
        "|---:|---|---:|---:|---:|",
    ]
    for row in speed_rows:
        lines.append(
            f"| {float(row['x_mm']):.0f} | {row['speed_gate']} | {float(row['max_ratio']):.6f} | {float(row['superluminal_fraction']):.6g} | {int(row['missing_rindex'])} |"
        )
    lines += [
        "",
        "## Flags",
        "",
    ]
    lines.extend([f"- `{flag}`" for flag in flags] or ["- none"])
    lines += [
        "",
        "## Rutas",
        "",
        "- `out/EXEC_28/scan11_run_inventory.csv`",
        "- `out/EXEC_28/scan11_npe_summary.csv`",
        "- `out/EXEC_28/scan11_light_sharing.csv`",
        "- `out/EXEC_28/scan11_endpoint_timing.csv`",
        "- `out/EXEC_28/scan11_combined_estimators.csv`",
        "- `out/EXEC_28/scan11_spatial.csv`",
        "- `out/EXEC_28/scan11_veff_fit.json`",
        "- `out/EXEC_28/real_anchor_comparison.json`",
        "- `out/EXEC_28/results_exec28.json`",
        "- `out/EXEC_28/png/`",
        "- raw ROOT/logs: `out/EXEC_28/raw/` (not committed)",
        "",
        "## Git",
        "",
        "- Branch: `exec28-endonly-scan11-weighted`",
        f"- Tag: `{tag}`",
        "- Commit: pending at report-writing time",
        "- Push printed, not executed:",
        "",
        "```bash",
        "git push origin exec28-endonly-scan11-weighted",
        "```",
        "",
        "## Proximo paso",
        "",
        str(results["next_step"]),
        "",
    ]
    path.write_text("\n".join(lines), encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("out/EXEC_28"))
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--events", type=int, default=200)
    parser.add_argument("--wall-clock-min", type=float, default=90.0)
    args = parser.parse_args()

    if args.events < 200:
        raise RuntimeError("EXEC_28 hard floor is 200 events per position")
    repo = args.repo.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    write_json(out_dir / "config_snapshot.json", CONFIG)
    tag_rows = subprocess.check_output(["git", "tag", "--list", "EXEC_28-pre-*", "--sort=-creatordate"], cwd=repo, text=True).splitlines()
    tag = tag_rows[0] if tag_rows else "unknown"

    deadline = time.monotonic() + args.wall_clock_min * 60.0
    inventory: list[dict[str, object]] = []
    partial = False
    for x in RUN_ORDER:
        name = f"scan11_x{int(x)}".replace("-", "m")
        config = RunConfig(name=name, x_mm=x, events=args.events)
        existing = out_dir / "raw" / name / "photon_hits.root"
        if existing.exists():
            item = {
                "name": name,
                "x_mm_nominal": x,
                "n_events_requested": args.events,
                "runtime_s": float((out_dir / "raw" / name / "runtime_s.txt").read_text().strip()) if (out_dir / "raw" / name / "runtime_s.txt").exists() else 0.0,
                "root_path": str(existing),
                "log_path": str(out_dir / "raw" / name / "run.log"),
                "audit_dir": str(out_dir / "raw" / name / "audit"),
                "returncode": 0,
                "timed_out": False,
                "completed": True,
            }
        else:
            remaining = deadline - time.monotonic()
            if remaining <= 60.0:
                partial = True
                break
            print(f"[EXEC_28] running x={x:g} events={args.events} remaining_min={remaining/60.0:.1f}", flush=True)
            item = run_config(config, repo, out_dir, args.threads, remaining)
            if item.get("timed_out") or item.get("returncode") != 0:
                partial = True
                inventory.append(item)
                break
        if item.get("completed"):
            series = event_series(pathlib.Path(str(item["root_path"])))
            item["n_events_completed"] = int(series["n_events"][0])
            item["gun_x_mm_median"] = float(series["gun_x_median"][0])
        inventory.append(item)

    inventory.sort(key=lambda row: float(row["x_mm_nominal"]))
    write_csv(out_dir / "scan11_run_inventory.csv", inventory)
    write_json(out_dir / "scan11_run_inventory.json", inventory)
    completed = [row for row in inventory if row.get("completed")]
    if len(completed) < 5:
        raise RuntimeError("ABORTADO-EN-T1: fewer than five positions completed")

    analysis = analyze_completed(inventory, out_dir)
    speed_rows = []
    for x in (-690.0, 0.0, 690.0):
        row = next((item for item in completed if abs(float(item["x_mm_nominal"]) - x) < 1.0e-9), None)
        if row is None:
            continue
        speed = speed_summary(pathlib.Path(str(row["audit_dir"])))
        speed_rows.append({"x_mm": x, **speed})
    write_csv(out_dir / "speed_audit_scan11.csv", speed_rows)
    write_json(out_dir / "speed_audit_scan11.json", speed_rows)
    speed_pass = bool(speed_rows) and all(row["speed_gate"] == "PASS" for row in speed_rows)
    if not speed_pass:
        raise RuntimeError("ABORTADO-EN-T7: speed/RINDEX gate failed")

    flags = []
    if partial or len(completed) < len(POSITIONS):
        flags.append("scan11_partial")
    npe_rows = analysis["npe"]
    if any(float(row["Npe_far_mean"]) < 5.0 for row in npe_rows):
        flags.append("far_end_starved")
    if any(float(row["Npe_near_mean"]) > 1000.0 for row in npe_rows):
        flags.append("near_end_high")
    flags += mirror_flags(npe_rows)
    if any(float(row["Npe_far_mean"]) < 5.0 for row in npe_rows):
        flags.append("spatial_fit_far_end_starved")
    edge_spatial = [row for row in analysis["spatial"] if abs(float(row["x_mm"])) > 500.0]
    if any(row["sigma_x_cm"] is None or float(row["sigma_x_cm"]) > 10.0 for row in edge_spatial):
        flags.append("sigma_x_edge_unreliable")
    flags = sorted(set(flags))

    combined = analysis["combined"]
    spatial = analysis["spatial"]
    center = min(combined, key=lambda row: abs(float(row["x_mm"])))
    center_spatial = min(spatial, key=lambda row: abs(float(row["x_mm"])))
    w = np.array([float(row["sigma_wmean_static_ps"]) for row in combined])
    g = np.array([float(row["sigma_GLS_residual_ps"]) for row in combined])
    sx = np.array([float(row["sigma_x_cm"]) for row in spatial if row["sigma_x_cm"] is not None])
    max_ratio = max(float(row["max_ratio"]) for row in speed_rows)
    max_super = max(float(row["superluminal_fraction"]) for row in speed_rows)
    missing = sum(int(row["missing_rindex"]) for row in speed_rows)
    results = {
        "exec_id": "EXEC_28",
        "verdict": "COMPLETADO-CON-FLAGS" if flags else "COMPLETADO",
        "config": {
            "scint_air_finish": CONFIG["scint_air_finish"],
            "sigma_alpha_rad": CONFIG["sigma_alpha_rad"],
            "reflector_R": CONFIG["reflector_R"],
            "surface_loss": CONFIG["surface_loss"],
        },
        "scan": {
            "n_positions_requested": 11,
            "n_positions_completed": len(completed),
            "events_per_position_target": args.events,
            "partial": partial or len(completed) < 11,
        },
        "timing": {
            "center_sigma_tavg_ps": center["sigma_tavg_ps"],
            "center_sigma_wmean_static_ps": center["sigma_wmean_static_ps"],
            "center_sigma_GLS_ps": center["sigma_GLS_residual_ps"],
            "median_sigma_wmean_static_ps": float(np.median(w)),
            "max_sigma_wmean_static_ps": float(np.max(w)),
            "median_sigma_GLS_ps": float(np.median(g)),
            "max_sigma_GLS_ps": float(np.max(g)),
        },
        "spatial": {
            "v_eff_cm_per_ns": analysis["fit"].get("v_eff_cm_per_ns"),
            "center_sigma_x_cm": center_spatial["sigma_x_cm"],
            "median_sigma_x_cm": float(np.median(sx)) if sx.size else None,
            "fit_ok": analysis["fit"].get("fit_ok"),
        },
        "speed_gate": {
            "pass": speed_pass,
            "max_ratio": max_ratio,
            "superluminal_fraction": max_super,
            "missing_rindex": missing,
        },
        "flags": flags,
        "next_step": (
            "EXEC_29 scan31 END-only with weighted/GLS estimators."
            if not partial and analysis["fit"].get("fit_ok")
            else "Complete remaining scan11 positions or review spatial far-end starvation before scan31."
        ),
    }
    write_json(out_dir / "results_exec28.json", results)
    write_report(repo / "EXEC_28_REPORT.md", tag, inventory, analysis, speed_rows, flags, results)
    print(json.dumps(clean(results), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
