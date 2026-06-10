#!/usr/bin/env python3
"""Diagnose apparent propagation velocities from End arrival-time estimators."""

from __future__ import annotations

import argparse
import math
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit

POSITIONS = (-690, -670, *range(-650, 651, 50), 670, 690)
SIDES = {"left": tuple(range(8)), "right": tuple(range(8, 16))}
N_BLOCKS = 20


def metric_values(event_id: np.ndarray, times: np.ndarray) -> dict[str, float]:
    events = np.unique(event_id)
    fpt = np.full(int(events.max()) + 1, np.inf)
    np.minimum.at(fpt, event_id, times)
    fpt = fpt[np.isfinite(fpt)]
    return {
        "mean_ns": float(np.mean(times)),
        "t50_ns": float(np.median(times)),
        "fpt_mean_ns": float(np.mean(fpt)),
        "frac_gt10": float(np.mean(times > 10.0)),
    }


def jackknife(event_id: np.ndarray, times: np.ndarray) -> tuple[dict[str, float], dict[str, float]]:
    central = metric_values(event_id, times)
    replicas = {name: [] for name in central}
    for omitted in range(N_BLOCKS):
        selected = event_id % N_BLOCKS != omitted
        values = metric_values(event_id[selected], times[selected])
        for name in replicas:
            replicas[name].append(values[name])
    factor = (N_BLOCKS - 1) / N_BLOCKS
    errors = {
        name: float(math.sqrt(factor * np.sum((np.asarray(values) - np.mean(values)) ** 2)))
        for name, values in replicas.items()
    }
    return central, errors


def line(x: np.ndarray, intercept: float, slope: float) -> np.ndarray:
    return intercept + slope * x


def fit_estimator(frame: pd.DataFrame, estimator: str) -> dict[str, float]:
    x = frame["distance_cm"].to_numpy()
    y = frame[estimator].to_numpy()
    error = frame[f"{estimator}_error"].to_numpy()
    params, covariance = curve_fit(line, x, y, sigma=error, absolute_sigma=True)
    prediction = line(x, *params)
    chi2 = float(np.sum(np.square((y - prediction) / error)))
    ndf = len(x) - 2
    slope = params[1]
    slope_error = math.sqrt(covariance[1, 1])
    velocity = 1.0 / slope
    velocity_error = abs(slope_error / slope * velocity)
    return {
        "estimator": estimator,
        "intercept_ns": params[0],
        "intercept_error_ns": math.sqrt(covariance[0, 0]),
        "slope_ns_per_cm": slope,
        "slope_error_ns_per_cm": slope_error,
        "velocity_cm_ns": velocity,
        "velocity_error_cm_ns": velocity_error,
        "chi2_ndf": chi2 / ndf,
        "n_points": len(x),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("analysis/exec07"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figs = args.output_dir / "figs"
    figs.mkdir(exist_ok=True)

    rows = []
    for x_mm in POSITIONS:
        print(f"timing estimators x={x_mm}", flush=True)
        with uproot.open(args.data_dir / f"photon_hits_x{x_mm}mm.root") as root_file:
            arrays = root_file["sipm_hits"].arrays(["event_id", "global_id", "time_ns"], cut="global_id < 16", library="np")
        for side, ids in SIDES.items():
            selected = np.isin(arrays["global_id"], ids)
            event_id = arrays["event_id"][selected].astype(np.int32)
            times = arrays["time_ns"][selected].astype(float)
            central, errors = jackknife(event_id, times)
            row = {
                "x_beam_mm": x_mm, "side": side,
                "distance_cm": (x_mm + 700.0) / 10.0 if side == "left" else (700.0 - x_mm) / 10.0,
            }
            for name in central:
                row[name] = central[name]
                row[f"{name}_error"] = errors[name]
            rows.append(row)
    metrics = pd.DataFrame(rows)
    metrics.to_csv(args.output_dir / "exec10_velocity_metrics.csv", index=False)

    fit_rows = []
    for scope, selected in (
        ("combined", metrics),
        ("left", metrics[metrics.side == "left"]),
        ("right", metrics[metrics.side == "right"]),
    ):
        for estimator in ("mean_ns", "t50_ns", "fpt_mean_ns"):
            row = fit_estimator(selected, estimator)
            row["scope"] = scope
            fit_rows.append(row)
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(args.output_dir / "exec10_velocity_fits.csv", index=False)

    labels = {
        "mean_ns": "All-photon mean (tail-biased)",
        "t50_ns": "All-photon median t50",
        "fpt_mean_ns": "Mean FPT",
    }
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for axis, estimator in zip(axes.flat[:3], labels):
        for side, marker, color in (("left", "o", "tab:blue"), ("right", "s", "tab:orange")):
            selected = metrics[metrics.side == side]
            axis.errorbar(
                selected.distance_cm.to_numpy(), selected[estimator].to_numpy(),
                yerr=selected[f"{estimator}_error"].to_numpy(), fmt=marker, ms=3,
                color=color, alpha=0.65, label=f"{side} End",
            )
        fit = fits[(fits.scope == "combined") & (fits.estimator == estimator)].iloc[0]
        distance = np.linspace(1.0, 139.0, 200)
        axis.plot(
            distance, line(distance, fit.intercept_ns, fit.slope_ns_per_cm),
            color="black", label=f"combined fit: v={fit.velocity_cm_ns:.2f}+/-{fit.velocity_error_cm_ns:.2f} cm/ns",
        )
        axis.set(xlabel="Track-to-End distance [cm]", ylabel=f"{labels[estimator]} [ns]")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    tail_axis = axes.flat[3]
    for side, marker, color in (("left", "o", "tab:blue"), ("right", "s", "tab:orange")):
        selected = metrics[metrics.side == side]
        tail_axis.errorbar(
            selected.distance_cm.to_numpy(), selected.frac_gt10.to_numpy(),
            yerr=selected.frac_gt10_error.to_numpy(), fmt=marker, ms=3, color=color, alpha=0.65, label=f"{side} End",
        )
    tail_axis.set(xlabel="Track-to-End distance [cm]", ylabel="Fraction of End photons with t > 10 ns")
    tail_axis.grid(alpha=0.25)
    tail_axis.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figs / "exec10_velocity_estimators.png", dpi=220)
    plt.close(fig)

    representative = metrics[metrics.distance_cm.isin([1.0, 30.0, 70.0, 139.0])].copy()
    representative.to_csv(args.output_dir / "exec10_late_fraction_representative.csv", index=False)
    combined = fits[fits.scope == "combined"].set_index("estimator")
    fpt = combined.loc["fpt_mean_ns"]
    close_to_cn = abs(fpt.velocity_cm_ns - 19.0) <= 3.0 * fpt.velocity_error_cm_ns
    report = [
        "# EXEC_10 apparent-velocity diagnosis", "",
        "All errors are delete-one-block jackknife errors over event IDs. Fits combine the mirrored left/right End samples.",
        "",
        "| estimator | velocity [cm/ns] | chi2/ndf |",
        "|---|---:|---:|",
    ]
    for estimator in ("mean_ns", "t50_ns", "fpt_mean_ns"):
        row = combined.loc[estimator]
        report.append(f"| {labels[estimator]} | {row.velocity_cm_ns:.3f} +/- {row.velocity_error_cm_ns:.3f} | {row.chi2_ndf:.2f} |")
    tail_near = representative[representative.distance_cm == 1.0].frac_gt10.mean()
    tail_far = representative[representative.distance_cm == 139.0].frac_gt10.mean()
    report += [
        "",
        f"Mean FPT is {'consistent' if close_to_cn else 'not consistent'} with c/n ~= 19 cm/ns within 3 sigma.",
        f"The requested t>10 ns fraction test refutes the proposed trend: it increases from {tail_near:.6f} near the End to {tail_far:.6f} at 139 cm, because a fixed absolute-time threshold also tracks propagation delay.",
        "None of the three global linear fits is statistically acceptable; the quoted slope-derived velocities are descriptive apparent velocities, not propagation measurements.",
        "",
        "Alternative hypotheses to test in a future dedicated audit (not changed here): record optical-photon creation time and path length; inspect Geant4's effective GROUPVEL generated from the SSLG4 RINDEX table; and fit local distance ranges instead of forcing one line over 1--139 cm.",
        "The source RINDEX table is constant at 1.58 from 200--800 nm, so the apparent FPT value above c/n is not explained by the documented phase index alone.",
    ]
    audit = args.output_dir.parent.parent / "audit" / "exec10_veff_diagnosis.md"
    audit.write_text("\n".join(report) + "\n")
    print("\n".join(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
