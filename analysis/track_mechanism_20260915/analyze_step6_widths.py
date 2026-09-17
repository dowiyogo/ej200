#!/usr/bin/env python3
"""Compute the non-SE_gap Step 6 width and L2 tables from the derived tree."""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

EXPECTED_SHA256 = "764c643e6ccfc3d0b96c7fc20d6753e86139233599ca19359c8d95e585844290"
MATERIALS = {0: "EJ-200", 1: "EJ-204", 2: "EJ-230"}
POSITIONS = (-650, -500, -200, 0, 200, 500, 650)
IQR_SCALE = 1.349
HISTOGRAM_BIN_NS = 0.004
LOCAL_HALF_WIDTH_NS = 0.320
GAUSSIAN_WINDOW_Q68 = 2.0


def sha256(path):
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gaussian_fit(values):
    median = float(np.median(values))
    q16, q84 = np.quantile(values, (0.16, 0.84))
    q68 = (q84 - q16) / 2.0
    fit_low = median - GAUSSIAN_WINDOW_Q68 * q68
    fit_high = median + GAUSSIAN_WINDOW_Q68 * q68
    hist_low = median - LOCAL_HALF_WIDTH_NS
    hist_high = median + LOCAL_HALF_WIDTH_NS
    edges = np.arange(hist_low, hist_high + HISTOGRAM_BIN_NS, HISTOGRAM_BIN_NS)
    counts, edges = np.histogram(values, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    selected = (centers >= fit_low) & (centers <= fit_high) & (counts > 0)
    if np.count_nonzero(selected) < 3:
        return {"gaussian_mean_ns": np.nan, "gaussian_sigma_ns": np.nan,
                "gaussian_chi2": np.nan, "gaussian_ndf": 0,
                "gaussian_chi2_ndf": np.nan, "gaussian_fit_reliable": False,
                "gaussian_status": "NOT_AVAILABLE"}
    x = centers[selected]
    y = counts[selected].astype(float)
    try:
        from scipy.optimize import curve_fit
        def model(x, amplitude, mean, sigma):
            return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
        initial = [float(y.max()), median, max(q68, HISTOGRAM_BIN_NS)]
        parameters, covariance = curve_fit(model, x, y, p0=initial,
                                            bounds=([0.0, hist_low, HISTOGRAM_BIN_NS / 4],
                                                    [np.inf, hist_high, 2.0]),
                                            maxfev=10000)
        expected = model(x, *parameters)
        chi2 = float(np.sum((y - expected) ** 2 / np.maximum(expected, 1.0)))
        ndf = max(0, len(x) - 3)
        chi2_ndf = chi2 / ndf if ndf else np.nan
        return {"gaussian_mean_ns": float(parameters[1]),
                "gaussian_sigma_ns": float(abs(parameters[2])),
                "gaussian_chi2": chi2, "gaussian_ndf": ndf,
                "gaussian_chi2_ndf": chi2_ndf,
                "gaussian_fit_reliable": bool(np.isfinite(chi2_ndf) and chi2_ndf <= 5.0),
                "gaussian_status": "PASS"}
    except Exception as error:
        return {"gaussian_mean_ns": np.nan, "gaussian_sigma_ns": np.nan,
                "gaussian_chi2": np.nan, "gaussian_ndf": 0,
                "gaussian_chi2_ndf": np.nan, "gaussian_fit_reliable": False,
                "gaussian_status": f"NOT_AVAILABLE: {type(error).__name__}"}


def width_row(material, x_mm, observable, values):
    values = np.asarray(values, dtype=float)
    rms = float(np.std(values, ddof=0))
    q25, q75 = np.quantile(values, (0.25, 0.75))
    row = {"material": material, "x_mm": int(x_mm), "observable": observable,
           "n_events": len(values), "mean_ns": float(np.mean(values)),
           "rms_ns": rms, "iqr_over_1p349_ns": float((q75 - q25) / IQR_SCALE)}
    row.update(gaussian_fit(values))
    row["gaussian_flag"] = "*" if (np.isfinite(row["gaussian_chi2_ndf"]) and
                                    row["gaussian_chi2_ndf"] > 5.0) else ""
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--derived", type=Path, required=True)
    parser.add_argument("--mixture", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    actual_sha = sha256(args.derived)
    if actual_sha != EXPECTED_SHA256:
        raise RuntimeError(f"derived tree SHA-256 mismatch: {actual_sha}")
    with uproot.open(args.derived) as root_file:
        arrays = root_file["derived_events"].arrays(
            ["material_code", "x_mm", "t_left_ns", "t_right_ns", "t0_ns"],
            library="np")
    rows = []
    l2_rows = []
    for code, material in MATERIALS.items():
        for x_mm in POSITIONS:
            mask = (arrays["material_code"] == code) & (arrays["x_mm"] == x_mm)
            values = {name: arrays[name][mask] for name in ("t_left_ns", "t_right_ns", "t0_ns")}
            for observable, data in (("tL", values["t_left_ns"]),
                                     ("tR", values["t_right_ns"]),
                                     ("T0", values["t0_ns"])):
                rows.append(width_row(material, x_mm, observable, data))
            rms_left = float(np.std(values["t_left_ns"], ddof=0))
            rms_right = float(np.std(values["t_right_ns"], ddof=0))
            rms_t0 = float(np.std(values["t0_ns"], ddof=0))
            bound = (rms_left + rms_right) / 2.0
            margin = rms_t0 - bound
            l2_rows.append({"material": material, "x_mm": int(x_mm),
                            "rms_t0_ns": rms_t0, "rms_end_average_ns": bound,
                            "margin_ns": margin, "se_gap_ns": "NOT_AVAILABLE",
                            "margin_flag": "CHECK" if margin >= -0.001 else ""})
    mixture = pd.read_csv(args.mixture)
    mixture_rows = mixture[mixture["gun_x_mm"].abs() == 650].copy()
    mixture_rows = mixture_rows[["material", "gun_x_mm", "sigma_mixture_ns",
                                 "corr_cherenkov_tprop_angle_penalty"]]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "t6_widths.csv", index=False, float_format="%.12g")
    pd.DataFrame(l2_rows).to_csv(args.output_dir / "t6_l2_guardrail.csv", index=False, float_format="%.12g")
    mixture_rows.to_csv(args.output_dir / "t6_mixture_and_correlation.csv", index=False, float_format="%.12g")
    metadata = {"derived_sha256": actual_sha, "histogram_bin_ns": HISTOGRAM_BIN_NS,
                "histogram_range_half_width_ns": LOCAL_HALF_WIDTH_NS,
                "gaussian_fit_window": "median +/- 2*q68, q68=(q84-q16)/2",
                "se_gap_ns": "NOT_AVAILABLE",
                "se_gap_note": "SE_gap delta-method by paired influence function is not implemented in this repository.",
                "sources": [str(args.derived.resolve()), str(args.mixture.resolve())]}
    (args.output_dir / "t6_widths.meta.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "width_rows": len(rows), "l2_rows": len(l2_rows),
                      "mixture_rows": len(mixture_rows), "derived_sha256": actual_sha}, indent=2))


if __name__ == "__main__":
    main()
