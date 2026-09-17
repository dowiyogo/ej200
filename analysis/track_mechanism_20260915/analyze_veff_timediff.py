#!/usr/bin/env python3
"""Tabulate time-difference propagation observables for EXEC_46 campaigns."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit

MATERIALS = {0: "EJ-200", 1: "EJ-204", 2: "EJ-230"}
POSITIONS = np.array([-650, -500, -200, 0, 200, 500, 650], dtype=float)
REFERENCE_MM_NS = 155.0
EXPECTED_HASHES = {
    "v2": "764c643e6ccfc3d0b96c7fc20d6753e86139233599ca19359c8d95e585844290",
    "old": "816c302c3ed37f7eb29855bcad7e734a1b6a77238f60a30361635b68377e341d",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_root_and_meta(path, frame, metadata):
    root_path = path.with_suffix(".root")
    numeric = {}
    for column in frame.columns:
        values = frame[column].to_numpy()
        if np.issubdtype(values.dtype, np.number):
            numeric[column] = values
    with uproot.recreate(root_path) as root_file:
        root_file.mktree("data", {name: value.dtype for name, value in numeric.items()})
        root_file["data"].extend(numeric)
    metadata = dict(metadata, csv=str(path.resolve()), root=str(root_path.resolve()),
                    root_sha256=sha256(root_path))
    path.with_suffix(".meta.json").write_text(json.dumps(metadata, indent=2) + "\n")


def weighted_fit(x, y, se, cubic=False):
    if cubic:
        design = np.column_stack([x, x ** 3])
    else:
        design = np.column_stack([np.ones(len(x)), x])
    w = 1.0 / (se * se)
    covariance = np.linalg.inv((design.T * w) @ design)
    beta = covariance @ (design.T @ (w * y))
    residual = y - design @ beta
    chi2 = float(np.sum((residual / se) ** 2))
    ndf = len(y) - len(beta)
    return beta, np.sqrt(np.diag(covariance)), residual, chi2, ndf


def load_tree(path, expected_hash):
    actual = sha256(path)
    if actual != expected_hash:
        raise RuntimeError(f"{path}: SHA-256 mismatch: {actual}")
    with uproot.open(path) as root_file:
        return root_file["derived_events"].arrays(
            ["material_code", "x_mm", "t_left_ns", "t_right_ns", "t0_ns"], library="np")


def read_campaign(path, label, expected_hash):
    arrays = load_tree(path, expected_hash)
    rows = []
    for code, material in MATERIALS.items():
        for x in POSITIONS.astype(int):
            mask = (arrays["material_code"] == code) & (arrays["x_mm"] == x)
            left = arrays["t_left_ns"][mask]
            right = arrays["t_right_ns"][mask]
            t0 = arrays["t0_ns"][mask]
            diff = left - right
            rows.append({"campaign": label, "material": material, "x_mm": x,
                         "mean_tdiff_ns": float(diff.mean()),
                         "sem_tdiff_ns": float(diff.std(ddof=1) / math.sqrt(len(diff))),
                         "mean_tL_ns": float(left.mean()),
                         "sem_tL_ns": float(left.std(ddof=1) / math.sqrt(len(left))),
                         "mean_tR_ns": float(right.mean()),
                         "sem_tR_ns": float(right.std(ddof=1) / math.sqrt(len(right))),
                         "mean_T0_ns": float(t0.mean()),
                         "sem_T0_ns": float(t0.std(ddof=1) / math.sqrt(len(t0))),
                         "n_events": len(diff)})
    return pd.DataFrame(rows)


def campaign_tables(cells):
    u2, u3, u4, u5 = [], [], [], []
    for (campaign, material), group in cells.groupby(["campaign", "material"], sort=True):
        group = group.sort_values("x_mm")
        x = group.x_mm.to_numpy(dtype=float)
        y = group.mean_tdiff_ns.to_numpy()
        se = group.sem_tdiff_ns.to_numpy()
        beta, err, residual, chi2, ndf = weighted_fit(x, y, se)
        beta_c, err_c, residual_c, chi2_c, ndf_c = weighted_fit(x, y, se, cubic=True)
        row = {"campaign": campaign, "material": material,
               "slope_tdiff_ns_per_mm": beta[1], "slope_error_ns_per_mm": err[1],
               "intercept_ns": beta[0], "intercept_error_ns": err[0],
               "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf,
               "chi2_flag": "*" if chi2 / ndf > 5 else ""}
        u2.append(row)
        for observable, ykey in (("tL", "mean_tL_ns"), ("tR", "mean_tR_ns")):
            beta_o, err_o, resid_o, chi2_o, ndf_o = weighted_fit(
                x, group[ykey].to_numpy(), group["sem_" + observable + "_ns"].to_numpy())
            velocity = 1.0 / abs(beta_o[1])
            velocity_error = err_o[1] / (beta_o[1] ** 2)
            u3.append({"campaign": campaign, "material": material, "observable": observable,
                       "slope_ns_per_mm": beta_o[1], "slope_error_ns_per_mm": err_o[1],
                       "intercept_ns": beta_o[0], "intercept_error_ns": err_o[0],
                       "chi2": chi2_o, "ndf": ndf_o, "chi2_ndf": chi2_o / ndf_o,
                       "chi2_flag": "*" if chi2_o / ndf_o > 5 else "",
                       "v_eff_one_end_mm_ns": velocity, "v_eff_error_mm_ns": velocity_error})
        slope = beta[1]
        slope_error = err[1]
        u3.append({"campaign": campaign, "material": material, "observable": "tdiff",
                   "slope_ns_per_mm": slope, "slope_error_ns_per_mm": slope_error,
                   "intercept_ns": beta[0], "intercept_error_ns": err[0],
                   "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf,
                   "chi2_flag": "*" if chi2 / ndf > 5 else "",
                   "v_eff_two_end_mm_ns": 2.0 / abs(slope),
                   "v_eff_error_mm_ns": 2.0 * slope_error / (slope ** 2)})
        for _, point in group.iterrows():
            predicted = beta[0] + beta[1] * point.x_mm
            u4.append({"campaign": campaign, "material": material, "x_mm": int(point.x_mm),
                       "mean_tdiff_ns": point.mean_tdiff_ns, "sem_tdiff_ns": point.sem_tdiff_ns,
                       "linear_prediction_ns": predicted,
                       "linear_residual_ps": 1000.0 * (point.mean_tdiff_ns - predicted),
                       "linear_chi2_ndf": chi2 / ndf, "linear_chi2_flag": "*" if chi2 / ndf > 5 else ""})
        local_velocities = []
        for first, second in zip(group.iloc[:-1].itertuples(), group.iloc[1:].itertuples()):
            dx = second.x_mm - first.x_mm
            ds = second.mean_tdiff_ns - first.mean_tdiff_ns
            local_slope = ds / dx
            local_v2 = 2.0 / abs(local_slope)
            local_v1 = 1.0 / abs(local_slope)
            local_velocities.append(local_v2)
            u4.append({"campaign": campaign, "material": material,
                       "x_low_mm": first.x_mm, "x_high_mm": second.x_mm,
                       "secant_slope_ns_per_mm": local_slope,
                       "v_eff_two_end_mm_ns": local_v2,
                       "v_eff_one_end_mm_ns": local_v1})
        left_fit = next(r for r in u3 if r["campaign"] == campaign and r["material"] == material and r["observable"] == "tL")
        right_fit = next(r for r in u3 if r["campaign"] == campaign and r["material"] == material and r["observable"] == "tR")
        u5.append({"campaign": campaign, "material": material,
                   "v_eff_two_end_mm_ns": 2.0 / abs(slope),
                   "v_eff_two_end_error_mm_ns": 2.0 * slope_error / slope ** 2,
               "v_eff_one_end_L_mm_ns": left_fit["v_eff_one_end_mm_ns"],
               "v_eff_one_end_L_error_mm_ns": left_fit["v_eff_error_mm_ns"],
               "v_eff_one_end_R_mm_ns": right_fit["v_eff_one_end_mm_ns"],
               "v_eff_one_end_R_error_mm_ns": right_fit["v_eff_error_mm_ns"],
                   "chi2_ndf_linear": chi2 / ndf,
                   "chi2_flag": "*" if chi2 / ndf > 5 else "",
                   "local_v_eff_two_end_min_mm_ns": min(local_velocities),
                   "local_v_eff_two_end_max_mm_ns": max(local_velocities),
                   "experimental_reference_mm_ns": REFERENCE_MM_NS})
        u4.append({"campaign": campaign, "material": material,
                   "cubic_a1_ns_per_mm": beta_c[0], "cubic_a1_error_ns_per_mm": err_c[0],
                   "cubic_a3_ns_per_mm3": beta_c[1], "cubic_a3_error_ns_per_mm3": err_c[1],
                   "cubic_chi2": chi2_c, "cubic_ndf": ndf_c,
                   "cubic_chi2_ndf": chi2_c / ndf_c,
                   "cubic_chi2_flag": "*" if chi2_c / ndf_c > 5 else ""})
    return pd.DataFrame(u2), pd.DataFrame(u3), pd.DataFrame(u4), pd.DataFrame(u5)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v2", type=Path, required=True)
    parser.add_argument("--old", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    v2 = read_campaign(args.v2, "v2", EXPECTED_HASHES["v2"])
    old = read_campaign(args.old, "old", EXPECTED_HASHES["old"])
    cells = pd.concat([v2, old], ignore_index=True)
    u1 = cells[["campaign", "material", "x_mm", "mean_tdiff_ns", "sem_tdiff_ns", "n_events"]]
    u2, u3, u4, u5 = campaign_tables(cells)
    for name, frame in (("u1_tdiff_cells", u1), ("u2_slope_fits", u2),
                        ("u3_veff", u3), ("u4_linearity", u4), ("u5_contrast", u5)):
        path = args.output_dir / f"{name}.csv"
        frame.to_csv(path, index=False, float_format="%.12g")
        with uproot.recreate(path.with_suffix(".root")) as root_file:
            numeric = {c: frame[c].to_numpy() for c in frame.columns if np.issubdtype(frame[c].dtype, np.number)}
            root_file.mktree("data", {c: v.dtype for c, v in numeric.items()})
            root_file["data"].extend(numeric)
        path.with_suffix(".meta.json").write_text(json.dumps({"csv": str(path.resolve()), "source_hashes": EXPECTED_HASHES, "experimental_reference_mm_ns": REFERENCE_MM_NS}, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "u1_rows": len(u1), "u2_rows": len(u2), "u3_rows": len(u3), "u4_rows": len(u4), "u5_rows": len(u5)}, indent=2))


if __name__ == "__main__":
    main()
