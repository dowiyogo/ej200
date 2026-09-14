#!/usr/bin/env python3
"""EXEC41 Block 1: corrected V2 analysis over the existing EXEC40 ROOT."""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import uproot
from scipy.stats import chi2

N = 2000
NBOOT = 500
SEED = 41091301
LOWER_BOUND = 0.226172
FACE_ORDER = ["+Z", "-Z", "-Y", "-X_SENSOR", "+X_SENSOR", "+Y_SENSOR",
              "OPEN_+Y_OR_+/-X_UNRESOLVED"]
POST_TO_FACE = {
    "AirGapZPlusPV": "+Z",
    "AirGapZMinusPV": "-Z",
    "AirGapYMinusPV": "-Y",
    "EndSiPMLeft_PV": "-X_SENSOR",
    "EndSiPMRight_PV": "+X_SENSOR",
    "TopSiPMPV": "+Y_SENSOR",
    "WorldPV": "OPEN_+Y_OR_+/-X_UNRESOLVED",
}

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")

def sidecars(out, stem, columns, metadata):
    data = {key: np.asarray(value) for key, value in columns.items()}
    with (out / f"{stem}.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(data)
        writer.writerows(zip(*data.values()))
    with uproot.recreate(out / f"{stem}.root") as root:
        root["data"] = data
    dump(out / f"{stem}.meta.json", {
        **metadata,
        "columns": list(data),
        "rows": len(next(iter(data.values()))),
        "csv_sha256": sha(out / f"{stem}.csv"),
        "root_sha256": sha(out / f"{stem}.root"),
    })

def bootstrap_ratio(numerator, denominator, weights):
    values = (weights @ numerator) / (weights @ denominator)
    return {
        "value": float(numerator.sum() / denominator.sum()),
        "se": float(values.std(ddof=1)),
        "numerator": int(numerator.sum()),
        "denominator": int(denominator.sum()),
    }

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--preregistration", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    weights = np.stack([
        np.bincount(rng.integers(0, N, N), minlength=N) for _ in range(NBOOT)
    ]).astype(float)
    face_counts = np.zeros((N, len(FACE_ORDER), 2), dtype=np.int64)
    hist = np.zeros((N, 20), dtype=np.int64)
    unmapped = {}
    selected_rows = 0
    tree = uproot.open(args.root)["first_bar_encounters"]
    for chunk in tree.iterate(
        ["event_id", "source", "post_volume", "outcome", "exiting_bar",
         "cos_incidence"], step_size="64 MB", library="np"
    ):
        event = chunk["event_id"].astype(np.int64)
        selected = (chunk["source"] == 1) & (chunk["exiting_bar"] == 1)
        event = event[selected]
        post = np.asarray(chunk["post_volume"])[selected]
        outcome = chunk["outcome"][selected]
        cosine = chunk["cos_incidence"][selected]
        selected_rows += len(event)
        for value in np.unique(post):
            name = str(value)
            if name not in POST_TO_FACE:
                unmapped[name] = unmapped.get(name, 0) + int(np.count_nonzero(post == value))
                continue
            face_index = FACE_ORDER.index(POST_TO_FACE[name])
            mask = post == value
            escaped = mask & (outcome == 2) & (
                name.startswith("AirGap") or name == "WorldPV")
            face_counts[:, face_index, 0] += np.bincount(event[mask], minlength=N)
            face_counts[:, face_index, 1] += np.bincount(event[escaped], minlength=N)
        valid = np.isfinite(cosine) & (cosine >= 0.0) & (cosine <= 1.0)
        bins = np.minimum((cosine[valid] * 20).astype(np.int64), 19)
        hist += np.bincount(20 * event[valid] + bins,
                            minlength=N * 20).reshape(N, 20)
    if unmapped:
        raise RuntimeError(f"Unmapped production volume names: {unmapped}")
    if selected_rows != int(face_counts[:, :, 0].sum()):
        raise RuntimeError("Face decomposition does not close")

    aggregate = bootstrap_ratio(face_counts[:, :, 1].sum(axis=1),
                                face_counts[:, :, 0].sum(axis=1), weights)
    aggregate["lower_3se"] = aggregate["value"] - 3 * aggregate["se"]
    aggregate["lower_bound"] = LOWER_BOUND
    aggregate["status"] = "PASS" if aggregate["lower_3se"] >= LOWER_BOUND else "FAIL"
    faces = []
    total_encounters = face_counts[:, :, 0].sum()
    total_escapes = face_counts[:, :, 1].sum()
    for index, face in enumerate(FACE_ORDER):
        entry = bootstrap_ratio(face_counts[:, index, 1],
                                face_counts[:, index, 0], weights)
        share_values = (weights @ face_counts[:, index, 0]) / (
            weights @ face_counts[:, :, 0].sum(axis=1))
        entry.update({
            "face": face,
            "encounter_share": float(face_counts[:, index, 0].sum() / total_encounters),
            "encounter_share_se": float(share_values.std(ddof=1)),
            "escape_contribution_to_all_encounters": float(
                face_counts[:, index, 1].sum() / total_encounters),
            "escape_share": float(face_counts[:, index, 1].sum() / total_escapes),
            "compatible_with_0p226172_3se": bool(
                abs(entry["value"] - LOWER_BOUND) <= 3 * entry["se"]),
        })
        faces.append(entry)

    nearest_indices = [FACE_ORDER.index("+Z"), FACE_ORDER.index("-Z")]
    nearest_escapes = int(face_counts[:, nearest_indices, 1].sum())
    other_escapes = int(total_escapes - nearest_escapes)
    mechanism = {
        "nearest_Z_escaped": nearest_escapes,
        "other_faces_escaped": other_escapes,
        "other_faces_share_of_escaped": float(other_escapes / total_escapes),
        "aggregate_excess_over_bound": aggregate["value"] - LOWER_BOUND,
        "interpretation": "Arithmetic decomposition only; unresolved WorldPV rows are not assigned to +Y or +/-X.",
    }

    edges = np.linspace(0.0, 1.0, 21)
    observed = hist.sum(axis=0) / hist.sum()
    bootstrap_hist = weights @ hist
    bootstrap_hist /= bootstrap_hist.sum(axis=1)[:, None]
    covariance = np.cov(bootstrap_hist, rowvar=False, ddof=1)
    rank = int(np.linalg.matrix_rank(covariance))
    def gof(expected):
        delta = observed - expected
        statistic = float(delta @ np.linalg.pinv(covariance) @ delta)
        return {"chi_square_event_bootstrap": statistic, "rank": rank,
                "p_value": float(chi2.sf(statistic, rank))}
    uniform = np.full(20, 0.05)
    flux = np.diff(edges ** 2)
    low_mu = {
        "range": [0.0, 0.2],
        "observed_probability": float(observed[:4].sum()),
        "uniform_probability": 0.2,
        "difference": float(observed[:4].sum() - 0.2),
        "bootstrap_se": float(bootstrap_hist[:, :4].sum(axis=1).std(ddof=1)),
    }
    angular = {
        "status": "DESCRIPTIVE_NO_PASS_FAIL",
        "uniform": gof(uniform),
        "old_flux_p_2mu": gof(flux),
        "low_mu": low_mu,
        "mean_cosine_binned": float(np.sum(observed * (edges[:-1] + edges[1:]) / 2)),
        "photons": int(hist.sum()),
    }
    old = {
        "escape_fraction": aggregate["value"],
        "escape_se": aggregate["se"],
        "expected_interval": [0.36, 0.40],
        "status": "FAIL_RETAINED_HYPOTHESIS_PHYSICALLY_INAPPROPRIATE",
        "p_2mu": angular["old_flux_p_2mu"],
    }
    metadata = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_root": str(args.root), "source_root_sha256": sha(args.root),
        "N": N, "seeds": [26092601, 8349041], "bootstrap_seed": SEED,
        "bootstrap_replicates": NBOOT, "bootstrap_unit": "generated event",
        "preregistration": str(args.preregistration),
        "preregistration_sha256": sha(args.preregistration),
        "script_sha256": sha(Path(__file__)),
    }
    sidecars(args.out, "v2_face_decomposition", {
        "face": [row["face"] for row in faces],
        "encounters": [row["denominator"] for row in faces],
        "encounter_share": [row["encounter_share"] for row in faces],
        "encounter_share_se": [row["encounter_share_se"] for row in faces],
        "escaped": [row["numerator"] for row in faces],
        "conditional_escape_fraction": [row["value"] for row in faces],
        "conditional_escape_se": [row["se"] for row in faces],
        "escape_contribution_to_all_encounters": [row["escape_contribution_to_all_encounters"] for row in faces],
        "escape_share": [row["escape_share"] for row in faces],
        "compatible_with_0p226172_3se": [row["compatible_with_0p226172_3se"] for row in faces],
    }, metadata)
    sidecars(args.out, "v2_uniform_cosine", {
        "bin_low": edges[:-1], "bin_high": edges[1:],
        "photons": hist.sum(axis=0), "observed_probability": observed,
        "probability_se": np.sqrt(np.diag(covariance)),
        "uniform_probability": uniform, "old_flux_probability": flux,
    }, metadata)
    sidecars(args.out, "v2_face_events", {
        "event_id": np.arange(N),
        **{f"{face}_encounters": face_counts[:, i, 0] for i, face in enumerate(FACE_ORDER)},
        **{f"{face}_escaped": face_counts[:, i, 1] for i, face in enumerate(FACE_ORDER)},
        **{f"cos_bin_{i:02}": hist[:, i] for i in range(20)},
    }, metadata)
    result = {
        "metadata": metadata, "H1": aggregate, "faces": faces,
        "mechanism": mechanism, "H2": angular, "old_EXEC40_hypothesis": old,
        "geometry_mm": {"source_x": 0, "source_y": 0, "source_z_range": [-5, 5],
                        "nearest_Z_distance_range": [0, 5],
                        "minus_Y_distance": 30, "plus_Y_plane_distance": 30,
                        "minus_X_distance": 700, "plus_X_distance": 700},
        "face_resolution_limitation": "BarPV->WorldPV cannot distinguish open +Y from +/-X without position or normal components.",
    }
    dump(args.out / "v2_results.json", result)
    print(json.dumps({"H1": aggregate, "H2": angular,
                      "mechanism": mechanism}, indent=2))

if __name__ == "__main__":
    main()
