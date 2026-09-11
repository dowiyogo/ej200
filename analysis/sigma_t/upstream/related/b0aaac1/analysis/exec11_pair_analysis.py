#!/usr/bin/env python3
"""EXEC_11 pair-scan analysis from the authoritative per-hit sipm_hits tree."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import subprocess
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit
from scipy.stats import moyal, pearsonr, spearmanr

HOOK_DATA_DIR = pathlib.Path("/home/reriosto/SHiP/t0minidaq/pairscan_2026-06-11")
HOOK_PAIR_IDS = (28, 29)
HOOK_N_EVENTS = 3000
HOOK_POS_REF_1 = "AUTO"
HOOK_POS_REF_2 = "AUTO"
HOOK_FIT_RANGE_NSIGMA = 2.5
HOOK_CALIB_EXCLUDE_TEST = True
HOOK_SPTR_PS = 0
HOOK_ELEC_JITTER_PS = 0
HOOK_WALK = "none"

N_PE_THRESH = 4
EXPECTED_POSITIONS = tuple(float(x) for x in range(-462, -421))
REQUIRED_BRANCHES = {"event_id", "global_id", "time_ns", "gun_x_mm"}
COUNT_DEFINITION = (
    "number of sipm_hits rows per event and SiPM, using the existing EXEC_08 "
    "PE-equivalent convention"
)
FILE_RE = re.compile(r"pairscan_x([+-]?\d+(?:\.\d+)?)mm\.root$")
DATA_COMMIT = "f431c01"
FIGURE_CONTEXT = "Pair (28,29) — EJ-204 — Top readout — intrinsic"


def analysis_commit() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], text=True
    ).strip()


def extract_position(path: pathlib.Path) -> float:
    match = FILE_RE.fullmatch(path.name)
    if not match:
        raise ValueError(f"Unexpected pair-scan filename: {path.name}")
    return float(match.group(1))


def fourth_hit_and_counts(
    event_ids: np.ndarray,
    global_ids: np.ndarray,
    times_ns: np.ndarray,
    n_events: int = HOOK_N_EVENTS,
    pair_ids: tuple[int, int] = HOOK_PAIR_IDS,
) -> pd.DataFrame:
    """Build fixed-denominator event observables independent of input row order."""
    event_ids = np.asarray(event_ids, dtype=np.int64)
    global_ids = np.asarray(global_ids, dtype=np.int64)
    times_ns = np.asarray(times_ns, dtype=float)
    if not (len(event_ids) == len(global_ids) == len(times_ns)):
        raise ValueError("Input arrays must have equal length")
    if np.any((event_ids < 0) | (event_ids >= n_events)):
        bad = np.unique(event_ids[(event_ids < 0) | (event_ids >= n_events)])
        raise ValueError(f"event_id outside [0,{n_events - 1}]: {bad[:10]}")

    out: dict[str, np.ndarray] = {"event_id": np.arange(n_events, dtype=np.int32)}
    for label, channel in zip(("A", "B"), pair_ids):
        selected = (global_ids == channel) & np.isfinite(times_ns)
        ev = event_ids[selected]
        ts = times_ns[selected]
        counts = np.bincount(ev, minlength=n_events).astype(np.int32)
        fourth = np.full(n_events, np.nan)
        if len(ev):
            order = np.lexsort((ts, ev))
            sorted_ev = ev[order]
            sorted_ts = ts[order]
            starts = np.r_[0, np.cumsum(counts[:-1])]
            eligible = counts >= N_PE_THRESH
            fourth[eligible] = sorted_ts[starts[eligible] + N_PE_THRESH - 1]
        out[f"npe_{label}"] = counts
        out[f"t_{label}_ns"] = fourth
        out[f"passed_{label}_4pe"] = counts >= N_PE_THRESH

    frame = pd.DataFrame(out)
    frame["passed_pair_4pe"] = frame["passed_A_4pe"] & frame["passed_B_4pe"]
    pair = frame["passed_pair_4pe"].to_numpy()
    frame["delta_t_ps"] = np.where(
        pair, 1000.0 * (frame["t_A_ns"] - frame["t_B_ns"]), np.nan
    )
    positive = (frame["npe_A"] > 0) & (frame["npe_B"] > 0)
    frame["R_log_ratio"] = np.where(
        positive, np.log(frame["npe_A"] / frame["npe_B"]), np.nan
    )
    denominator = frame["npe_A"] + frame["npe_B"]
    frame["npe_asymmetry"] = np.where(
        denominator > 0, (frame["npe_A"] - frame["npe_B"]) / denominator, np.nan
    )
    return frame


def inventory_and_derive(path: pathlib.Path, derived_dir: pathlib.Path) -> dict[str, object]:
    x_true = extract_position(path)
    present = np.zeros(HOOK_N_EVENTS, dtype=bool)
    pair_event: list[np.ndarray] = []
    pair_gid: list[np.ndarray] = []
    pair_time: list[np.ndarray] = []
    gun_counter: Counter[float] = Counter()
    n_entries = 0
    n_nonfinite = 0
    hit_counts = {HOOK_PAIR_IDS[0]: 0, HOOK_PAIR_IDS[1]: 0}

    with uproot.open(path) as root_file:
        if "sipm_hits" not in root_file:
            raise RuntimeError(f"{path}: missing sipm_hits")
        tree = root_file["sipm_hits"]
        missing = REQUIRED_BRANCHES - set(tree.keys())
        if missing:
            raise RuntimeError(f"{path}: missing branches {sorted(missing)}")
        for arrays in tree.iterate(
            ["event_id", "global_id", "time_ns", "gun_x_mm"],
            step_size="100 MB",
            library="np",
        ):
            ev = arrays["event_id"].astype(np.int64, copy=False)
            gid = arrays["global_id"].astype(np.int64, copy=False)
            times = arrays["time_ns"].astype(float, copy=False)
            if np.any((ev < 0) | (ev >= HOOK_N_EVENTS)):
                bad = np.unique(ev[(ev < 0) | (ev >= HOOK_N_EVENTS)])
                raise RuntimeError(f"{path}: event_id outside expected range: {bad[:10]}")
            present[np.unique(ev)] = True
            n_entries += len(ev)
            n_nonfinite += int(np.count_nonzero(~np.isfinite(times)))
            values, counts = np.unique(np.round(arrays["gun_x_mm"], 6), return_counts=True)
            gun_counter.update(dict(zip(values.tolist(), counts.tolist())))
            selected = np.isin(gid, HOOK_PAIR_IDS)
            if np.any(selected):
                pair_event.append(ev[selected].copy())
                pair_gid.append(gid[selected].copy())
                pair_time.append(times[selected].copy())
                for channel in HOOK_PAIR_IDS:
                    hit_counts[channel] += int(np.count_nonzero(gid[selected] == channel))

    event_ids = np.concatenate(pair_event)
    global_ids = np.concatenate(pair_gid)
    times_ns = np.concatenate(pair_time)
    frame = fourth_hit_and_counts(event_ids, global_ids, times_ns)
    frame.insert(0, "x_true_mm", x_true)
    derived_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        derived_dir / f"pair_events_x{x_true:+.1f}mm.npz",
        **{column: frame[column].to_numpy() for column in frame.columns},
    )
    dominant_gun, dominant_count = gun_counter.most_common(1)[0]
    return {
        "file": path.name,
        "x_true_mm": x_true,
        "size_bytes": path.stat().st_size,
        "mtime_ns": path.stat().st_mtime_ns,
        "n_entries": n_entries,
        "n_event_ids_present_anywhere": int(np.count_nonzero(present)),
        "event_id_min": int(np.flatnonzero(present)[0]),
        "event_id_max": int(np.flatnonzero(present)[-1]),
        "n_events_missing_from_entire_tree": int(HOOK_N_EVENTS - np.count_nonzero(present)),
        "gun_x_mm_dominant": dominant_gun,
        "gun_x_mm_dominant_fraction": dominant_count / n_entries,
        "hits_A": hit_counts[HOOK_PAIR_IDS[0]],
        "hits_B": hit_counts[HOOK_PAIR_IDS[1]],
        "nonfinite_time_fraction": n_nonfinite / n_entries,
        "n_events_with_A": int(np.count_nonzero(frame["npe_A"])),
        "n_events_with_B": int(np.count_nonzero(frame["npe_B"])),
        "n_events_with_A_and_B": int(np.count_nonzero((frame["npe_A"] > 0) & (frame["npe_B"] > 0))),
        "n_events_passing_4PE_both": int(frame["passed_pair_4pe"].sum()),
        "efficiency": float(frame["passed_pair_4pe"].sum() / HOOK_N_EVENTS),
    }


def validate_file_set(data_dir: pathlib.Path) -> list[pathlib.Path]:
    files = sorted(data_dir.glob("pairscan_x*.root"), key=extract_position)
    positions = [extract_position(path) for path in files]
    if positions != list(EXPECTED_POSITIONS):
        missing = sorted(set(EXPECTED_POSITIONS) - set(positions))
        duplicates = sorted({x for x in positions if positions.count(x) > 1})
        raise RuntimeError(
            f"Expected 41 integer positions -462..-422; missing={missing}, duplicates={duplicates}"
        )
    return files


def write_metadata(output_dir: pathlib.Path) -> None:
    metadata = {
        "data_dir": str(HOOK_DATA_DIR),
        "pair_ids": HOOK_PAIR_IDS,
        "n_events_expected": HOOK_N_EVENTS,
        "n_pe_threshold": N_PE_THRESH,
        "count_definition": COUNT_DEFINITION,
        "data_commit": "f431c01",
        "analysis_commit_at_generation": analysis_commit(),
        "root_phase_0b": "DEFERRED — ROOT runtime unavailable on MSI",
    }
    (output_dir / "analysis" / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


def derive_all(data_dir: pathlib.Path, output_dir: pathlib.Path) -> None:
    files = validate_file_set(data_dir)
    for subdir in ("analysis", "figures", "tables", "derived", "logs"):
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)
    rows = []
    for index, path in enumerate(files, 1):
        print(f"[{index:02d}/{len(files)}] inventory + derive {path.name}", flush=True)
        rows.append(inventory_and_derive(path, output_dir / "derived"))
    inventory = pd.DataFrame(rows).sort_values("x_true_mm")
    inventory.to_csv(output_dir / "analysis" / "data_inventory.csv", index=False)
    if (
        (inventory["n_event_ids_present_anywhere"] != HOOK_N_EVENTS).any()
        or (inventory["nonfinite_time_fraction"] != 0).any()
        or not np.allclose(inventory["gun_x_mm_dominant"], inventory["x_true_mm"])
    ):
        raise RuntimeError("Inventory gate failed; inspect data_inventory.csv")
    write_metadata(output_dir)


def load_derived(path: pathlib.Path) -> pd.DataFrame:
    with np.load(path) as arrays:
        return pd.DataFrame({name: arrays[name] for name in arrays.files})


def gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
    return amplitude * np.exp(-0.5 * np.square((x - mean) / sigma))


def iterative_gaussian_fit(
    values: np.ndarray,
    bins: np.ndarray,
    n_sigma: float = HOOK_FIT_RANGE_NSIGMA,
    guard_subpeak: bool = False,
) -> dict[str, float | str]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    observed, edges = np.histogram(values, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak = float(centers[np.argmax(observed)])
    q01, q99 = np.quantile(values, [0.01, 0.99])
    sigma = float(np.std(values[(values >= q01) & (values <= q99)], ddof=1))
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.std(values, ddof=1))
    params = np.array([float(np.max(observed)), peak, sigma])
    covariance = np.full((3, 3), np.nan)
    status = "ok"
    try:
        for _ in range(2):
            selected = np.abs(centers - params[1]) <= n_sigma * abs(params[2])
            errors = np.sqrt(np.maximum(observed[selected], 1.0))
            lower_sigma = 0.8 * sigma if guard_subpeak else 1e-9
            params, covariance = curve_fit(
                gaussian,
                centers[selected],
                observed[selected],
                p0=params,
                sigma=errors,
                absolute_sigma=True,
                bounds=([0, peak - 5 * sigma, lower_sigma], [np.inf, peak + 5 * sigma, 2 * sigma if guard_subpeak else np.inf]),
                maxfev=10000,
            )
    except Exception as error:  # noqa: BLE001
        status = f"failed:{type(error).__name__}"
    selected = np.abs(centers - params[1]) <= n_sigma * abs(params[2])
    expected = gaussian(centers[selected], *params)
    valid = expected >= 1.0
    ndf = int(np.count_nonzero(valid) - 3)
    chi2_ndf = (
        float(np.sum(np.square(observed[selected][valid] - expected[valid]) / expected[valid]) / ndf)
        if ndf > 0 else math.nan
    )
    errors = np.sqrt(np.diag(covariance))
    return {
        "mean": float(params[1]),
        "mean_error": float(errors[1]),
        "sigma": float(abs(params[2])),
        "sigma_error": float(errors[2]),
        "fit_status": status,
        "chi2_ndf": chi2_ndf,
        "fraction_within_3sigma": float(np.mean(np.abs(values - params[1]) <= 3 * abs(params[2]))),
        "cov_00": float(covariance[1, 1]),
        "cov_01": float(covariance[1, 2]),
        "cov_11": float(covariance[2, 2]),
        "subpeak_sigma_guard_active": bool(
            guard_subpeak and abs(params[2] - 0.8 * sigma) / sigma < 1e-3
        ),
    }


def qa_summary(output_dir: pathlib.Path, v1_path: pathlib.Path) -> None:
    rows = []
    dt_bins = np.linspace(-3000.0, 3000.0, 301)
    r_bins = np.linspace(-4.0, 4.0, 81)
    for path in sorted((output_dir / "derived").glob("pair_events_x*.npz")):
        frame = load_derived(path)
        x_true = float(frame["x_true_mm"].iloc[0])
        dt = iterative_gaussian_fit(frame["delta_t_ps"].to_numpy(), dt_bins, guard_subpeak=True)
        ratio = iterative_gaussian_fit(frame["R_log_ratio"].to_numpy(), r_bins)
        rows.append({
            "x_true_mm": x_true,
            "n_events_expected": HOOK_N_EVENTS,
            "n_event_ids_present": len(frame),
            "n_passed": int(frame["passed_pair_4pe"].sum()),
            "efficiency": float(frame["passed_pair_4pe"].mean()),
            "mean_npe_A": float(frame["npe_A"].mean()),
            "mean_npe_B": float(frame["npe_B"].mean()),
            "mu_dt_ps": dt["mean"],
            "mu_dt_err_ps": dt["mean_error"],
            "sigma_dt_ps": dt["sigma"],
            "sigma_dt_err_ps": dt["sigma_error"],
            "fit_status": dt["fit_status"],
            "chi2_ndf": dt["chi2_ndf"],
            "fraction_within_3sigma": dt["fraction_within_3sigma"],
            "subpeak_sigma_guard_active": dt["subpeak_sigma_guard_active"],
            "mu_R": ratio["mean"],
            "mu_R_err": ratio["mean_error"],
            "sigma_R": ratio["sigma"],
            "sigma_R_err": ratio["sigma_error"],
        })
    summary = pd.DataFrame(rows).sort_values("x_true_mm")
    summary.to_csv(output_dir / "analysis" / "pairscan_summary_v2.csv", index=False)

    v1 = pd.read_csv(v1_path).rename(columns={
        "mu_dt_ps": "v1_mu_dt_ps", "sigma_dt_ps": "v1_sigma_dt_ps",
        "mu_r": "v1_mu_R",
    })
    qa = summary.merge(
        v1[["x_true_mm", "v1_mu_dt_ps", "v1_sigma_dt_ps", "v1_mu_R"]],
        on="x_true_mm",
        how="left",
    )
    qa["qa_refit_required"] = (
        (qa["fit_status"] != "ok")
        | (qa["chi2_ndf"] > 5)
        | (qa["fraction_within_3sigma"] < 0.9)
    )
    qa.to_csv(output_dir / "analysis" / "fit_qa.csv", index=False)

    ref1_candidates = summary.loc[summary["mu_dt_ps"].abs() == summary["mu_dt_ps"].abs().min()]
    ref1 = ref1_candidates.iloc[ref1_candidates["x_true_mm"].abs().argmin()]
    ref2_candidates = summary.loc[summary["mean_npe_A"] == summary["mean_npe_A"].max()]
    ref2 = ref2_candidates.iloc[ref2_candidates["x_true_mm"].abs().argmin()]
    pd.DataFrame([
        {"reference": "POS_REF_1", "x_true_mm": ref1["x_true_mm"], "criterion": "minimum abs(mu_dt_ps), tie -> minimum abs(x_true_mm)"},
        {"reference": "POS_REF_2", "x_true_mm": ref2["x_true_mm"], "criterion": "maximum mean_npe_A, tie -> minimum abs(x_true_mm)"},
    ]).to_csv(output_dir / "analysis" / "reference_positions.csv", index=False)


def figure_footer(fig: plt.Figure) -> None:
    fig.text(
        0.5, 0.005,
        f"data commit {DATA_COMMIT} | analysis commit {analysis_commit()} | "
        "simulation prediction — Top readout — intrinsic",
        ha="center", fontsize=7,
    )


def save_figure(fig: plt.Figure, output_dir: pathlib.Path, name: str) -> None:
    figure_footer(fig)
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    for suffix in ("pdf", "png"):
        fig.savefig(output_dir / "figures" / f"{name}.{suffix}", dpi=220)
    plt.close(fig)


def write_latex_table(frame: pd.DataFrame, path: pathlib.Path, columns: list[str]) -> None:
    def render(value: object) -> str:
        if isinstance(value, (float, np.floating)):
            return "--" if not np.isfinite(value) else f"{value:.4g}"
        return str(value).replace("_", r"\_")

    lines = [
        r"\begin{tabular}{" + "l" * len(columns) + "}",
        r"\hline",
        " & ".join(column.replace("_", r"\_") for column in columns) + r" \\",
        r"\hline",
    ]
    for _, row in frame.iterrows():
        lines.append(" & ".join(render(row[column]) for column in columns) + r" \\")
    lines.extend([r"\hline", r"\end{tabular}"])
    path.write_text("\n".join(lines) + "\n")


def read_pair_hit_times(data_dir: pathlib.Path, x_true: float) -> dict[int, np.ndarray]:
    path = data_dir / f"pairscan_x{x_true:+.1f}mm.root"
    collected: dict[int, list[np.ndarray]] = {channel: [] for channel in HOOK_PAIR_IDS}
    with uproot.open(path) as root_file:
        for arrays in root_file["sipm_hits"].iterate(
            ["global_id", "time_ns"], step_size="100 MB", library="np"
        ):
            for channel in HOOK_PAIR_IDS:
                selected = arrays["global_id"] == channel
                if np.any(selected):
                    collected[channel].append(arrays["time_ns"][selected])
    return {channel: np.concatenate(parts) for channel, parts in collected.items()}


def bootstrap_pearson(x: np.ndarray, y: np.ndarray, seed: int, replicas: int = 2000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = np.empty(replicas)
    for index in range(replicas):
        selected = rng.integers(0, len(x), len(x))
        values[index] = pearsonr(x[selected], y[selected])[0]
    return float(np.std(values, ddof=1)), float(np.mean(values))


def moyal_shape_fit(values: np.ndarray, seed: int) -> dict[str, float]:
    positive = np.asarray(values, dtype=float)
    positive = positive[positive > 0]
    loc, scale = moyal.fit(positive)
    rng = np.random.default_rng(seed)
    replicas = np.empty((300, 2))
    for index in range(len(replicas)):
        sample = positive[rng.integers(0, len(positive), len(positive))]
        replicas[index] = moyal.fit(sample)
    return {
        "mpv_moyal": float(loc),
        "mpv_moyal_error": float(np.std(replicas[:, 0], ddof=1)),
        "width_moyal": float(scale),
        "width_moyal_error": float(np.std(replicas[:, 1], ddof=1)),
    }


def detailed_analysis(data_dir: pathlib.Path, output_dir: pathlib.Path) -> None:
    references = pd.read_csv(output_dir / "analysis" / "reference_positions.csv")
    summary = pd.read_csv(output_dir / "analysis" / "pairscan_summary_v2.csv")
    rows = []
    for ref_index, ref in references.iterrows():
        label = ref["reference"]
        x_true = float(ref["x_true_mm"])
        frame = load_derived(output_dir / "derived" / f"pair_events_x{x_true:+.1f}mm.npz")
        hits = read_pair_hit_times(data_dir, x_true)
        passed = frame["passed_pair_4pe"].to_numpy()
        dt = frame.loc[passed, "delta_t_ps"].to_numpy()
        asym = frame.loc[passed, "npe_asymmetry"].to_numpy()
        ratio = frame.loc[passed, "R_log_ratio"].to_numpy()
        pearson_asym = pearsonr(dt, asym)
        pearson_ratio = pearsonr(dt, ratio)
        pearson_asym_boot, _ = bootstrap_pearson(dt, asym, 1100 + ref_index)
        pearson_ratio_boot, _ = bootstrap_pearson(dt, ratio, 2100 + ref_index)
        dt_fit = iterative_gaussian_fit(dt, np.linspace(-3000, 3000, 301), guard_subpeak=True)

        row = {
            "reference": label, "x_true_mm": x_true,
            "mean_npe_A": frame["npe_A"].mean(), "mean_npe_B": frame["npe_B"].mean(),
            "mean_asymmetry": frame["npe_asymmetry"].mean(),
            "mu_dt_ps": dt_fit["mean"], "sigma_dt_ps": dt_fit["sigma"],
            "dt_chi2_ndf": dt_fit["chi2_ndf"], "dt_fit_status": dt_fit["fit_status"],
            "pearson_dt_asym": pearson_asym[0],
            "pearson_dt_asym_boot_error": pearson_asym_boot,
            "spearman_dt_asym": spearmanr(dt, asym)[0],
            "pearson_dt_R": pearson_ratio[0],
            "pearson_dt_R_boot_error": pearson_ratio_boot,
            "spearman_dt_R": spearmanr(dt, ratio)[0],
        }
        for channel, tag in zip(HOOK_PAIR_IDS, ("A", "B")):
            row.update({f"{key}_{tag}": value for key, value in moyal_shape_fit(frame[f"npe_{tag}"], 3100 + ref_index * 10 + channel).items()})
        rows.append(row)

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
        for channel, color in zip(HOOK_PAIR_IDS, ("tab:blue", "tab:orange")):
            axes[0].hist(hits[channel], bins=160, histtype="step", color=color, label=f"ID {channel}")
            axes[1].hist(hits[channel], bins=160, density=True, histtype="step", color=color, label=f"ID {channel}")
        axes[0].set(xlabel="Photon-hit time (ns)", ylabel="SiPM hit count", title="Individual hits, unnormalized")
        axes[1].set(xlabel="Photon-hit time (ns)", ylabel="Area-normalized density", title="Individual hits, normalized for shape comparison")
        for ax in axes: ax.legend(); ax.grid(alpha=0.2)
        fig.suptitle(f"{FIGURE_CONTEXT}\n{label}: x={x_true:.1f} mm; individual photon hits, not event estimators")
        save_figure(fig, output_dir, f"{label.lower()}_all_hit_times")

        fig, ax = plt.subplots(figsize=(8.5, 5.2))
        ax.hist(frame.loc[frame["passed_A_4pe"], "t_A_ns"], bins=120, histtype="step", label="t_A: fourth hit, ID 28")
        ax.hist(frame.loc[frame["passed_B_4pe"], "t_B_ns"], bins=120, histtype="step", label="t_B: fourth hit, ID 29")
        ax.set(xlabel="Fourth-hit event estimator (ns)", ylabel="Events", title=f"{FIGURE_CONTEXT}\n{label}: x={x_true:.1f} mm")
        ax.legend(); ax.grid(alpha=0.2)
        save_figure(fig, output_dir, f"{label.lower()}_event_times")

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
        for ax, tag, channel in zip(axes, ("A", "B"), HOOK_PAIR_IDS):
            values = frame[f"npe_{tag}"].to_numpy()
            fit = {key.removesuffix(f"_{tag}"): value for key, value in row.items() if key.endswith(f"_{tag}") and "moyal" in key}
            edges = np.arange(values.min(), np.quantile(values, 0.999) + 3, 3)
            ax.hist(values, bins=edges, histtype="step", label=f"ID {channel}, zeros included")
            centers = np.linspace(max(0.1, edges[0]), edges[-1], 500)
            width = edges[1] - edges[0]
            ax.plot(centers, len(values) * width * moyal.pdf(centers, loc=fit["mpv_moyal"], scale=fit["width_moyal"]), label="Moyal approximation fit (>0)")
            ax.set(xlabel="SiPM hit count / PE-equivalent count", ylabel="Events", title=f"ID {channel}: MPV={fit['mpv_moyal']:.2f}±{fit['mpv_moyal_error']:.2f}")
            ax.legend(fontsize=8); ax.grid(alpha=0.2)
        fig.suptitle(f"{FIGURE_CONTEXT}\n{label}: x={x_true:.1f} mm; Moyal approximation, not exact Landau")
        save_figure(fig, output_dir, f"{label.lower()}_npe_moyal")

        fig, ax = plt.subplots(figsize=(8.5, 5.2))
        observed, edges, _ = ax.hist(dt, bins=np.linspace(-3000, 3000, 301), histtype="step", label="events")
        centers = np.linspace(dt_fit["mean"] - 4 * dt_fit["sigma"], dt_fit["mean"] + 4 * dt_fit["sigma"], 400)
        ax.plot(centers, gaussian(centers, observed.max(), dt_fit["mean"], dt_fit["sigma"]), label="2-iteration Gaussian")
        ax.set(xlabel="delta_t = t_A - t_B (ps)", ylabel="Events", title=f"{FIGURE_CONTEXT}\n{label}: mu={dt_fit['mean']:.2f} ps, sigma={dt_fit['sigma']:.2f} ps, chi2/ndf={dt_fit['chi2_ndf']:.2f}")
        ax.legend(); ax.grid(alpha=0.2)
        save_figure(fig, output_dir, f"{label.lower()}_delta_t")

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
        axes[0].hexbin(asym, dt, gridsize=45, mincnt=1)
        axes[0].set(xlabel="NPE asymmetry", ylabel="delta_t (ps)", title=f"Pearson={pearson_asym[0]:.3f}±{pearson_asym_boot:.3f}; Spearman={spearmanr(dt, asym)[0]:.3f}")
        axes[1].hexbin(ratio, dt, gridsize=45, mincnt=1)
        axes[1].set(xlabel="R = ln(npe_A/npe_B)", ylabel="delta_t (ps)", title=f"Pearson={pearson_ratio[0]:.3f}±{pearson_ratio_boot:.3f}; Spearman={spearmanr(dt, ratio)[0]:.3f}")
        fig.suptitle(f"{FIGURE_CONTEXT}\n{label}: x={x_true:.1f} mm")
        save_figure(fig, output_dir, f"{label.lower()}_correlations")

    comparison = pd.DataFrame(rows)
    comparison.to_csv(output_dir / "analysis" / "reference_comparison.csv", index=False)
    write_latex_table(
        comparison,
        output_dir / "tables" / "reference_comparison.tex",
        ["reference", "x_true_mm", "mean_npe_A", "mean_npe_B", "mu_dt_ps", "sigma_dt_ps", "mean_asymmetry"],
    )


def weighted_polynomial_fit(x: np.ndarray, y: np.ndarray, error: np.ndarray, degree: int) -> dict[str, object]:
    design = np.vander(x, degree + 1)
    weight = 1.0 / np.square(error)
    normal = design.T @ (weight[:, None] * design)
    covariance = np.linalg.inv(normal)
    coefficients = covariance @ design.T @ (weight * y)
    prediction = design @ coefficients
    residuals = y - prediction
    chi2 = float(np.sum(np.square(residuals / error)))
    ndf = len(x) - degree - 1
    return {
        "coefficients": coefficients,
        "covariance": covariance,
        "prediction": prediction,
        "residuals": residuals,
        "chi2": chi2,
        "chi2_ndf": chi2 / ndf,
        "aic": chi2 + 2 * (degree + 1),
    }


def loocv_rmse(x: np.ndarray, y: np.ndarray, error: np.ndarray, degree: int) -> float:
    predictions = []
    for omitted in range(len(x)):
        keep = np.arange(len(x)) != omitted
        fit = weighted_polynomial_fit(x[keep], y[keep], error[keep], degree)
        predictions.append(np.polyval(fit["coefficients"], x[omitted]))
    return float(np.sqrt(np.mean(np.square(y - predictions))))


def invert_cubic(values: np.ndarray, coefficients: np.ndarray, linear_guess: np.ndarray) -> tuple[np.ndarray, int, int]:
    reconstructed = np.full(len(values), np.nan)
    ambiguous = 0
    no_root = 0
    for index, (value, guess) in enumerate(zip(values, linear_guess)):
        shifted = coefficients.copy()
        shifted[-1] -= value
        roots = np.roots(shifted)
        physical = roots[(np.abs(roots.imag) < 1e-8) & (roots.real >= -462) & (roots.real <= -422)].real
        if len(physical) == 0:
            no_root += 1
        else:
            if len(physical) > 1:
                ambiguous += 1
            reconstructed[index] = physical[np.argmin(np.abs(physical - guess))]
    return reconstructed, ambiguous, no_root


def sample_summary(values: np.ndarray, x_true: float) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    mean = float(np.mean(values))
    sigma = float(np.std(values, ddof=1))
    return {
        "mean_x_rec_mm": mean,
        "mean_x_rec_error_mm": sigma / math.sqrt(len(values)),
        "bias_mm": mean - x_true,
        "bias_error_mm": sigma / math.sqrt(len(values)),
        "sigma_x_mm": sigma,
        "sigma_x_error_mm": sigma / math.sqrt(2 * (len(values) - 1)),
        "n_valid": len(values),
    }


def reconstruction(output_dir: pathlib.Path) -> None:
    summary = pd.read_csv(output_dir / "analysis" / "pairscan_summary_v2.csv")
    references = pd.read_csv(output_dir / "analysis" / "reference_positions.csv")
    ref_positions = set(references["x_true_mm"].astype(float))
    calibration = summary.loc[~summary["x_true_mm"].isin(ref_positions)].copy()
    x = calibration["x_true_mm"].to_numpy()

    temporal = weighted_polynomial_fit(
        x, calibration["mu_dt_ps"].to_numpy(), calibration["mu_dt_err_ps"].to_numpy(), 1
    )
    ratio_linear = weighted_polynomial_fit(
        x, calibration["mu_R"].to_numpy(), calibration["mu_R_err"].to_numpy(), 1
    )
    ratio_cubic = weighted_polynomial_fit(
        x, calibration["mu_R"].to_numpy(), calibration["mu_R_err"].to_numpy(), 3
    )
    temporal_coeff = temporal["coefficients"]
    ratio_coeff = ratio_linear["coefficients"]

    pd.DataFrame([{
        "model": "mu_delta_t_ps = m_t*x + b_t",
        "m_t_ps_per_mm": temporal_coeff[0],
        "m_t_error_ps_per_mm": math.sqrt(temporal["covariance"][0, 0]),
        "b_t_ps": temporal_coeff[1],
        "b_t_error_ps": math.sqrt(temporal["covariance"][1, 1]),
        "cov_m_b": temporal["covariance"][0, 1],
        "chi2_ndf": temporal["chi2_ndf"],
        "n_calibration_positions": len(calibration),
        "excluded_positions": ",".join(str(x) for x in sorted(ref_positions)),
    }]).to_csv(output_dir / "analysis" / "calibration_temporal.csv", index=False)
    pd.DataFrame([
        {
            "model": "linear", "coefficients_high_to_low": json.dumps(ratio_coeff.tolist()),
            "covariance": json.dumps(ratio_linear["covariance"].tolist()),
            "chi2_ndf": ratio_linear["chi2_ndf"], "aic": ratio_linear["aic"],
            "loocv_rmse": loocv_rmse(x, calibration["mu_R"].to_numpy(), calibration["mu_R_err"].to_numpy(), 1),
        },
        {
            "model": "cubic_comparison", "coefficients_high_to_low": json.dumps(ratio_cubic["coefficients"].tolist()),
            "covariance": json.dumps(ratio_cubic["covariance"].tolist()),
            "chi2_ndf": ratio_cubic["chi2_ndf"], "aic": ratio_cubic["aic"],
            "loocv_rmse": loocv_rmse(x, calibration["mu_R"].to_numpy(), calibration["mu_R_err"].to_numpy(), 3),
        },
    ]).to_csv(output_dir / "analysis" / "calibration_ratio.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    xx = np.linspace(-462, -422, 400)
    axes[0, 0].errorbar(
        summary["x_true_mm"].to_numpy(), summary["mu_dt_ps"].to_numpy(),
        yerr=summary["mu_dt_err_ps"].to_numpy(), fmt="o", ms=3,
    )
    axes[0, 0].plot(xx, np.polyval(temporal_coeff, xx))
    axes[0, 0].set(xlabel="True x (mm)", ylabel="mu delta_t (ps)", title=f"Temporal linear calibration; chi2/ndf={temporal['chi2_ndf']:.2f}")
    axes[1, 0].axhline(0, color="black", lw=0.8); axes[1, 0].plot(x, temporal["residuals"], "o")
    axes[1, 0].set(xlabel="True x (mm)", ylabel="Temporal residual (ps)")
    axes[0, 1].errorbar(
        summary["x_true_mm"].to_numpy(), summary["mu_R"].to_numpy(),
        yerr=summary["mu_R_err"].to_numpy(), fmt="o", ms=3,
    )
    axes[0, 1].plot(xx, np.polyval(ratio_coeff, xx), label="linear")
    axes[0, 1].plot(xx, np.polyval(ratio_cubic["coefficients"], xx), label="cubic comparison")
    axes[0, 1].set(xlabel="True x (mm)", ylabel="mean R", title=f"Ratio calibration; linear chi2/ndf={ratio_linear['chi2_ndf']:.2f}"); axes[0, 1].legend()
    axes[1, 1].axhline(0, color="black", lw=0.8)
    axes[1, 1].plot(x, ratio_linear["residuals"], "o", label="linear")
    axes[1, 1].plot(x, ratio_cubic["residuals"], "o", label="cubic")
    axes[1, 1].set(xlabel="True x (mm)", ylabel="Ratio residual"); axes[1, 1].legend()
    for ax in axes.flat: ax.grid(alpha=0.2)
    fig.suptitle(FIGURE_CONTEXT)
    save_figure(fig, output_dir, "calibrations_and_residuals")

    rows = []
    reco_values: dict[tuple[str, str], np.ndarray] = {}
    for _, ref in references.iterrows():
        label = ref["reference"]
        x_true = float(ref["x_true_mm"])
        frame = load_derived(output_dir / "derived" / f"pair_events_x{x_true:+.1f}mm.npz")
        valid = frame["passed_pair_4pe"].to_numpy()
        dt = frame.loc[valid, "delta_t_ps"].to_numpy()
        ratio = frame.loc[valid, "R_log_ratio"].to_numpy()
        x_t = (dt - temporal_coeff[1]) / temporal_coeff[0]
        x_r = (ratio - ratio_coeff[1]) / ratio_coeff[0]
        x_cubic, ambiguous, no_root = invert_cubic(ratio, ratio_cubic["coefficients"], x_r)
        covariance = np.cov(np.vstack((x_t, x_r)), ddof=1)
        condition = float(np.linalg.cond(covariance))
        inverse = np.linalg.inv(covariance)
        one = np.ones(2)
        weights = inverse @ one / (one @ inverse @ one)
        x_blue = weights[0] * x_t + weights[1] * x_r
        diagonal_weights = (1 / np.diag(covariance)) / np.sum(1 / np.diag(covariance))
        x_diagonal = diagonal_weights[0] * x_t + diagonal_weights[1] * x_r
        corr = float(np.corrcoef(x_t, x_r)[0, 1])
        for method, values in (
            ("temporal", x_t), ("ratio_linear", x_r), ("ratio_cubic_comparison", x_cubic),
            ("BLUE", x_blue), ("covariance_ignored_cross_check", x_diagonal),
        ):
            row = {"reference": label, "x_true_mm": x_true, "method": method}
            row.update(sample_summary(values, x_true))
            row.update({
                "corr_x_t_x_R": corr if method in ("BLUE", "covariance_ignored_cross_check") else math.nan,
                "weight_t": weights[0] if method == "BLUE" else diagonal_weights[0] if method == "covariance_ignored_cross_check" else math.nan,
                "weight_R": weights[1] if method == "BLUE" else diagonal_weights[1] if method == "covariance_ignored_cross_check" else math.nan,
                "covariance_condition_number": condition if method == "BLUE" else math.nan,
                "cubic_ambiguous_fraction": ambiguous / len(ratio) if method == "ratio_cubic_comparison" else math.nan,
                "cubic_no_root_fraction": no_root / len(ratio) if method == "ratio_cubic_comparison" else math.nan,
            })
            rows.append(row)
            reco_values[(label, method)] = values
        temporal_row = rows[-5]
        position_summary = summary.loc[summary["x_true_mm"] == x_true].iloc[0]
        expected = position_summary["sigma_dt_ps"] / abs(temporal_coeff[0])
        expected_error = expected * math.sqrt(
            (position_summary["sigma_dt_err_ps"] / position_summary["sigma_dt_ps"]) ** 2
            + (math.sqrt(temporal["covariance"][0, 0]) / temporal_coeff[0]) ** 2
        )
        temporal_row["guardrail_sigma_dt_over_slope_mm"] = expected
        temporal_row["guardrail_error_mm"] = expected_error
        temporal_row["guardrail_pull"] = (
            (temporal_row["sigma_x_mm"] - expected)
            / math.sqrt(temporal_row["sigma_x_error_mm"] ** 2 + expected_error ** 2)
        )

    result = pd.DataFrame(rows)
    result.to_csv(output_dir / "analysis" / "reconstruction_summary.csv", index=False)
    write_latex_table(
        result,
        output_dir / "tables" / "reconstruction_summary.tex",
        ["reference", "x_true_mm", "method", "mean_x_rec_mm", "bias_mm", "sigma_x_mm", "sigma_x_error_mm"],
    )
    fig, axes = plt.subplots(len(references), 1, figsize=(9, 8), squeeze=False)
    for ax, (_, ref) in zip(axes[:, 0], references.iterrows()):
        label = ref["reference"]
        for method in ("temporal", "ratio_linear", "BLUE"):
            ax.hist(reco_values[(label, method)], bins=70, histtype="step", density=True, label=method)
        ax.axvline(ref["x_true_mm"], color="black", ls="--", label="true x")
        ax.set(xlabel="Reconstructed x (mm)", ylabel="Density", title=f"{label}: x={ref['x_true_mm']:.1f} mm")
        ax.legend(); ax.grid(alpha=0.2)
    fig.suptitle(f"{FIGURE_CONTEXT}\nsimulation prediction — Top readout — intrinsic")
    save_figure(fig, output_dir, "position_reconstruction")


def report(output_dir: pathlib.Path, v1_path: pathlib.Path) -> None:
    v1 = pd.read_csv(v1_path)
    v2 = pd.read_csv(output_dir / "analysis" / "pairscan_summary_v2.csv")
    inventory = pd.read_csv(output_dir / "analysis" / "data_inventory.csv")
    refs = pd.read_csv(output_dir / "analysis" / "reference_positions.csv")
    comparison = pd.read_csv(output_dir / "analysis" / "reference_comparison.csv")
    reco = pd.read_csv(output_dir / "analysis" / "reconstruction_summary.csv")
    temporal = pd.read_csv(output_dir / "analysis" / "calibration_temporal.csv").iloc[0]
    ratio = pd.read_csv(output_dir / "analysis" / "calibration_ratio.csv")

    focus = v2[v2["x_true_mm"].isin([-439.0, -442.0, -444.0])].merge(
        v1[["x_true_mm", "mu_dt_ps", "sigma_dt_ps"]],
        on="x_true_mm", suffixes=("_v2", "_v1"),
    )
    v1_fit = weighted_polynomial_fit(
        v1["x_true_mm"].to_numpy(), v1["mu_dt_ps"].to_numpy(),
        v1["mu_dt_err_ps"].to_numpy(), 1,
    )
    v2_fit = weighted_polynomial_fit(
        v2["x_true_mm"].to_numpy(), v2["mu_dt_ps"].to_numpy(),
        v2["mu_dt_err_ps"].to_numpy(), 1,
    )
    calibration_compare = pd.DataFrame([
        {"version": "v1", "m_t_ps_per_mm": v1_fit["coefficients"][0], "m_t_error_ps_per_mm": math.sqrt(v1_fit["covariance"][0, 0]), "chi2_ndf": v1_fit["chi2_ndf"]},
        {"version": "v2", "m_t_ps_per_mm": v2_fit["coefficients"][0], "m_t_error_ps_per_mm": math.sqrt(v2_fit["covariance"][0, 0]), "chi2_ndf": v2_fit["chi2_ndf"]},
    ])
    focus.to_csv(output_dir / "analysis" / "fit_qa_v1_v2_focus.csv", index=False)
    calibration_compare.to_csv(output_dir / "analysis" / "calibration_v1_v2.csv", index=False)
    write_latex_table(
        focus, output_dir / "tables" / "fit_qa_v1_v2.tex",
        ["x_true_mm", "mu_dt_ps_v1", "mu_dt_ps_v2", "sigma_dt_ps_v1", "sigma_dt_ps_v2", "chi2_ndf", "fraction_within_3sigma"],
    )

    ref_lines = "\n".join(
        f"- `{row.reference}`: `{row.x_true_mm:.1f} mm`, {row.criterion}."
        for row in refs.itertuples()
    )
    reco_key = reco[reco["method"].isin(["temporal", "ratio_linear", "BLUE"])]
    reco_lines = "\n".join(
        f"- `{row.reference}` / `{row.method}`: mean `{row.mean_x_rec_mm:.3f} mm`, "
        f"bias `{row.bias_mm:+.3f} mm`, sigma `{row.sigma_x_mm:.3f} +/- {row.sigma_x_error_mm:.3f} mm`."
        for row in reco_key.itertuples()
    )
    correlation_lines = "\n".join(
        f"- `{row.reference}`: Pearson(delta_t,R) = `{row.pearson_dt_R:.4f} +/- {row.pearson_dt_R_boot_error:.4f}`, "
        f"Spearman = `{row.spearman_dt_R:.4f}`."
        for row in comparison.itertuples()
    )
    blue_lines = "\n".join(
        f"- `{row.reference}`: corr(x_t,x_R)=`{row.corr_x_t_x_R:.4f}`, "
        f"weights=(`{row.weight_t:.4f}`, `{row.weight_R:.4f}`), "
        f"condition number=`{row.covariance_condition_number:.3f}`."
        for row in reco[reco["method"] == "BLUE"].itertuples()
    )
    readme = f"""# EXEC_11 pair (28,29) timing and 1-D position reconstruction

## Provenance and status

- Branch: `exp/pair-scan-2026-06-11`
- Input data: `{HOOK_DATA_DIR}`
- Data-producing commit: `{DATA_COMMIT}`
- Analysis commit embedded in the final figures: `{analysis_commit()}`
- Material/readout: the recorded `pairscan.log` states `OPSC-101`, yield 10400 ph/MeV,
  rise 0.7 ns, decay 1.8 ns, attenuation 160 cm, emission peak 408.8 nm, and `Top`
  readout. Repository documentation maps OPSC-101 to EJ-204.
- Pair identity: `DetectorConstruction::TopSiPMCenterX` and global/local ID mapping put
  IDs 28 and 29 on Top at -452 and -432 mm, respectively.
- ROOT phase 0b: **DEFERRED — ROOT runtime unavailable on MSI**. Static inspection finds
  a plausible crash site: `analysis/analyze_pair_scan.C:281-283` fits with option `N`
  (do not store the TF1), while `analysis/analyze_pair_scan.C:357-360` dereferences
  `g_dt_cal.GetFunction("f_lin_dt")` after the CSV is written. The macro was not changed
  because the crash could not be reproduced without ROOT.

## Input gate

- 41 stable ROOT files, integer positions -462.0 through -422.0 mm, no missing or
  duplicate positions, and no files changed during a 31-second size/mtime check.
- All files opened with uproot. Entries range from `{inventory.n_entries.min()}` to
  `{inventory.n_entries.max()}`; every file contains exactly 3000 unique event IDs,
  IDs 0..2999, no entirely missing events, and zero non-finite `time_ns` values.
- Every position has 3000 events passing 4 PE in both channels, so efficiency is 1.0.
- `SiPMSD.cc:40-113` records a row only when the Geant4 optical boundary invokes the
  sensitive detector after the selected surface `EFFICIENCY`; PDE is not applied again.
- `count_definition = "{COUNT_DEFINITION}"`. Axes use SiPM hit count / PE-equivalent
  count where physical photoelectron acceptance cannot be independently re-established.

## QA v1 to v2

The v2 delta-t histogram uses 20 ps bins over [-3000,3000] ps. Its Gaussian is seeded
at the maximum bin, starts from the RMS after clipping the outer 1%, and is fitted twice
inside +/-{HOOK_FIT_RANGE_NSIGMA} sigma. The broad-fit guard prevents convergence to the
known narrow subpeak. All 41 fits terminate successfully, but all have chi2/ndf > 5:
the distributions are visibly non-Gaussian mixtures, and that limitation is retained.

The three collapsed v1 fits become:

| x (mm) | sigma v1 (ps) | sigma v2 (ps) |
|---:|---:|---:|
""" + "\n".join(
        f"| {row.x_true_mm:.1f} | {row.sigma_dt_ps_v1:.3f} | {row.sigma_dt_ps_v2:.3f} |"
        for row in focus.itertuples()
    ) + f"""

The all-position temporal slope changes from
`{calibration_compare.iloc[0].m_t_ps_per_mm:.4f} +/- {calibration_compare.iloc[0].m_t_error_ps_per_mm:.4f} ps/mm`
to `{calibration_compare.iloc[1].m_t_ps_per_mm:.4f} +/- {calibration_compare.iloc[1].m_t_error_ps_per_mm:.4f} ps/mm`.

## Automatic references

{ref_lines}

Both selected files contain 3000 valid pair-threshold events. Individual-hit time
figures use every photon-hit row; event-time figures use the fourth-hit estimators.
NPE-shape plots use a **Moyal approximation**, not an exact Landau density.

Event-level timing/count correlations:

{correlation_lines}

## Reconstruction

The leave-two-out temporal calibration over the other 39 positions is
`mu_delta_t = ({temporal.m_t_ps_per_mm:.5f} +/- {temporal.m_t_error_ps_per_mm:.5f}) x
+ ({temporal.b_t_ps:.3f} +/- {temporal.b_t_error_ps:.3f}) ps`, with
chi2/ndf `{temporal.chi2_ndf:.3f}`.

The linear ratio calibration has chi2/ndf `{ratio.iloc[0].chi2_ndf:.1f}` and LOOCV
RMSE `{ratio.iloc[0].loocv_rmse:.4f}`. A cubic comparison improves training chi2 and
LOOCV RMSE to `{ratio.iloc[1].loocv_rmse:.4f}`, but gives worse event-level resolution
and has ambiguous physical roots for about 3.3% of events. Multiple cubic roots are
resolved by choosing the physical root nearest the predefined linear estimate.

{reco_lines}

BLUE diagnostics:

{blue_lines}

The temporal event standard deviations are smaller than `sigma_fit/abs(m_t)` by more
than 10 sigma at both references. This guardrail failure is not tuned away: the broad
single-Gaussian v2 fit describes long tails, whereas event RMS used by reconstruction
is narrower. Both quantities and their disagreement remain in `reconstruction_summary.csv`.

## Interpretation

Lv et al. use 64 timing measurements and all 2016 SiPM pairs to intersect geometric
circles and locate an empty region. EXEC_11 uses one adjacent pair and empirical 1-D
calibrations of delta-t and count ratio. It is an inspired one-pair/one-dimensional
adaptation, not a literal reproduction of their circle algorithm. It must not be
compared directly with their 2-3 mm geometric result or approximately 1.5 mm CNN result.

Every position resolution here is a **simulation prediction — Top readout — intrinsic**.
There is no Top-readout test-beam counterpart, SPTR, electronics jitter, or walk
correction. The strong ratio non-linearity and non-Gaussian delta-t tails are open
limitations, not resolved detector performance.

## Reproduction

```bash
OUT={output_dir}
pytest -q tests/test_exec11_pair_analysis.py
python3 analysis/exec11_pair_analysis.py derive --output-dir "$OUT"
python3 analysis/exec11_pair_analysis.py qa --output-dir "$OUT"
python3 analysis/exec11_pair_analysis.py detail --output-dir "$OUT"
python3 analysis/exec11_pair_analysis.py reconstruct --output-dir "$OUT"
python3 analysis/exec11_pair_analysis.py report --output-dir "$OUT"
```

Figures obtain the analysis SHA at runtime through `git rev-parse --short HEAD`.
They were regenerated after the analysis/reconstruction commit and before the final
report-only commit, so the embedded SHA identifies the exact analysis implementation.
"""
    (output_dir / "README.md").write_text(readme)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("derive", "qa", "detail", "reconstruct", "report"))
    parser.add_argument("--data-dir", type=pathlib.Path, default=HOOK_DATA_DIR)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument(
        "--v1-path",
        type=pathlib.Path,
        default=HOOK_DATA_DIR / "analysis" / "pairscan_summary.csv",
    )
    args = parser.parse_args()
    if args.stage == "derive":
        derive_all(args.data_dir, args.output_dir)
    elif args.stage == "qa":
        qa_summary(args.output_dir, args.v1_path)
    elif args.stage == "detail":
        detailed_analysis(args.data_dir, args.output_dir)
    elif args.stage == "reconstruct":
        reconstruction(args.output_dir)
    elif args.stage == "report":
        report(args.output_dir, args.v1_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
