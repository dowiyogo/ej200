#!/usr/bin/env python3
"""EXEC_10 N_pe fluctuation analysis with a Moyal-Gaussian Landau proxy."""

from __future__ import annotations

import argparse
import math
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.stats import moyal, poisson

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from exec07.common import EXPECTED_POSITIONS_MM, N_CHANNELS, N_EVENTS  # noqa: E402

KEY_POSITIONS = (-690, -650, -400, 0, 400, 650, 690)
N_BLOCKS = 20


def counts_matrix(path: pathlib.Path) -> np.ndarray:
    with uproot.open(path) as root_file:
        arrays = root_file["sipm_hits"].arrays(["event_id", "global_id"], library="np")
    combined = arrays["event_id"].astype(np.int64) * N_CHANNELS + arrays["global_id"]
    return np.bincount(combined, minlength=N_EVENTS * N_CHANNELS).reshape(N_EVENTS, N_CHANNELS)


def fano_jackknife(counts: np.ndarray) -> tuple[float, float, float, float]:
    mean = float(np.mean(counts))
    fano = float(np.var(counts) / mean) if mean else math.nan
    replicas = []
    means = []
    event_ids = np.arange(len(counts))
    for omitted in range(N_BLOCKS):
        selected = counts[event_ids % N_BLOCKS != omitted]
        replica_mean = float(np.mean(selected))
        means.append(replica_mean)
        replicas.append(float(np.var(selected) / replica_mean) if replica_mean else math.nan)
    factor = (N_BLOCKS - 1) / N_BLOCKS
    mean_error = math.sqrt(factor * np.sum((np.asarray(means) - np.mean(means)) ** 2))
    fano_error = math.sqrt(factor * np.sum((np.asarray(replicas) - np.mean(replicas)) ** 2))
    return mean, mean_error, fano, fano_error


def poisson_chi2_ndf(counts: np.ndarray) -> float:
    mean = float(np.mean(counts))
    support = np.arange(int(np.min(counts)), int(np.max(counts)) + 1)
    observed = np.bincount(counts.astype(int) - support[0], minlength=len(support)).astype(float)
    expected = poisson.pmf(support, mean) * len(counts)
    valid = expected >= 5.0
    ndf = np.count_nonzero(valid) - 1
    return float(np.sum((observed[valid] - expected[valid]) ** 2 / expected[valid]) / ndf) if ndf > 0 else math.nan


def moyal_gauss_pdf(
    x: np.ndarray, amplitude: float, mpv: float, width: float, sigma_gauss: float,
) -> np.ndarray:
    step = 0.25
    low = min(float(np.min(x)) - 8.0 * (width + sigma_gauss), 0.0)
    high = float(np.max(x)) + 12.0 * (width + sigma_gauss)
    grid = np.arange(low, high + step, step)
    density = moyal.pdf(grid, loc=mpv, scale=width)
    convolved = gaussian_filter1d(density, sigma=max(sigma_gauss / step, 1e-3), mode="constant")
    return amplitude * np.interp(x, grid, convolved, left=0.0, right=0.0)


def fit_moyal_gauss(counts: np.ndarray) -> tuple[dict[str, float], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    q01, q25, q50, q75, q995 = np.quantile(counts, [0.01, 0.25, 0.5, 0.75, 0.995])
    bin_width = max(1.0, (q75 - q25) / 18.0)
    edges = np.arange(max(0.0, q01 - 2.0 * bin_width), q995 + 3.0 * bin_width, bin_width)
    observed, edges = np.histogram(counts, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    errors = np.sqrt(np.maximum(observed, 1.0))
    width0 = max((q75 - q25) / 2.0, 1.0)
    initial = (len(counts) * bin_width, q50 - 0.4 * width0, width0, max(math.sqrt(q50), 1.0))
    upper_scale = max(float(np.max(counts)), 20.0)
    params, covariance = curve_fit(
        moyal_gauss_pdf,
        centers,
        observed,
        p0=initial,
        sigma=errors,
        absolute_sigma=True,
        bounds=(
            (0.0, 0.0, 0.1, 0.1),
            (np.inf, upper_scale * 1.2, upper_scale, upper_scale),
        ),
        maxfev=3000,
    )
    expected = moyal_gauss_pdf(centers, *params)
    valid = expected >= 5.0
    ndf = np.count_nonzero(valid) - len(params)
    chi2_ndf = (
        float(np.sum((observed[valid] - expected[valid]) ** 2 / expected[valid]) / ndf)
        if ndf > 0 else math.nan
    )
    errors_fit = np.sqrt(np.diag(covariance))
    result = {
        "mpv": params[1],
        "mpv_error": errors_fit[1],
        "width": params[2],
        "width_error": errors_fit[2],
        "sigma_gauss": params[3],
        "sigma_gauss_error": errors_fit[3],
        "chi2_ndf": chi2_ndf,
    }
    return result, (edges, observed, expected)


def fit_fano(frame: pd.DataFrame) -> dict[str, float]:
    mean = frame["mean"].to_numpy()
    fano = frame["fano"].to_numpy()
    error = frame["fano_error"].to_numpy()
    weight = 1.0 / np.square(error)
    c = float(np.sum(weight * mean * (fano - 1.0)) / np.sum(weight * mean * mean))
    c_error = float(1.0 / math.sqrt(np.sum(weight * mean * mean)))
    residual = (fano - (1.0 + c * mean)) / error
    chi2 = float(np.sum(residual * residual))
    ndf = len(frame) - 1
    return {"c": c, "c_error": c_error, "chi2_ndf": chi2 / ndf, "n_points": len(frame)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("analysis/exec07"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figs = args.output_dir / "figs"
    figs.mkdir(exist_ok=True)

    all_rows = []
    key_counts: dict[tuple[int, int], np.ndarray] = {}
    for x_mm in EXPECTED_POSITIONS_MM:
        print(f"Fano x={x_mm}", flush=True)
        matrix = counts_matrix(args.data_dir / f"photon_hits_x{x_mm}mm.root")
        for channel in range(N_CHANNELS):
            mean, mean_error, fano, fano_error = fano_jackknife(matrix[:, channel])
            all_rows.append({
                "x_beam_mm": x_mm, "channel": channel, "mean": mean,
                "mean_error": mean_error, "fano": fano, "fano_error": fano_error,
            })
            if x_mm in KEY_POSITIONS:
                key_counts[(x_mm, channel)] = matrix[:, channel]
    fano_frame = pd.DataFrame(all_rows)
    fano_fit = fit_fano(fano_frame[fano_frame["mean"] > 0.5])
    fano_frame.to_csv(args.output_dir / "exec10_fano_by_channel.csv", index=False)
    pd.DataFrame([fano_fit]).to_csv(args.output_dir / "exec10_fano_fit.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    plotted = fano_frame[fano_frame["mean"] > 0.5]
    ax.errorbar(
        plotted["mean"].to_numpy(), plotted["fano"].to_numpy(),
        xerr=plotted["mean_error"].to_numpy(), yerr=plotted["fano_error"].to_numpy(),
        fmt=".", ms=2, alpha=0.22, lw=0.3, label="86 channels x 31 positions",
    )
    x_fit = np.geomspace(plotted["mean"].min(), plotted["mean"].max(), 300)
    ax.plot(
        x_fit, 1.0 + fano_fit["c"] * x_fit, color="black", lw=2,
        label=(
            f"F=1+c< Npe >, c={fano_fit['c']:.5f}+/-{fano_fit['c_error']:.5f}, "
            f"chi2/ndf={fano_fit['chi2_ndf']:.1f}"
        ),
    )
    ax.axhline(1.0, color="tab:red", ls="--", label="Poisson F=1")
    ax.set(xscale="log", yscale="log", xlabel="Mean N_pe per event", ylabel="F = var(N_pe) / mean(N_pe)")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figs / "exec10_fano_vs_mean.png", dpi=220)
    plt.close(fig)

    fit_rows = []
    example = None
    for x_mm in KEY_POSITIONS:
        for channel in range(N_CHANNELS):
            counts = key_counts[(x_mm, channel)]
            mean, mean_error, fano, fano_error = fano_jackknife(counts)
            regime = "poisson-only" if mean < 20.0 else "intermediate-both" if mean <= 50.0 else "moyal-gauss"
            row = {
                "x_beam_mm": x_mm, "channel": channel, "mean": mean,
                "mean_error": mean_error, "fano": fano, "fano_error": fano_error,
                "fit_regime": regime, "poisson_chi2_ndf": poisson_chi2_ndf(counts),
                "mpv": math.nan, "mpv_error": math.nan, "width": math.nan,
                "width_error": math.nan, "sigma_gauss": math.nan,
                "sigma_gauss_error": math.nan, "chi2_ndf": math.nan, "fit_status": "not-requested",
            }
            if mean >= 20.0:
                try:
                    result, curve = fit_moyal_gauss(counts)
                    row.update(result)
                    row["fit_status"] = "ok"
                    if x_mm == 0 and channel == 51:
                        example = (counts, curve, row.copy())
                except Exception as error:  # noqa: BLE001
                    row["fit_status"] = f"failed:{type(error).__name__}"
            fit_rows.append(row)
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(args.output_dir / "exec10_landau_mpv.csv", index=False)

    if example is None:
        candidates = fits[(fits.fit_status == "ok") & (fits["mean"] > 50)].sort_values("mean", ascending=False)
        selected = candidates.iloc[0]
        counts = key_counts[(int(selected.x_beam_mm), int(selected.channel))]
        _, curve = fit_moyal_gauss(counts)
        example = (counts, curve, selected.to_dict())
    counts, (edges, observed, expected), row = example
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.stairs(observed, edges, label="data", color="black")
    ax.plot(centers, expected, lw=2, color="tab:blue", label="Moyal (Landau proxy) convolved with Gaussian")
    ax.set(
        xlabel="Detected N_pe per event", ylabel="Events per N_pe bin",
        title=f"x={int(row['x_beam_mm'])} mm, channel {int(row['channel'])}: chi2/ndf={row['chi2_ndf']:.2f}",
    )
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figs / "exec10_landau_fit_example.png", dpi=220)
    plt.close(fig)

    failed = fits[(fits["mean"] >= 20.0) & (fits.fit_status != "ok")]
    sigma_at_bound = fits[(fits.fit_status == "ok") & (fits.sigma_gauss < 0.11)]
    report = [
        "# EXEC_10 Landau/Moyal N_pe audit", "",
        "The fitted high-N_pe model is scipy.stats.moyal convolved numerically with a Gaussian.",
        "Moyal is a practical analytic approximation to Landau, not the exact ROOT Landau density.",
        "",
        f"- Global F fit: c={fano_fit['c']:.7f} +/- {fano_fit['c_error']:.7f}; chi2/ndf={fano_fit['chi2_ndf']:.2f} over {fano_fit['n_points']} channel-position points.",
        f"- Moyal-Gaussian fits requested for mean >=20 N_pe: {len(fits[fits['mean'] >= 20.0])}; failed: {len(failed)}.",
        f"- High-N_pe fits (mean >50) with chi2/ndf >5: {len(fits[(fits['mean'] > 50.0) & (fits.chi2_ndf > 5.0)])}.",
        f"- Fits with sigma_Gauss at the 0.1 N_pe lower bound: {len(sigma_at_bound)}; sigma_Gauss is not identifiable there, while MPV and Landau width remain reported.",
        "",
        "The CSV is intended as future input for comparison with test-beam Landau MPVs through the pending ToT(NPE) calibration.",
    ]
    (args.output_dir.parent.parent / "audit" / "exec10_landau_analysis.md").write_text("\n".join(report) + "\n")
    print("\n".join(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
