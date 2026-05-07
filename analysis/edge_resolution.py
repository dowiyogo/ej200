#!/usr/bin/env python3
"""
edge_resolution.py — Resolución temporal y light yield en la región de borde.

Lee photon_hits_run*.root del scan_edge.mac y produce:
  - sigma_t [ps] vs x ∈ [600, 700] mm (paso 10 mm) — end SiPMs y top SiPMs
  - <N_pe> por evento vs x
  - Fracción de eventos sin hits en uno o ambos extremos
  - Distribución del FPT en x = 650, 670, 690, 700 mm (cuatro paneles)
  - Localización del 'breakdown point' (x donde sigma_t crece >50% sobre el centro)

Outputs:
  edge_resolution.pdf
  edge_lightyield.pdf
  edge_fpt_panels.pdf
  edge_summary.csv  (columnas: x_mm, sigma_end_ps, sigma_top_ps, npe_left, npe_right, npe_top, frac_dead_events)
"""

import argparse
import glob
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from resolution_vs_x_fixed import fit_fpt_distribution, gauss, load_event_level_data

FACE_LEFT = 0
FACE_RIGHT = 1
FACE_TOP = 2


def default_root_files():
    files = sorted(glob.glob("photon_hits_run*.root"))
    if files:
        return files
    for name in ("photon_hits_merged.root", "photon_hits.root"):
        if pathlib.Path(name).exists():
            return [name]
    return []


def load_edge_event_data(root_files, step_size="100 MB"):
    chunks = []
    event_offset = 0
    for path in root_files:
        df = load_event_level_data(path, step_size=step_size)
        if "event_id" in df.columns and len(df) > 0:
            df = df.copy()
            df["event_id"] = df["event_id"].astype(np.int64) + event_offset
            event_offset = int(df["event_id"].max()) + 1
        chunks.append(df)
    if not chunks:
        raise RuntimeError("No ROOT input files found.")
    return pd.concat(chunks, ignore_index=True)


def fpt_for_faces(event_face, faces):
    sub = event_face[event_face["face_type"].isin(faces)]
    if len(sub) == 0:
        return pd.DataFrame(columns=["event_id", "gun_x_mm", "fpt_ns"])
    return (sub.groupby(["event_id", "gun_x_mm"], sort=False)
               .agg(fpt_ns=("fpt_ns", "min"))
               .reset_index())


def resolution_vs_x(event_face, faces):
    fpt = fpt_for_faces(event_face, faces)
    rows = []
    for x_val, grp in fpt.groupby("gun_x_mm", sort=True):
        mu, sigma, sigma_err = fit_fpt_distribution(grp["fpt_ns"].to_numpy())
        rows.append({
            "x_mm": float(x_val),
            "mu_ns": mu,
            "sigma_ns": sigma,
            "sigma_err_ns": sigma_err,
            "n_events_with_hits": len(grp),
        })
    return pd.DataFrame(rows)


def lightyield_vs_x(event_face):
    event_counts = (event_face.pivot_table(index=["event_id", "gun_x_mm"],
                                           columns="face_type",
                                           values="n_ph",
                                           aggfunc="sum",
                                           fill_value=0)
                              .reset_index())
    event_counts.columns.name = None
    for face in (FACE_LEFT, FACE_RIGHT, FACE_TOP):
        if face not in event_counts.columns:
            event_counts[face] = 0
    event_counts = event_counts.rename(columns={
        FACE_LEFT: "npe_left",
        FACE_RIGHT: "npe_right",
        FACE_TOP: "npe_top",
    })
    prof = (event_counts.groupby("gun_x_mm")[["npe_left", "npe_right", "npe_top"]]
                        .mean()
                        .reset_index()
                        .rename(columns={"gun_x_mm": "x_mm"}))
    return prof, event_counts


def dead_fraction_by_x(event_counts, expected_events_per_position, expected_x=None):
    rows = []
    if expected_x is None:
        x_values = sorted(event_counts["gun_x_mm"].unique())
    else:
        x_values = list(expected_x)

    grouped = dict(tuple(event_counts.groupby("gun_x_mm", sort=True))) if len(event_counts) else {}
    for x_val in x_values:
        grp = grouped.get(x_val, pd.DataFrame())
        n_missing = expected_events_per_position - len(grp)
        dead_known = 0
        if len(grp) > 0:
            dead_known = ((grp["npe_left"] <= 0) | (grp["npe_right"] <= 0)).sum()
        dead = max(0, n_missing) + int(dead_known)
        rows.append({
            "x_mm": float(x_val),
            "frac_dead_events": dead / float(expected_events_per_position),
        })
    return pd.DataFrame(rows)


def find_breakdown(res_end):
    ok = res_end["sigma_ns"].notna() & (res_end["sigma_ns"] > 0)
    df = res_end[ok].sort_values("x_mm")
    if len(df) == 0:
        return np.nan, np.nan
    centre = df.iloc[0]["sigma_ns"]
    threshold = 1.5 * centre
    hit = df[df["sigma_ns"] > threshold]
    if len(hit) == 0:
        return np.nan, threshold
    return float(hit.iloc[0]["x_mm"]), float(threshold)


def plot_resolution(res_end, res_top, breakdown_x, threshold_ns, out_pdf):
    fig, ax = plt.subplots(figsize=(9, 5))
    for label, df, color, marker in [
        ("End SiPMs", res_end, "#542788", "o"),
        ("Top SiPMs", res_top, "#4dac26", "s"),
    ]:
        if len(df) == 0:
            continue
        ok = df["sigma_ns"].notna() & (df["sigma_ns"] > 0)
        sub = df[ok]
        ax.errorbar(sub["x_mm"].to_numpy(), (sub["sigma_ns"] * 1e3).to_numpy(),
                    yerr=(sub["sigma_err_ns"] * 1e3).to_numpy(),
                    marker=marker, lw=1.8, ms=5, capsize=3,
                    color=color, label=label)
    if not np.isnan(breakdown_x):
        ax.axvline(breakdown_x, color="crimson", ls="--", lw=1.2,
                   label=f"Breakdown x = {breakdown_x:.0f} mm")
    if not np.isnan(threshold_ns):
        ax.axhline(threshold_ns * 1e3, color="crimson", ls=":", lw=1.0,
                   label="+50% over first edge point")
    ax.set_xlabel("Muon x [mm]")
    ax.set_ylabel(r"Timing resolution $\sigma_t$ [ps]")
    ax.set_title("Edge timing resolution")
    ax.grid(alpha=0.3)
    ax.legend(framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"  -> {out_pdf}")


def plot_lightyield(summary, out_pdf):
    fig, ax = plt.subplots(figsize=(9, 5))
    styles = {
        "npe_left": ("End left", "#2166ac", "^"),
        "npe_right": ("End right", "#d6604d", "v"),
        "npe_top": ("Top", "#4dac26", "s"),
    }
    for col, (label, color, marker) in styles.items():
        ax.plot(summary["x_mm"].to_numpy(), summary[col].to_numpy(), marker=marker, lw=1.7,
                color=color, label=label)
    ax2 = ax.twinx()
    ax2.plot(summary["x_mm"].to_numpy(), summary["frac_dead_events"].to_numpy(), color="black",
             marker="x", lw=1.3, ls="--", label="Dead event fraction")
    ax.set_xlabel("Muon x [mm]")
    ax.set_ylabel(r"$\langle N_{pe}\rangle$ per event")
    ax2.set_ylabel("Fraction with missing end hits")
    ax.set_title("Edge light yield and dead-event fraction")
    ax.grid(alpha=0.3)
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"  -> {out_pdf}")


def plot_fpt_panels(event_face, x_targets, out_pdf):
    fpt = fpt_for_faces(event_face, [FACE_LEFT, FACE_RIGHT])
    if len(fpt) == 0:
        return
    x_vals = sorted(fpt["gun_x_mm"].unique())
    chosen = [min(x_vals, key=lambda v: abs(v - target)) for target in x_targets]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharey=False)
    for ax, x_val in zip(axes.flat, chosen):
        times = fpt.loc[fpt["gun_x_mm"] == x_val, "fpt_ns"].to_numpy()
        if len(times) == 0:
            ax.set_title(f"x = {x_val:.0f} mm")
            continue
        mu, sigma, sigma_err = fit_fpt_distribution(times, n_bins=30)
        lo, hi = np.percentile(times, [1, 99])
        margin = max(0.5 * (hi - lo), 0.5)
        counts, edges, _ = ax.hist(times, bins=30, range=(lo - margin, hi + margin),
                                   color="#808080", alpha=0.7)
        if not np.isnan(sigma):
            xx = np.linspace(lo - margin, hi + margin, 300)
            binw = edges[1] - edges[0]
            ax.plot(xx, gauss(xx, mu, sigma, len(times) * binw),
                    color="crimson", lw=1.8,
                    label=fr"$\sigma$={sigma*1e3:.0f}$\pm${sigma_err*1e3:.0f} ps")
        ax.set_title(f"x = {x_val:.0f} mm")
        ax.set_xlabel("FPT [ns]")
        ax.set_ylabel("Events / bin")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9)
    fig.suptitle("End-SiPM first photon time near the edge")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"  -> {out_pdf}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_files", nargs="*", help="ROOT files (default: photon_hits_run*.root)")
    parser.add_argument("--label", default="edge", help="Label stored in CSV metadata")
    parser.add_argument("--events-per-position", type=int, default=500)
    parser.add_argument("--x-min", type=float, default=600.0,
                        help="Expected scan minimum x [mm], used to keep fully dead positions in the CSV")
    parser.add_argument("--x-max", type=float, default=700.0,
                        help="Expected scan maximum x [mm], used to keep fully dead positions in the CSV")
    parser.add_argument("--x-step", type=float, default=10.0,
                        help="Expected scan x step [mm], used to keep fully dead positions in the CSV")
    parser.add_argument("--step-size", default="100 MB")
    args = parser.parse_args()

    root_files = args.root_files or default_root_files()
    expected_x = np.arange(args.x_min, args.x_max + 0.5 * args.x_step, args.x_step)
    if not root_files:
        print("WARNING: no ROOT files found; writing all expected positions as fully dead.")
        summary = pd.DataFrame({
            "x_mm": expected_x,
            "sigma_end_ps": np.nan,
            "sigma_top_ps": np.nan,
            "npe_left": 0.0,
            "npe_right": 0.0,
            "npe_top": 0.0,
            "frac_dead_events": 1.0,
        })
        summary.to_csv("edge_summary.csv", index=False)
        print("  -> edge_summary.csv")
        return

    print(f"Loading {len(root_files)} ROOT file(s)")
    event_face = load_edge_event_data(root_files, step_size=args.step_size)

    res_end = resolution_vs_x(event_face, [FACE_LEFT, FACE_RIGHT])
    res_top = resolution_vs_x(event_face, [FACE_TOP])
    light, event_counts = lightyield_vs_x(event_face)
    dead = dead_fraction_by_x(event_counts, args.events_per_position, expected_x=expected_x)

    observed_x = sorted(event_face["gun_x_mm"].unique())
    x_values = sorted(set(float(x) for x in expected_x).union(float(x) for x in observed_x))
    summary = pd.DataFrame({"x_mm": x_values})
    for df, cols in [
        (res_end.rename(columns={"sigma_ns": "sigma_end_ns"}), ["x_mm", "sigma_end_ns"]),
        (res_top.rename(columns={"sigma_ns": "sigma_top_ns"}), ["x_mm", "sigma_top_ns"]),
        (light, ["x_mm", "npe_left", "npe_right", "npe_top"]),
        (dead, ["x_mm", "frac_dead_events"]),
    ]:
        summary = summary.merge(df[cols], on="x_mm", how="left")
    summary["sigma_end_ps"] = summary["sigma_end_ns"] * 1e3
    summary["sigma_top_ps"] = summary["sigma_top_ns"] * 1e3
    for col in ["npe_left", "npe_right", "npe_top"]:
        summary[col] = summary[col].fillna(0.0)
    summary["frac_dead_events"] = summary["frac_dead_events"].fillna(1.0)
    summary["label"] = args.label

    breakdown_x, threshold_ns = find_breakdown(res_end)
    print(f"Breakdown point: {breakdown_x if not np.isnan(breakdown_x) else 'not found'}")

    plot_resolution(res_end, res_top, breakdown_x, threshold_ns, "edge_resolution.pdf")
    plot_lightyield(summary, "edge_lightyield.pdf")
    plot_fpt_panels(event_face, [650, 670, 690, 700], "edge_fpt_panels.pdf")

    out_cols = ["x_mm", "sigma_end_ps", "sigma_top_ps",
                "npe_left", "npe_right", "npe_top", "frac_dead_events"]
    summary[out_cols].to_csv("edge_summary.csv", index=False)
    print("  -> edge_summary.csv")


if __name__ == "__main__":
    main()
