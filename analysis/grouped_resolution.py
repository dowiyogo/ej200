#!/usr/bin/env python3
"""
grouped_resolution.py — sigma_t vs x para diferentes esquemas de agrupación
de canales SiPM (Sum-of-N), emulando la lógica analógica del FastIC+.

Modelos de agrupación implementados:
  - Single channel (referencia): tiempo = primer fotón del SiPM individual
  - Sum-of-4 end:    grupos {0..3, 4..7, 8..11, 12..15}      -> 4 grupos
  - Sum-of-8 end:    grupos {0..7, 8..15}                    -> 2 grupos
  - Sum-of-4 top:    grupos consecutivos en IDs 16..85       -> 18 grupos
  - Sum-of-5 top:    grupos consecutivos en IDs 16..85       -> 14 grupos
  - Full-end (OR):   un solo grupo de los 16 end SiPMs       -> 1 grupo
  - Full-top (OR):   un solo grupo de los 70 top SiPMs       -> 1 grupo

Trigger: el k-ésimo fotón temprano del grupo (k = threshold en p.e.).
Default k=4. Usa k=1 para reproducir FPT.

NOTA DE MODELACIÓN
==================
Este script emula la lógica Sum-of-N del FastIC+ (SOP-3.2 SHiP) a nivel
de tiempos discretos de fotones. Es una APROXIMACIÓN PRIMERA al
comportamiento real, válida bajo los siguientes supuestos:

  1. El pulso del SiPM individual tiene rise-time mucho menor que la
     separación inter-fotón media (válido para nuestros pulsos de ~ns con
     decay 2.1 ns y ~10^4 fotones generados por MIP).
  2. La suma analógica preserva el ordenamiento temporal de los frentes
     leading; en realidad un pulso desplazado puede caer dentro del rise
     del anterior y mover el cruce de threshold (efecto walk).
  3. El threshold se define en fotoelectrones equivalentes; en el ASIC
     real es un voltaje y depende de la ganancia (OV, T).

Para una emulación más fiel se requiere sintetizar formas de onda por
cada SiPM (convolución del tren de tiempos con un kernel SPE) y aplicar
un constant-fraction discriminator sobre la suma. Eso queda para una
fase futura del código (ver issue #TBD).

Uso:
    python analysis/grouped_resolution.py photon_hits.root
    python analysis/grouped_resolution.py photon_hits.root --threshold 4
    python analysis/grouped_resolution.py photon_hits.root --threshold 4 --groupings sum4_end sum8_end full_end

Outputs:
    grouped_resolution_thr<k>.pdf
    grouped_resolution_thr<k>.csv
    grouped_resolution_compare.pdf
"""

import argparse
import glob
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

from resolution_vs_x_fixed import fit_fpt_distribution

TOP_IDS = tuple(range(16, 86))


def chunks(ids, size):
    return [tuple(ids[i:i + size]) for i in range(0, len(ids), size)]


GROUPINGS = {
    "single_end_left":  [(i,) for i in range(0, 8)],
    "single_end_right": [(i,) for i in range(8, 16)],
    "single_top":       [(i,) for i in TOP_IDS],
    "sum4_end":         [tuple(range(0, 4)), tuple(range(4, 8)),
                         tuple(range(8, 12)), tuple(range(12, 16))],
    "sum8_end":         [tuple(range(0, 8)), tuple(range(8, 16))],
    "sum4_top":         chunks(TOP_IDS, 4),
    "sum5_top":         chunks(TOP_IDS, 5),
    "full_end":         [tuple(range(0, 16))],
    "full_top":         [TOP_IDS],
}

COLORS = {
    "single_end_left": "#2166ac",
    "single_end_right": "#d6604d",
    "single_top": "#4dac26",
    "sum4_end": "#542788",
    "sum8_end": "#b2182b",
    "sum4_top": "#1b7837",
    "sum5_top": "#a6dba0",
    "full_end": "#000000",
    "full_top": "#666666",
}


def find_default_root_files(args_files):
    if args_files:
        return args_files
    for pattern in ("photon_hits_merged.root", "photon_hits.root", "photon_hits_run*.root"):
        files = sorted(glob.glob(pattern))
        if files:
            return files
    return []


def load_data_full(root_files, step_size="100 MB"):
    chunks = []
    event_offset = 0
    for path in root_files:
        with uproot.open(path) as f:
            tree = f["sipm_hits"]
            available = set(tree.keys())
            cols = ["event_id", "global_id", "time_ns"]
            if "gun_x_mm" in available:
                cols.append("gun_x_mm")
            max_event = -1
            for chunk in tree.iterate(cols, library="pd", step_size=step_size):
                if len(chunk) == 0:
                    continue
                if "gun_x_mm" not in chunk.columns:
                    chunk["gun_x_mm"] = 0.0
                max_event = max(max_event, int(chunk["event_id"].max()))
                chunk = chunk[["event_id", "gun_x_mm", "global_id", "time_ns"]].copy()
                chunk["event_id"] = chunk["event_id"].astype(np.int64) + event_offset
                chunks.append(chunk)
            event_offset += max_event + 1 if max_event >= 0 else 0
    if not chunks:
        raise RuntimeError("No photon hits found in input ROOT files.")
    df = pd.concat(chunks, ignore_index=True)
    print(f"Loaded {len(df):,} photon hits from {len(root_files)} ROOT file(s)")
    return df


def kth_photon_time(times_ns: np.ndarray, k: int) -> float:
    """Returns the k-th smallest time (1-indexed). Returns NaN if fewer than k photons."""
    if len(times_ns) < k:
        return np.nan
    return float(np.partition(times_ns, k - 1)[k - 1])


def cluster_times_by_event(df, group_ids, threshold_pe):
    """
    Para cada evento, calcula el tiempo del cluster (unión de SiPMs en group_ids)
    como el threshold-ésimo fotón más temprano.
    Retorna DataFrame: event_id, gun_x_mm, t_cluster_ns
    """
    sub = df[df["global_id"].isin(group_ids)]
    if len(sub) == 0:
        return pd.DataFrame(columns=["event_id", "gun_x_mm", "t_cluster_ns"])
    return (sub.groupby(["event_id", "gun_x_mm"])
              .agg(t_cluster_ns=("time_ns", lambda v: kth_photon_time(v.values, threshold_pe)))
              .reset_index()
              .dropna())


def resolution_for_grouping(df, grouping_name, threshold_pe):
    """
    Para un grouping (lista de tuplas de global_id), calcula sigma_t vs x.

    Estrategia: para cada evento, toma el cluster con el TIEMPO MÁS TEMPRANO
    (modelo OR de los grupos del lado relevante). Esto emula que la electrónica
    elige el grupo que disparó primero.
    """
    groups = GROUPINGS[grouping_name]
    cluster_dfs = [cluster_times_by_event(df, g, threshold_pe) for g in groups]
    cluster_dfs = [d for d in cluster_dfs if len(d) > 0]
    if not cluster_dfs:
        return pd.DataFrame()
    all_clusters = pd.concat(cluster_dfs, ignore_index=True)
    earliest = (all_clusters.groupby(["event_id", "gun_x_mm"])
                            .agg(t_trigger=("t_cluster_ns", "min"))
                            .reset_index())

    rows = []
    for x_val, grp in earliest.groupby("gun_x_mm"):
        mu, sigma, sigma_err = fit_fpt_distribution(grp["t_trigger"].values)
        rows.append(dict(x_mm=x_val, sigma_ns=sigma, sigma_err=sigma_err,
                         mu_ns=mu, n_events=len(grp)))
    return pd.DataFrame(rows).sort_values("x_mm").reset_index(drop=True)


def describe_point(df, x_target):
    if len(df) == 0:
        return np.nan
    sub = df[df["sigma_ns"].notna() & (df["sigma_ns"] > 0)]
    if len(sub) == 0:
        return np.nan
    idx = (sub["x_mm"] - x_target).abs().idxmin()
    return float(sub.loc[idx, "sigma_ns"] * 1e3)


def plot_grouped_resolution(results, threshold_pe, out_pdf):
    fig, ax = plt.subplots(figsize=(11, 5.8))
    for name, df in results.items():
        if df is None or len(df) == 0:
            continue
        ok = df["sigma_ns"].notna() & (df["sigma_ns"] > 0)
        sub = df[ok]
        if len(sub) == 0:
            continue
        central = sub.loc[sub["x_mm"].abs() < 500, "sigma_ns"]
        central_ps = central.mean() * 1e3 if len(central) else np.nan
        sigma650_ps = describe_point(sub, 650.0)
        label = f"{name} (central {central_ps:.0f} ps, x=650 {sigma650_ps:.0f} ps)"
        ax.errorbar(sub["x_mm"].to_numpy(), (sub["sigma_ns"] * 1e3).to_numpy(),
                    yerr=(sub["sigma_err"] * 1e3).to_numpy(),
                    marker="o", ms=4, lw=1.5, capsize=2,
                    color=COLORS.get(name), label=label)

    ax.axhline(50, color="black", ls="--", lw=1.0, label="SHiP goal 50 ps")
    ax.axhline(100, color="black", ls=":", lw=1.2, label="Betancourt et al. 2020 target 100 ps")
    ax.set_xlabel("Muon x [mm]")
    ax.set_ylabel(r"Timing resolution $\sigma_t$ [ps]")
    ax.set_title(f"Grouped SiPM timing resolution (threshold = {threshold_pe} p.e.)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"  -> {out_pdf}")


def plot_compare(results, out_pdf):
    keys = ["single_end_left", "sum4_end", "sum8_end", "full_end",
            "single_top", "sum4_top", "full_top"]
    subset = {k: results[k] for k in keys if k in results and len(results[k]) > 0}
    if not subset:
        return
    plot_grouped_resolution(subset, "compare", out_pdf)


def save_csv(results, out_csv):
    rows = []
    for name, df in results.items():
        if df is None or len(df) == 0:
            continue
        out = df.copy()
        out["grouping"] = name
        rows.append(out[["grouping", "x_mm", "mu_ns", "sigma_ns", "sigma_err", "n_events"]])
    if rows:
        pd.concat(rows, ignore_index=True).to_csv(out_csv, index=False)
        print(f"  -> {out_csv}")


def print_summary(results, threshold_pe):
    print(f"\nSummary for threshold = {threshold_pe} p.e.")
    print(f"{'grouping':<18} {'central [ps]':>14} {'x=650 [ps]':>12} {'n x bins':>8}")
    for name, df in results.items():
        if df is None or len(df) == 0:
            continue
        ok = df["sigma_ns"].notna() & (df["sigma_ns"] > 0)
        sub = df[ok]
        central = sub.loc[sub["x_mm"].abs() < 500, "sigma_ns"]
        central_ps = central.mean() * 1e3 if len(central) else np.nan
        sigma650_ps = describe_point(sub, 650.0)
        print(f"{name:<18} {central_ps:>14.1f} {sigma650_ps:>12.1f} {len(sub):>8}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_files", default=None, nargs="*",
                        help="ROOT file(s), default: photon_hits_merged.root, photon_hits.root, or photon_hits_run*.root")
    parser.add_argument("--threshold", type=int, default=4,
                        help="Threshold in photoelectrons (k-th photon trigger). Default: 4")
    parser.add_argument("--groupings", nargs="+",
                        default=["single_end_left", "sum4_end", "sum8_end", "full_end",
                                 "single_top", "sum4_top", "full_top"])
    parser.add_argument("--out-prefix", default="grouped_resolution")
    parser.add_argument("--step-size", default="100 MB")
    args = parser.parse_args()

    if args.threshold < 1:
        print("ERROR: --threshold must be >= 1")
        sys.exit(1)

    root_files = find_default_root_files(args.root_files)
    if not root_files:
        print("ERROR: no ROOT files found.")
        sys.exit(1)

    df = load_data_full(root_files, step_size=args.step_size)

    results = {}
    for grouping in args.groupings:
        if grouping not in GROUPINGS:
            print(f"  [skip] unknown grouping: {grouping}")
            continue
        print(f"  Computing {grouping} (threshold={args.threshold} p.e.)...")
        results[grouping] = resolution_for_grouping(df, grouping, args.threshold)

    plot_grouped_resolution(results, args.threshold,
                            f"{args.out_prefix}_thr{args.threshold}.pdf")
    save_csv(results, f"{args.out_prefix}_thr{args.threshold}.csv")
    plot_compare(results, "grouped_resolution_compare.pdf")
    print_summary(results, args.threshold)


if __name__ == "__main__":
    main()
