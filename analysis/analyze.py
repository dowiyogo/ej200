#!/usr/bin/env python3
"""
analyze.py — análisis de photon_hits.root generado por ej200_v2
Requiere: uproot, numpy, matplotlib, pandas

Uso:
    python analysis/analyze.py                     # lee photon_hits.root
    python analysis/analyze.py mi_archivo.root
"""

import sys
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import uproot
import pandas as pd

# ── Constantes de geometría (deben coincidir con DetectorConstruction) ────────
N_END_SIPMS  = 8    # por lado
N_TOP_SIPMS  = 20
BAR_HALF_X   = 700  # mm
TOP_POSITIONS_MM = np.linspace(-665, 665, N_TOP_SIPMS)  # centros x de top SiPMs

FACE_NAMES = {0: "End left (−X)", 1: "End right (+X)", 2: "Top (+Y face)"}
COLORS     = {0: "#2166ac",        1: "#d6604d",        2: "#4dac26"}
SENSOR_NAMES = {0: "Broadcom NUV-MT 14M", 1: "Hamamatsu S13360"}


def compute_dcfd_time(times_ns: np.ndarray, fraction: float = 0.14) -> float:
    """dCFD sobre timestamps de fotoelectrones con interpolacion lineal."""
    times = np.sort(np.asarray(times_ns, dtype=float))
    n = len(times)
    if n == 0:
        return np.nan

    threshold = fraction * n
    counts = np.arange(1, n + 1, dtype=float)
    idx = np.searchsorted(counts, threshold, side="left")

    if idx <= 0:
        return float(times[0])
    if np.isclose(counts[idx], threshold):
        return float(times[idx])

    t0, t1 = times[idx - 1], times[idx]
    y0, y1 = counts[idx - 1], counts[idx]
    if np.isclose(t1, t0):
        return float(t1)
    return float(t0 + (threshold - y0) * (t1 - t0) / (y1 - y0))


def build_dcfd_table(df: pd.DataFrame, electronics_sigma_ps: float = 30.0,
                     fraction: float = 0.14, seed: int = 12345) -> pd.DataFrame:
    """Reconstruye el timestamp por evento y cara usando dCFD al 14%."""
    rows = []
    rng = np.random.default_rng(seed)

    group_cols = ["event_id", "face_type"]
    if "gun_x_mm" in df.columns:
        group_cols.append("gun_x_mm")

    for keys, grp in df.groupby(group_cols):
        if len(group_cols) == 3:
            event_id, face_type, gun_x_mm = keys
        else:
            event_id, face_type = keys
            gun_x_mm = 0.0

        dcfd_time = compute_dcfd_time(grp["time_ns"].to_numpy(), fraction=fraction)
        if np.isnan(dcfd_time):
            continue

        elec_jitter_ns = rng.normal(0.0, electronics_sigma_ps * 1e-3)
        rows.append({
            "event_id": event_id,
            "face_type": face_type,
            "gun_x_mm": gun_x_mm,
            "n_pe": len(grp),
            "dcfd_time_ns": dcfd_time,
            "electronics_jitter_ns": elec_jitter_ns,
            "measured_time_ns": dcfd_time + elec_jitter_ns,
        })

    return pd.DataFrame(rows)

# ─────────────────────────────────────────────────────────────────────────────
def load_data(root_file: str) -> pd.DataFrame:
    with uproot.open(root_file) as f:
        tree = f["sipm_hits"]
        # Columnas base — siempre presentes
        base_cols = ["event_id", "face_type", "global_id", "local_id",
                     "time_ns", "time_raw_ns", "sptr_jitter_ns",
                     "energy_eV", "wl_nm", "pde",
                     "x_mm", "y_mm", "z_mm"]
        # gun_x_mm se agrego en v3; toleramos ntuples viejos sin ella
        available = tree.keys()
        extra_cols = [
            col for col in
            ["gun_x_mm", "sensor_model_id", "sensor_sptr_sigma_ns", "electronics_sigma_ns"]
            if col in available
        ]
        cols = [col for col in base_cols if col in available] + extra_cols
        df = tree.arrays(cols, library="pd")
        if "gun_x_mm" not in df.columns:
            df["gun_x_mm"] = 0.0   # placeholder para compatibilidad
        if "time_raw_ns" not in df.columns:
            df["time_raw_ns"] = df["time_ns"]
        if "sptr_jitter_ns" not in df.columns:
            df["sptr_jitter_ns"] = 0.0
    print(f"Loaded {len(df):,} detected photons from '{root_file}'")
    return df


def photons_per_sipm(df: pd.DataFrame) -> pd.DataFrame:
    """Total detected photons per SiPM (global_id)."""
    return (df.groupby(["global_id", "face_type", "local_id"])
              .size()
              .reset_index(name="n_photons"))


def plot_photons_per_sipm(df: pd.DataFrame, out: str = "photons_per_sipm.pdf"):
    """Bar chart: detected photons vs SiPM ID, coloured by face type."""
    pps = photons_per_sipm(df)

    fig, ax = plt.subplots(figsize=(14, 5))
    for face, grp in pps.groupby("face_type"):
        ax.bar(grp["global_id"], grp["n_photons"],
               color=COLORS[face], label=FACE_NAMES[face], alpha=0.85)

    ax.set_xlabel("Global SiPM ID", fontsize=12)
    ax.set_ylabel("Detected photons", fontsize=12)
    ax.set_title("Detected photons per SiPM", fontsize=13)
    ax.legend()
    ax.grid(axis="y", alpha=0.4)

    # Separator lines between face groups
    ax.axvline(N_END_SIPMS - 0.5,     color="gray", ls="--", lw=0.8)
    ax.axvline(2*N_END_SIPMS - 0.5,   color="gray", ls="--", lw=0.8)

    ax.set_xticks(pps["global_id"].values)
    ax.tick_params(axis="x", labelsize=7)

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print(f"  → {out}")
    plt.close(fig)


def plot_top_sipm_profile(df: pd.DataFrame, out: str = "top_sipm_profile.pdf"):
    """
    Longitudinal profile: photons detected by each top SiPM vs its x-position.
    This is the key plot to validate that ALL top SiPMs collect light.
    """
    top = df[df["face_type"] == 2].copy()
    pps = (top.groupby("global_id").size().reset_index(name="n_photons"))
    pps["local_id"] = pps["global_id"] - 2*N_END_SIPMS
    pps["x_mm"]     = TOP_POSITIONS_MM[pps["local_id"].values]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(pps["x_mm"], pps["n_photons"], width=60, color=COLORS[2], alpha=0.85)
    ax.set_xlabel("SiPM x position [mm]", fontsize=12)
    ax.set_ylabel("Detected photons", fontsize=12)
    ax.set_title("Top SiPM longitudinal profile (full bar length)", fontsize=13)
    ax.set_xlim(-BAR_HALF_X - 20, BAR_HALF_X + 20)
    ax.axvline(0, color="k", ls=":", lw=0.8, label="Bar centre")
    ax.grid(axis="y", alpha=0.4)
    ax.legend()

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print(f"  → {out}")
    plt.close(fig)


def plot_arrival_time(df: pd.DataFrame, out: str = "arrival_time.pdf"):
    """Arrival-time distribution per face type with SPTR-smearing included."""
    fig, ax = plt.subplots(figsize=(9, 4))
    tmax = df["time_ns"].quantile(0.995)

    for face in [0, 1, 2]:
        sub = df[df["face_type"] == face]["time_ns"]
        if len(sub) == 0:
            continue
        ax.hist(sub, bins=80, range=(0, tmax),
                color=COLORS[face], label=FACE_NAMES[face],
                histtype="step", linewidth=1.5)

    ax.set_xlabel("Arrival time [ns]", fontsize=12)
    ax.set_ylabel("Counts / bin", fontsize=12)
    ax.set_title("Optical photon arrival-time spectrum", fontsize=13)
    ax.legend()
    ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print(f"  → {out}")
    plt.close(fig)


def plot_dcfd_time(df: pd.DataFrame, out: str = "dcfd_time.pdf"):
    """Distribucion del tiempo reconstruido por dCFD al 14%."""
    dcfd = build_dcfd_table(df)
    if dcfd.empty:
        print("  [skip] dcfd_time: no hay eventos con fotoelectrones")
        return

    fig, ax = plt.subplots(figsize=(9, 4))
    for face in [0, 1, 2]:
        sub = dcfd[dcfd["face_type"] == face]["measured_time_ns"]
        if len(sub) == 0:
            continue
        ax.hist(sub, bins=60, color=COLORS[face], label=FACE_NAMES[face],
                histtype="step", linewidth=1.5)

    ax.set_xlabel("dCFD time + electronics jitter [ns]", fontsize=12)
    ax.set_ylabel("Events / bin", fontsize=12)
    ax.set_title("Muon timestamp reconstructed with 14% dCFD", fontsize=13)
    ax.legend()
    ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print(f"  → {out}")
    plt.close(fig)


def plot_wavelength(df: pd.DataFrame, out: str = "wavelength.pdf"):
    """Wavelength spectrum of detected photons (weighted by PDE already applied)."""
    fig, ax = plt.subplots(figsize=(9, 4))

    for face in [0, 1, 2]:
        sub = df[df["face_type"] == face]["wl_nm"]
        if len(sub) == 0:
            continue
        ax.hist(sub, bins=60, range=(350, 600),
                color=COLORS[face], label=FACE_NAMES[face],
                histtype="stepfilled", alpha=0.5)

    ax.set_xlabel("Wavelength [nm]", fontsize=12)
    ax.set_ylabel("Detected photons / bin", fontsize=12)
    ax.set_title("Wavelength spectrum of detected photons", fontsize=13)
    ax.legend()
    ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print(f"  → {out}")
    plt.close(fig)


def plot_end_asymmetry(df: pd.DataFrame, out: str = "end_asymmetry.pdf"):
    """
    End-SiPM asymmetry as a function of muon x-position.
    Useful for position reconstruction using the ratio (L-R)/(L+R).

    This plot requires events from multiple x positions (scan.mac).
    If all events are at x=0, the plot will show a flat distribution.
    """
    end = df[df["face_type"].isin([0, 1])].copy()
    if len(end) == 0:
        print("  [skip] end_asymmetry: no end-SiPM hits")
        return

    # Group by event and face, count photons
    evts = end.groupby(["event_id", "face_type"]).size().unstack(fill_value=0)
    evts.columns.name = None
    evts = evts.rename(columns={0: "left", 1: "right"})
    total = evts["left"] + evts["right"]
    evts = evts[total > 0]
    asym = (evts["left"] - evts["right"]) / (evts["left"] + evts["right"])

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(asym, bins=50, range=(-1, 1), color="#7b2d8b", alpha=0.8)
    ax.set_xlabel("Asymmetry (L−R)/(L+R)", fontsize=12)
    ax.set_ylabel("Events / bin", fontsize=12)
    ax.set_title("End-SiPM charge asymmetry", fontsize=13)
    ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print(f"  → {out}")
    plt.close(fig)


def print_summary(df: pd.DataFrame):
    n_evt = df["event_id"].nunique()
    total = len(df)
    print(f"\n{'='*55}")
    print(f"  Total detected photons : {total:>10,}")
    print(f"  Events with ≥1 photon  : {n_evt:>10,}")
    print(f"  Mean photons / event   : {total/n_evt:>10.1f}")
    print("-"*55)
    for face, grp in df.groupby("face_type"):
        print(f"  {FACE_NAMES[face]:<22}: {len(grp):>8,}  "
              f"({100*len(grp)/total:.1f}%)")
    if "sensor_model_id" in df.columns:
        sensor_id = int(df["sensor_model_id"].mode().iloc[0])
        sptr_sigma = df.get("sensor_sptr_sigma_ns", pd.Series([np.nan])).iloc[0]
        elec_sigma = df.get("electronics_sigma_ns", pd.Series([np.nan])).iloc[0]
        print("-"*55)
        print(f"  Sensor model          : {SENSOR_NAMES.get(sensor_id, sensor_id)}")
        print(f"  SPTR sigma            : {sptr_sigma*1e3:.1f} ps")
        print(f"  Electronics sigma     : {elec_sigma*1e3:.1f} ps")
    print("="*55 + "\n")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    root_file = sys.argv[1] if len(sys.argv) > 1 else "photon_hits.root"

    if not pathlib.Path(root_file).exists():
        print(f"ERROR: file '{root_file}' not found.")
        print("Run the simulation first: ./ej200_bar_sim -m macros/run.mac")
        sys.exit(1)

    df = load_data(root_file)
    print_summary(df)

    print("Generating plots…")
    plot_photons_per_sipm   (df)
    plot_top_sipm_profile   (df)
    plot_arrival_time       (df)
    plot_dcfd_time          (df)
    plot_wavelength         (df)
    plot_end_asymmetry      (df)
    print("Done.")
