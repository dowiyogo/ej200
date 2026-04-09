
#!/usr/bin/env python3
"""
resolution_vs_x_fixed.py — versión robusta en memoria para scans combinados.

Cambios principales respecto a resolution_vs_x.py:
  - NO carga todo el TTree en pandas de una sola vez.
  - Lee por chunks con uproot.iterate.
  - Reduce el ntuple a nivel de evento/cara:
        (event_id, gun_x_mm, face_type) -> n_ph, fpt_ns
  - Luego construye todas las figuras desde ese resumen compacto.
  - Agrupa por (event_id, gun_x_mm), no sólo por event_id.

Uso:
    python resolution_vs_x_fixed.py
    python resolution_vs_x_fixed.py photon_hits_merged.root
    python resolution_vs_x_fixed.py --out resolution_vs_x.pdf
"""

import sys
import argparse
from typing import Optional
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.optimize import curve_fit
import uproot

# ── Constantes de geometria (deben coincidir con DetectorConstruction.hh) ────
N_END_SIPMS = 8          # por lado
N_TOP_SIPMS = 20
BAR_HALF_X  = 700.0      # mm

COLORS = {
    "end_left":  "#2166ac",
    "end_right": "#d6604d",
    "end_both":  "#542788",
    "top":       "#4dac26",
}


def gauss(x, mu, sigma, A):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_fpt_distribution(times_ns, n_bins=40):
    """
    Ajusta una Gaussiana a la distribucion de FPT.
    Retorna (mu [ns], sigma [ns], sigma_err [ns]).
    Retorna (nan, nan, nan) si hay muy pocos eventos.
    """
    times_ns = np.asarray(times_ns, dtype=float)
    if len(times_ns) < 10:
        return np.nan, np.nan, np.nan

    lo, hi = np.percentile(times_ns, [1, 99])
    margin = max(0.5 * (hi - lo), 0.5)
    lo -= margin
    hi += margin

    counts, edges = np.histogram(times_ns, bins=n_bins, range=(lo, hi))
    centres = 0.5 * (edges[:-1] + edges[1:])

    mu0 = float(np.mean(times_ns))
    sigma0 = float(np.std(times_ns))
    if sigma0 < 1e-6:
        sigma0 = 0.1
    A0 = float(max(counts.max(), 1))

    try:
        popt, pcov = curve_fit(
            gauss, centres, counts,
            p0=[mu0, sigma0, A0],
            bounds=([lo, 0, 0], [hi, hi - lo, 10 * A0]),
            maxfev=5000
        )
        perr = np.sqrt(np.diag(pcov))
        mu_fit, sigma_fit, _ = popt
        sigma_err = perr[1]
        return mu_fit, abs(sigma_fit), sigma_err
    except Exception:
        # Fallback robusto
        sigma = float(np.std(times_ns))
        return mu0, sigma, sigma / np.sqrt(len(times_ns))


def find_default_root_file(user_arg: Optional[str]) -> str:
    """
    Si el usuario no pasa archivo, probamos primero el merged y luego el clásico.
    """
    if user_arg:
        return user_arg

    candidates = [
        "photon_hits_merged.root",
        "photon_hits.root",
    ]
    for c in candidates:
        if pathlib.Path(c).exists():
            return c
    return candidates[0]


def load_event_level_data(root_file: str, step_size="100 MB") -> pd.DataFrame:
    """
    Lee el ntuple por chunks y lo reduce a un DataFrame compacto con columnas:
        event_id, gun_x_mm, face_type, n_ph, fpt_ns

    Esto evita cargar todo el TTree completo a memoria.
    """
    with uproot.open(root_file) as f:
        tree = f["sipm_hits"]
        available = set(tree.keys())
        cols = ["event_id", "face_type", "time_ns"]
        if "gun_x_mm" in available:
            cols.append("gun_x_mm")

        reduced_chunks = []
        n_rows_total = 0

        for chunk in tree.iterate(cols, library="pd", step_size=step_size):
            if "gun_x_mm" not in chunk.columns:
                chunk["gun_x_mm"] = 0.0

            n_rows_total += len(chunk)

            # Reducir a nivel de evento/cara dentro del chunk
            grp = (chunk.groupby(["event_id", "gun_x_mm", "face_type"], sort=False)
                        .agg(n_ph=("time_ns", "size"),
                             fpt_ns=("time_ns", "min"))
                        .reset_index())

            reduced_chunks.append(grp)

    if not reduced_chunks:
        raise RuntimeError("No se encontraron entradas en 'sipm_hits'.")

    # Unir resúmenes parciales y volver a reducir por si un mismo evento quedó partido
    event_face = pd.concat(reduced_chunks, ignore_index=True)
    event_face = (event_face.groupby(["event_id", "gun_x_mm", "face_type"], sort=False)
                           .agg(n_ph=("n_ph", "sum"),
                                fpt_ns=("fpt_ns", "min"))
                           .reset_index())

    print(f"Cargados {n_rows_total:,} fotones detectados desde '{root_file}'")
    print(f"Reducido a {len(event_face):,} filas evento/cara")
    return event_face


def first_photon_time_from_event_face(event_face: pd.DataFrame, face_mask) -> pd.DataFrame:
    """
    Construye FPT por evento para la seleccion indicada.

    Retorna columnas: event_id, gun_x_mm, fpt_ns
    """
    sub = event_face.loc[face_mask, ["event_id", "gun_x_mm", "fpt_ns"]].copy()
    if len(sub) == 0:
        return pd.DataFrame(columns=["event_id", "gun_x_mm", "fpt_ns"])

    # Primer fotón entre las caras seleccionadas
    fpt = (sub.groupby(["event_id", "gun_x_mm"], sort=False)
              .agg(fpt_ns=("fpt_ns", "min"))
              .reset_index())
    return fpt


def compute_resolution_vs_x(event_face: pd.DataFrame, face_label: str, face_mask):
    fpt = first_photon_time_from_event_face(event_face, face_mask)
    if len(fpt) == 0:
        print(f"  [!] {face_label}: sin hits, salteando.")
        return pd.DataFrame()

    rows = []
    for x_val, grp in fpt.groupby("gun_x_mm", sort=True):
        times = grp["fpt_ns"].to_numpy()
        mu, sigma, sigma_err = fit_fpt_distribution(times)
        rows.append({
            "x_mm":       x_val,
            "sigma_ns":   sigma,
            "sigma_err":  sigma_err,
            "mu_ns":      mu,
            "n_events":   len(times),
        })

    return pd.DataFrame(rows).sort_values("x_mm").reset_index(drop=True)


def plot_resolution_vs_x(results: dict, out_pdf: str):
    fig, ax = plt.subplots(figsize=(11, 5.5))

    style = {
        "end_both": dict(color=COLORS["end_both"], marker="o",  ms=6,  lw=1.8, label="End SiPMs (L+R combinados)"),
        "top":      dict(color=COLORS["top"],      marker="s",  ms=6,  lw=1.8, label="Top SiPMs"),
        "end_left": dict(color=COLORS["end_left"], marker="^",  ms=5,  lw=1.2,
                         ls="--", alpha=0.6, label="End izquierdo (−X)"),
        "end_right":dict(color=COLORS["end_right"],marker="v",  ms=5,  lw=1.2,
                         ls="--", alpha=0.6, label="End derecho (+X)"),
    }

    for key, df_res in results.items():
        if df_res is None or len(df_res) == 0:
            continue
        ok = df_res["sigma_ns"].notna() & (df_res["sigma_ns"] > 0)
        sub = df_res[ok]
        if len(sub) == 0:
            continue

        x = sub["x_mm"].to_numpy()
        sigma = sub["sigma_ns"].to_numpy() * 1e3
        err = sub["sigma_err"].to_numpy() * 1e3

        st = style.get(key, {})
        ax.errorbar(x, sigma, yerr=err, capsize=3, capthick=1,
                    **{k: v for k, v in st.items()})

    ax.axhline(100, color="black", ls=":", lw=1.5, label="Objetivo TD SHiP (100 ps)")
    ax.set_xlabel("Posicion longitudinal del muon $x$ [mm]", fontsize=13)
    ax.set_ylabel(r"Resolucion temporal $\sigma_t$ [ps]", fontsize=13)
    ax.set_title("Resolucion temporal vs posicion longitudinal — EJ-200, SHiP TD",
                 fontsize=13, pad=10)

    ax.set_xlim(-BAR_HALF_X - 30, BAR_HALF_X + 30)
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.grid(which="major", alpha=0.35)
    ax.grid(which="minor", alpha=0.12)
    ax.legend(fontsize=10, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"  → {out_pdf}")
    plt.close(fig)


def plot_fpt_examples(event_face: pd.DataFrame, x_positions_mm, out_pdf: str,
                      face_mask, face_label: str):
    fpt = first_photon_time_from_event_face(event_face, face_mask)
    if len(fpt) == 0:
        return

    x_vals = sorted(fpt["gun_x_mm"].unique())
    chosen = []
    for target in x_positions_mm:
        best = min(x_vals, key=lambda v: abs(v - target))
        if best not in chosen:
            chosen.append(best)

    n = len(chosen)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4.5), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, x_val in zip(axes, chosen):
        times = fpt.loc[fpt["gun_x_mm"] == x_val, "fpt_ns"].to_numpy()
        mu, sigma, sigma_err = fit_fpt_distribution(times, n_bins=30)

        lo, hi = np.percentile(times, [1, 99])
        margin = max(0.5 * (hi - lo), 0.5)
        hist_counts, hist_edges, _ = ax.hist(
            times, bins=30, range=(lo - margin, hi + margin),
            color="#888", alpha=0.7, label="datos"
        )

        if not np.isnan(sigma):
            x_fit = np.linspace(lo - margin, hi + margin, 300)
            binw = hist_edges[1] - hist_edges[0]
            A_est = len(times) * binw
            ax.plot(
                x_fit, gauss(x_fit, mu, sigma, A_est),
                color="crimson", lw=2,
                label=fr"$\mu$={mu:.2f} ns" + "\n"
                      + fr"$\sigma$={sigma*1e3:.0f}±{sigma_err*1e3:.0f} ps"
            )

        ax.set_title(f"{face_label}\n$x$ = {x_val:.0f} mm", fontsize=11)
        ax.set_xlabel("FPT [ns]", fontsize=11)
        ax.set_ylabel("Eventos / bin" if ax == axes[0] else "", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.35)

    fig.suptitle("Distribucion del primer foton detectado (FPT)", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"  → {out_pdf}")
    plt.close(fig)


def plot_asymmetry_vs_x(event_face: pd.DataFrame, out_pdf: str):
    end = event_face[event_face["face_type"].isin([0, 1])].copy()
    if len(end) == 0:
        print("  [skip] asymmetry_vs_x: sin hits en end SiPMs")
        return

    evts = (end.pivot_table(index=["event_id", "gun_x_mm"],
                            columns="face_type",
                            values="n_ph",
                            aggfunc="sum",
                            fill_value=0)
               .reset_index())
    evts.columns.name = None
    if 0 not in evts.columns:
        evts[0] = 0
    if 1 not in evts.columns:
        evts[1] = 0
    evts = evts.rename(columns={0: "left", 1: "right"})

    total = evts["left"] + evts["right"]
    evts = evts[total > 0].copy()
    evts["asym"] = (evts["left"] - evts["right"]) / (evts["left"] + evts["right"])

    prof = evts.groupby("gun_x_mm")["asym"].agg(["mean", "std", "count"]).reset_index()
    prof["sem"] = prof["std"] / np.sqrt(prof["count"].clip(lower=1))

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.errorbar(prof["gun_x_mm"].to_numpy(), prof["mean"].to_numpy(), yerr=prof["sem"].to_numpy(),
                color=COLORS["end_both"], marker="o", ms=6, lw=1.8, capsize=3,
                label=r"$\langle A \rangle \pm \sigma_{\langle A \rangle}$")

    if prof["gun_x_mm"].nunique() >= 2:
        p = np.polyfit(prof["gun_x_mm"].to_numpy(), prof["mean"].to_numpy(), 1)
        x_fit = np.linspace(-BAR_HALF_X, BAR_HALF_X, 300)
        ax.plot(x_fit, np.polyval(p, x_fit), "k--", lw=1.5,
                label=f"Ajuste lineal: slope = {p[0]*1000:.2f} m$^{{-1}}$")

    ax.set_xlabel("Posicion del muon $x$ [mm]", fontsize=12)
    ax.set_ylabel("Asimetria $(N_L - N_R)/(N_L + N_R)$", fontsize=12)
    ax.set_title("Asimetria de carga end-SiPMs vs posicion longitudinal", fontsize=13)
    ax.set_xlim(-BAR_HALF_X - 30, BAR_HALF_X + 30)
    ax.axhline(0, color="gray", lw=0.8, ls=":")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"  → {out_pdf}")
    plt.close(fig)


def plot_n_photons_vs_x(event_face: pd.DataFrame, out_pdf: str):
    face_names = {0: "End izq. (−X)", 1: "End der. (+X)", 2: "Top (+Y)"}
    face_colors = {0: COLORS["end_left"], 1: COLORS["end_right"], 2: COLORS["top"]}

    prof = (event_face.groupby(["gun_x_mm", "face_type"])["n_ph"]
                     .agg(["mean", "sem"])
                     .reset_index())

    fig, ax = plt.subplots(figsize=(11, 4.5))
    for face, grp in prof.groupby("face_type"):
        ax.errorbar(grp["gun_x_mm"].to_numpy(), grp["mean"].to_numpy(), yerr=grp["sem"].to_numpy(),
                    color=face_colors.get(face, "black"), marker="o", ms=5, lw=1.6, capsize=2,
                    label=face_names.get(face, f"face {face}"))

    ax.set_xlabel("Posicion del muon $x$ [mm]", fontsize=12)
    ax.set_ylabel("Fotones detectados por evento", fontsize=12)
    ax.set_title("Fotones detectados vs posicion longitudinal", fontsize=13)
    ax.set_xlim(-BAR_HALF_X - 30, BAR_HALF_X + 30)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"  → {out_pdf}")
    plt.close(fig)


def print_summary_table(results: dict):
    print("\n" + "=" * 65)
    print(f"  {'Configuracion':<22} {'sigma_t medio [ps]':>18}  {'max-min [ps]':>12}")
    print("-" * 65)
    for key, df_res in results.items():
        if df_res is None or len(df_res) == 0:
            continue
        ok = df_res["sigma_ns"].notna() & (df_res["sigma_ns"] > 0)
        s = df_res.loc[ok, "sigma_ns"].to_numpy() * 1e3
        if len(s) == 0:
            continue
        label = {
            "end_left":  "End izquierdo",
            "end_right": "End derecho",
            "end_both":  "End L+R (OR)",
            "top":       "Top SiPMs",
        }.get(key, key)
        print(f"  {label:<22} {np.mean(s):>18.1f}  {np.max(s)-np.min(s):>12.1f}")
    print("=" * 65 + "\n")


def save_summary_csv(results: dict, out_csv: str = "resolution_vs_x_summary.csv"):
    rows = []
    for label, df_res in results.items():
        if df_res is None or len(df_res) == 0:
            continue
        sub = df_res.copy()
        sub["config"] = label
        rows.append(sub[["config", "x_mm", "mu_ns", "sigma_ns", "sigma_err", "n_events"]])

    if not rows:
        return

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(out_csv, index=False)
    print(f"  → {out_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Resolucion temporal sigma_t vs posicion x — EJ-200 TD SHiP")
    parser.add_argument("root_file", nargs="?", default=None,
                        help="Archivo ROOT con ntuple sipm_hits (default: photon_hits_merged.root o photon_hits.root)")
    parser.add_argument("--out", default="resolution_vs_x.pdf",
                        help="PDF de salida para la figura principal (default: resolution_vs_x.pdf)")
    parser.add_argument("--step-size", default="100 MB",
                        help="Tamano del chunk para uproot.iterate (default: 100 MB)")
    args = parser.parse_args()

    root_file = find_default_root_file(args.root_file)
    if not pathlib.Path(root_file).exists():
        print(f"ERROR: '{root_file}' no encontrado.")
        print("Prueba con: python merge_runs.py")
        sys.exit(1)

    event_face = load_event_level_data(root_file, step_size=args.step_size)

    n_pos = event_face["gun_x_mm"].nunique()
    print(f"Posiciones de muon encontradas: {n_pos}  "
          f"({event_face['gun_x_mm'].min():.0f} a {event_face['gun_x_mm'].max():.0f} mm)")
    if n_pos < 3:
        print("AVISO: se encontraron < 3 posiciones distintas.")
        print("Para la figura de tesis usa scan.mac (21 posiciones).")

    mask_left  = event_face["face_type"] == 0
    mask_right = event_face["face_type"] == 1
    mask_end   = event_face["face_type"].isin([0, 1])
    mask_top   = event_face["face_type"] == 2

    print("\nCalculando FPT y ajustando Gaussianas…")
    results = {
        "end_left":  compute_resolution_vs_x(event_face, "End izquierdo", mask_left),
        "end_right": compute_resolution_vs_x(event_face, "End derecho", mask_right),
        "end_both":  compute_resolution_vs_x(event_face, "End L+R OR", mask_end),
        "top":       compute_resolution_vs_x(event_face, "Top SiPMs", mask_top),
    }

    print("\nGenerando figuras…")
    plot_resolution_vs_x(results, args.out)

    x_examples = [-BAR_HALF_X * 0.9, 0.0, BAR_HALF_X * 0.9]
    plot_fpt_examples(event_face, x_examples, "fpt_dist_end.pdf", mask_end, "End SiPMs")
    plot_fpt_examples(event_face, x_examples, "fpt_dist_top.pdf", mask_top, "Top SiPMs")
    plot_asymmetry_vs_x(event_face, "asymmetry_vs_x.pdf")
    plot_n_photons_vs_x(event_face, "n_photons_vs_x.pdf")
    save_summary_csv(results)

    print_summary_table(results)
    print("Listo.")


if __name__ == "__main__":
    main()
