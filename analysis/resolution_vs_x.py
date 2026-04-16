#!/usr/bin/env python3
"""
resolution_vs_x.py — Resolucion temporal sigma_t vs posicion longitudinal x

Figura central de la tesis: compara end SiPMs vs top SiPMs.

Requiere el ntuple sipm_hits producido por scan.mac (21 posiciones, 200
eventos cada una).  El campo gun_x_mm identifica la posicion del muon.

Uso:
    python analysis/resolution_vs_x.py                    # photon_hits.root
    python analysis/resolution_vs_x.py mi_scan.root
    python analysis/resolution_vs_x.py scan.root --out tesis_fig4.pdf
"""

import sys
import argparse
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.optimize import curve_fit
from scipy.stats import norm as sp_norm
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


def compute_dcfd_time(times_ns, fraction=0.14):
    """
    Reconstruye el tiempo de cruce dCFD a partir del tren temporal de PEs.
    La amplitud del pulso se aproxima como N_PE y el umbral es fraction*N_PE.
    """
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

# ── Helpers ───────────────────────────────────────────────────────────────────

def load_data(root_file: str) -> pd.DataFrame:
    with uproot.open(root_file) as f:
        tree = f["sipm_hits"]
        df = tree.arrays(
            ["event_id", "face_type", "global_id", "local_id",
             "time_ns", "energy_eV", "wl_nm", "pde",
             "x_mm", "y_mm", "z_mm", "gun_x_mm"],
            library="pd"
        )
    print(f"Cargados {len(df):,} fotones detectados desde '{root_file}'")
    return df


def dcfd_time_per_event(df: pd.DataFrame, face_mask,
                        fraction: float = 0.14,
                        electronics_sigma_ps: float = 30.0,
                        seed: int = 12345) -> pd.DataFrame:
    """
    Calcula el tiempo dCFD por evento para los SiPMs seleccionados.
    Al tiempo dCFD se le suma un jitter FastIC+ gaussiano de 30 ps RMS.

    Retorna DataFrame con columnas: event_id, gun_x_mm, dcfd_time_ns
    """
    sub = df[face_mask][["event_id", "gun_x_mm", "time_ns"]].copy()
    rows = []
    rng = np.random.default_rng(seed)

    for event_id, grp in sub.groupby("event_id"):
        t_dcfd = compute_dcfd_time(grp["time_ns"].to_numpy(), fraction=fraction)
        if np.isnan(t_dcfd):
            continue
        rows.append({
            "event_id": event_id,
            "gun_x_mm": grp["gun_x_mm"].iloc[0],
            "dcfd_time_ns": t_dcfd + rng.normal(0.0, electronics_sigma_ps * 1e-3),
        })

    return pd.DataFrame(rows)


def gauss(x, mu, sigma, A):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_fpt_distribution(times_ns, n_bins=40):
    """
    Ajusta una Gaussiana a la distribucion de FPT.
    Retorna (mu [ns], sigma [ns], sigma_err [ns]).
    Retorna (nan, nan, nan) si el ajuste falla o hay <10 eventos.
    """
    if len(times_ns) < 10:
        return np.nan, np.nan, np.nan

    lo, hi = np.percentile(times_ns, [1, 99])
    margin  = max(0.5 * (hi - lo), 0.5)   # minimo 0.5 ns de rango
    lo -= margin
    hi += margin

    counts, edges = np.histogram(times_ns, bins=n_bins, range=(lo, hi))
    centres = 0.5 * (edges[:-1] + edges[1:])

    # Estimados iniciales: media muestral, desv. est. muestral, maximo del hist.
    mu0    = float(np.mean(times_ns))
    sigma0 = float(np.std(times_ns))
    if sigma0 < 1e-6:
        sigma0 = 0.1
    A0     = float(counts.max())

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
        # Si scipy falla, devolver sigma RMS (estimador robusto)
        return mu0, float(np.std(times_ns)), float(np.std(times_ns) / np.sqrt(len(times_ns)))


# ── Analisis principal ────────────────────────────────────────────────────────

def compute_resolution_vs_x(df: pd.DataFrame, face_label: str, face_mask):
    """
    Para cada posicion gun_x_mm, calcula el tiempo dCFD de todos los eventos y
    ajusta una Gaussiana.  Retorna DataFrame con:
        x_mm, sigma_ns, sigma_err_ns, mu_ns, n_events
    """
    t_evt = dcfd_time_per_event(df, face_mask)
    if len(t_evt) == 0:
        print(f"  [!] {face_label}: sin hits, salteando.")
        return pd.DataFrame()

    rows = []
    for x_val, grp in t_evt.groupby("gun_x_mm"):
        times = grp["dcfd_time_ns"].values
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
    """
    Figura principal: sigma_t [ps] vs x [mm] para end SiPMs y top SiPMs.
    """
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
        ok  = df_res["sigma_ns"].notna() & (df_res["sigma_ns"] > 0)
        sub = df_res[ok]
        if len(sub) == 0:
            continue

        x     = sub["x_mm"].values
        sigma = sub["sigma_ns"].values * 1e3          # ns → ps
        err   = sub["sigma_err"].values * 1e3

        st = style.get(key, {})
        ax.errorbar(x, sigma, yerr=err,
                    capsize=3, capthick=1,
                    **{k: v for k, v in st.items()})

    # Linea de referencia: objetivo de diseno del TD de SHiP = 100 ps
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


def plot_fpt_examples(df: pd.DataFrame, x_positions_mm, out_pdf: str,
                      face_mask, face_label: str):
    """
    Distribucion del tiempo dCFD (con ajuste Gaussiano) para tres posiciones
    representativas: centro, cuarto de barra, extremo.
    """
    t_evt = dcfd_time_per_event(df, face_mask)
    if len(t_evt) == 0:
        return

    x_vals = sorted(t_evt["gun_x_mm"].unique())
    # Seleccionar las tres posiciones mas cercanas a lo pedido
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
        times = t_evt[t_evt["gun_x_mm"] == x_val]["dcfd_time_ns"].values
        mu, sigma, sigma_err = fit_fpt_distribution(times, n_bins=30)

        lo, hi = np.percentile(times, [1, 99])
        margin = max(0.5 * (hi - lo), 0.5)
        ax.hist(times, bins=30, range=(lo - margin, hi + margin),
                color="#888", alpha=0.7, label="datos")

        if not np.isnan(sigma):
            x_fit = np.linspace(lo - margin, hi + margin, 300)
            A_est = len(times) * (hi - lo + 2 * margin) / 30
            ax.plot(x_fit, gauss(x_fit, mu, sigma, A_est),
                    color="crimson", lw=2,
                    label=fr"$\mu$={mu:.2f} ns" + "\n"
                          + fr"$\sigma$={sigma*1e3:.0f}±{sigma_err*1e3:.0f} ps")

        ax.set_title(f"{face_label}\n$x$ = {x_val:.0f} mm", fontsize=11)
        ax.set_xlabel("dCFD time [ns]", fontsize=11)
        ax.set_ylabel("Eventos / bin" if ax == axes[0] else "", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.35)

    fig.suptitle("Distribucion del timestamp dCFD reconstruido", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"  → {out_pdf}")
    plt.close(fig)


def plot_asymmetry_vs_x(df: pd.DataFrame, out_pdf: str):
    """
    Asimetria de carga <(N_L - N_R)/(N_L + N_R)> vs posicion x.
    Muestra la linealidad del observable de reconstruccion de posicion.
    """
    end = df[df["face_type"].isin([0, 1])].copy()
    if len(end) == 0:
        print("  [skip] asymmetry_vs_x: sin hits en end SiPMs")
        return

    # Contar fotones por evento y cara
    evts = (end.groupby(["event_id", "gun_x_mm", "face_type"])
               .size()
               .unstack(fill_value=0)
               .reset_index())
    evts.columns.name = None
    evts = evts.rename(columns={0: "left", 1: "right"})

    total = evts["left"] + evts["right"]
    evts  = evts[total > 0].copy()
    evts["asym"] = (evts["left"] - evts["right"]) / (evts["left"] + evts["right"])

    # Media y desv. est. por posicion
    prof = evts.groupby("gun_x_mm")["asym"].agg(["mean", "std", "count"]).reset_index()
    prof["sem"] = prof["std"] / np.sqrt(prof["count"])

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.errorbar(prof["gun_x_mm"].to_numpy(), prof["mean"].to_numpy(), yerr=prof["sem"].to_numpy(),
                color=COLORS["end_both"], marker="o", ms=6, lw=1.8, capsize=3,
                label=r"$\langle A \rangle \pm \sigma_{\langle A \rangle}$")

    # Ajuste lineal (requiere al menos 2 posiciones distintas)
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


def plot_n_photons_vs_x(df: pd.DataFrame, out_pdf: str):
    """
    Numero medio de fotones detectados por evento vs posicion x,
    separado por face_type.  Observable de uniformidad.
    """
    face_names = {0: "End izq. (−X)", 1: "End der. (+X)", 2: "Top (+Y)"}
    face_colors = {0: COLORS["end_left"], 1: COLORS["end_right"], 2: COLORS["top"]}

    # Fotones por evento y cara
    per_evt = (df.groupby(["event_id", "gun_x_mm", "face_type"])
                 .size()
                 .reset_index(name="n_ph"))

    # Promedio por posicion y cara
    prof = (per_evt.groupby(["gun_x_mm", "face_type"])["n_ph"]
                   .agg(["mean", "sem"])
                   .reset_index())

    fig, ax = plt.subplots(figsize=(11, 4.5))
    for face, grp in prof.groupby("face_type"):
        ax.errorbar(grp["gun_x_mm"].to_numpy(), grp["mean"].to_numpy(), yerr=grp["sem"].to_numpy(),
                    color=face_colors[face], marker="o", ms=5, lw=1.6, capsize=2,
                    label=face_names[face])

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
    """Imprime tabla de sigma_t promedio y variacion para la comparacion de tesis."""
    print("\n" + "=" * 65)
    print(f"  {'Configuracion':<22} {'sigma_t medio [ps]':>18}  {'max-min [ps]':>12}")
    print("-" * 65)
    for key, df_res in results.items():
        if df_res is None or len(df_res) == 0:
            continue
        ok = df_res["sigma_ns"].notna() & (df_res["sigma_ns"] > 0)
        s  = df_res.loc[ok, "sigma_ns"].values * 1e3  # ps
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


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Resolucion temporal sigma_t vs posicion x — EJ-200 TD SHiP")
    parser.add_argument("root_file", nargs="?", default="photon_hits.root",
                        help="Archivo ROOT con ntuple sipm_hits (default: photon_hits.root)")
    parser.add_argument("--out", default="resolution_vs_x.pdf",
                        help="PDF de salida para la figura principal (default: resolution_vs_x.pdf)")
    args = parser.parse_args()

    if not pathlib.Path(args.root_file).exists():
        print(f"ERROR: '{args.root_file}' no encontrado.")
        print("Primero corre la simulacion:  ./ej200_bar_sim -m macros/scan.mac")
        sys.exit(1)

    df = load_data(args.root_file)

    # Verificar que gun_x_mm esta presente (requiere ej200_v3+)
    if "gun_x_mm" not in df.columns:
        print("ERROR: el ntuple no contiene 'gun_x_mm'.")
        print("Recompila con ej200_v3 y vuelve a correr el scan.")
        sys.exit(1)

    n_pos = df["gun_x_mm"].nunique()
    print(f"Posiciones de muon encontradas: {n_pos}  "
          f"({df['gun_x_mm'].min():.0f} a {df['gun_x_mm'].max():.0f} mm)")

    if n_pos < 3:
        print("AVISO: se encontraron < 3 posiciones distintas.")
        print("Para la figura de tesis usa scan.mac (21 posiciones).")

    # ── Mascaras de face_type ────────────────────────────────────────────────
    mask_left  = df["face_type"] == 0
    mask_right = df["face_type"] == 1
    mask_end   = df["face_type"].isin([0, 1])  # OR: primer foton en cualquier end
    mask_top   = df["face_type"] == 2

    # ── Calcular sigma_t vs x ────────────────────────────────────────────────
    print("\nCalculando tiempos dCFD y ajustando Gaussianas…")
    results = {}
    results["end_left"]  = compute_resolution_vs_x(df, "End izquierdo", mask_left)
    results["end_right"] = compute_resolution_vs_x(df, "End derecho",   mask_right)
    results["end_both"]  = compute_resolution_vs_x(df, "End L+R OR",    mask_end)
    results["top"]       = compute_resolution_vs_x(df, "Top SiPMs",     mask_top)

    # ── Figura principal ─────────────────────────────────────────────────────
    print("\nGenerando figuras…")
    plot_resolution_vs_x(results, args.out)

    # ── Distribuciones FPT en tres posiciones representativas ────────────────
    x_examples = [-BAR_HALF_X * 0.9, 0.0, BAR_HALF_X * 0.9]  # extremo, centro, extremo

    plot_fpt_examples(df, x_examples,
                      out_pdf="fpt_dist_end.pdf",
                      face_mask=mask_end,
                      face_label="End SiPMs")

    plot_fpt_examples(df, x_examples,
                      out_pdf="fpt_dist_top.pdf",
                      face_mask=mask_top,
                      face_label="Top SiPMs")

    # ── Asimetria vs x ───────────────────────────────────────────────────────
    plot_asymmetry_vs_x(df, out_pdf="asymmetry_vs_x.pdf")

    # ── Fotones por evento vs x ───────────────────────────────────────────────
    plot_n_photons_vs_x(df, out_pdf="n_photons_vs_x.pdf")

    # ── Resumen en pantalla ───────────────────────────────────────────────────
    print_summary_table(results)
    print("Listo.")


if __name__ == "__main__":
    main()
