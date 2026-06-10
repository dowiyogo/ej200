#!/usr/bin/env python3
"""Physical fits, QA plots, and generated Beamer tables for CODEX_EXEC_09."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/analysis_sigma_vs_x_2026-06-10"
OUT = BASE / "deep_dive"


def chi2_ndf(y, prediction, error, npar):
    return float(np.sum(((y - prediction) / error) ** 2) / (len(y) - npar))


def save_attenuation(scan, rows):
    selected = scan[np.abs(scan.posicion_mm) <= 500].copy()
    x = np.concatenate([selected.posicion_mm, selected.posicion_mm])
    side = np.concatenate([np.zeros(len(selected)), np.ones(len(selected))])
    y = np.concatenate([selected.NPE_L, selected.NPE_R])
    error = np.sqrt(np.maximum(y, 1.0) / 10000.0)

    def model(inputs, n0_left, n0_right, lambda_mm):
        xx, ss = inputs
        distance = np.where(ss == 0, 700.0 + xx, 700.0 - xx)
        n0 = np.where(ss == 0, n0_left, n0_right)
        return n0 * np.exp(-distance / lambda_mm)

    popt, pcov = curve_fit(model, (x, side), y, sigma=error, absolute_sigma=True,
                           p0=[2500, 2500, 1200], maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    prediction = model((x, side), *popt)
    quality = chi2_ndf(y, prediction, error, 3)
    rows.append(("attenuation", "lambda_eff_mm", popt[2], perr[2], quality))

    grid = np.linspace(-500, 500, 300)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.errorbar(selected.posicion_mm.to_numpy(), selected.NPE_L.to_numpy(),
                yerr=np.sqrt(selected.NPE_L.to_numpy() / 10000),
                fmt="o", label="End-L")
    ax.errorbar(selected.posicion_mm.to_numpy(), selected.NPE_R.to_numpy(),
                yerr=np.sqrt(selected.NPE_R.to_numpy() / 10000),
                fmt="s", label="End-R")
    ax.plot(grid, model((grid, np.zeros_like(grid)), *popt), label="Fit End-L")
    ax.plot(grid, model((grid, np.ones_like(grid)), *popt), label="Fit End-R")
    ax.set(title="Simulación Geant4 — End readout: atenuación efectiva",
           xlabel="Posición x [mm]", ylabel="NPE medio por evento")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.text(0.03, 0.04,
            rf"$\lambda_{{eff}}={popt[2]/10:.1f}\pm{perr[2]/10:.1f}$ cm"
            "\n" rf"$\lambda_{{bulk}}=160$ cm"
            "\n" rf"$\chi^2/ndf={quality:.1f}$",
            transform=ax.transAxes, bbox=dict(facecolor="white", alpha=0.9))
    fig.tight_layout()
    fig.savefig(OUT / "attenuation_effective.png", dpi=160)
    plt.close(fig)
    return popt[2], perr[2], quality


def save_photo_scaling(scan, rows):
    npe_far = np.minimum(scan.NPE_L.to_numpy(), scan.NPE_R.to_numpy())
    sigma = scan.sigma_single_ps.to_numpy()
    error = scan.err_sigma_single_ps.to_numpy()

    def model(npe, p0, p1):
        return np.sqrt(p0 * p0 + p1 * p1 / npe)

    popt, pcov = curve_fit(model, npe_far, sigma, sigma=error, absolute_sigma=True,
                           p0=[70, 900], bounds=(0, np.inf), maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    quality = chi2_ndf(sigma, model(npe_far, *popt), error, 2)
    rows.extend([
        ("photo_scaling", "p0_ps", popt[0], perr[0], quality),
        ("photo_scaling", "p1_ps_sqrt_npe", popt[1], perr[1], quality),
    ])

    z = 1.0 / np.sqrt(npe_far)
    linear, linear_cov = np.polyfit(z, sigma, 1, w=1 / error, cov=True)
    grid_npe = np.linspace(npe_far.min(), npe_far.max(), 300)
    grid_z = np.linspace(z.min(), z.max(), 300)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].errorbar(npe_far, sigma, yerr=error, fmt="o")
    axes[0].plot(grid_npe, model(grid_npe, *popt), color="red")
    axes[0].set(xlabel=r"$NPE_{far}$", ylabel=r"$\sigma_{single}$ [ps]",
                title="Modelo cuadrático")
    axes[1].errorbar(z, sigma, yerr=error, fmt="o")
    axes[1].plot(grid_z, np.polyval(linear, grid_z), color="red")
    axes[1].set(xlabel=r"$1/\sqrt{NPE_{far}}$", ylabel=r"$\sigma_{single}$ [ps]",
                title="Panel lineal equivalente")
    for ax in axes:
        ax.grid(True, alpha=0.3)
    axes[0].text(0.04, 0.68,
                 rf"$p_0={popt[0]:.1f}\pm{perr[0]:.1f}$ ps"
                 "\n" rf"$p_1={popt[1]:.0f}\pm{perr[1]:.0f}$ ps$\sqrt{{NPE}}$"
                 "\n" rf"$\chi^2/ndf={quality:.1f}$",
                 transform=axes[0].transAxes, bbox=dict(facecolor="white", alpha=0.9))
    fig.suptitle("Simulación Geant4 — End readout: escalamiento fotoestadístico")
    fig.tight_layout()
    fig.savefig(OUT / "photo_scaling.png", dpi=160)
    plt.close(fig)
    return popt, perr, quality


def save_velocity(metrics, rows):
    x = metrics.position_mm.to_numpy()
    mean = metrics.fit_mean_ns.to_numpy()
    error = metrics.fit_sigma_lr_ps.to_numpy() / 1000.0 / np.sqrt(metrics.n_eff.to_numpy())
    coeff, cov = np.polyfit(x, mean, 1, w=1 / error, cov=True)
    slope, offset = coeff
    slope_err, offset_err = np.sqrt(np.diag(cov))
    velocity = 2.0 / slope
    velocity_err = abs(2.0 * slope_err / slope**2)
    quality = chi2_ndf(mean, np.polyval(coeff, x), error, 2)
    rows.append(("velocity", "v_eff_mm_ns", velocity, velocity_err, quality))

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.errorbar(x, mean, yerr=error, fmt="o", label=r"$\mu(\Delta T_{LR})$")
    grid = np.linspace(x.min(), x.max(), 300)
    ax.plot(grid, np.polyval(coeff, grid), color="red", label="Fit lineal")
    ax.set(title="Simulación Geant4 — End readout: velocidad efectiva",
           xlabel="Posición x [mm]", ylabel=r"$\mu(\Delta T_{LR})$ [ns]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.text(0.04, 0.72,
            rf"$v_{{eff}}={velocity:.1f}\pm{velocity_err:.1f}$ mm/ns"
            "\n" rf"$c/n\simeq190$ mm/ns"
            "\n" rf"$\chi^2/ndf={quality:.1f}$",
            transform=ax.transAxes, bbox=dict(facecolor="white", alpha=0.9))
    fig.tight_layout()
    fig.savefig(OUT / "velocity_effective.png", dpi=160)
    plt.close(fig)
    return velocity, velocity_err, quality


def save_gaussianity(scan):
    ratio = scan.fwhm_ps / (2.355 * scan.sigma_LR_ps)
    table = pd.DataFrame({"position_mm": scan.posicion_mm, "fwhm_over_2p355sigma": ratio})
    table.to_csv(OUT / "gaussianity.csv", index=False)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(scan.posicion_mm.to_numpy(), ratio.to_numpy(), "o-")
    ax.axhspan(0.95, 1.05, color="green", alpha=0.18, label="Banda ±5%")
    ax.axhline(1.0, color="black", linestyle="--")
    ax.set(title="Simulación Geant4 — End readout: gaussianidad",
           xlabel="Posición x [mm]", ylabel=r"FWHM / $(2.355\,\sigma_{gauss})$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "gaussianity_qa.png", dpi=160)
    plt.close(fig)
    return ratio


def save_ablation():
    data = pd.read_csv(OUT / "ablation_results.csv")
    labels = ["Etapa 1\nopen", "Open +\nfit espejo", "Wrapped +\nfit Etapa 1",
              "Wrapped +\nfit espejo", "Wrapped +\np25", "Wrapped +\np50"]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    bars = ax.bar(labels, data.sigma_single_ps, yerr=data.sigma_single_err_ps,
                  capsize=4, color=["#777777", "#999999", "#2b8cbe", "#0868ac",
                                    "#fdae6b", "#e6550d"])
    ax.axhline(88, color="red", linestyle="--", label="Test beam reportado")
    ax.set(title="Simulación Geant4 — End readout: ablación en x=0",
           ylabel=r"$\sigma_{single}$ [ps]")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    for bar, value in zip(bars, data.sigma_single_ps):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 12, f"{value:.1f}",
                ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "ablation_x0.png", dpi=160)
    plt.close(fig)


def generated_tables(scan, metrics, ratio, physical_rows):
    selected_x = [-690, -500, -350, -200, 0, 200, 350, 500, 690]
    qa = scan.merge(metrics[["position_mm", "n_eff"]],
                    left_on="posicion_mm", right_on="position_mm")
    qa["ratio"] = ratio
    qa.to_csv(OUT / "qa_metrics.csv", index=False)
    qa = qa[qa.position_mm.isin(selected_x)]
    with (OUT / "qa_table.tex").open("w") as out:
        out.write("\\begin{tabular}{rrrrrrr}\\toprule\n")
        out.write("$x$ & $N_{eff}$ & $\\sigma_{single}$ & $\\chi^2_g/ndf$ & "
                  "FWHM/$2.355\\sigma$ & MPV & $\\chi^2_L/ndf$\\\\\\midrule\n")
        for _, r in qa.iterrows():
            out.write(f"{r.position_mm:.0f} & {r.n_eff:.0f} & "
                      f"${r.sigma_single_ps:.1f}\\pm{r.err_sigma_single_ps:.1f}$ & "
                      f"{r.chi2_ndf_gaus:.2f} & {r['ratio']:.2f} & "
                      f"{r.MPV_NPE:.1f} & {r.chi2_ndf_landau:.2f}\\\\\n")
        out.write("\\bottomrule\\end{tabular}\n")

    ablation = pd.read_csv(OUT / "ablation_results.csv")
    with (OUT / "ablation_table.tex").open("w") as out:
        out.write("\\begin{tabular}{lrr}\\toprule\nVariante & $\\sigma_{single}$ [ps] & $N_{eff}$\\\\\\midrule\n")
        names = ["Etapa 1 histórica reproducida", "Open + fit espejo",
                 "Wrapped + fit Etapa 1", "Wrapped + fit espejo",
                 "Wrapped + espejo + p25", "Wrapped + espejo + p50"]
        for name, (_, r) in zip(names, ablation.iterrows()):
            out.write(f"{name} & ${r.sigma_single_ps:.1f}\\pm{r.sigma_single_err_ps:.1f}$ "
                      f"& {r.n_eff:.0f}\\\\\n")
        out.write("\\bottomrule\\end{tabular}\n")

    trace = [
        ("SPR", "$\\tau_r=0.5$ ns, $\\tau_f=5$ ns", "igual"),
        ("Umbral LE", "4 PE absoluto", "igual"),
        ("Trigger End", "primer SUM4 por extremo", "igual"),
        ("Corrección walk", "no", "no"),
        ("Corte energía", "no", "NPE $\\geq$p25/p50"),
        ("Fit", "mediana$\\pm8$MAD; 4 it. $\\pm2\\sigma$", "pico + ventana local 4 ns"),
        ("Selección", "ambos extremos disparan", "igual"),
        ("Dataset/geometría", "Top abierta", "End wrapped, $+Y$ envuelta"),
    ]
    with (OUT / "traceability_table.tex").open("w") as out:
        out.write("\\begin{tabular}{p{0.19\\textwidth}p{0.35\\textwidth}p{0.35\\textwidth}}\\toprule\n")
        out.write("Elemento & Etapa 1 provisional & Espejo TB EXEC\\_08\\\\\\midrule\n")
        for row in trace:
            out.write(" & ".join(row) + "\\\\\n")
        out.write("\\bottomrule\\end{tabular}\n")

    fits = pd.DataFrame(physical_rows, columns=["fit", "parameter", "value", "error", "chi2_ndf"])
    fits.to_csv(OUT / "fit_results.csv", index=False)
    with (OUT / "fit_summary.tex").open("w") as out:
        values = {r.parameter: r for r in fits.itertuples()}
        ablation = pd.read_csv(OUT / "ablation_results.csv").set_index("variant")
        out.write(f"\\newcommand{{\\LambdaEffCm}}{{{values['lambda_eff_mm'].value/10:.3f}}}\n")
        out.write(f"\\newcommand{{\\LambdaEffErrCm}}{{{values['lambda_eff_mm'].error/10:.3f}}}\n")
        out.write("\\newcommand{\\LambdaBulkCm}{160}\n")
        out.write(f"\\newcommand{{\\AttenChi}}{{{values['lambda_eff_mm'].chi2_ndf:.1f}}}\n")
        out.write(f"\\newcommand{{\\PhotoPzero}}{{{values['p0_ps'].value:.1f}}}\n")
        out.write(f"\\newcommand{{\\PhotoPzeroErr}}{{{values['p0_ps'].error:.1f}}}\n")
        out.write(f"\\newcommand{{\\PhotoPone}}{{{values['p1_ps_sqrt_npe'].value:.0f}}}\n")
        out.write(f"\\newcommand{{\\PhotoChi}}{{{values['p0_ps'].chi2_ndf:.1f}}}\n")
        out.write(f"\\newcommand{{\\VelocityEff}}{{{values['v_eff_mm_ns'].value:.1f}}}\n")
        out.write(f"\\newcommand{{\\VelocityEffErr}}{{{values['v_eff_mm_ns'].error:.1f}}}\n")
        out.write(f"\\newcommand{{\\VelocityChi}}{{{values['v_eff_mm_ns'].chi2_ndf:.1f}}}\n")
        out.write(f"\\newcommand{{\\LightSpeedOverN}}{{{300/1.58:.1f}}}\n")
        out.write(f"\\newcommand{{\\AblationOldStage}}{{{ablation.loc['historical_stage1_reproduced','sigma_single_ps']:.1f}}}\n")
        out.write(f"\\newcommand{{\\AblationWrappedStage}}{{{ablation.loc['wrapped_geometry_stage1_fit','sigma_single_ps']:.1f}}}\n")
        out.write(f"\\newcommand{{\\AblationWrappedMirror}}{{{ablation.loc['wrapped_geometry_mirror_fit','sigma_single_ps']:.1f}}}\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    scan = pd.read_csv(BASE / "sigma_vs_x_End.csv")
    metrics = pd.read_csv(OUT / "event_metrics.csv")
    physical_rows = []
    save_attenuation(scan, physical_rows)
    save_photo_scaling(scan, physical_rows)
    save_velocity(metrics, physical_rows)
    ratio = save_gaussianity(scan)
    save_ablation()
    generated_tables(scan, metrics, ratio, physical_rows)


if __name__ == "__main__":
    main()
