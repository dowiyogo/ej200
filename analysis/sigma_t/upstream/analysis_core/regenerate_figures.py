#!/usr/bin/env python3.12
"""
regenerate_figures.py — CP3.5: Regenerate EJ-230 figures with filtering rules.

Reads existing CSVs and ROOT sidecars from CP3 (no fits recomputed).
Applies:
  Rule 1: omit positions with efficiency < 5%
  Rule 2: mark TOP_NEAREST N=20 with asterisk in legend

Overwrites existing PNG/PDF in outputs/.
"""

import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─── config ────────────────────────────────────────────────────────────────────
BASE     = Path("/home/reriosto/SHiP/analysis_core")
OUTDIR   = BASE / "outputs"
CFG_PATH = BASE / "config" / "exec16_config.yaml"

cfg       = yaml.safe_load(CFG_PATH.read_text())
EFF_FLOOR = cfg.get("EFFICIENCY_FLOOR", 0.05)
CHI2_WARN = cfg.get("CHI2_NDF_WARN", 3.0)
CORE_FRAC = cfg.get("CORE_NOT_GAUSSIAN_FRACTION", 0.30)
REPR_POS  = cfg.get("REPR_POSITIONS_MM", [-690, 0, 690])
N_VALUES  = cfg.get("N_VALUES", [4, 20])

# ── material parameters set dynamically in main() ──────────────────────────
MAT_SLUG  = "ej230_endtop"   # default; overridden by --material arg
MAT_LABEL = "EJ-230"
BRANCH    = "feat/ej230-sslg4"
SHA       = "ca2f1c3"
N_EV      = 300
N_POS     = 31

# groups and their display-friendly slugs (CSV filename part)
GROUPS = cfg.get("GROUPS", [])
SLUG_MAP = {g: g.lower() for g in GROUPS}    # e.g. "TOP_SUM4" → "top_sum4"

# ─── helpers ──────────────────────────────────────────────────────────────────

def _info(ax, group, N, asterisk=False):
    """Small metadata box in the lower-right corner of a subplot."""
    lbl = f"* single-SiPM high-N" if asterisk else ""
    ax.text(
        0.99, 0.01,
        f"{MAT_LABEL} | {BRANCH} @ {SHA}\n"
        f"{N_EV} ev/pos | {N_POS} pos | N={N} | EXEC_16 {lbl}",
        transform=ax.transAxes, fontsize=6, ha="right", va="bottom",
        family="monospace", bbox=dict(boxstyle="round,pad=0.2", alpha=0.1),
    )


def _save(fig, name: str, exts=("png", "pdf")):
    """Save figure as both PNG and PDF, overwriting existing files."""
    paths = []
    for ext in exts:
        p = str(OUTDIR / f"{name}.{ext}")
        fig.savefig(p, dpi=150, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


def _nearest_valid_x(df_filt: pd.DataFrame, target_x: int) -> int | None:
    """Return x_gun_mm in df_filt closest to target_x, or None if df_filt empty."""
    if df_filt.empty:
        return None
    xs = df_filt["x_gun_mm"].values
    idx = int(np.argmin(np.abs(xs - target_x)))   # index of minimum distance
    return int(xs[idx])


def _asterisk_note(group: str, N: int) -> str:
    """Return '* single-SiPM, high-N regime' note if applicable, else ''."""
    return ("  (* single-SiPM, high-N regime)"
            if group == "TOP_NEAREST" and N == 20 else "")


# ─── Type A: resolution vs position ───────────────────────────────────────────

def plot_resolution(df: pd.DataFrame, group: str, slug: str) -> list:
    """
    Two-panel: N=4 (left), N=20 (right).
    Only plot points with efficiency >= EFF_FLOOR.
    TOP_NEAREST N=20 gets asterisk in legend.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"{MAT_LABEL} EndTop — {group} | σ_fit(x)",
                 fontsize=12, fontweight="bold")

    colors = {4: "#2166AC", 20: "#D6604D"}   # blue / red

    for ax, N in zip(axes, N_VALUES):
        dN   = df[(df["N"] == N) & (df["efficiency"] >= EFF_FLOOR)].copy()
        col  = colors.get(N, "gray")
        star = _asterisk_note(group, N)
        is_ast = bool(star)

        if len(dN) < 3:
            # not enough data after filtering — annotate and leave empty
            ax.text(0.5, 0.5,
                    f"N={N}: only {len(dN)} positions\nwith eff ≥ {EFF_FLOOR:.0%}",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=10, color="gray")
            ax.set_title(f"N={N}  [insufficient data]", fontsize=10)
        else:
            # combine TF1 error and bootstrap error → take the larger
            yerr = np.maximum(
                dN["sigma_fit_err_ps"].fillna(0).values,
                dN["bootstrap_err_ps"].fillna(0).values,
            )
            ax.errorbar(
                dN["x_gun_mm"].values,
                dN["sigma_fit_ps"].values,
                yerr=yerr,
                fmt="o-", color=col, ms=5, lw=1.5, capsize=3, elinewidth=1,
                label=f"N={N}{' *' if is_ast else ''}  σ_fit(x)",
            )
            if is_ast:
                # extra text in legend to explain the asterisk
                ax.plot([], [], " ", label="* single-SiPM, high-N regime")

            ax.set_title(f"N={N}{star}", fontsize=10)

        ax.set_xlabel("x_gun [mm]", fontsize=10)
        ax.set_ylabel("σ_fit [ps]", fontsize=10)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper left", fontsize=8)
        _info(ax, group, N, asterisk=is_ast)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return _save(fig, f"{MAT_SLUG}_{slug}_resolution_N4_N20")


# ─── Type B: χ²/NDF diagnostic ────────────────────────────────────────────────

def plot_chi2ndf(df: pd.DataFrame, group: str, slug: str) -> tuple:
    """
    Two-panel: N=4 (left), N=20 (right).
    Only plot positions with eff >= EFF_FLOOR.
    Returns (saved_files, updated_flags) where updated_flags is a dict
    {N: bool} indicating whether the group remains flagged after filtering.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"{MAT_LABEL} EndTop — {group} | χ²/NDF(x)  [fit quality diagnostic]",
                 fontsize=12, fontweight="bold")

    colors = {4: "#2166AC", 20: "#D6604D"}
    updated_flags = {}

    for ax, N in zip(axes, N_VALUES):
        dN    = df[(df["N"] == N) & (df["efficiency"] >= EFF_FLOOR)].copy()
        col   = colors.get(N, "gray")
        valid = dN.dropna(subset=["chi2_ndf"])

        if len(valid) >= 3:
            # stepped line: each position is a horizontal step
            ax.step(valid["x_gun_mm"].values, valid["chi2_ndf"].values,
                    where="mid", color=col, lw=1.5, label=f"χ²/NDF (N={N})")

            # highlight points above warning threshold with vertical red marker
            above = valid[valid["chi2_ndf"] > CHI2_WARN]
            if not above.empty:
                ax.scatter(above["x_gun_mm"].values, above["chi2_ndf"].values,
                           marker="*", color="red", s=80, zorder=5,
                           label=f"χ²/NDF > {CHI2_WARN:.0f}")

            # recalculate flagging fraction after filter
            n_above = len(above)
            frac    = n_above / len(valid) if len(valid) > 0 else 0.0
            updated_flags[N] = frac > CORE_FRAC   # True if still core_not_gaussian

        else:
            ax.text(0.5, 0.5,
                    f"N={N}: {len(valid)} valid points after eff filter",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=10, color="gray")
            updated_flags[N] = False

        # reference lines
        ax.axhline(1.0,      color="green", ls="-",  lw=1, alpha=0.5,
                   label="χ²/NDF = 1 (ideal)")
        ax.axhline(CHI2_WARN, color="red",  ls="--", lw=1.2,
                   label=f"warn = {CHI2_WARN:.0f}")

        ax.set_title(f"N={N}", fontsize=10)
        ax.set_xlabel("x_gun [mm]", fontsize=10)
        ax.set_ylabel("χ²/NDF", fontsize=10)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=8)
        _info(ax, group, N)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    saved = _save(fig, f"{MAT_SLUG}_{slug}_chi2ndf_N4_N20")
    return saved, updated_flags


# ─── Type C: overlay histograms + fits ────────────────────────────────────────

def _read_th1(rfile: uproot.ReadOnlyFile, name: str):
    """
    Read a TH1F from an open uproot file. Returns (centers, counts) or (None, None).
    """
    try:
        h = rfile[name]
        edges   = h.axis().edges()
        centers = 0.5 * (edges[:-1] + edges[1:])   # bin centres in ns
        counts  = h.values().astype(float)
        return centers, counts
    except Exception:
        return None, None


def plot_overlay_for_N(df: pd.DataFrame, group: str, slug: str,
                       rfile: uproot.ReadOnlyFile, N: int) -> list:
    """
    One figure per N: 3 subplots showing histograms + Gaussian core overlay.

    For each representative position in REPR_POS, if that position has
    eff >= EFF_FLOOR, use it; otherwise substitute the nearest valid position.
    """
    is_ast = (group == "TOP_NEAREST" and N == 20)

    # build filtered DataFrame for this N
    dN_filt = df[(df["N"] == N) & (df["efficiency"] >= EFF_FLOOR)].copy()

    if len(dN_filt) < 2:
        # not enough valid positions for an overlay — skip
        print(f"    skip overlay {group} N={N}: only {len(dN_filt)} valid positions")
        return []

    # determine which x positions to show (fallback to nearest valid if needed)
    show_pos = []
    for target_x in REPR_POS:
        # check if target_x itself is valid
        row_exact = dN_filt[dN_filt["x_gun_mm"] == target_x]
        if not row_exact.empty:
            show_pos.append(int(target_x))   # use the requested position
        else:
            # fall back to nearest valid position
            nearest = _nearest_valid_x(dN_filt, target_x)
            if nearest is not None:
                show_pos.append(nearest)

    # deduplicate while preserving order
    show_pos = list(dict.fromkeys(show_pos))

    if not show_pos:
        return []

    ncols = len(show_pos)
    fig, axes = plt.subplots(1, ncols, figsize=(5.5 * ncols, 5))
    if ncols == 1:
        axes = [axes]

    fig.suptitle(
        f"{MAT_LABEL} EndTop — {group} | histogram + Gaussian core fit  [N={N}]"
        + (" *" if is_ast else ""),
        fontsize=12, fontweight="bold"
    )

    fit_k = cfg.get("FIT_WINDOW_SIGMAS", 2.0)   # window = peak ± k*σ_MAD

    for ax, x_mm in zip(axes, show_pos):
        row = dN_filt[dN_filt["x_gun_mm"] == x_mm]
        if row.empty:
            ax.set_title(f"x = {x_mm:+d} mm — no data", fontsize=9)
            continue
        row = row.iloc[0]   # single row

        sf_ps  = row["sigma_fit_ps"]       # ns → ps (already in ps in CSV)
        sf_err = row["sigma_fit_err_ps"]
        mu_ns  = row["mu_fit_ns"]
        sm_ps  = row["sigma_mad_ps"]
        ss_ps  = row["sigma_std_ps"]
        chi2   = row["chi2_ndf"]
        eff    = row["efficiency"]

        # get histogram from ROOT sidecar
        x_tag  = f"x{'m' if x_mm < 0 else 'p'}{abs(x_mm)}mm"
        h_name = f"h_N{N}_{x_tag}"
        centers, counts = _read_th1(rfile, h_name)

        if centers is None:
            ax.text(0.5, 0.5, "TH1 not found", ha="center", va="center",
                    transform=ax.transAxes, color="gray")
            ax.set_title(f"x={x_mm:+d}mm N={N}", fontsize=9)
            continue

        # ── plot histogram (convert ns to ps for display) ──────────────────
        edges_ps = np.concatenate([[centers[0] - (centers[1]-centers[0])/2],
                                   centers + (centers[1]-centers[0])/2]) * 1000
        ax.stairs(counts, edges_ps, fill=True, color="#4878CF", alpha=0.6,
                  label=f"data  n={int(counts.sum())}")

        # ── fit window shading ─────────────────────────────────────────────
        if not (math.isnan(mu_ns) or math.isnan(sm_ps)):
            win_lo_ps = (mu_ns - fit_k * sm_ps / 1000) * 1000   # sm_ps is already ps
            win_hi_ps = (mu_ns + fit_k * sm_ps / 1000) * 1000
            ax.axvspan(win_lo_ps, win_hi_ps, alpha=0.12, color="red", zorder=1,
                       label=f"window [μ ± {fit_k:.0f}·σ_MAD]")

        # ── Gaussian fit curve (reconstructed from CSV scalars) ────────────
        if not (math.isnan(sf_ps) or math.isnan(mu_ns)):
            sig_ns = sf_ps / 1000           # ps → ns
            # amplitude: scale Gaussian so it matches the histogram peak near mu
            bin_width_ns = float(centers[1] - centers[0]) if len(centers) > 1 else 0.1
            A = counts.sum() * bin_width_ns / (sig_ns * np.sqrt(2 * np.pi))

            t_plot_ns  = np.linspace(centers.min(), centers.max(), 400)
            t_plot_ps  = t_plot_ns * 1000
            y_gauss    = A * np.exp(-0.5 * ((t_plot_ns - mu_ns) / sig_ns) ** 2)

            ax.plot(t_plot_ps, y_gauss, color="red", lw=2.0, zorder=4,
                    label=(f"Gauss fit:\n"
                           f"  μ = {mu_ns*1000:.1f} ± {row['mu_fit_err_ns']*1000:.1f} ps\n"
                           f"  σ = {sf_ps:.1f} ± {sf_err:.1f} ps\n"
                           f"  χ²/NDF = {chi2:.2f}\n"
                           f"  eff = {eff:.2f}" +
                           ("\n  * single-SiPM, high-N" if is_ast else "")))

            # vertical reference lines (all in ps)
            ax.axvline(mu_ns*1000,               color="red",    ls="--", lw=1.2, alpha=0.7)
            ax.axvline((mu_ns + sig_ns)*1000,    color="red",    ls=":",  lw=1.0)
            ax.axvline((mu_ns - sig_ns)*1000,    color="red",    ls=":",  lw=1.0,
                       label=f"±σ_fit = {sf_ps:.0f} ps")
            if not math.isnan(sm_ps):
                ax.axvline((mu_ns + sm_ps/1000)*1000, color="orange", ls="-.", lw=1.0)
                ax.axvline((mu_ns - sm_ps/1000)*1000, color="orange", ls="-.", lw=1.0,
                           label=f"σ_MAD = {sm_ps:.0f} ps")
            if not math.isnan(ss_ps):
                ax.axvline((mu_ns + ss_ps/1000)*1000, color="green",  ls=":",  lw=0.8, alpha=0.5)
                ax.axvline((mu_ns - ss_ps/1000)*1000, color="green",  ls=":",  lw=0.8, alpha=0.5,
                           label=f"σ_std = {ss_ps:.0f} ps")

        sub_title = f"{group} @ x = {x_mm:+d} mm | N={N}"
        if x_mm not in REPR_POS:
            sub_title += f"\n(nearest to {min(REPR_POS, key=lambda p: abs(p-x_mm)):+d} mm)"
        ax.set_title(sub_title, fontsize=9)
        ax.set_xlabel("t_N [ps]", fontsize=9)
        ax.set_ylabel("events / bin", fontsize=9)
        ax.legend(loc="upper left", fontsize=7, framealpha=0.85)
        ax.grid(True, alpha=0.2)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return _save(fig, f"{MAT_SLUG}_{slug}_overlay_pos_m690_0_p690_N{N}")


# ─── main loop ────────────────────────────────────────────────────────────────

def main():
    import argparse as _ap
    _parser = _ap.ArgumentParser(description="Regenerate figures with eff-filter")
    _parser.add_argument("--material", default="ej230",
                         choices=["ej230", "ej204"],
                         help="Material to process (ej230 or ej204)")
    _args = _parser.parse_args()

    # ── set module-level material variables from CLI argument ────────────────
    global MAT_SLUG, MAT_LABEL, BRANCH, SHA
    _MAT_CFG = {
        "ej230": ("ej230_endtop", "EJ-230", "feat/ej230-sslg4", "ca2f1c3"),
        "ej204": ("ej204_endtop", "EJ-204", "feat/endtop-sslg4", "a0368c4"),
    }
    MAT_SLUG, MAT_LABEL, BRANCH, SHA = _MAT_CFG[_args.material]

    all_saved = []
    updated_flag_status = {}    # {(group, N): bool}
    n_removed_by_group  = {}    # {group: {N: int}}

    print(f"\n{'='*60}")
    print(f"  Figure regeneration — {MAT_LABEL}  eff≥5% filter")
    print(f"{'='*60}\n")

    for group in GROUPS:
        slug     = SLUG_MAP[group]
        csv_path = OUTDIR / f"{MAT_SLUG}_{slug}_results.csv"
        root_path= OUTDIR / f"{MAT_SLUG}_{slug}_sidecars.root"

        if not csv_path.exists():
            print(f"  WARN: CSV not found for {group}, skipping")
            continue

        df = pd.read_csv(csv_path)
        n_removed_by_group[group] = {}

        print(f"  {group} …", end=" ", flush=True)

        # count filtered points per N
        for N in N_VALUES:
            dN_all  = df[df["N"] == N]
            dN_filt = dN_all[dN_all["efficiency"] >= EFF_FLOOR]
            n_removed_by_group[group][N] = len(dN_all) - len(dN_filt)

        # Type A: resolution vs position
        saved_a = plot_resolution(df, group, slug)
        all_saved.extend(saved_a)

        # Type B: χ²/NDF diagnostic (returns updated flags)
        saved_b, upd_flags = plot_chi2ndf(df, group, slug)
        all_saved.extend(saved_b)
        for N, flagged in upd_flags.items():
            updated_flag_status[(group, N)] = flagged

        # Type C: overlay — needs ROOT sidecar for TH1
        if root_path.exists():
            with uproot.open(str(root_path)) as rfile:
                for N in N_VALUES:
                    saved_c = plot_overlay_for_N(df, group, slug, rfile, N)
                    all_saved.extend(saved_c)
        else:
            print(f"(no ROOT sidecar found for {group})", end=" ")

        print("done")

    # ─── summary ──────────────────────────────────────────────────────────────
    still_flagged = [(g, N) for (g, N), v in updated_flag_status.items() if v]

    lines = [
        "═" * 65,
        f"Figure Regeneration Summary — {MAT_LABEL} Figures with eff≥5% Filter",
        "═" * 65,
        "",
        "Filters applied:",
        f"  - Omit positions with efficiency < {EFF_FLOOR:.0%}",
        "  - Mark TOP_NEAREST N=20 with asterisk (*) in legend",
        "",
        f"Groups processed: {len(GROUPS)}",
        f"Figures regenerated: {len(all_saved)} files",
        "",
        "Points removed per group (eff < 5%):",
    ]
    for g in GROUPS:
        for N in N_VALUES:
            nr = n_removed_by_group.get(g, {}).get(N, 0)
            if nr > 0:
                lines.append(f"  {g:20s} N={N}: {nr} / {N_POS} points removed")
    lines += [
        "",
        f"Groups still flagged 'core_not_gaussian' after filter: "
        f"{[f'{g}/N={N}' for g,N in still_flagged] or 'none'}",
        "",
        "Updated figures:",
        f"  {str(OUTDIR)}/{MAT_SLUG}_*.{{png,pdf}}",
        "",
        "Ready for CP4 (EJ-204 pipeline).",
        "═" * 65,
    ]
    summary_text = "\n".join(lines)
    print("\n" + summary_text)

    sp = OUTDIR / "cp3_5_summary_ej230_figures_updated.txt"
    sp.write_text(summary_text)
    print(f"\n  Summary: {sp}")
    print(f"\n  Total files: {len(all_saved)} PNG/PDF")
    for p in all_saved[:6]:
        print(f"    {p}")
    if len(all_saved) > 6:
        print(f"    … ({len(all_saved)-6} more)")


if __name__ == "__main__":
    main()
