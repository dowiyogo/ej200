#!/usr/bin/env python3.12
"""
sidecar.py — Write the three-file sidecar for each (material, group) pair.

Three files per group:
    .root      — TH1F histograms + TF1 fit curves + TGraphErrors(σ_fit vs x)
    .csv       — all numeric results, 31 rows × 2 N values per group
    .meta.json — provenance: input SHA-256, fit config, branch, timestamp

Rules:
    - NEVER overwrite an existing .root (abort on collision).
    - meta.json is written last so its existence signals a complete run.
"""

import json
import hashlib
import datetime
import math
from pathlib import Path

import numpy as np
import pandas as pd
import ROOT


# ─── collision guard ──────────────────────────────────────────────────────────

def _check_no_overwrite(path: str) -> None:
    """Abort if file already exists — protects against silent data loss."""
    if Path(path).exists():
        raise SystemExit(
            f"[SIDECAR ABORT] Would overwrite {path}. "
            "Delete or rename the existing file before re-running."
        )


# ─── ROOT sidecar ─────────────────────────────────────────────────────────────

def write_root_sidecar(results: dict, group_slug: str, material_slug: str,
                        output_dir: str) -> str:
    """
    Write one ROOT file per group with histograms, fits, and σ-vs-x graphs.

    Parameters
    ----------
    results      : dict  {N: [{result_dict, 'x_mm': int, 't_N': ndarray}, ...]}
    group_slug   : str   e.g. "top_sum4"
    material_slug: str   e.g. "ej230_endtop"
    output_dir   : str

    Returns
    -------
    str — absolute path of the written .root file
    """
    outpath = str(Path(output_dir) / f"{material_slug}_{group_slug}_sidecars.root")
    _check_no_overwrite(outpath)

    tfile = ROOT.TFile(outpath, "RECREATE")   # new file — _check_no_overwrite guarantees it's new

    for N, rows in sorted(results.items()):
        # ── TGraphErrors: σ_fit [ps] vs x_gun [mm] for this N ─────────────────
        valid = [r for r in rows
                 if r.get("flag") == "ok" and not math.isnan(r.get("sigma_fit", float("nan")))]

        if valid:
            graph = ROOT.TGraphErrors(len(valid))   # preallocate with correct size
            graph.SetName(f"g_sigma_N{N}")
            graph.SetTitle(
                f"#sigma_{{fit}} vs x_{{gun}} | N={N};"
                f"x_{{gun}} [mm];#sigma_{{fit}} [ps]"
            )
            for i, r in enumerate(valid):
                sig_ps  = r["sigma_fit"] * 1000      # ns → ps
                # reported error = max(TF1 error, bootstrap error)
                err_ps  = max(r.get("sigma_fit_err", 0)*1000,
                              r.get("bootstrap_err",  0)*1000)
                graph.SetPoint(i, float(r["x_mm"]), sig_ps)      # (x, y)
                graph.SetPointError(i, 0.0, err_ps)               # (ex=0, ey)
            graph.Write()

        # ── per-position TH1F and TF1 ─────────────────────────────────────────
        for row in rows:
            x_mm  = row["x_mm"]
            x_tag = f"x{'m' if x_mm < 0 else 'p'}{abs(x_mm)}mm"  # "xm690mm" / "xp0mm"

            h = row.get("h_root")    # TH1F from fit_engine (None if fit failed)
            f = row.get("f_root")    # TF1   from fit_engine

            if h is not None:
                h.SetName(f"h_N{N}_{x_tag}")
                h.SetTitle(
                    f"{group_slug.upper()} N={N} x={x_mm} mm;"
                    f"t_{{N}} [ns];events / bin"
                )
                h.Write()       # serialize TH1F to the open TFile

            if f is not None:
                f.SetName(f"f_N{N}_{x_tag}")
                f.SetTitle(f"Gauss core fit | N={N} x={x_mm} mm")
                f.Write()       # serialize TF1 to the open TFile

    tfile.Close()   # flush and release
    print(f"    → ROOT sidecar: {outpath}")
    return outpath


# ─── CSV sidecar ──────────────────────────────────────────────────────────────

def write_csv_sidecar(results: dict, group_slug: str, material_slug: str,
                       output_dir: str) -> str:
    """
    Write a CSV of all numeric fit results (one row per position per N).

    Columns: x_gun_mm, group, N, mu_fit_ns, mu_fit_err_ns,
             sigma_fit_ns, sigma_fit_err_ns, sigma_fit_ps, sigma_fit_err_ps,
             sigma_mad_ps, sigma_std_ps, bootstrap_err_ps,
             chi2_ndf, fit_status, efficiency, n_events, flag.
    """
    outpath = str(Path(output_dir) / f"{material_slug}_{group_slug}_results.csv")

    out_rows = []
    for N, rows in sorted(results.items()):
        for row in rows:
            nan = float("nan")
            sf  = row.get("sigma_fit",     nan)     # ns
            sfe = row.get("sigma_fit_err", nan)     # ns
            sm  = row.get("sigma_mad",     nan)     # ns
            ss  = row.get("sigma_std",     nan)     # ns
            be  = row.get("bootstrap_err", nan)     # ns

            def _ps(v):
                """Convert ns to ps, returning nan if input is nan."""
                return v * 1000 if not math.isnan(v) else nan

            out_rows.append({
                "x_gun_mm":         row.get("x_mm",        nan),
                "group":            group_slug,
                "N":                N,
                "mu_fit_ns":        row.get("mu_fit",       nan),
                "mu_fit_err_ns":    row.get("mu_fit_err",   nan),
                "sigma_fit_ns":     sf,
                "sigma_fit_err_ns": sfe,
                "sigma_fit_ps":     _ps(sf),
                "sigma_fit_err_ps": _ps(sfe),
                "sigma_mad_ps":     _ps(sm),
                "sigma_std_ps":     _ps(ss),
                "bootstrap_err_ps": _ps(be),
                "chi2_ndf":         row.get("chi2_ndf",     nan),
                "fit_status":       row.get("fit_status",   -1),
                "efficiency":       row.get("efficiency",   nan),
                "n_events":         row.get("n_events",     0),
                "flag":             row.get("flag",         "unknown"),
            })

    df = pd.DataFrame(out_rows)
    df.to_csv(outpath, index=False, float_format="%.6g")
    print(f"    → CSV sidecar:  {outpath}  ({len(df)} rows)")
    return outpath


# ─── JSON metadata sidecar ────────────────────────────────────────────────────

def write_metadata_json(results: dict, group_slug: str, material_slug: str,
                         cfg: dict, input_files: list, runtime_sha: str,
                         geant4_version: str, output_dir: str) -> str:
    """
    Write a JSON provenance record including SHA-256 of each input ROOT file.

    Parameters
    ----------
    results        : dict
    group_slug     : str
    material_slug  : str
    cfg            : dict — exec16_config.yaml
    input_files    : list — absolute paths to the 31 input ROOT files for this material
    runtime_sha    : str  — git SHA verified at runtime
    geant4_version : str
    output_dir     : str
    """
    outpath = str(Path(output_dir) / f"{material_slug}_{group_slug}_metadata.json")

    # compute SHA-256 for each input file: allows tracing outputs to exact simulation run
    input_sha256 = {}
    for fpath in sorted(str(f) for f in input_files):
        h = hashlib.sha256()
        with open(fpath, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)         # stream in 64 kB chunks to control RAM
        input_sha256[Path(fpath).name] = h.hexdigest()

    # determine material key in MATERIALS dict
    mat_key = "EJ-230" if "ej230" in material_slug else "EJ-204"
    mat_cfg = cfg.get("MATERIALS", {}).get(mat_key, {})

    n_valid = sum(1 for N, rows in results.items()
                  for row in rows if row.get("flag") == "ok")
    n_total = sum(len(rows) for rows in results.values())

    meta = {
        "exec_tag":          cfg.get("EXEC_TAG", "EXEC_16"),
        "material":          material_slug,
        "group":             group_slug,
        "branch":            mat_cfg.get("branch",       "unknown"),
        "sha_declared":      mat_cfg.get("expected_sha", "unknown"),
        "sha_runtime":       runtime_sha,
        "n_positions":       31,
        "scan_positions_mm": cfg.get("SCAN_POSITIONS_MM", []),
        "n_values":          cfg.get("N_VALUES", [4, 20]),
        "n_valid_fits":      n_valid,
        "n_total_fits":      n_total,
        "fit_config": {
            "fit_window_sigmas":  cfg.get("FIT_WINDOW_SIGMAS",  2.0),
            "min_events_for_fit": cfg.get("MIN_EVENTS_FOR_FIT", 30),
            "binning_strategy":   cfg.get("BINNING_STRATEGY",   "sqrt_n"),
            "n_bootstrap":        cfg.get("N_BOOTSTRAP",        300),
            "random_seed":        cfg.get("RANDOM_SEED",        20260618),
            "fit_options":        cfg.get("FIT_OPTIONS",        "R Q S 0"),
            "chi2_ndf_warn":      cfg.get("CHI2_NDF_WARN",      3.0),
        },
        "input_root_base":   cfg.get("INPUT_ROOT_BASE", ""),
        "ttree_name":        cfg.get("TTREE_NAME", "sipm_hits"),
        "input_sha256":      input_sha256,
        "geant4_version":    geant4_version,
        "timestamp":         datetime.datetime.now().isoformat(),
    }

    Path(outpath).write_text(json.dumps(meta, indent=2))
    print(f"    → JSON metadata:{outpath}")
    return outpath
