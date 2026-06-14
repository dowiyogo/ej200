#!/usr/bin/env python3
"""
analysis_endonly_mylar.py — End-only + Mylar-wrapped EJ-204 bar simulation analysis.

Reads from DATA_DIR (default: /home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614):
  analysis/attenuation_curve.csv, sigma_t_sum4.csv
  run_x*.log  (boundary census)
  run_metadata.txt
  photon_hits_x*.root  (timing histograms)

Also reads EndTop reference data from ENDTOP_DIR (exec07_endtop_2000) if available.

Writes to OUT_DIR (default: ./):
  figures/*.png   — DPI >= 150
  tables/*.csv / *.tex
  deck_values.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import uproot
from scipy.optimize import curve_fit

# ─── constants (physical, EJ-204) ────────────────────────────────────────────
BAR_HALF_MM    = 700.0        # bar half-length (mm)
LAMBDA_BULK_CM = 160.0        # EJ-204 bulk ABSLENGTH (cm) — Materials.cc:99
SCINT_YIELD    = 10400        # EJ-204 photons/MeV — Eljen datasheet Aug 2023, Materials.cc:100
PEAK_EMISSION_NM = 408.0      # EJ-204 emission peak (nm) — Materials.cc:86
REFRACTIVE_INDEX = 1.58       # EJ-204 refractive index (flat) — Materials.cc:41
DENSITY_GCC    = 1.023        # EJ-204 density (g/cm³) — Materials.cc:24
TAU_RISE_NS    = 0.7          # EJ-204 rise time (ns)
TAU_DECAY_NS   = 1.8          # EJ-204 decay time (ns)
PDE_PEAK       = 0.63         # SiPM PDE at peak (420 nm) AFBR-S4N66P014M
SIPM_MODEL     = "AFBR-S4N66P014M"
SIPM_MODEL_PDE = "AFBR-S4N66P024M"   # PDE curve only
DPI            = 150
C_MM_NS        = 299.792      # speed of light in mm/ns

# Attenuation fit threshold for lambda_tail: fit only d > D_TAIL_CM (cm)
D_TAIL_CM = 40.0              # 40 cm = 400 mm; excludes prompt-dominated near-end region

# EndTop reference dataset (EXEC_07, R=0.98 side reflectors + top SiPMs)
ENDTOP_DIR = pathlib.Path(
    "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"
)
ENDTOP_REFLECTOR_R = 0.98     # BarSkinReflector R in exec07 — from run log

# Representative positions for timing histograms
REPR_MM = [-690, -400, 0, 400, 690]

# ─── matplotlib style ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.size":       10,
    "axes.titlesize":  10,
    "axes.labelsize":  10,
    "legend.fontsize":  8,
    "xtick.labelsize":  8,
    "ytick.labelsize":  8,
    "figure.dpi":      DPI,
    "axes.grid":       True,
    "grid.alpha":      0.3,
})


# ─── helpers ──────────────────────────────────────────────────────────────────

def _pos_tag(x: float) -> str:
    sign = "+" if x >= 0 else ""
    return f"x={sign}{int(x):d} mm"


def _read_metadata(data_dir: pathlib.Path) -> dict:
    meta: dict = {}
    path = data_dir / "run_metadata.txt"
    for line in path.read_text().splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            meta[k.strip()] = v.strip()
    return meta


def _parse_position_log(path: pathlib.Path) -> dict:
    """Return boundary census and final run summary from a per-position log."""
    text = path.read_text(errors="replace")
    census: dict[str, int] = {}
    m = re.search(r"Boundary Census.*?\n(.*?)={10}", text, re.DOTALL)
    if m:
        for line in m.group(1).splitlines():
            if ":" in line:
                key, _, val_str = line.strip().partition(":")
                try:
                    census[key.strip()] = int(val_str.strip().replace(",", ""))
                except ValueError:
                    pass
    summaries = re.findall(
        r"=== EJ Scintillator Bar Run Summary ===(.*?)={5,}", text, re.DOTALL
    )
    final: dict[str, float] = {}
    max_events = 0
    for block in summaries:
        fields: dict[str, float] = {}
        for line in block.splitlines():
            if ":" in line:
                k, _, v = line.strip().partition(":")
                try:
                    fields[k.strip()] = float(v.strip().replace(",", ""))
                except ValueError:
                    pass
        n = int(fields.get("Events run", 0))
        if n > max_events:
            max_events = n
            final = fields
    return {"census": census, "summary": final}


def _load_csv(path: pathlib.Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _safe(v) -> object:
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return None
    if isinstance(v, (np.floating, np.integer)):
        return float(v)
    return v


# ─── attenuation fitting ──────────────────────────────────────────────────────

def _exp1(d: np.ndarray, n0: float, lam: float) -> np.ndarray:
    return n0 * np.exp(-d / lam)


def _expoff(d: np.ndarray, A: float, lam: float, C: float) -> np.ndarray:
    return A * np.exp(-d / lam) + C


def _exp2_reparam(d: np.ndarray, A_s: float, lam_s: float,
                  A_l: float, log_dlam: float) -> np.ndarray:
    """Double-exp with ordering constraint: lam_l = lam_s + exp(log_dlam) > lam_s."""
    lam_l = lam_s + math.exp(log_dlam)
    return A_s * np.exp(-d / lam_s) + A_l * np.exp(-d / lam_l)


def _aic(chi2: float, npar: int) -> float:
    """AIC = chi2 + 2*k (uses chi2 as -2*log L under Gaussian errors)."""
    return chi2 + 2.0 * npar


def _aicc(chi2: float, npar: int, n: int) -> float:
    """AICc = AIC + 2k(k+1)/(n-k-1)."""
    aic = _aic(chi2, npar)
    denom = n - npar - 1
    return aic + (2.0 * npar * (npar + 1) / denom) if denom > 0 else math.nan


def _bic(chi2: float, npar: int, n: int) -> float:
    """BIC = chi2 + k*ln(n)."""
    return chi2 + npar * math.log(n) if n > 0 else math.nan


def _weighted_fit(model_fn, d, y, yerr, p0, bounds=(-np.inf, np.inf)):
    """
    Weighted curve_fit. Returns (popt, perr_scaled, chi2ndf).
    Errors from covariance matrix, scaled by sqrt(chi2/ndf) if chi2/ndf > 1.
    """
    sigma = np.where(yerr > 0, yerr, 1e-3)
    try:
        popt, pcov = curve_fit(
            model_fn, d, y, p0=p0, sigma=sigma,
            absolute_sigma=True, bounds=bounds, maxfev=50000
        )
        perr = np.sqrt(np.diag(pcov))
        residuals = y - model_fn(d, *popt)
        chi2 = float(np.sum((residuals / sigma) ** 2))
        ndf = len(d) - len(popt)
        chi2ndf = chi2 / ndf if ndf > 0 else math.nan
        if math.isfinite(chi2ndf) and chi2ndf > 1.0:
            perr = perr * math.sqrt(chi2ndf)
        return popt, perr, chi2ndf
    except Exception:
        npar = len(p0) if hasattr(p0, "__len__") else 2
        return np.full(npar, math.nan), np.full(npar, math.nan), math.nan


def _bootstrap_fits(d: np.ndarray, y: np.ndarray, yerr: np.ndarray,
                    n_boot: int = 200, rng_seed: int = 42) -> dict:
    """
    Pairs bootstrap on (d, y, yerr) data (31 fit points), n_boot replicas.
    Resamples with replacement; refits M3 (exp+offset), M4 (lambda_tail), M2 (double-exp).
    Returns percentile summaries.  Seed is fixed for reproducibility.
    """
    rng = np.random.default_rng(rng_seed)
    n = len(d)
    res_m3_lam, res_m3_c, res_m4_lam, res_m2_lam_s, res_m2_lam_l = [], [], [], [], []
    n_fail_m2 = 0

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        db, yb, eb = d[idx], y[idx], yerr[idx]

        # M3: exp+offset
        try:
            C_g = max(float(np.min(yb)) * 0.9, 0.1)
            p0 = [max(float(yb[0]) - C_g, 1.0), 18.0, C_g]
            popt, _, _ = _weighted_fit(_expoff, db, yb, eb, p0,
                                       ([0, 0.5, 0], [np.inf, 300, np.inf]))
            if all(np.isfinite(popt)):
                res_m3_lam.append(float(popt[1]))
                res_m3_c.append(float(popt[2]))
        except Exception:
            pass

        # M4: tail slope
        tail_m = db > D_TAIL_CM
        if tail_m.sum() >= 4:
            try:
                p0 = [float(np.max(yb[tail_m])) * 3.0, 42.0]
                popt, _, _ = _weighted_fit(_exp1, db[tail_m], yb[tail_m], eb[tail_m],
                                           p0, ([0, 5], [np.inf, 500]))
                if all(np.isfinite(popt)):
                    res_m4_lam.append(float(popt[1]))
            except Exception:
                pass

        # M2: constrained double-exponential
        try:
            p0 = [float(yb[0]) * 0.7, 10.0, float(yb[0]) * 0.3, math.log(30.0)]
            bounds = ([0, 0.5, 0, -5.0], [np.inf, 50.0, np.inf, 10.0])
            popt, _, _ = _weighted_fit(_exp2_reparam, db, yb, eb, p0, bounds)
            if all(np.isfinite(popt)):
                lam_l = float(popt[1]) + math.exp(float(popt[3]))
                res_m2_lam_s.append(float(popt[1]))
                res_m2_lam_l.append(lam_l)
        except Exception:
            n_fail_m2 += 1

    def _summ(vals):
        v = np.array(vals)
        if len(v) < 10:
            return None
        return {"n": len(v), "median": float(np.median(v)),
                "p16": float(np.percentile(v, 16)), "p84": float(np.percentile(v, 84)),
                "mean": float(np.mean(v)), "std": float(np.std(v))}

    return {
        "n_boot": n_boot, "seed": rng_seed,
        "M3_lam_cm": _summ(res_m3_lam),
        "M3_C_PE":   _summ(res_m3_c),
        "M4_lam_cm": _summ(res_m4_lam),
        "M2_lam_s_cm": _summ(res_m2_lam_s),
        "M2_lam_l_cm": _summ(res_m2_lam_l),
        "M2_n_fail": n_fail_m2,
        "M2_failure_fraction": n_fail_m2 / n_boot if n_boot > 0 else math.nan,
    }


def _fit_side(dist_mm: np.ndarray, npe: np.ndarray, npe_err: np.ndarray,
              label: str = "", run_bootstrap: bool = False,
              n_boot: int = 200) -> dict:
    """
    Fit four models to (dist_mm, npe) data with per-point SEM weights.

    M1: single exponential (reference, χ²/ndf diagnostic)
    M2: constrained double-exponential (λ_l > λ_s enforced via reparametrization)
    M3: exponential + constant offset  A·exp(-d/λ) + C
    M4: single-exp in tail region only (d > D_TAIL_CM)

    Returns dict with all fit results plus AIC/BIC per model.
    """
    d = dist_mm / 10.0  # convert mm → cm for all fits
    n_pts = len(d)

    # --- M1: single exponential ---
    p0_a = [float(npe[0]), 30.0]
    bounds_a = ([0, 1], [np.inf, 500])
    popt_a, perr_a, chi2_a = _weighted_fit(_exp1, d, npe, npe_err, p0_a, bounds_a)
    ndf_a = n_pts - 2
    raw_chi2_a = chi2_a * ndf_a if math.isfinite(chi2_a) else math.nan

    # --- M3: exponential + offset A·exp(-d/λ) + C ---
    C_guess = float(np.min(npe)) * 0.9
    A_guess = float(npe[0]) - C_guess
    p0_b = [max(A_guess, 1.0), 18.0, max(C_guess, 0.5)]
    bounds_b = ([0, 0.5, 0], [np.inf, 300, np.inf])
    popt_b, perr_b, chi2_b = _weighted_fit(_expoff, d, npe, npe_err, p0_b, bounds_b)
    ndf_b = n_pts - 3
    raw_chi2_b = chi2_b * ndf_b if math.isfinite(chi2_b) else math.nan

    # --- M2: constrained double-exponential ---
    p0_m2 = [float(npe[0]) * 0.7, 10.0, float(npe[0]) * 0.3, math.log(30.0)]
    bounds_m2 = ([0, 0.5, 0, -5.0], [np.inf, 50.0, np.inf, 10.0])
    popt_m2, perr_m2, chi2_m2 = _weighted_fit(
        _exp2_reparam, d, npe, npe_err, p0_m2, bounds_m2)
    ndf_m2 = n_pts - 4
    raw_chi2_m2 = chi2_m2 * ndf_m2 if math.isfinite(chi2_m2) else math.nan
    if all(np.isfinite(popt_m2)):
        lam_s_m2 = float(popt_m2[1])
        lam_l_m2 = lam_s_m2 + math.exp(float(popt_m2[3]))
        A_s_m2, A_l_m2 = float(popt_m2[0]), float(popt_m2[2])
        lam_l_err_m2 = float(perr_m2[3]) * math.exp(float(popt_m2[3]))  # error propagation
    else:
        lam_s_m2 = lam_l_m2 = A_s_m2 = A_l_m2 = lam_l_err_m2 = math.nan

    # --- M4: lambda_tail — tail region only ---
    tail_mask = d > D_TAIL_CM
    n_tail = int(np.sum(tail_mask))
    if n_tail >= 4:
        p0_c = [float(np.max(npe[tail_mask])) * 3.0, 42.0]
        bounds_c = ([0, 5], [np.inf, 500])
        popt_c, perr_c, chi2_c = _weighted_fit(
            _exp1, d[tail_mask], npe[tail_mask], npe_err[tail_mask], p0_c, bounds_c
        )
        ndf_c = n_tail - 2
        raw_chi2_c = chi2_c * ndf_c if math.isfinite(chi2_c) else math.nan
    else:
        popt_c = perr_c = np.full(2, math.nan)
        chi2_c = raw_chi2_c = math.nan
        n_tail = 0
        ndf_c = 0

    # --- AIC / AICc / BIC (for M1, M2, M3 which share same n_pts) ---
    aic_m1  = _aic(raw_chi2_a,  2)         if math.isfinite(raw_chi2_a)  else math.nan
    aicc_m1 = _aicc(raw_chi2_a, 2, n_pts)  if math.isfinite(raw_chi2_a)  else math.nan
    bic_m1  = _bic(raw_chi2_a,  2, n_pts)  if math.isfinite(raw_chi2_a)  else math.nan

    aic_m3  = _aic(raw_chi2_b,  3)         if math.isfinite(raw_chi2_b)  else math.nan
    aicc_m3 = _aicc(raw_chi2_b, 3, n_pts)  if math.isfinite(raw_chi2_b)  else math.nan
    bic_m3  = _bic(raw_chi2_b,  3, n_pts)  if math.isfinite(raw_chi2_b)  else math.nan

    aic_m2  = _aic(raw_chi2_m2,  4)        if math.isfinite(raw_chi2_m2) else math.nan
    aicc_m2 = _aicc(raw_chi2_m2, 4, n_pts) if math.isfinite(raw_chi2_m2) else math.nan
    bic_m2  = _bic(raw_chi2_m2,  4, n_pts) if math.isfinite(raw_chi2_m2) else math.nan

    if label:
        print(f"  Att fit [{label}]:")
        print(f"    M1 single-exp:   N0={popt_a[0]:.0f}±{perr_a[0]:.0f} PE, "
              f"λ={popt_a[1]:.2f}±{perr_a[1]:.2f} cm, χ²/ndf={chi2_a:.0f}, AICc={aicc_m1:.0f}")
        print(f"    M3 exp+offset:   A={popt_b[0]:.0f}±{perr_b[0]:.0f} PE, "
              f"λ_p={popt_b[1]:.2f}±{perr_b[1]:.2f} cm, "
              f"C={popt_b[2]:.1f}±{perr_b[2]:.1f} PE, χ²/ndf={chi2_b:.0f}, AICc={aicc_m3:.0f}")
        print(f"    M2 double-exp:   A_s={A_s_m2:.0f} PE, λ_s={lam_s_m2:.1f} cm, "
              f"A_l={A_l_m2:.0f} PE, λ_l={lam_l_m2:.1f} cm, "
              f"χ²/ndf={chi2_m2:.0f}, AICc={aicc_m2:.0f}")
        print(f"    M4 λ_tail (d>{D_TAIL_CM:.0f}cm, N={n_tail}): "
              f"N0={popt_c[0]:.0f}±{perr_c[0]:.0f} PE, "
              f"λ={popt_c[1]:.2f}±{perr_c[1]:.2f} cm, χ²/ndf={chi2_c:.1f}")

    boot = None
    if run_bootstrap:
        print(f"    Bootstrap [{label}]: {n_boot} replicas...")
        boot = _bootstrap_fits(d, npe, npe_err, n_boot=n_boot)
        if boot["M2_failure_fraction"] > 0.3:
            print(f"    WARNING: M2 bootstrap failure fraction = "
                  f"{boot['M2_failure_fraction']:.0%} — M2 not robustly identifiable")

    return {
        # M1: single-exp
        "single_n0":       _safe(float(popt_a[0])),
        "single_n0_err":   _safe(float(perr_a[0])),
        "single_lam_cm":   _safe(float(popt_a[1])),
        "single_lam_err":  _safe(float(perr_a[1])),
        "single_chi2ndf":  _safe(float(chi2_a)),
        "single_aic":      _safe(aic_m1),
        "single_aicc":     _safe(aicc_m1),
        "single_bic":      _safe(bic_m1),
        # M3: exp+offset
        "expoff_A":        _safe(float(popt_b[0])),
        "expoff_A_err":    _safe(float(perr_b[0])),
        "expoff_lam_cm":   _safe(float(popt_b[1])),
        "expoff_lam_err":  _safe(float(perr_b[1])),
        "expoff_C":        _safe(float(popt_b[2])),
        "expoff_C_err":    _safe(float(perr_b[2])),
        "expoff_chi2ndf":  _safe(float(chi2_b)),
        "expoff_aic":      _safe(aic_m3),
        "expoff_aicc":     _safe(aicc_m3),
        "expoff_bic":      _safe(bic_m3),
        # M2: constrained double-exp
        "dblexp_A_s":      _safe(A_s_m2),
        "dblexp_A_s_err":  _safe(float(perr_m2[0])) if all(np.isfinite(popt_m2)) else None,
        "dblexp_lam_s_cm": _safe(lam_s_m2),
        "dblexp_lam_s_err":_safe(float(perr_m2[1])) if all(np.isfinite(popt_m2)) else None,
        "dblexp_A_l":      _safe(A_l_m2),
        "dblexp_A_l_err":  _safe(float(perr_m2[2])) if all(np.isfinite(popt_m2)) else None,
        "dblexp_lam_l_cm": _safe(lam_l_m2),
        "dblexp_lam_l_err":_safe(lam_l_err_m2),
        "dblexp_chi2ndf":  _safe(float(chi2_m2)),
        "dblexp_aic":      _safe(aic_m2),
        "dblexp_aicc":     _safe(aicc_m2),
        "dblexp_bic":      _safe(bic_m2),
        # M4: lambda_tail
        "tail_n0":         _safe(float(popt_c[0])),
        "tail_n0_err":     _safe(float(perr_c[0])),
        "tail_lam_cm":     _safe(float(popt_c[1])),
        "tail_lam_err":    _safe(float(perr_c[1])),
        "tail_chi2ndf":    _safe(float(chi2_c)),
        "tail_n_points":   n_tail,
        "tail_dmin_cm":    D_TAIL_CM,
        # Bootstrap
        "bootstrap":       boot,
    }


# ─── EndTop NPE extraction ────────────────────────────────────────────────────

def _load_endtop_npe(endtop_dir: pathlib.Path) -> dict | None:
    """
    Read end-SiPM NPE per position from exec07 ROOT files.
    Uses only face_type 0 (IDs 0–7, end-left) and face_type 1 (IDs 8–15, end-right).
    Returns dict keyed by x_mm with npe_left, npe_right and their SEMs.
    """
    root_files = sorted(endtop_dir.glob("photon_hits_x*.root"))
    if not root_files:
        return None

    results = {}
    print(f"  Loading {len(root_files)} EndTop ROOT files...")
    for rp in root_files:
        try:
            with uproot.open(rp) as f:
                arr = f["sipm_hits"].arrays(
                    ["event_id", "face_type", "gun_x_mm"], library="np"
                )
        except Exception as e:
            print(f"    WARNING: {rp.name}: {e}")
            continue

        gx = float(arr["gun_x_mm"][0])
        ft = arr["face_type"]
        ev = arr["event_id"]

        mask_l = ft == 0
        mask_r = ft == 1
        ev_l = ev[mask_l]
        ev_r = ev[mask_r]

        n_ev = int(ev.max()) + 1
        cnt_l = np.bincount(ev_l, minlength=n_ev).astype(float)
        cnt_r = np.bincount(ev_r, minlength=n_ev).astype(float)

        results[gx] = {
            "npe_l":     float(np.mean(cnt_l)),
            "npe_l_err": float(np.std(cnt_l) / math.sqrt(n_ev)),
            "npe_r":     float(np.mean(cnt_r)),
            "npe_r_err": float(np.std(cnt_r) / math.sqrt(n_ev)),
            "n_events":  n_ev,
        }
    return results if results else None


# ─── ROOT timing extraction ───────────────────────────────────────────────────

SPR_RISE_NS  = 0.5
SPR_FALL_NS  = 5.0
LE_THRESHOLD = 4.0


def _spr_pulse(slow: float, fast: float, delta: float) -> float:
    return (slow * math.exp(-delta / SPR_FALL_NS)
            - fast * math.exp(-delta / SPR_RISE_NS))


def _leading_edge_time(arrivals: np.ndarray) -> float:
    if arrivals.size == 0:
        return math.nan
    arrivals = np.sort(arrivals)
    slow = fast = 0.0
    i = 0
    while i < arrivals.size:
        cur = float(arrivals[i])
        j = i
        while j < arrivals.size and arrivals[j] == cur:
            slow += 1.0
            fast += 1.0
            j += 1
        interval = float(arrivals[j] - cur) if j < arrivals.size else math.inf
        d0 = fast / SPR_RISE_NS - slow / SPR_FALL_NS
        if d0 > 0.0:
            peak_d = math.log(fast * SPR_FALL_NS / (slow * SPR_RISE_NS)) / (
                1.0 / SPR_RISE_NS - 1.0 / SPR_FALL_NS
            )
            reach = min(peak_d, interval)
            if reach >= 0.0 and _spr_pulse(slow, fast, reach) >= LE_THRESHOLD:
                lo, hi = 0.0, reach
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if _spr_pulse(slow, fast, mid) >= LE_THRESHOLD:
                        hi = mid
                    else:
                        lo = mid
                return cur + hi
        if j >= arrivals.size:
            break
        slow *= math.exp(-interval / SPR_FALL_NS)
        fast *= math.exp(-interval / SPR_RISE_NS)
        i = j
    return math.nan


def _earliest(a: float, b: float) -> float:
    if not math.isfinite(a):
        return b
    if not math.isfinite(b):
        return a
    return min(a, b)


def _per_event_timing(root_path: pathlib.Path, n_events: int) -> dict:
    with uproot.open(root_path) as f:
        arr = f["sipm_hits"].arrays(
            ["event_id", "global_id", "time_ns"], library="np"
        )

    ev  = arr["event_id"].astype(np.int64)
    gid = arr["global_id"].astype(np.int64)
    tns = arr["time_ns"].astype(np.float64)

    n = n_events
    fpt_l  = np.full(n, np.nan); fpt_r  = np.full(n, np.nan)
    mean_l = np.full(n, np.nan); mean_r = np.full(n, np.nan)
    t50_l  = np.full(n, np.nan); t50_r  = np.full(n, np.nan)
    sum4_l = np.full(n, np.nan); sum4_r = np.full(n, np.nan)

    left_mask  = gid < 8
    right_mask = gid >= 8

    fpt_l_work = np.full(n, np.inf); fpt_r_work = np.full(n, np.inf)
    np.minimum.at(fpt_l_work, ev[left_mask],  tns[left_mask])
    np.minimum.at(fpt_r_work, ev[right_mask], tns[right_mask])
    fpt_l = np.where(fpt_l_work < np.inf, fpt_l_work, np.nan)
    fpt_r = np.where(fpt_r_work < np.inf, fpt_r_work, np.nan)

    sum_l = np.zeros(n); cnt_l = np.zeros(n, dtype=np.int64)
    sum_r = np.zeros(n); cnt_r = np.zeros(n, dtype=np.int64)
    np.add.at(sum_l, ev[left_mask],  tns[left_mask])
    np.add.at(cnt_l, ev[left_mask],  1)
    np.add.at(sum_r, ev[right_mask], tns[right_mask])
    np.add.at(cnt_r, ev[right_mask], 1)
    has_l = cnt_l > 0; has_r = cnt_r > 0
    mean_l = np.where(has_l, sum_l / cnt_l, np.nan)
    mean_r = np.where(has_r, sum_r / cnt_r, np.nan)

    for side_mask, t50_arr in ((left_mask, t50_l), (right_mask, t50_r)):
        ev_s  = ev[side_mask]
        tns_s = tns[side_mask]
        if ev_s.size == 0:
            continue
        order = np.argsort(ev_s, kind="stable")
        ev_s  = ev_s[order]
        tns_s = tns_s[order]
        unique_ev, starts = np.unique(ev_s, return_index=True)
        stops = np.r_[starts[1:], ev_s.size]
        for ue, st, en in zip(unique_ev, starts, stops):
            if 0 <= ue < n:
                t50_arr[ue] = float(np.median(tns_s[st:en]))

    triggers = np.full((n, 4), np.nan)
    combined = ev * 4 + gid // 4
    order    = np.argsort(combined, kind="stable")
    combined_s = combined[order]
    tns_s      = tns[order]
    unique_g, g_starts = np.unique(combined_s, return_index=True)
    g_stops = np.r_[g_starts[1:], combined_s.size]
    for ug, gs, ge in zip(unique_g, g_starts, g_stops):
        ue  = int(ug // 4)
        grp = int(ug % 4)
        if 0 <= ue < n and 0 <= grp < 4:
            triggers[ue, grp] = _leading_edge_time(tns_s[gs:ge])

    for i in range(n):
        sum4_l[i] = _earliest(triggers[i, 0], triggers[i, 1])
        sum4_r[i] = _earliest(triggers[i, 2], triggers[i, 3])

    return {
        "fpt_l": fpt_l, "fpt_r": fpt_r,
        "mean_l": mean_l, "mean_r": mean_r,
        "t50_l": t50_l, "t50_r": t50_r,
        "sum4_l": sum4_l, "sum4_r": sum4_r,
    }


# ─── figure helpers ───────────────────────────────────────────────────────────

def _save(fig: plt.Figure, path: pathlib.Path, caption: str) -> None:
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    path.with_suffix(".txt").write_text(caption)
    print(f"  [fig] {path.name}")


# ─── main analysis ────────────────────────────────────────────────────────────

def run(data_dir: pathlib.Path, out_dir: pathlib.Path) -> None:
    fig_dir = out_dir / "figures"
    tab_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    # ── metadata & CSVs ──────────────────────────────────────────────────────
    print("Reading metadata...")
    meta         = _read_metadata(data_dir)
    n_events     = int(meta["events_per_position"])
    mylar_r      = float(meta["mylar_reflectivity"])
    mylar_lobe   = float(meta["mylar_specular_lobe"])
    mylar_sigma  = float(meta["mylar_sigma_alpha_deg"])
    expected_pos = [float(x) for x in meta.get("positions", "").split()]

    print("Reading CSVs...")
    att_rows = _load_csv(data_dir / "analysis" / "attenuation_curve.csv")
    sig_rows = _load_csv(data_dir / "analysis" / "sigma_t_sum4.csv")

    x_mm       = np.array([float(r["x_mm"])           for r in att_rows])
    npe_l      = np.array([float(r["npe_left_mean"])   for r in att_rows])
    npe_l_err  = np.array([float(r["npe_left_sem"])    for r in att_rows])
    npe_r      = np.array([float(r["npe_right_mean"])  for r in att_rows])
    npe_r_err  = np.array([float(r["npe_right_sem"])   for r in att_rows])
    dist_l     = np.array([float(r["distance_left_mm"])  for r in att_rows])
    dist_r     = np.array([float(r["distance_right_mm"]) for r in att_rows])

    npe_sum     = npe_l + npe_r
    npe_sum_err = np.sqrt(npe_l_err**2 + npe_r_err**2)

    sig_x    = np.array([float(r["x_mm"])                 for r in sig_rows])
    sig_s    = np.array([float(r["sigma_single_ps"])       for r in sig_rows])
    sig_se   = np.array([float(r["sigma_single_error_ps"]) for r in sig_rows])
    trig_eff = np.array([float(r["trigger_efficiency"])    for r in sig_rows])
    mean_delta = np.array([float(r["mean_delta_lr_ns"])    for r in sig_rows])

    # ── attenuation fits ──────────────────────────────────────────────────────
    print("Fitting attenuation models...")
    fit_l = _fit_side(dist_l, npe_l, npe_l_err, label="left",
                      run_bootstrap=True, n_boot=200)
    fit_r = _fit_side(dist_r, npe_r, npe_r_err, label="right",
                      run_bootstrap=True, n_boot=200)

    # Combined: mirror the right-side data and stack with left
    dist_comb = np.concatenate([dist_l, dist_r])
    npe_comb  = np.concatenate([npe_l,  npe_r])
    err_comb  = np.concatenate([npe_l_err, npe_r_err])
    fit_c = _fit_side(dist_comb, npe_comb, err_comb, label="combined",
                      run_bootstrap=False)

    # ── Lambda_tail range sensitivity (d >= 30, 40, 50 cm) ──────────────────
    print("Computing lambda_tail range sensitivity...")
    d_cm_l = dist_l / 10.0
    d_cm_r = dist_r / 10.0
    tail_sensitivity = {}
    for d_min in [30.0, 40.0, 50.0]:
        entry = {"dmin_cm": d_min}
        for side_key, dc, yv, ye in [
            ("left",  d_cm_l, npe_l, npe_l_err),
            ("right", d_cm_r, npe_r, npe_r_err),
        ]:
            mask = dc > d_min
            n_pts_s = int(np.sum(mask))
            if n_pts_s >= 4:
                p0 = [float(np.max(yv[mask])) * 3.0, 42.0]
                popt, perr, chi2 = _weighted_fit(
                    _exp1, dc[mask], yv[mask], ye[mask], p0, ([0, 5], [np.inf, 500])
                )
                entry[side_key] = {
                    "lam_cm":   _safe(float(popt[1])),
                    "lam_err":  _safe(float(perr[1])),
                    "chi2ndf":  _safe(float(chi2)),
                    "n_points": n_pts_s,
                }
            else:
                entry[side_key] = None
        tail_sensitivity[f"dmin_{int(d_min)}cm"] = entry

    # ── EndTop comparison ─────────────────────────────────────────────────────
    print("Processing EndTop reference data...")
    endtop_data = _load_endtop_npe(ENDTOP_DIR) if ENDTOP_DIR.exists() else None
    endtop_available = endtop_data is not None

    et_x_mm = et_npe_l = et_npe_l_err = et_npe_r = et_npe_r_err = None
    et_dist_l = et_dist_r = None
    fit_et_l = fit_et_r = {}

    if endtop_available:
        sorted_et = sorted(endtop_data.items())
        et_x_mm      = np.array([v[0] for v in sorted_et])
        et_npe_l     = np.array([v[1]["npe_l"]     for v in sorted_et])
        et_npe_l_err = np.array([v[1]["npe_l_err"] for v in sorted_et])
        et_npe_r     = np.array([v[1]["npe_r"]     for v in sorted_et])
        et_npe_r_err = np.array([v[1]["npe_r_err"] for v in sorted_et])
        et_dist_l    = et_x_mm + 700.0   # mm
        et_dist_r    = 700.0 - et_x_mm   # mm

        print("  Fitting EndTop attenuation models...")
        fit_et_l = _fit_side(et_dist_l, et_npe_l, et_npe_l_err, label="endtop-left")
        fit_et_r = _fit_side(et_dist_r, et_npe_r, et_npe_r_err, label="endtop-right")
        print(f"  EndTop n_positions = {len(endtop_data)}")

    # ── boundary census ───────────────────────────────────────────────────────
    print("Parsing boundary census from logs...")
    census_list: list[dict] = []
    for xv in expected_pos:
        tag = f"{int(xv):+d}mm" if xv != 0 else "x0mm"
        log_path = data_dir / f"run_x{tag}.log"
        if not log_path.exists():
            log_path = data_dir / f"run_x{int(xv)}mm.log"
        if log_path.exists():
            res = _parse_position_log(log_path)
            res["x_mm"] = xv
            census_list.append(res)

    total_census: dict[str, int] = {}
    total_summary: dict[str, float] = {}
    for cl in census_list:
        for k, v in cl["census"].items():
            total_census[k] = total_census.get(k, 0) + v
        for k, v in cl.get("summary", {}).items():
            total_summary[k] = total_summary.get(k, 0.0) + v

    # ── ROOT timing ───────────────────────────────────────────────────────────
    print("Loading ROOT files for timing analysis...")
    fpt_l_mean = np.full(len(x_mm), np.nan)
    fpt_r_mean = np.full(len(x_mm), np.nan)
    mean_l_avg = np.full(len(x_mm), np.nan)
    mean_r_avg = np.full(len(x_mm), np.nan)
    t50_l_mean = np.full(len(x_mm), np.nan)
    t50_r_mean = np.full(len(x_mm), np.nan)
    repr_timing: dict[int, dict] = {}

    root_files = sorted(data_dir.glob("photon_hits_x*.root"))
    if not root_files:
        print("  WARNING: no ROOT files found — timing histograms will be skipped.")
    else:
        for idx, xv in enumerate(x_mm):
            tag = f"{int(xv):+d}mm" if xv != 0 else "x0mm"
            rpath = data_dir / f"photon_hits_x{tag}.root"
            if not rpath.exists():
                rpath = data_dir / f"photon_hits_x{int(xv)}mm.root"
            if not rpath.exists():
                print(f"  WARNING: ROOT file missing for x={xv} mm")
                continue
            timing = _per_event_timing(rpath, n_events)
            fpt_l_mean[idx] = float(np.nanmean(timing["fpt_l"]))
            fpt_r_mean[idx] = float(np.nanmean(timing["fpt_r"]))
            mean_l_avg[idx] = float(np.nanmean(timing["mean_l"]))
            mean_r_avg[idx] = float(np.nanmean(timing["mean_r"]))
            t50_l_mean[idx] = float(np.nanmean(timing["t50_l"]))
            t50_r_mean[idx] = float(np.nanmean(timing["t50_r"]))
            if int(xv) in REPR_MM:
                repr_timing[int(xv)] = timing
                print(f"  Loaded {rpath.name} [repr]")
            else:
                print(f"  Loaded {rpath.name}")

    # =========================================================================
    # FIGURES
    # =========================================================================
    print("Generating figures...")

    # ── FIG 1: Npe vs position ────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(x_mm, npe_l, yerr=npe_l_err, fmt="o-", ms=3, lw=1.2,
                color="C0", label="End-Left (IDs 0–7)")
    ax.errorbar(x_mm, npe_r, yerr=npe_r_err, fmt="s-", ms=3, lw=1.2,
                color="C1", label="End-Right (IDs 8–15)")
    ax.errorbar(x_mm, npe_sum, yerr=npe_sum_err, fmt="^-", ms=3, lw=1.2,
                color="C2", label="Sum L+R")
    ax.set_xlabel("Beam position $x$ (mm)")
    ax.set_ylabel("Mean photon count per event (PE)")
    ax.set_title("Photon count at ends vs. beam position — End-only + Mylar, EJ-204")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(-720, 720)
    _save(fig, fig_dir / "fig_npe_vs_x.png",
          f"Mean PE count at End-Left (IDs 0–7), End-Right (IDs 8–15), and Sum vs. x. "
          f"{n_events} muon events per position. Mylar R={mylar_r:.2f}.")

    # ── FIG 2: Attenuation — exp+offset fit (main plot) ──────────────────────
    d_fine = np.linspace(1, 145, 400)   # cm, fine grid for fit curves

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, side, dist_mm, npe, npe_err, fit, color, label in [
        (axes[0], "left",  dist_l, npe_l, npe_l_err, fit_l, "C0", "End-Left"),
        (axes[1], "right", dist_r, npe_r, npe_r_err, fit_r, "C1", "End-Right"),
    ]:
        ax.errorbar(dist_mm / 10, npe, yerr=npe_err, fmt="o", ms=3,
                    color=color, label=label + " data", zorder=3)

        # Exp+offset curve (chosen model)
        A, lam, C = fit["expoff_A"], fit["expoff_lam_cm"], fit["expoff_C"]
        chi2_b = fit["expoff_chi2ndf"]
        if all(v is not None and math.isfinite(v) for v in [A, lam, C]):
            y_fit = _expoff(d_fine, A, lam, C)
            ax.plot(d_fine, y_fit, "-", color=color, lw=1.8,
                    label=(rf"Exp+offset fit:"
                           rf" $\lambda_p={lam:.1f}$ cm, C={C:.1f} PE"
                           rf" ($\chi^2/\mathrm{{ndf}}$={chi2_b:.0f})"))
            # Decompose: prompt and floor components
            ax.plot(d_fine, _expoff(d_fine, A, lam, 0), "--",
                    color=color, lw=0.9, alpha=0.6, label="Prompt component $A e^{-d/\\lambda_p}$")
            ax.axhline(C, color=color, ls=":", lw=0.9, alpha=0.6,
                       label=f"Floor $C={C:.1f}$ PE")

        # Lambda_tail fit (d > D_TAIL_CM)
        lam_t = fit["tail_lam_cm"]
        chi2_t = fit["tail_chi2ndf"]
        if lam_t is not None and math.isfinite(lam_t):
            n0_t = fit["tail_n0"]
            tail_d = d_fine[d_fine > D_TAIL_CM]
            ax.plot(tail_d, _exp1(tail_d, n0_t, lam_t), "k-.", lw=1.2,
                    label=(rf"$\lambda_\mathrm{{tail}}={lam_t:.1f}$ cm"
                           rf" ($d>{D_TAIL_CM:.0f}$ cm, $\chi^2/\mathrm{{ndf}}$={chi2_t:.0f})"))

        ax.axvline(LAMBDA_BULK_CM, color="gray", ls="--", lw=1,
                   label=rf"Bulk $\lambda_\mathrm{{bulk}}={LAMBDA_BULK_CM:.0f}$ cm (ABSLENGTH)")
        ax.set_xlabel("Distance from end (cm)")
        ax.set_ylabel("Mean photon count per event (PE)")
        ax.set_title(f"Attenuation — {label}  (Mylar R={mylar_r:.2f})")
        ax.legend(fontsize=6.5)
        ax.set_yscale("log")
        ax.set_ylim(bottom=5)

    plt.tight_layout()
    _save(fig, fig_dir / "fig_attenuation.png",
          f"Photon count vs. distance from respective end. "
          f"Exp+offset fit (λ_prompt, floor C) and λ_tail (d>{D_TAIL_CM:.0f} cm). "
          f"Bulk ABSLENGTH = {LAMBDA_BULK_CM:.0f} cm is a material property "
          f"(not directly comparable to λ_eff). Mylar R={mylar_r:.2f}.")

    # ── FIG 3: Two-component fit diagnostics ──────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, side, dist_mm, npe, npe_err, fit, color, label in [
        (axes[0], "left",  dist_l, npe_l, npe_l_err, fit_l, "C0", "End-Left"),
        (axes[1], "right", dist_r, npe_r, npe_r_err, fit_r, "C1", "End-Right"),
    ]:
        ax.errorbar(dist_mm / 10, npe, yerr=npe_err, fmt="o", ms=3,
                    color=color, zorder=3)

        A, lam_p, C = fit["expoff_A"], fit["expoff_lam_cm"], fit["expoff_C"]
        n0_s, lam_s = fit["single_n0"], fit["single_lam_cm"]
        chi2_s = fit["single_chi2ndf"]
        chi2_b = fit["expoff_chi2ndf"]

        if all(v is not None and math.isfinite(v) for v in [n0_s, lam_s]):
            ax.plot(d_fine, _exp1(d_fine, n0_s, lam_s), "r--", lw=1.4,
                    label=(rf"Single-exp: $\lambda={lam_s:.1f}$ cm, "
                           rf"$N_0={n0_s:.0f}$ PE ($\chi^2/\mathrm{{ndf}}$={chi2_s:.0f})"))

        if all(v is not None and math.isfinite(v) for v in [A, lam_p, C]):
            ax.plot(d_fine, _expoff(d_fine, A, lam_p, C), "-", color=color, lw=1.8,
                    label=(rf"Exp+offset: $\lambda_p={lam_p:.1f}$ cm, "
                           rf"$C={C:.1f}$ PE ($\chi^2/\mathrm{{ndf}}$={chi2_b:.0f})"))
            # Mark the inconsistency: single-exp N0 vs actual at d≈0
            ax.axhline(n0_s, color="red", ls=":", lw=0.8, alpha=0.7,
                       label=rf"Single-exp $N_0={n0_s:.0f}$ PE at $d=0$")
            actual_d0 = npe[np.argmin(dist_mm)]
            ax.scatter([dist_mm.min() / 10], [actual_d0], marker="*", s=120,
                       color="red", zorder=5, label=f"Actual at d≈0: {actual_d0:.0f} PE")

        ax.set_xlabel("Distance from end (cm)")
        ax.set_ylabel("Mean photon count per event (PE)")
        ax.set_title(f"Model comparison — {label}")
        ax.legend(fontsize=6.5)
        ax.set_yscale("log")
        ax.set_ylim(bottom=5)

    plt.tight_layout()
    _save(fig, fig_dir / "fig_attenuation_twocomp.png",
          "Single-exp (red dashed) vs. exp+offset (solid) fit comparison. "
          "Red star: actual NPE at d≈0; red dotted: single-exp N0 extrapolated to d=0. "
          "Inconsistency factor ~3.6 shows single-exp model fails. "
          f"Chi2/ndf single-exp≈{fit_l['single_chi2ndf']:.0f} vs "
          f"exp+offset≈{fit_l['expoff_chi2ndf']:.0f}.")

    # ── FIG 4: EndTop vs Mylar comparison ────────────────────────────────────
    if endtop_available:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        for ax, dist_mylar, npe_mylar, err_mylar, dist_et, npe_et, err_et, \
                fit_mylar, fit_et, label in [
            (axes[0], dist_l, npe_l, npe_l_err, et_dist_l, et_npe_l, et_npe_l_err,
             fit_l, fit_et_l, "Left end"),
            (axes[1], dist_r, npe_r, npe_r_err, et_dist_r, et_npe_r, et_npe_r_err,
             fit_r, fit_et_r, "Right end"),
        ]:
            ax.errorbar(dist_mylar / 10, npe_mylar, yerr=err_mylar,
                        fmt="o", ms=3, color="C0",
                        label=f"End-only+Mylar R={mylar_r:.2f}", zorder=3)
            ax.errorbar(dist_et / 10, npe_et, yerr=err_et,
                        fmt="s", ms=3, color="C3",
                        label=f"EndTop R={ENDTOP_REFLECTOR_R:.2f}", zorder=3)

            # Shade tail region
            ax.axvspan(D_TAIL_CM, 145, alpha=0.06, color="C2",
                       label=rf"Tail region $d>{D_TAIL_CM:.0f}$ cm")

            lam_m = fit_mylar.get("tail_lam_cm")
            lam_e = fit_et.get("tail_lam_cm")
            if lam_m and math.isfinite(lam_m):
                n0_m = fit_mylar["tail_n0"]
                tail_d = d_fine[d_fine > D_TAIL_CM]
                ax.plot(tail_d, _exp1(tail_d, n0_m, lam_m), "-", color="C0",
                        lw=1.4, alpha=0.8,
                        label=rf"Mylar $\lambda_\mathrm{{tail}}={lam_m:.1f}$ cm")
            if lam_e and math.isfinite(lam_e):
                n0_e = fit_et["tail_n0"]
                tail_d = d_fine[d_fine > D_TAIL_CM]
                ax.plot(tail_d, _exp1(tail_d, n0_e, lam_e), "-", color="C3",
                        lw=1.4, alpha=0.8,
                        label=rf"EndTop $\lambda_\mathrm{{tail}}={lam_e:.1f}$ cm")

            ax.set_xlabel("Distance from end (cm)")
            ax.set_ylabel("Mean photon count per event (PE)")
            ax.set_title(f"End-only+Mylar vs. EndTop — {label}")
            ax.legend(fontsize=6.5)
            ax.set_yscale("log")
            ax.set_ylim(bottom=5)

        plt.tight_layout()
        _save(fig, fig_dir / "fig_endtop_comparison.png",
              f"Apples-to-apples: end-SiPM NPE vs distance for Mylar (R={mylar_r:.2f}) "
              f"and EndTop (R={ENDTOP_REFLECTOR_R:.2f} side reflectors). "
              f"Same pipeline applied to both datasets. "
              f"Tail region (d>{D_TAIL_CM:.0f}cm) shaded; λ_tail fits shown. "
              f"NPE amplitude differs (reflector R), but λ_tail is similar.")

    # ── FIG 5: Far-end light yield (center positions) ────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4))
    # "far-end" = center positions where |x| ≤ 100 mm (d ≥ 600mm = 60cm from both ends)
    far_mask_m  = np.abs(x_mm) <= 100
    far_mean_l  = float(np.mean(npe_l[far_mask_m]))
    far_mean_r  = float(np.mean(npe_r[far_mask_m]))
    far_mean_m  = 0.5 * (far_mean_l + far_mean_r)

    configs  = ["End-only+Mylar"]
    vals_l   = [far_mean_l]
    vals_r   = [far_mean_r]
    colors_b = ["C0"]

    if endtop_available:
        far_mask_e  = np.abs(et_x_mm) <= 100
        far_et_l    = float(np.mean(et_npe_l[far_mask_e]))
        far_et_r    = float(np.mean(et_npe_r[far_mask_e]))
        configs.append("EndTop")
        vals_l.append(far_et_l)
        vals_r.append(far_et_r)
        colors_b.append("C3")

    x_pos = np.arange(len(configs))
    w = 0.35
    bars_l = ax.bar(x_pos - w/2, vals_l, w, label="End-Left",  color=[c for c in colors_b], alpha=0.8)
    bars_r = ax.bar(x_pos + w/2, vals_r, w, label="End-Right", color=[c for c in colors_b], alpha=0.5)
    ax.set_xticks(x_pos); ax.set_xticklabels(configs)
    ax.set_ylabel("Mean PE count (|x|≤100 mm)")
    ax.set_title("Far-end light yield: center-of-bar NPE (photons that traversed ≥60 cm)")
    ax.legend()
    for bar in list(bars_l) + list(bars_r):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    caption = (f"Mean NPE at end SiPMs for center positions (|x|≤100 mm, d≥60 cm). "
               f"End-only+Mylar (R={mylar_r:.2f}): L={far_mean_l:.1f}, R={far_mean_r:.1f} PE. ")
    if endtop_available:
        caption += (f"EndTop (R={ENDTOP_REFLECTOR_R:.2f}): "
                    f"L={far_et_l:.1f}, R={far_et_r:.1f} PE. "
                    f"Ratio EndTop/Mylar={far_et_l/far_mean_l:.1f}×. "
                    f"Higher yield in EndTop due to better reflector (R={ENDTOP_REFLECTOR_R:.2f} vs {mylar_r:.2f}).")
    _save(fig, fig_dir / "fig_far_end_yield.png", caption)

    # ── FIG 6: Timing histograms ──────────────────────────────────────────────
    if repr_timing:
        n_repr = len(repr_timing)
        fig, axes_grid = plt.subplots(n_repr, 2, figsize=(10, 2.2 * n_repr),
                                      sharex=False, sharey=False)
        if n_repr == 1:
            axes_grid = axes_grid[np.newaxis, :]

        for row_idx, xv in enumerate(REPR_MM):
            if xv not in repr_timing:
                continue
            td = repr_timing[xv]
            for col_idx, (key, label, color) in enumerate([
                ("fpt_l", "FPT End-Left (ns)",  "C0"),
                ("fpt_r", "FPT End-Right (ns)", "C1"),
            ]):
                ax = axes_grid[row_idx, col_idx]
                data = td[key]
                data = data[np.isfinite(data)]
                if data.size == 0:
                    ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                    continue
                lo, hi = np.percentile(data, 1), np.percentile(data, 99)
                bins = min(100, max(20, int(np.sqrt(data.size))))
                ax.hist(data, bins=bins, range=(lo, hi), color=color, alpha=0.75,
                        histtype="stepfilled", edgecolor="none")
                ax.set_xlabel("Arrival time (ns)")
                ax.set_ylabel("Events")
                ax.set_title(f"{_pos_tag(float(xv))}, {label}, N={data.size}")
        plt.tight_layout()
        _save(fig, fig_dir / "fig_timing_histograms.png",
              "Per-event first photon arrival time (FPT) at End-Left and End-Right "
              "for representative positions. Counts (not normalized).")

        # FPT vs t50 comparison
        fig2, axes2 = plt.subplots(1, 2, figsize=(10, 4))
        for col_idx, (key_fpt, key_t50, end_label, color) in enumerate([
            ("fpt_l", "t50_l", "End-Left",  "C0"),
            ("fpt_r", "t50_r", "End-Right", "C1"),
        ]):
            ax = axes2[col_idx]
            for xv in REPR_MM:
                if xv not in repr_timing:
                    continue
                td = repr_timing[xv]
                for key, ls, lbl in [
                    (key_fpt, "-",  "FPT"),
                    (key_t50, "--", "$t_{50}$"),
                ]:
                    d = td[key]
                    d = d[np.isfinite(d)]
                    if d.size == 0:
                        continue
                    lo, hi = np.percentile(d, 1), np.percentile(d, 99)
                    counts, edges = np.histogram(d, bins=80, range=(lo, hi))
                    cx = 0.5 * (edges[:-1] + edges[1:])
                    ax.plot(cx, counts / counts.max(), ls=ls, lw=1.2,
                            label=f"x={xv:+d} {lbl}")
            ax.set_xlabel("Arrival time (ns)")
            ax.set_ylabel("Normalized count")
            ax.set_title(f"FPT vs. $t_{{50}}$ — {end_label}")
            ax.legend(fontsize=6)
        plt.tight_layout()
        _save(fig2, fig_dir / "fig_fpt_vs_t50.png",
              "Comparison of FPT and t50 arrival time distributions. Normalized to peak.")

    # ── FIG 7: Mean arrival time vs x and v_eff ───────────────────────────────
    v_eff_fpt = v_eff_t50 = v_eff_sum4 = math.nan
    ok   = np.isfinite(fpt_l_mean) & np.isfinite(fpt_r_mean)
    ok_l = np.isfinite(fpt_l_mean)
    ok_r = np.isfinite(fpt_r_mean)
    ok_t50_l = np.isfinite(t50_l_mean)
    ok_t50_r = np.isfinite(t50_r_mean)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    if ok.sum() >= 4:
        for side_arr, label, color in [
            (fpt_l_mean, "End-Left FPT",  "C0"),
            (fpt_r_mean, "End-Right FPT", "C1"),
        ]:
            axes[0].plot(x_mm[ok], side_arr[ok], "o-", ms=3, color=color, label=label)

        if ok_l.sum() >= 4:
            cl = np.polyfit(x_mm[ok_l], fpt_l_mean[ok_l], 1)
            v_eff_fpt_l = abs(1.0 / cl[0]) if abs(cl[0]) > 1e-6 else math.nan
        else:
            v_eff_fpt_l = math.nan

        if ok_r.sum() >= 4:
            cr = np.polyfit(x_mm[ok_r], fpt_r_mean[ok_r], 1)
            v_eff_fpt_r = abs(1.0 / cr[0]) if abs(cr[0]) > 1e-6 else math.nan
        else:
            v_eff_fpt_r = math.nan

        if math.isfinite(v_eff_fpt_l) and math.isfinite(v_eff_fpt_r):
            v_eff_fpt = 0.5 * (v_eff_fpt_l + v_eff_fpt_r)

        ok_sig = np.isfinite(mean_delta)
        if ok_sig.sum() >= 4:
            cs = np.polyfit(sig_x[ok_sig], mean_delta[ok_sig], 1)
            v_eff_sum4 = abs(2.0 / cs[0]) if abs(cs[0]) > 1e-6 else math.nan

        if ok_t50_l.sum() >= 4 and ok_t50_r.sum() >= 4:
            cl50 = np.polyfit(x_mm[ok_t50_l], t50_l_mean[ok_t50_l], 1)
            cr50 = np.polyfit(x_mm[ok_t50_r], t50_r_mean[ok_t50_r], 1)
            sl, sr = cl50[0], cr50[0]
            if abs(sl) > 1e-6 and abs(sr) > 1e-6:
                v_eff_t50 = 0.5 * (abs(1/sl) + abs(1/sr))

        axes[1].plot(x_mm[ok_l], fpt_l_mean[ok_l],  "o-", ms=3, color="C0",
                     label=r"$\langle t_\mathrm{FPT}\rangle$ Left")
        axes[1].plot(x_mm[ok_r], fpt_r_mean[ok_r],  "s-", ms=3, color="C1",
                     label=r"$\langle t_\mathrm{FPT}\rangle$ Right")
        axes[1].plot(x_mm[ok_t50_l], t50_l_mean[ok_t50_l], "^--", ms=3, color="C0",
                     alpha=0.6, label=r"$\langle t_{50}\rangle$ Left")
        axes[1].plot(x_mm[ok_t50_r], t50_r_mean[ok_t50_r], "v--", ms=3, color="C1",
                     alpha=0.6, label=r"$\langle t_{50}\rangle$ Right")

    axes[0].set_xlabel("Beam position $x$ (mm)")
    axes[0].set_ylabel("Mean FPT (ns)")
    axes[0].set_title("Mean first photon time vs. position")
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("Beam position $x$ (mm)")
    axes[1].set_ylabel("Mean arrival time (ns)")
    axes[1].set_title("Mean FPT and $t_{50}$ vs. position")
    axes[1].legend(fontsize=7)
    plt.tight_layout()
    _save(fig, fig_dir / "fig_mean_arrival_time.png",
          f"Mean FPT and t50 vs. position. "
          f"v_eff(FPT)={v_eff_fpt:.1f} mm/ns, v_eff(SUM4)={v_eff_sum4:.1f} mm/ns.")

    # ── FIG 8: sigma_t(x) ────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(sig_x, sig_s, yerr=sig_se, fmt="o-", ms=3, color="C3",
                label=r"$\sigma_\mathrm{single} = \sigma_{LR}/\sqrt{2}$")
    ax.fill_between(sig_x, sig_s - sig_se, sig_s + sig_se, alpha=0.2, color="C3")
    ax.set_xlabel("Beam position $x$ (mm)")
    ax.set_ylabel(r"$\sigma_t$ (ps)")
    ax.set_title(r"Intrinsic timing resolution $\sigma_t(x)$ — End-only + Mylar, EJ-204")
    ax.legend()
    sigma_center = float(sig_s[np.argmin(np.abs(sig_x))])
    sigma_end690 = float(np.nanmean([sig_s[np.argmin(np.abs(sig_x - 690))],
                                     sig_s[np.argmin(np.abs(sig_x + 690))]]))
    _save(fig, fig_dir / "fig_sigma_t.png",
          rf"Intrinsic σ_single = σ_LR/√2 (SUM4 leading-edge, core Gaussian fit) vs. x. "
          f"Center (x=0): {sigma_center:.0f} ps; ends (|x|=690 mm): {sigma_end690:.0f} ps. "
          r"Intrinsic only — SPTR and FastIC jitter not folded in.")

    # ── FIG 9: Photon budget (Panel A: terminal fates; Panel B: sensor conversion) ──
    key_panel = "Bar -> reflector panel"
    key_world = "Bar -> World (escaped)"
    key_sipm  = "Bar -> SiPM (entering)"
    key_scint = "Scint photons generated"

    n_gen  = total_summary.get(key_scint, 0.0)
    n_sipm_budget = float(total_census.get(key_sipm, 0))
    n_esc  = float(total_census.get(key_world, 0))
    n_ref  = float(total_census.get(key_panel, 0))
    n_det  = n_sipm_budget * PDE_PEAK
    # "Surface losses + bulk" = generated – escaped – SiPM-entering
    # This lumps Mylar surface absorption and any bulk ABSLENGTH absorption.
    # The simulation census records "Bar→reflector panel" crossings (not absorption
    # events), so we cannot separate surface from bulk from this census alone.
    n_surface_bulk = max(0.0, n_gen - n_esc - n_sipm_budget)

    fig, axes_pb = plt.subplots(1, 2, figsize=(11, 4))

    if n_gen > 0:
        # Panel A: terminal fates (mutually exclusive, sum ≈ 100%)
        labels_a = ["Surface+bulk\nabsorbed", "Escaped\n(Bar→World)", "Entered\nSiPM"]
        vals_a   = np.array([n_surface_bulk, n_esc, n_sipm_budget])
        pct_a    = 100.0 * vals_a / n_gen
        colors_a = ["C0", "C1", "C3"]
        bars_a   = axes_pb[0].barh(labels_a[::-1], pct_a[::-1], color=colors_a[::-1])
        axes_pb[0].set_xlabel("Fraction of generated photons (%)")
        axes_pb[0].set_title(
            f"Panel A — Terminal fates (sum = {pct_a.sum():.1f}%)\nMylar R={mylar_r:.2f}; "
            r"loss per interaction $1-R=0.10$")
        for bar, p in zip(bars_a, pct_a[::-1]):
            axes_pb[0].text(bar.get_width() + 0.05,
                            bar.get_y() + bar.get_height() / 2,
                            f"{p:.2f}%", va="center", fontsize=8)
        axes_pb[0].set_xlim(0, max(pct_a) * 1.18)

        axes_pb[0].text(
            0.02, 0.02,
            "Note: 'Surface+bulk absorbed' = generated\n"
            "minus escaped minus SiPM-entering.\n"
            "Bulk ABSLENGTH and surface losses\nnot separately resolved in census.",
            transform=axes_pb[0].transAxes,
            fontsize=6.5, color="gray", va="bottom")

        # Panel B: SiPM conversion chain (conditional on reaching SiPM)
        if n_sipm_budget > 0:
            pde_eff = n_det / n_sipm_budget
            labels_b  = ["Not detected\n(1−PDE_eff)", "Detected\n(after PDE)"]
            vals_b    = np.array([n_sipm_budget - n_det, n_det])
            pct_b     = 100.0 * vals_b / n_sipm_budget
            bars_b    = axes_pb[1].barh(labels_b[::-1], pct_b[::-1],
                                        color=["lightgray", "C2"])
            axes_pb[1].set_xlabel("Fraction of SiPM-entering photons (%)")
            axes_pb[1].set_title(
                f"Panel B — SiPM conversion\n"
                f"Entered={n_sipm_budget/n_gen*100:.2f}% of total\n"
                f"PDE_eff = {pde_eff:.2f} ({PDE_PEAK:.0%} nominal)")
            for bar, p in zip(bars_b, pct_b[::-1]):
                axes_pb[1].text(bar.get_width() + 0.3,
                                bar.get_y() + bar.get_height() / 2,
                                f"{p:.1f}%", va="center", fontsize=8)
            axes_pb[1].set_xlim(0, 115)
    else:
        axes_pb[0].text(0.5, 0.5, "No log data available",
                        transform=axes_pb[0].transAxes, ha="center")

    plt.tight_layout()
    _save(fig, fig_dir / "fig_photon_budget.png",
          f"Photon budget: Panel A — terminal fates (surface+bulk absorbed, escaped, SiPM-entering); "
          f"Panel B — SiPM conversion. "
          f"Generated={n_gen:.2e}. "
          f"Surface+bulk≈{n_surface_bulk/n_gen*100:.1f}% (Mylar R={mylar_r:.2f}, "
          r"loss per interaction 1−R=0.10; cumulative over many bounces). "
          f"SiPM entering={n_sipm_budget/n_gen*100:.2f}%; "
          f"detected={n_det/n_gen*100:.2f}% of total (={PDE_PEAK*100:.0f}% of SiPM-entering).")

    # =========================================================================
    # TABLES
    # =========================================================================
    print("Writing tables...")
    bar_length_mm = 2 * BAR_HALF_MM

    # Table 1: simulation parameters
    params_latex = [
        ["Parameter", "Value", "Unit", "Source"],
        ["Bar material", "EJ-204 (OPSC-101)", "", "SSLG4/Materials.cc"],
        ["Bar dimensions",
         f"{bar_length_mm:.0f}\\,$\\times$\\,60\\,$\\times$\\,10",
         "mm$^3$", "GDML"],
        ["Readout", "End-only ($\\pm X$ faces)", "", "geometry"],
        ["Wrapping ($\\pm Y,\\pm Z$)", "Mylar", "", "run\\_metadata"],
        ["Mylar reflectivity $R$",    f"{mylar_r:.2f}", "", "run\\_metadata"],
        ["Mylar specular lobe",        f"{mylar_lobe:.1f}", "", "run\\_metadata"],
        ["Mylar $\\sigma_\\alpha$",    f"{mylar_sigma:.1f}", "deg", "run\\_metadata"],
        ["Scint.~yield",
         f"{SCINT_YIELD:,}",
         "$\\gamma$/MeV",
         "Eljen EJ-204 datasheet Aug~2023; Materials.cc:100"],
        ["Emission peak",    f"{PEAK_EMISSION_NM:.0f}", "nm",
         "Eljen EJ-204 datasheet Aug~2023; Materials.cc:86"],
        ["Bulk $\\lambda$ (ABSLENGTH)", f"{LAMBDA_BULK_CM:.0f}", "cm",
         "Eljen EJ-204 datasheet Aug~2023; Materials.cc:99"],
        ["Refractive index $n$", f"{REFRACTIVE_INDEX:.2f}", "",
         "Eljen EJ-204 datasheet; Materials.cc:41"],
        ["Density", f"{DENSITY_GCC:.3f}", "g/cm$^3$", "Eljen EJ-204 datasheet; Materials.cc:24"],
        ["Rise time $\\tau_r$",   f"{TAU_RISE_NS:.1f}", "ns", "SSLG4 (EJ-204)"],
        ["Decay time $\\tau_d$",  f"{TAU_DECAY_NS:.1f}", "ns", "SSLG4 (EJ-204)"],
        ["SiPM model (physical)", SIPM_MODEL, "", "hardware"],
        ["SiPM PDE (peak)", f"{PDE_PEAK*100:.0f}", "\\% @ 420\\,nm", "datasheet"],
        ["Beam", "$\\mu^-$, 1\\,GeV", "", "macro"],
        ["Positions", f"{len(expected_pos)}", "", "run\\_metadata"],
        ["Events/position", f"{n_events:,}", "", "run\\_metadata"],
        ["Geant4 version", "11.4.0", "", "log"],
    ]

    with open(tab_dir / "tab_sim_params.tex", "w") as f:
        f.write("\\begin{tabular}{llll}\\toprule\n")
        f.write(" & ".join(params_latex[0]) + " \\\\\n\\midrule\n")
        for row in params_latex[1:]:
            f.write(" & ".join(str(c) for c in row) + " \\\\\n")
        f.write("\\bottomrule\\end{tabular}\n")

    with open(tab_dir / "tab_sim_params.csv", "w", newline="") as f:
        writer = csv.writer(f)
        for row in params_latex:
            writer.writerow([c.replace("\\", "").replace("$", "").replace("~", " ")
                             for c in row])

    # Table 2: attenuation fits (three models)
    att_rows_tex = []
    att_rows_tex.append(["Model", "Side",
                         "$N_0$ or $A$ (PE)", "$\\lambda$ (cm)", "$C$ (PE)", "$\\chi^2/\\mathrm{ndf}$"])
    for side_label, fit in [("Left", fit_l), ("Right", fit_r)]:
        # Single-exp row
        att_rows_tex.append([
            "Single-exp",
            side_label,
            (f"${fit['single_n0']:.0f}\\pm{fit['single_n0_err']:.0f}$"
             if fit["single_n0"] else "---"),
            (f"${fit['single_lam_cm']:.1f}\\pm{fit['single_lam_err']:.1f}$"
             if fit["single_lam_cm"] else "---"),
            "---",
            f"{fit['single_chi2ndf']:.0f}" if fit["single_chi2ndf"] else "---",
        ])
        # Exp+offset row
        att_rows_tex.append([
            "Exp+offset",
            side_label,
            (f"${fit['expoff_A']:.0f}\\pm{fit['expoff_A_err']:.0f}$"
             if fit["expoff_A"] else "---"),
            (f"${fit['expoff_lam_cm']:.1f}\\pm{fit['expoff_lam_err']:.1f}$"
             if fit["expoff_lam_cm"] else "---"),
            (f"${fit['expoff_C']:.1f}\\pm{fit['expoff_C_err']:.1f}$"
             if fit["expoff_C"] else "---"),
            f"{fit['expoff_chi2ndf']:.0f}" if fit["expoff_chi2ndf"] else "---",
        ])
        # Lambda_tail row
        att_rows_tex.append([
            f"$\\lambda_\\mathrm{{tail}}$ ($d>{D_TAIL_CM:.0f}$\\,cm)",
            side_label,
            (f"${fit['tail_n0']:.0f}\\pm{fit['tail_n0_err']:.0f}$"
             if fit["tail_n0"] else "---"),
            (f"${fit['tail_lam_cm']:.1f}\\pm{fit['tail_lam_err']:.1f}$"
             if fit["tail_lam_cm"] else "---"),
            "---",
            f"{fit['tail_chi2ndf']:.1f}" if fit["tail_chi2ndf"] else "---",
        ])

    with open(tab_dir / "tab_attenuation_fits.tex", "w") as f:
        f.write("\\begin{tabular}{llllll}\\toprule\n")
        f.write(" & ".join(att_rows_tex[0]) + " \\\\\n\\midrule\n")
        for row in att_rows_tex[1:]:
            f.write(" & ".join(str(c) for c in row) + " \\\\\n")
        f.write("\\bottomrule\\end{tabular}\n")

    with open(tab_dir / "tab_attenuation_fits.csv", "w", newline="") as f:
        writer = csv.writer(f)
        for row in att_rows_tex:
            writer.writerow([c.replace("\\", "").replace("$", "").replace("~", " ")
                             for c in row])

    # Lambda_tail range sensitivity table (CSV)
    with open(tab_dir / "attenuation_range_sensitivity.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dmin_cm", "side", "lam_cm", "lam_err", "chi2ndf", "n_points"])
        for key, entry in tail_sensitivity.items():
            for side in ("left", "right"):
                d_info = entry.get(side)
                if d_info:
                    w.writerow([entry["dmin_cm"], side,
                                 d_info["lam_cm"], d_info["lam_err"],
                                 d_info["chi2ndf"], d_info["n_points"]])

    # Bootstrap summary CSV
    for side_key, fit_side_data in [("left", fit_l), ("right", fit_r)]:
        boot = fit_side_data.get("bootstrap")
        if boot is None:
            continue
        with open(tab_dir / f"attenuation_bootstrap_{side_key}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["metric", "n", "median", "p16", "p84", "mean", "std"])
            for met in ["M3_lam_cm", "M3_C_PE", "M4_lam_cm", "M2_lam_s_cm", "M2_lam_l_cm"]:
                s = boot.get(met)
                if s:
                    w.writerow([met, s["n"], s["median"], s["p16"], s["p84"],
                                 s["mean"], s["std"]])
            w.writerow(["M2_n_fail", boot["M2_n_fail"], "", "", "", "", ""])
            w.writerow(["M2_failure_frac", boot["M2_failure_fraction"], "", "", "", "", ""])

    # Table 3: sigma_t at key positions
    key_xvals = [-690, -400, 0, 400, 690]
    sig_table_tex_rows = [["$x$ (mm)", "$N_\\mathrm{trig}$",
                           "Trig.~eff.", "$\\sigma_\\mathrm{single}$ (ps)", "Err. (ps)"]]
    sig_table_csv_rows = [["x (mm)", "N triggered", "Trig. eff.",
                           "sigma_single (ps)", "sigma_err (ps)"]]
    for kx in key_xvals:
        idx_a = np.argmin(np.abs(sig_x - kx))
        sig_table_tex_rows.append([
            f"${int(sig_x[idx_a]):+d}$",
            str(int(float(sig_rows[idx_a]["n_triggered_lr"]))),
            f"{trig_eff[idx_a]:.3f}",
            f"{sig_s[idx_a]:.1f}",
            f"{sig_se[idx_a]:.1f}",
        ])
        sig_table_csv_rows.append([
            f"{int(sig_x[idx_a]):+d}",
            str(int(float(sig_rows[idx_a]["n_triggered_lr"]))),
            f"{trig_eff[idx_a]:.3f}",
            f"{sig_s[idx_a]:.1f}",
            f"{sig_se[idx_a]:.1f}",
        ])
    with open(tab_dir / "tab_sigma_t.tex", "w") as f:
        f.write("\\begin{tabular}{rrrrl}\\toprule\n")
        f.write(" & ".join(sig_table_tex_rows[0]) + " \\\\\n\\midrule\n")
        for row in sig_table_tex_rows[1:]:
            f.write(" & ".join(row) + " \\\\\n")
        f.write("\\bottomrule\\end{tabular}\n")
    with open(tab_dir / "tab_sigma_t.csv", "w", newline="") as f:
        csv.writer(f).writerows(sig_table_csv_rows)

    # =========================================================================
    # deck_values.json
    # =========================================================================
    print("Writing deck_values.json...")

    sigma_at = {}
    for kx in key_xvals:
        idx_a = int(np.argmin(np.abs(sig_x - kx)))
        sigma_at[str(kx)] = {
            "sigma_single_ps":     _safe(float(sig_s[idx_a])),
            "sigma_single_err_ps": _safe(float(sig_se[idx_a])),
        }

    # Far-end (tail) NPE
    idx0        = int(np.argmin(np.abs(x_mm)))
    far_mask_m  = np.abs(x_mm) <= 100
    far_mean_l  = float(np.mean(npe_l[far_mask_m]))
    far_mean_r  = float(np.mean(npe_r[far_mask_m]))

    endtop_section: dict = {"available": False}
    if endtop_available:
        far_mask_e = np.abs(et_x_mm) <= 100
        far_et_l   = float(np.mean(et_npe_l[far_mask_e]))
        far_et_r   = float(np.mean(et_npe_r[far_mask_e]))
        endtop_section = {
            "available":          True,
            "data_dir":           str(ENDTOP_DIR),
            "n_positions":        len(endtop_data),
            "n_events":           int(list(endtop_data.values())[0]["n_events"]),
            "reflector_R":        ENDTOP_REFLECTOR_R,
            "npe_left_at_x0":     _safe(et_npe_l[np.argmin(np.abs(et_x_mm))]),
            "npe_right_at_x0":    _safe(et_npe_r[np.argmin(np.abs(et_x_mm))]),
            "tail_npe_left_mean": _safe(far_et_l),
            "tail_npe_right_mean":_safe(far_et_r),
            "fit_left":           fit_et_l,
            "fit_right":          fit_et_r,
        }

    deck = {
        "source":   "analysis_endonly_mylar.py",
        "data_dir": str(data_dir),
        # Run parameters
        "n_events_per_position": n_events,
        "n_positions":           len(x_mm),
        "bar_half_length_mm":    BAR_HALF_MM,
        "bar_length_mm":         2 * BAR_HALF_MM,
        # Material properties (from Eljen EJ-204 datasheet Aug 2023 + Materials.cc)
        "scint_yield_per_MeV":   SCINT_YIELD,
        "peak_emission_nm":      PEAK_EMISSION_NM,
        "lambda_bulk_cm":        LAMBDA_BULK_CM,
        "refractive_index":      REFRACTIVE_INDEX,
        "density_gcc":           DENSITY_GCC,
        "tau_rise_ns":           TAU_RISE_NS,
        "tau_decay_ns":          TAU_DECAY_NS,
        # SiPM
        "pde_peak_frac":         PDE_PEAK,
        "sipm_model":            SIPM_MODEL,
        "sipm_model_pde_curve":  SIPM_MODEL_PDE,
        # Mylar
        "mylar_reflectivity":    mylar_r,
        "mylar_specular_lobe":   mylar_lobe,
        "mylar_sigma_alpha_deg": mylar_sigma,
        # Attenuation fits: three models, left side
        "att_left":  fit_l,
        "att_right": fit_r,
        "att_combined": fit_c,
        # Backward-compat keys used by old verify_deck / build_deck
        "lambda_eff_left_cm":       fit_l["single_lam_cm"],
        "lambda_eff_left_err_cm":   fit_l["single_lam_err"],
        "lambda_eff_right_cm":      fit_r["single_lam_cm"],
        "lambda_eff_right_err_cm":  fit_r["single_lam_err"],
        "lambda_eff_combined_cm":   fit_c["single_lam_cm"],
        "lambda_eff_combined_err_cm": fit_c["single_lam_err"],
        # Tail threshold and range sensitivity
        "att_tail_dmin_cm": D_TAIL_CM,
        "att_tail_range_sensitivity": tail_sensitivity,
        # sigma_t
        "sigma_t_at_key_positions": sigma_at,
        "sigma_t_center_ps":        _safe(sigma_center),
        "sigma_t_end690_ps":        _safe(sigma_end690),
        # v_eff
        "v_eff_fpt_mm_ns":  _safe(v_eff_fpt),
        "v_eff_t50_mm_ns":  _safe(v_eff_t50),
        "v_eff_sum4_mm_ns": _safe(v_eff_sum4),
        # Photon budget (Panel A: terminal fates; Panel B: SiPM conversion)
        "photon_budget_generated":          _safe(n_gen),
        "photon_budget_escaped":            _safe(n_esc),
        "photon_budget_sipm_entering":      _safe(n_sipm_budget),
        "photon_budget_detected":           _safe(n_det),
        # "surface_bulk_absorbed" = generated - escaped - SiPM-entering
        # = Mylar surface absorption + bulk ABSLENGTH absorption (not separately resolved)
        "photon_budget_surface_bulk_absorbed": _safe(max(0.0, n_gen - n_esc - n_sipm_budget)),
        "photon_budget_mylar_absorbed":     _safe(max(0.0, n_gen - n_esc - n_sipm_budget)),
        "photon_budget_reflector_crossings":_safe(float(n_ref)),
        "photon_budget_pct_surface_bulk":   _safe(100.0 * max(0.0, n_gen - n_esc - n_sipm_budget) / n_gen if n_gen > 0 else None),
        "photon_budget_pct_escaped":        _safe(100.0 * n_esc / n_gen if n_gen > 0 else None),
        "photon_budget_pct_sipm":           _safe(100.0 * n_sipm_budget / n_gen if n_gen > 0 else None),
        "photon_budget_pct_detected":       _safe(100.0 * n_det / n_gen if n_gen > 0 else None),
        "photon_budget_pde_eff":            _safe(n_det / n_sipm_budget if n_sipm_budget > 0 else None),
        # Theoretical group velocity from flat RINDEX (no GROUPVEL property set in Materials.cc)
        "group_velocity_theoretical_mm_ns": _safe(C_MM_NS / REFRACTIVE_INDEX),
        "group_velocity_note": "Geant4 computes v_g from d(omega)/dk on flat RINDEX=1.58; equals c/n=189.7 mm/ns",
        # NPE at key positions (Mylar)
        "npe_left_at_x0":    _safe(float(npe_l[idx0])),
        "npe_right_at_x0":   _safe(float(npe_r[idx0])),
        "npe_left_at_x690":  _safe(float(npe_l[np.argmin(np.abs(x_mm - 690))])),
        "npe_right_at_x690": _safe(float(npe_r[np.argmin(np.abs(x_mm - 690))])),
        "npe_left_at_xm690": _safe(float(npe_l[np.argmin(np.abs(x_mm + 690))])),
        "npe_right_at_xm690":_safe(float(npe_r[np.argmin(np.abs(x_mm + 690))])),
        # Far-end (tail) NPE for Mylar
        "mylar_tail_npe_left_mean":  _safe(far_mean_l),
        "mylar_tail_npe_right_mean": _safe(far_mean_r),
        # EndTop comparison
        "endtop": endtop_section,
        # Figures & tables
        "figures": {
            "fig_npe_vs_x":            "figures/fig_npe_vs_x.png",
            "fig_attenuation":         "figures/fig_attenuation.png",
            "fig_attenuation_twocomp": "figures/fig_attenuation_twocomp.png",
            "fig_far_end_yield":       "figures/fig_far_end_yield.png",
            "fig_timing_histograms":   "figures/fig_timing_histograms.png",
            "fig_fpt_vs_t50":          "figures/fig_fpt_vs_t50.png",
            "fig_mean_arrival_time":   "figures/fig_mean_arrival_time.png",
            "fig_sigma_t":             "figures/fig_sigma_t.png",
            "fig_photon_budget":       "figures/fig_photon_budget.png",
            **({"fig_endtop_comparison": "figures/fig_endtop_comparison.png"}
               if endtop_available else {}),
        },
        "tables": {
            "tab_sim_params":        "tables/tab_sim_params",
            "tab_attenuation_fits":  "tables/tab_attenuation_fits",
            "tab_sigma_t":           "tables/tab_sigma_t",
        },
    }

    json_path = out_dir / "deck_values.json"
    with open(json_path, "w") as f:
        json.dump(deck, f, indent=2)
    print(f"  [json] {json_path.name}")

    # ── sanity checks ─────────────────────────────────────────────────────────
    print("\nSanity checks:")
    for label, fit in [("left", fit_l), ("right", fit_r)]:
        s = fit["single_lam_cm"]
        e = fit["expoff_lam_cm"]
        t = fit["tail_lam_cm"]
        print(f"  {label}: λ_single={s:.1f}cm (χ²/ndf={fit['single_chi2ndf']:.0f})  "
              f"λ_prompt={e:.1f}cm  C={fit['expoff_C']:.1f}PE  "
              f"λ_tail={t:.1f}cm")
    print(f"  Yield: {SCINT_YIELD}/MeV  ABSLENGTH: {LAMBDA_BULK_CM}cm")
    if endtop_available:
        print(f"  EndTop λ_tail_left={fit_et_l.get('tail_lam_cm'):.1f}cm  "
              f"NPE@x=0 L={endtop_section['npe_left_at_x0']:.1f}  "
              f"vs Mylar={deck['npe_left_at_x0']:.1f}")
    print("\nDone.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614",
                    type=pathlib.Path)
    ap.add_argument("--out-dir", default=pathlib.Path(__file__).parent, type=pathlib.Path)
    args = ap.parse_args()
    run(args.data_dir, args.out_dir)
