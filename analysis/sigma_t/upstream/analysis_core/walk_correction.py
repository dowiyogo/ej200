#!/usr/bin/env python3.12
"""
walk_correction.py — HOOK_WALK: leading-edge time-walk correction (EXEC_18).

Reutilizable entre campañas. Corrige el slewing de umbral N en el stream de
PE de un estimador dado. Variable de slewing: NPE total en el stream.

Physics:
    For leading-edge triggering with threshold N=1, events with more PE
    have earlier arrival times (negative correlation t_raw ↔ NPE).
    Walk model: t = α + β/√s  (monotonically decreasing with NPE s).
    Correction: t_corr = t_raw − (walk(s) − walk(s_ref))
    anchored at s_ref = median(s) so the mean is preserved.

Validation: |corr(t_corr, s)| ≈ 0 AND σ_corr ≤ σ_raw.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator


# ── Walk model: α + β/√s ─────────────────────────────────────────────────────

def _walk_model(s: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    return alpha + beta / np.sqrt(np.maximum(s, 1e-6))


# ── Robust median curve in bins of s ─────────────────────────────────────────

def _median_curve(t: np.ndarray, s: np.ndarray, n_bins: int = 20,
                  min_in_bin: int = 30) -> tuple:
    """
    Bin t by s (percentile-based equal-count bins), return median per bin.

    Returns (s_centers, t_medians) for bins with ≥ min_in_bin events.
    """
    valid = np.isfinite(t) & np.isfinite(s) & (s > 0)
    tv, sv = t[valid], s[valid]
    if len(tv) < min_in_bin * 3:
        return np.array([]), np.array([])
    pct = np.linspace(0, 100, n_bins + 1)
    edges = np.unique(np.percentile(sv, pct))
    if len(edges) < 3:
        return np.array([]), np.array([])
    centers, medians = [], []
    for i in range(len(edges) - 1):
        mask = (sv >= edges[i]) & (sv < edges[i + 1])
        if mask.sum() >= min_in_bin:
            centers.append(float(np.median(sv[mask])))
            medians.append(float(np.median(tv[mask])))
    return np.array(centers), np.array(medians)


# ── Fit walk curve ────────────────────────────────────────────────────────────

def fit_walk(t: np.ndarray, s: np.ndarray, n_bins: int = 20) -> dict:
    """
    Fit walk curve to t vs s.

    Tries α + β/√s (parametric). If residual correlation |r| > 0.15 after
    parametric fit, falls back to a monotone PchipInterpolator spline on the
    median curve.

    Returns dict with keys:
        form:       "parametric" or "spline"
        popt:       (alpha, beta) for parametric; None for spline
        spline:     PchipInterpolator or None
        s_ref:      reference NPE (median of s)
        walk_at_ref: walk(s_ref) for anchoring
        n_points:   number of (s,t) pairs used
        r_raw:      Pearson r(t, s) before correction
        converged:  bool
    """
    valid = np.isfinite(t) & np.isfinite(s) & (s > 0)
    tv, sv = t[valid], s[valid]
    s_ref = float(np.median(sv))
    r_raw = float(np.corrcoef(tv, sv)[0, 1]) if len(tv) > 5 else 0.0

    centers, medians = _median_curve(tv, sv, n_bins=n_bins)
    result = dict(popt=None, spline=None, s_ref=s_ref,
                  walk_at_ref=float(np.median(tv)),  # default: no correction
                  n_points=int(valid.sum()), r_raw=r_raw,
                  converged=False, form="none")

    if len(centers) < 4:
        return result

    # Try parametric fit
    try:
        # Initial guess: alpha = median(t), beta > 0 if r_raw < 0 (more PE → earlier t)
        beta_sign = -np.sign(r_raw) if abs(r_raw) > 0.05 else 1.0
        p0 = [float(np.median(medians)), beta_sign * float(np.std(medians)) * np.sqrt(float(np.median(centers)))]
        popt, _ = curve_fit(_walk_model, centers, medians, p0=p0, maxfev=3000)
        walk_at_ref = float(_walk_model(np.array([s_ref]), *popt)[0])
        result.update(form="parametric", popt=popt,
                      walk_at_ref=walk_at_ref, converged=True)
    except Exception:
        pass

    # Check closure with parametric; fall back to spline if needed
    if result["converged"]:
        t_corr_test = tv - (_walk_model(sv, *result["popt"]) - result["walk_at_ref"])
        r_corr = float(np.corrcoef(t_corr_test, sv)[0, 1])
        if abs(r_corr) > 0.15:
            result["converged"] = False  # parametric insufficient

    if not result["converged"] and len(centers) >= 4:
        # Monotone spline
        try:
            spl = PchipInterpolator(centers, medians, extrapolate=True)
            walk_at_ref = float(spl(s_ref))
            result.update(form="spline", spline=spl,
                          walk_at_ref=walk_at_ref, converged=True)
        except Exception:
            pass

    return result


# ── Apply correction ──────────────────────────────────────────────────────────

def apply_correction(t: np.ndarray, s: np.ndarray, fit: dict) -> np.ndarray:
    """
    Apply walk correction: t_corr = t − (walk(s) − walk(s_ref)).

    Events with nan in s or t get nan in t_corr.
    """
    t_corr = t.copy()
    valid = np.isfinite(t) & np.isfinite(s) & (s > 0)
    if not fit["converged"]:
        return t_corr   # no correction applied

    if fit["form"] == "parametric":
        walk_s = _walk_model(s[valid], *fit["popt"])
    elif fit["form"] == "spline":
        walk_s = fit["spline"](s[valid])
    else:
        return t_corr

    t_corr[valid] = t[valid] - (walk_s - fit["walk_at_ref"])
    return t_corr


# ── Validate correction ────────────────────────────────────────────────────────

def validate_correction(t_raw: np.ndarray, t_corr: np.ndarray,
                        s: np.ndarray, sigma_raw_ns: float,
                        sigma_corr_ns: float) -> dict:
    """
    Check that correction is valid:
      1. |corr(t_corr, s)| < 0.15  (residual correlation removed)
      2. sigma_corr <= sigma_raw    (timing improves or holds)

    Returns status dict with 'ok' flag, correlation before/after, delta_sigma.
    """
    valid = np.isfinite(t_raw) & np.isfinite(t_corr) & np.isfinite(s) & (s > 0)
    r_raw  = float(np.corrcoef(t_raw[valid],  s[valid])[0, 1]) if valid.sum() > 5 else 0.0
    r_corr = float(np.corrcoef(t_corr[valid], s[valid])[0, 1]) if valid.sum() > 5 else 0.0
    delta_sigma = sigma_corr_ns - sigma_raw_ns  # negative = improvement

    corr_ok  = abs(r_corr) < 0.15
    sigma_ok = sigma_corr_ns <= sigma_raw_ns * 1.02  # allow 2% rounding

    return dict(r_raw=r_raw, r_corr=r_corr, delta_sigma_ps=(delta_sigma * 1000),
                corr_ok=corr_ok, sigma_ok=sigma_ok,
                ok=(corr_ok and sigma_ok),
                status="ok" if (corr_ok and sigma_ok) else "walk_failed")
