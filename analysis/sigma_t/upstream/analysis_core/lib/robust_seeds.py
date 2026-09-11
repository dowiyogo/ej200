#!/usr/bin/env python3.12
"""
robust_seeds.py — Robust seed estimators for the Gaussian fit in EXEC_16.

Purpose
-------
Provides the MAD×1.4826 scaled sigma estimator (René's authoritative robust σ),
the histogram-mode peak finder, and a master gather_seeds() that packages all
seeds needed by fit_engine.py to initialize the PyROOT TF1 Gaussian fit.

Why not just use np.std as the final metric?
    np.std weights residuals by x², so a handful of tail events (delayed
    reflections, slow scintillation photons visible at ~1.4–1.6 ns in the
    TOP histograms) inflate it.  MAD weighs by |x|; it is insensitive to
    outliers beyond the scale of the core.  The factor 1.4826 makes it
    converge to σ for a Gaussian population.  Here σ_MAD serves as SEED
    only — the reported metric is always σ_fit from the core Gaussian fit.

Peak finder choice: histogram mode
    We bin the data and return the center of the highest-count bin.
    KDE was rejected because bandwidth selection introduces its own bias.
    Raw argmax of individual event times is noise-sensitive.
    The histogram mode is deterministic and resolution-controlled by n_bins.
"""

import numpy as np


# ─── binning rules ────────────────────────────────────────────────────────────

def compute_n_bins_sqrt(n: int) -> int:
    """
    Number of bins by the sqrt-N rule.

    Parameters
    ----------
    n : int  — number of data points.

    Returns
    -------
    int  — n_bins, floored at 10 to prevent degenerate TH1.

    Notes
    -----
    For N=300 events this gives ~17 bins — coarse but appropriate for a seed
    stage.  The Freedman-Diaconis alternative is available via BINNING_STRATEGY.
    """
    return max(10, int(np.ceil(np.sqrt(n))))  # floor at 10 → TH1 never has < 10 bins


def compute_n_bins_fd(v: np.ndarray) -> int:
    """
    Number of bins by the Freedman-Diaconis rule: width = 2*IQR / n^(1/3).

    Parameters
    ----------
    v : np.ndarray  — 1-D array of time values (ns).

    Returns
    -------
    int  — n_bins clamped to [10, 200].
    """
    n = len(v)
    q25, q75 = np.percentile(v, [25, 75])   # IQR is robust to outliers
    iqr = q75 - q25
    if iqr < 1e-9:          # all values are identical → degenerate; fall back
        return 20
    bw = 2.0 * iqr / (n ** (1.0 / 3.0))   # FD bin width formula
    rng = float(v.max() - v.min())
    n_bins = int(np.ceil(rng / bw)) if bw > 0 else 20
    return max(10, min(200, n_bins))         # clamp: never fewer than 10, never more than 200


# ─── MAD×1.4826 robust sigma ─────────────────────────────────────────────────

def compute_mad_sigma(v: np.ndarray) -> tuple:
    """
    Compute the MAD-scaled robust sigma and the sample median.

    σ_MAD = 1.4826 × median( |v − median(v)| )

    For a Gaussian distribution, MAD → σ/1.4826 as N→∞, so multiplying
    by 1.4826 gives an asymptotically unbiased estimator of σ that is
    robust to up to ~50% contamination in the tails.

    Parameters
    ----------
    v : np.ndarray
        1-D array of t_N values (ns), non-NaN, length >= 1.

    Returns
    -------
    sigma_mad : float
        Robust sigma in ns.  This is René's authoritative seed.
    median_v : float
        Median of v in ns.  Used as μ cross-check for the fit seed.

    Raises
    ------
    ValueError
        If v is empty.
    """
    if len(v) == 0:
        raise ValueError("compute_mad_sigma: cannot process an empty array")

    median_v = float(np.median(v))               # median: robust central location
    mad = float(np.median(np.abs(v - median_v))) # MAD: median of |deviations from median|
    sigma_mad = 1.4826 * mad                     # scale factor for Gaussian consistency

    # Floor at 0.1 ps expressed in ns: prevents the fit window from collapsing
    # if the distribution is extremely narrow (e.g., all events at the same time)
    sigma_mad = max(sigma_mad, 1e-4)

    return sigma_mad, median_v


# ─── histogram peak finder ────────────────────────────────────────────────────

def find_histogram_peak(v: np.ndarray, n_bins: int) -> tuple:
    """
    Find the mode of the distribution via histogram binning.

    Returns the center of the highest-count bin and its count.
    These are the (μ, amplitude) seeds for the PyROOT TF1("gaus").

    Parameters
    ----------
    v : np.ndarray
        1-D array of t_N values (ns).
    n_bins : int
        Number of histogram bins.

    Returns
    -------
    peak_center : float
        Center of the highest-count bin (ns).
    peak_count : float
        Count in that bin.  Amplitude seed for TF1.
    """
    counts, edges = np.histogram(v, bins=n_bins)  # bin the data
    peak_bin = int(np.argmax(counts))              # index of the tallest bin
    # center of the tallest bin = average of its left and right edges
    peak_center = float(0.5 * (edges[peak_bin] + edges[peak_bin + 1]))
    peak_count = float(counts[peak_bin])           # height = amplitude seed
    return peak_center, peak_count


# ─── master seed function ─────────────────────────────────────────────────────

def gather_seeds(v: np.ndarray, binning_strategy: str = "sqrt_n") -> dict:
    """
    Package all seeds needed to initialize the Gaussian TF1 fit.

    Parameters
    ----------
    v : np.ndarray
        1-D array of valid (non-NaN) t_N values (ns).
    binning_strategy : str
        "sqrt_n"  → √N rule (default, appropriate for N~30–300)
        "fd"      → Freedman-Diaconis rule (better for N>1000)

    Returns
    -------
    dict
        peak       : float — histogram-mode center (ns), used as μ seed
        peak_count : float — amplitude seed for TF1("gaus")
        sigma_seed : float — σ_MAD in ns (René's authoritative seed)
        median     : float — median (ns), sanity cross-check for peak
        n_bins     : int   — bins chosen by the rule (stored in meta.json)
        sigma_std  : float — np.std(v, ddof=1) (ns), stored for comparison ONLY

    Notes
    -----
    sigma_seed = σ_MAD is used ONLY to:
      (a) define the fit window [peak ± k*σ_MAD]
      (b) initialize the TF1 σ parameter
    The pipeline reports σ_fit from the Gaussian fit, never σ_MAD or σ_std
    as the final physics result.  sigma_std is saved in the CSV so the
    improvement from std → fit can be audited in the thesis.
    """
    n = len(v)

    # select binning rule from config hook
    if binning_strategy == "fd":
        n_bins = compute_n_bins_fd(v)        # Freedman-Diaconis for large N
    else:
        n_bins = compute_n_bins_sqrt(n)       # default √N rule

    # robust estimators
    sigma_mad, median_v = compute_mad_sigma(v)   # σ_MAD and median

    # histogram-mode peak: more stable than raw event maximum
    peak, peak_count = find_histogram_peak(v, n_bins)

    # raw std kept for comparison in sidecar CSV; NOT the reported metric
    sigma_std = float(np.std(v, ddof=1)) if n > 1 else float("nan")

    return {
        "peak":        peak,         # μ seed for TF1
        "peak_count":  peak_count,   # amplitude seed for TF1
        "sigma_seed":  sigma_mad,    # σ seed = σ_MAD (authoritative robust estimator)
        "median":      median_v,     # sanity check: should be close to peak for unimodal dist
        "n_bins":      n_bins,       # stored in meta.json for reproducibility
        "sigma_std":   sigma_std,    # comparison only — shows how much std overestimates core σ
    }
