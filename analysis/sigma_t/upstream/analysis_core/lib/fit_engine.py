#!/usr/bin/env python3.12
"""
fit_engine.py — PyROOT Gaussian fit on the timing-distribution core (EXEC_16).

Design
------
Each call to fit_core_gaussian() processes one 1-D array of t_N values
(one event = one entry) for a given (material, position, group, N) tuple.

Main fit:   PyROOT TF1("gaus") on a windowed TH1F.
Bootstrap:  scipy.optimize.curve_fit on histogram bins (much faster than
            300 ROOT fits per distribution).

Why windowed fit?
    The t_N distributions for TOP SiPMs have a Gaussian core corresponding
    to the arrival-time spread of prompt photons, plus an exponential/flat
    tail from reflected, scattered, and slow-component photons.  The tail
    is physically real but does not represent the timing precision of a fast
    trigger.  Fitting in [peak ± k*σ_MAD] isolates the core.

PyROOT object lifecycle:
    TH1F and TF1 objects are created with SetDirectory(0) so ROOT's global
    TDirectory does not own them; they live in Python's heap and are returned
    to the caller.  The caller (timing_fit_pipeline.py) decides when to write
    them to the sidecar .root file.
"""

import math
import numpy as np
import scipy.optimize as opt

# import ROOT here; the shebang python3.12 guarantees compatibility
import ROOT
ROOT.gROOT.SetBatch(True)          # no display: running in WSL2 without X
ROOT.gErrorIgnoreLevel = ROOT.kWarning  # suppress ROOT Info/Print messages

from .robust_seeds import gather_seeds   # MAD seeds from sibling module


# ─── Gaussian function for scipy bootstrap ────────────────────────────────────

def _gaussian(x: np.ndarray, A: float, mu: float, sigma: float) -> np.ndarray:
    """
    Unnormalised Gaussian: A * exp(-(x-mu)^2 / (2*sigma^2)).

    Used by scipy.optimize.curve_fit during bootstrap resampling.
    scipy is ~50x faster than PyROOT for the bootstrap stage because
    it avoids ROOT object construction overhead.

    Parameters
    ----------
    x : np.ndarray — bin centres (ns)
    A : float      — amplitude (counts)
    mu : float     — mean (ns)
    sigma : float  — standard deviation (ns)

    Returns
    -------
    np.ndarray — Gaussian values at x
    """
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ─── nan result factory ───────────────────────────────────────────────────────

def _nan_result(flag: str, n_events: int = 0) -> dict:
    """
    Return a fully-formed result dict populated with NaN for all fit quantities.

    Having a consistent schema regardless of success/failure means the caller
    can always write a CSV row without special-casing.

    Parameters
    ----------
    flag : str     — reason code ("insufficient_stats", "fit_failed", etc.)
    n_events : int — number of events that were available (for the CSV)

    Returns
    -------
    dict — see fit_core_gaussian() Returns section for schema.
    """
    nan = float("nan")
    return dict(
        sigma_fit=nan,        # core Gaussian σ from TF1 (ns)
        sigma_fit_err=nan,    # TF1 parameter error on σ (ns)
        mu_fit=nan,           # core Gaussian μ from TF1 (ns)
        mu_fit_err=nan,       # TF1 parameter error on μ (ns)
        chi2_ndf=nan,         # goodness-of-fit: χ²/NDF
        fit_status=-1,        # TFitResult.Status(); -1 = not attempted
        sigma_mad=nan,        # σ_MAD seed (ns) — reference
        sigma_std=nan,        # raw np.std (ns) — comparison only
        bootstrap_err=nan,    # empirical σ uncertainty from resampling (ns)
        n_events=n_events,
        efficiency=nan,       # set by caller (depends on total events)
        flag=flag,            # human-readable reason code
        h_root=None,          # TH1F (not created on failure)
        f_root=None,          # TF1  (not created on failure)
        fit_window=(nan, nan),
        n_bins=0,
    )


# ─── bootstrap using scipy ────────────────────────────────────────────────────

def _bootstrap_sigma(
    v: np.ndarray,
    n_bins: int,
    win_lo: float,
    win_hi: float,
    seeds: dict,
    n_boot: int,
    rng: np.random.Generator,
) -> float:
    """
    Estimate the standard error of σ_fit by parametric bootstrap.

    For each resample:
      1. Draw N samples with replacement from v.
      2. Histogram into n_bins within [win_lo, win_hi].
      3. Fit a Gaussian with scipy.optimize.curve_fit.
      4. Collect the fitted σ value.

    The standard deviation of the collected σ values is the bootstrap
    standard error — an empirical error that makes no assumption about
    the shape of the distribution.

    scipy is used instead of PyROOT for speed: 300 ROOT TH1/TF1 fits
    per distribution × 400 distributions would take several minutes.
    scipy on histogrammed bins runs ~50x faster with consistent results.

    Parameters
    ----------
    v : np.ndarray
        Full array of t_N values (ns).
    n_bins : int
        Number of bins for the bootstrap histogram (same as main fit).
    win_lo, win_hi : float
        Fit window boundaries (ns).
    seeds : dict
        Output of gather_seeds(); provides (peak, peak_count, sigma_seed).
    n_boot : int
        Number of bootstrap resamples (from config N_BOOTSTRAP).
    rng : np.random.Generator
        Seeded RNG for reproducibility.

    Returns
    -------
    float
        Standard error of σ_fit (ns).  NaN if < 10 successful resamples.
    """
    n = len(v)
    sigma_boots = []

    for _ in range(n_boot):
        # resample with replacement
        v_b = rng.choice(v, size=n, replace=True)

        # histogram only the window region — same strategy as the main fit
        counts_b, edges_b = np.histogram(v_b, bins=n_bins, range=(win_lo, win_hi))
        centers_b = 0.5 * (edges_b[:-1] + edges_b[1:])  # bin centres

        ok = counts_b > 0   # skip empty bins to avoid divide-by-zero in curve_fit
        if ok.sum() < 4:    # need at least 4 points to constrain 3-parameter Gaussian
            continue

        try:
            # initial parameter guesses: reuse seeds from main fit
            p0 = [seeds["peak_count"], seeds["peak"], seeds["sigma_seed"]]
            popt, _ = opt.curve_fit(
                _gaussian,
                centers_b[ok],
                counts_b[ok].astype(float),
                p0=p0,
                maxfev=200,    # fast tolerance; bootstrap only needs approximate fits
            )
            sigma_b = abs(popt[2])  # take absolute value: σ must be positive

            # reject unphysical bootstrap results: σ beyond 10× seed is a failed fit
            if sigma_b < 10 * seeds["sigma_seed"]:
                sigma_boots.append(sigma_b)
        except (RuntimeError, ValueError):
            # curve_fit did not converge for this resample → skip quietly
            continue

    if len(sigma_boots) < 10:          # too few successes → bootstrap is unreliable
        return float("nan")
    return float(np.std(sigma_boots, ddof=1))   # std of bootstrap σ values = standard error


# ─── main fit function ────────────────────────────────────────────────────────

def fit_core_gaussian(
    v: np.ndarray,
    cfg: dict,
    name_prefix: str = "fit",
) -> dict:
    """
    Fit a Gaussian to the core of a t_N timing distribution using PyROOT.

    This replaces np.std() as the resolution estimator.  The key change:
    instead of weighting all events equally (as std does), the fit restricts
    to the Gaussian core [peak ± k*σ_MAD], where tail events (slow photons,
    reflections) contribute negligible weight.

    Fit flow
    --------
    1. Statistics gate: if len(v) < MIN_EVENTS_FOR_FIT → return NaN + flag.
    2. Gather seeds: peak (μ), σ_MAD (σ), peak_count (amplitude) from robust_seeds.
    3. Fit window: [peak - k*σ_MAD, peak + k*σ_MAD] where k = FIT_WINDOW_SIGMAS.
    4. Build TH1F with n_bins, fill all events.
    5. Fit TF1("gaus") restricted to the window with options from config.
    6. Extract σ_fit, μ_fit, χ²/NDF, fit_status from TFitResult.
    7. Bootstrap with scipy to get empirical σ error.
    8. Return dict with scalars + ROOT objects (for sidecar).

    Parameters
    ----------
    v : np.ndarray
        1-D array of t_N values (ns) for one (position, group, N) cell.
        Pre-filtered: NaN already removed by the caller.
    cfg : dict
        Loaded exec16_config.yaml as a dict.
    name_prefix : str
        Prefix for ROOT object names (must be unique within one ROOT session).
        Convention: "{material}_{group}_x{x}mm_N{N}" e.g. "EJ230_TOP_SUM4_x0mm_N4"

    Returns
    -------
    dict with keys:
        sigma_fit     : float — core Gaussian σ (ns)
        sigma_fit_err : float — TF1 parameter error on σ (ns)
        mu_fit        : float — core Gaussian μ (ns)
        mu_fit_err    : float — TF1 parameter error on μ (ns)
        chi2_ndf      : float — χ²/NDF (goodness of fit)
        fit_status    : int   — 0 = converged; non-zero = problem
        sigma_mad     : float — σ_MAD seed (ns); reference for diagnostics
        sigma_std     : float — raw np.std (ns); comparison only, NOT the metric
        bootstrap_err : float — empirical std error of σ_fit from resampling (ns)
        n_events      : int   — number of events used
        efficiency    : float — NaN here; set by caller using n_events / n_total
        flag          : str   — "ok", "insufficient_stats", "fit_failed", etc.
        h_root        : ROOT.TH1F — histogram object for sidecar (caller manages)
        f_root        : ROOT.TF1  — fitted Gaussian for sidecar (caller manages)
        fit_window    : tuple — (win_lo, win_hi) in ns
        n_bins        : int   — number of bins used

    Raises
    ------
    Does not raise; failures are signalled via the `flag` field.
    """

    # ── 1. Statistics gate ────────────────────────────────────────────────────
    n = len(v)
    min_ev = cfg["MIN_EVENTS_FOR_FIT"]   # from config; typically 30
    if n < min_ev:
        # too few events to build a meaningful histogram + fit
        return _nan_result("insufficient_stats", n_events=n)

    # ── 2. Gather seeds from robust estimators ────────────────────────────────
    binning_strategy = cfg.get("BINNING_STRATEGY", "sqrt_n")
    seeds = gather_seeds(v, binning_strategy)    # returns peak, sigma_seed, median, n_bins, sigma_std
    peak        = seeds["peak"]                  # μ seed (histogram mode)
    sigma_seed  = seeds["sigma_seed"]            # σ_MAD: the authoritative seed
    peak_count  = seeds["peak_count"]            # amplitude seed
    n_bins      = seeds["n_bins"]                # chosen by √N or FD rule

    # ── 3. Fit window ─────────────────────────────────────────────────────────
    k        = cfg["FIT_WINDOW_SIGMAS"]          # typically 2.0 σ
    win_lo   = peak - k * sigma_seed             # lower bound of core window
    win_hi   = peak + k * sigma_seed             # upper bound of core window

    # sanity: window must be at least one bin wide
    if win_hi - win_lo < 1e-6:
        return _nan_result("degenerate_window", n_events=n)

    # ── 4. Build TH1F ─────────────────────────────────────────────────────────
    # Histogram covers the FULL data range (not just the window) so the tail
    # is visible in diagnostic overlay plots.  The fit is restricted to the window.
    v_min = float(v.min()) - 0.5 * (v.max() - v.min()) / n_bins  # small left margin
    v_max = float(v.max()) + 0.5 * (v.max() - v.min()) / n_bins  # small right margin

    h_name = f"h_{name_prefix}"  # ROOT name: must be unique in the session
    h = ROOT.TH1F(h_name, "", n_bins, v_min, v_max)
    h.SetDirectory(0)            # disown from ROOT global directory → Python manages lifetime

    for t in v:                  # fill the TH1 event-by-event
        h.Fill(float(t))

    # ── 5. Build TF1 and set initial parameters ───────────────────────────────
    f_name = f"f_{name_prefix}"
    # "gaus" is ROOT's built-in unnormalised Gaussian: p0 * exp(-0.5*((x-p1)/p2)^2)
    f = ROOT.TF1(f_name, "gaus", win_lo, win_hi)
    f.SetParameters(peak_count, peak, sigma_seed)  # (amplitude, μ, σ) seeds

    # ── 6. Fit ────────────────────────────────────────────────────────────────
    fit_options = cfg["FIT_OPTIONS"]   # "R Q S 0" from config
    # "R" restricts fit to [win_lo, win_hi]
    # "Q" quiet (no printout)
    # "S" return TFitResult object
    # "0" do not draw the function on the histogram automatically
    fit_result = h.Fit(f_name, fit_options, "", win_lo, win_hi)

    # ── 7. Extract results ────────────────────────────────────────────────────
    # GetParameter index: 0=amplitude, 1=μ, 2=σ  (for "gaus" TF1)
    sigma_fit     = abs(float(f.GetParameter(2)))    # abs: σ must be positive
    sigma_fit_err = abs(float(f.GetParError(2)))
    mu_fit        = float(f.GetParameter(1))
    mu_fit_err    = float(f.GetParError(1))

    # χ²/NDF: goodness of fit; ~1 means Gaussian is a good model for the core
    ndf = f.GetNDF()
    chi2_ndf = float(f.GetChisquare() / ndf) if ndf > 0 else float("nan")

    # fit_status: 0 = MINUIT converged, non-zero = problem
    # TFitResult.IsValid() checks that the error matrix is valid
    try:
        fit_status = int(fit_result.Status())
    except Exception:
        fit_status = -1   # TFitResult not available (shouldn't happen with "S")

    # flag converged/failed for QA
    if sigma_fit <= 0 or math.isnan(sigma_fit) or math.isinf(sigma_fit):
        return _nan_result("fit_failed_sigma_invalid", n_events=n)

    # ── 8. Bootstrap ─────────────────────────────────────────────────────────
    rng = np.random.default_rng(int(cfg["RANDOM_SEED"]))
    n_boot = int(cfg["N_BOOTSTRAP"])
    bootstrap_err = _bootstrap_sigma(v, n_bins, win_lo, win_hi, seeds, n_boot, rng)

    # ── 9. Return ─────────────────────────────────────────────────────────────
    return dict(
        sigma_fit     = sigma_fit,            # ← THE reported physics result (ns)
        sigma_fit_err = sigma_fit_err,        # TF1 error (ns)
        mu_fit        = mu_fit,               # fitted peak position (ns)
        mu_fit_err    = mu_fit_err,
        chi2_ndf      = chi2_ndf,             # goodness-of-fit; should be ~1 for Gaussian core
        fit_status    = fit_status,
        sigma_mad     = seeds["sigma_seed"],  # σ_MAD seed (ns) — reference, not the metric
        sigma_std     = seeds["sigma_std"],   # raw std (ns) — stored for comparison, not the metric
        bootstrap_err = bootstrap_err,        # empirical std error of σ_fit from resampling (ns)
        n_events      = n,
        efficiency    = float("nan"),         # caller sets this: n_events / n_total
        flag          = "ok",
        h_root        = h,                    # TH1F: full distribution, for sidecar + overlay
        f_root        = f,                    # TF1:  fitted Gaussian, for sidecar + overlay
        fit_window    = (win_lo, win_hi),     # stored in meta.json
        n_bins        = n_bins,
    )
