"""EXEC35 quantile widths. Imported function bodies are loaded, never mains."""
import ast
from pathlib import Path
import numpy as np
from scipy.stats import skew, kurtosis

PACKAGE = Path(__file__).resolve().parents[1]
CORE_SOURCE = PACKAGE/'upstream/related/b0aaac1/analysis/exec12t_timing_threshold_analysis.py'
IQR_SOURCE = PACKAGE/'upstream/related/420addf/analysis/optim/phase_ab.py'


def function_only(path, name, namespace):
    node = next(n for n in ast.parse(path.read_text()).body
                if isinstance(n, ast.FunctionDef) and n.name == name)
    module = ast.Module(body=[node], type_ignores=[])
    exec(compile(module, str(path), 'exec'), namespace)
    return namespace[name]


class ExactCentralQuantiles:
    """Explicit adapter: legacy 16/84 levels become requested 15.865/84.135."""
    def __getattr__(self, name):
        return getattr(np, name)

    def quantile(self, values, levels, *args, **kwargs):
        if np.array_equal(levels, [.16, .84]):
            levels = [.15865, .84135]
        return np.quantile(values, levels, *args, **kwargs)


legacy_robust = function_only(CORE_SOURCE, 'robust', dict(np=np, skew=skew, kurtosis=kurtosis))
exact_robust = function_only(CORE_SOURCE, 'robust', dict(np=ExactCentralQuantiles(), skew=skew, kurtosis=kurtosis))
rsigma = function_only(IQR_SOURCE, 'rsigma', dict(np=np))
NAMES = ['sigma_core', 'sigma_IQR', 'fwhm_over_2_355', 'q15_865', 'q50',
         'q84_135', 'q2_5', 'q97_5', 'asymmetry']


def histogram_edges(values):
    v = np.asarray(values)[np.isfinite(values)]
    if not len(v) or v.max() == v.min():
        return np.array([0., 1.])
    return np.linspace(v.min(), v.max(), int(np.ceil(np.sqrt(len(v)))) + 1)


def interpolated_fwhm(counts, edges):
    counts = np.asarray(counts, float)
    if not len(counts) or counts.max() <= 0:
        return float('nan')
    centers = (edges[:-1] + edges[1:]) / 2
    width = edges[1] - edges[0]
    x = np.r_[centers[0] - width, centers, centers[-1] + width]
    y = np.r_[0., counts, 0.]
    peak = int(np.argmax(y)); half = y[peak] / 2
    left = peak
    while left > 0 and y[left] > half:
        left -= 1
    right = peak
    while right < len(y)-1 and y[right] > half:
        right += 1
    if left == peak or right == peak:
        return float('nan')
    a = x[left] + (half-y[left])/(y[left+1]-y[left])*(x[left+1]-x[left])
    b = x[right-1] + (half-y[right-1])/(y[right]-y[right-1])*(x[right]-x[right-1])
    return float(b-a)


def estimates(values, edges):
    v = np.asarray(values); v = v[np.isfinite(v)]
    if not len(v):
        return np.full(len(NAMES), np.nan)
    core = exact_robust(v)['sigma_core']
    qlo, med, qhi, q025, q975 = np.quantile(v, [.15865, .5, .84135, .025, .975])
    fwhm = interpolated_fwhm(np.histogram(v, edges)[0], edges)
    return np.array([core, rsigma(v), fwhm/2.355, qlo, med, qhi, q025, q975,
                     qhi + qlo - 2*med])


def measure(events, replicas=500, seed=35091201, callback=None):
    """One timestamp/NaN per generated event; bootstrap the events, not bins."""
    events = np.asarray(events, float)
    finite = np.isfinite(events)
    edges = histogram_edges(events)
    point = estimates(events, edges)
    rng = np.random.default_rng(seed)
    samples = np.empty((replicas, len(NAMES)))
    extra = []
    for i in range(replicas):
        drawn = events[rng.integers(0, len(events), len(events))]
        samples[i] = estimates(drawn, edges)
        if callback:
            extra.append(callback(drawn[np.isfinite(drawn)], i))
    row = dict(n_eff=int(finite.sum()), N_generated_partition=len(events),
               efficiency=float(finite.mean()), discarded=int((~finite).sum()))
    for j, name in enumerate(NAMES):
        valid = samples[:, j][np.isfinite(samples[:, j])]
        row[name] = float(point[j])
        row[name+'_se'] = float(np.std(valid, ddof=1)) if len(valid)>1 else np.nan
        row[name+'_bootstrap_fraction'] = len(valid)/replicas
        row[name+'_ci_low'], row[name+'_ci_high'] = (np.quantile(valid, [.025, .975]).tolist()
                                                   if len(valid) else [np.nan, np.nan])
    return row, samples, edges, extra


def fixed_timestamps(event, gid, time_ns, channels, n_generated=10000):
    """Exact merged-stream order statistics with fixed IDs, no ranking/argmin."""
    selected = np.isin(gid, channels)
    ev, ts = event[selected], time_ns[selected]
    counts = np.bincount(ev, minlength=n_generated)
    order = np.lexsort((ts, ev)); ts = ts[order]
    starts = np.r_[0, np.cumsum(counts[:-1])]
    result = np.full((n_generated, 20), np.nan)
    for n in range(1, 21):
        eligible = counts >= n
        result[eligible, n-1] = ts[starts[eligible]+n-1]
    return result
