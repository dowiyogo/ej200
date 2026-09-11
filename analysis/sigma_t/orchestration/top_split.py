"""New EXEC34 orchestration; all imported estimator primitives stay unchanged.

Fit seeds AND histogram axes are learned in TRAIN. EVAL only estimates its
distribution's Gaussian parameters, with the frozen estimator construction.
The scoped dependency injection does not edit upstream files or functions.
ROOT fits run sequentially; this module is deliberately not thread safe.
"""
from contextlib import contextmanager
from pathlib import Path
import sys
import math

import numpy as np

CORE = Path(__file__).resolve().parents[1]/'upstream/analysis_core'
sys.path.insert(0, str(CORE))
from timing_fit_pipeline import best_gids_by_count, compute_tN
from lib import fit_engine
from lib.robust_seeds import gather_seeds

CONFIG = dict(FIT_WINDOW_SIGMAS=2.0, FIT_OPTIONS='R Q S 0',
              MIN_EVENTS_FOR_FIT=30, N_BOOTSTRAP=300, RANDOM_SEED=20260618,
              BINNING_STRATEGY='sqrt_n', EFFICIENCY_FLOOR=0.05)
N_VALUES = list(range(1, 21))


def generated_ids(n_generated, parity=None):
    ids = np.arange(n_generated, dtype=np.int64)
    return ids if parity is None else ids[ids % 2 == parity]


def select_channels(arr, parity):
    mask = arr['face_type'] == 2
    if parity is not None:
        mask &= arr['event_id'] % 2 == parity
    return best_gids_by_count(mask, arr['global_id'], 4).astype(int).tolist()


def timestamps(arr, channels, parity, n_generated, index):
    mask = (arr['face_type'] == 2) & np.isin(arr['global_id'], channels)
    if parity is not None:
        mask &= arr['event_id'] % 2 == parity
    return compute_tN(arr['event_id'], arr['time_ns'], mask, index,
                      generated_ids(n_generated, parity))


def fit_state(values, cfg):
    if len(values) < cfg['MIN_EVENTS_FOR_FIT']:
        return None
    seeds = gather_seeds(values, cfg['BINNING_STRATEGY'])
    padding = 0.5*(values.max()-values.min())/seeds['n_bins']
    return dict(seeds=seeds, histogram_lo=float(values.min()-padding),
                histogram_hi=float(values.max()+padding),
                random_seed=cfg['RANDOM_SEED'])


class FrozenROOT:
    def __init__(self, root, state):
        self.root, self.state = root, state

    def __getattr__(self, name):
        return getattr(self.root, name)

    def TH1F(self, name, title, bins, lo, hi):
        return self.root.TH1F(name, title, self.state['seeds']['n_bins'],
                             self.state['histogram_lo'], self.state['histogram_hi'])


@contextmanager
def frozen_fit_construction(state):
    original_seeds, original_root = fit_engine.gather_seeds, fit_engine.ROOT
    if state is not None:
        fit_engine.gather_seeds = lambda *args, **kwargs: dict(state['seeds'])
        fit_engine.ROOT = FrozenROOT(original_root, state)
    try:
        yield
    finally:
        fit_engine.gather_seeds, fit_engine.ROOT = original_seeds, original_root


def fit(values, denominator, cfg, name, state, root_file=None):
    if state is None:
        # No valid TRAIN construction exists. Do not derive one from EVAL.
        result = fit_engine._nan_result('no_valid_training_fit_state', len(values))
    else:
        with frozen_fit_construction(state):
            result = fit_engine.fit_core_gaussian(values, cfg, name)
    histogram, function = result.pop('h_root', None), result.pop('f_root', None)
    result['ndf'] = int(function.GetNDF()) if function else 0
    result['chi2'] = float(function.GetChisquare()) if function else float('nan')
    if root_file is not None:
        root_file.cd()
        if histogram: histogram.Write()
        if function: function.Write()
    result['n_eff'] = len(values)
    result['N_generated_partition'] = denominator
    result['efficiency'] = len(values)/denominator if denominator else 0.0
    result['discarded_fewer_than_N'] = denominator-len(values)
    result['sigma_ps'] = 1000*result['sigma_fit']
    result['fit_error_ps'] = 1000*result['sigma_fit_err']
    result['bootstrap_error_ps'] = 1000*result['bootstrap_err']
    finite_errors = [x for x in (result['fit_error_ps'], result['bootstrap_error_ps']) if math.isfinite(x)]
    result['uncertainty_ps'] = max(finite_errors) if finite_errors else float('nan')
    result['uncertainty_definition'] = 'max(ROOT covariance error, fixed-estimator bootstrap error); conditional on frozen selection'
    result['fit_valid'] = bool(result['fit_status'] == 0 and math.isfinite(result['sigma_ps']) and result['sigma_ps'] > 0)
    result['measurement_label'] = 'INTRINSIC timing resolution — electronics not included'
    if state:
        lo = state['seeds']['peak']-cfg['FIT_WINDOW_SIGMAS']*state['seeds']['sigma_seed']
        hi = state['seeds']['peak']+cfg['FIT_WINDOW_SIGMAS']*state['seeds']['sigma_seed']
        result['n_in_fit_window'] = int(np.count_nonzero((values >= lo) & (values <= hi)))
        result['n_outside_fit_window'] = len(values)-result['n_in_fit_window']
        result['histogram_underflow'] = int(np.count_nonzero(values < state['histogram_lo']))
        result['histogram_overflow'] = int(np.count_nonzero(values >= state['histogram_hi']))
    return result


def winner(curve):
    valid = [r for r in curve if r['fit_valid']]
    # No efficiency optimization or unrequested quality threshold in selection.
    # EFFICIENCY_FLOOR is evaluated independently by G-P.3.
    return min(valid, key=lambda r: (r['sigma_ps'], r['N'])) if valid else None


def learn(arr, n_generated, train_parity, cfg=None, root_file=None):
    cfg = dict(CONFIG if cfg is None else cfg)
    channels = select_channels(arr, train_parity)
    curve, states = [], {}
    denominator = len(generated_ids(n_generated, train_parity))
    prefix = 'all' if train_parity is None else f'train{train_parity}'
    for index in N_VALUES:
        values = timestamps(arr, channels, train_parity, n_generated, index)
        state = fit_state(values, cfg)
        states[index] = state
        row = fit(values, denominator, cfg, f'{prefix}_N{index}', state, root_file)
        row.update(N=index, sample='ALL' if train_parity is None else 'TRAIN')
        curve.append(row)
    selected = winner(curve)
    return dict(channels=channels, train_parity=train_parity, cfg=cfg,
                winner_N=selected['N'] if selected else None, states=states,
                train_curve=curve,
                frozen_from='ALL (biased comparator)' if train_parity is None else f'event_id % 2 == {train_parity}',
                frozen_items=['TOP channel IDs and ranking', 'winning photon index N',
                              'per-N peak/amplitude/MAD seeds', 'per-N histogram axes/bin count',
                              'per-N fit windows', 'bootstrap RNG seed and fit configuration'],
                calibration='none: no walk, SPTR or electronics')


def evaluate(arr, n_generated, model, root_file=None):
    parity = 1-model['train_parity']
    denominator = len(generated_ids(n_generated, parity))
    curve = []
    for index in N_VALUES:
        values = timestamps(arr, model['channels'], parity, n_generated, index)
        row = fit(values, denominator, model['cfg'],
                  f'eval{parity}_learned{model["train_parity"]}_N{index}',
                  model['states'][index], root_file)
        row.update(N=index, sample='EVAL')
        curve.append(row)
    selected = next((r for r in curve if r['N'] == model['winner_N']), None)
    return dict(curve=curve, primary=selected, eval_parity=parity,
                selection_rule='Apply TRAIN winner; EVAL curve is diagnostic only; no EVAL argmin')
