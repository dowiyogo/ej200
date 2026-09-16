#!/usr/bin/env python3
"""Analyze the authorized I2 BC-404 validation after its ROOT is complete."""

import hashlib
import json
import math
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import uproot

from build_step4_pairs import analyze_cell
from exec46_schema import HC_EV_NM, geant4_group_velocity_table


BASE = Path(__file__).resolve().parent
RUN = Path('/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650')
OUTPUT = BASE/'i2_bc404_validation'
BASELINE_PAIRS = BASE/'step4/step4_event_pairs.root'
EVENTS = 10_000
BOOTSTRAPS = 2_000
SEED = 460612
COMMAND = ('env PYTHONPATH=analysis/track_mechanism_20260915 python3 '
           'analysis/track_mechanism_20260915/analyze_i2_bc404_validation.py')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def folded(values):
    return np.minimum(values, 180.0-values)


def bootstrap_width(values, rng):
    estimates = np.empty(BOOTSTRAPS)
    for start in range(0, BOOTSTRAPS, 50):
        stop = min(start+50, BOOTSTRAPS)
        sample = values[rng.integers(0, len(values), size=(stop-start, len(values)))]
        q = np.quantile(sample, (.05, .95), axis=1)
        estimates[start:stop] = q[1]-q[0]
    return float(estimates.std(ddof=1))


def bootstrap_ratio(distance, time, rng):
    estimates = np.empty(BOOTSTRAPS)
    for start in range(0, BOOTSTRAPS, 50):
        stop = min(start+50, BOOTSTRAPS)
        sample = rng.integers(0, len(time), size=(stop-start, len(time)))
        estimates[start:stop] = distance[sample].mean(axis=1)/time[sample].mean(axis=1)
    return float(estimates.std(ddof=1))


def load_new():
    root_path = RUN/'photon_hits_run000.root'
    require(root_path.is_file() and (RUN/'.DONE.json').is_file(), 'I2 output incomplete')
    scratch = OUTPUT/'scratch'
    scratch.mkdir(parents=True, exist_ok=True)
    cell = {'cell_id': 'I2_EJ204_xm650', 'root_path': str(root_path),
            'x_mm': -650, 'material': 'EJ-204'}
    path = Path(analyze_cell((cell, 0, scratch)))
    with np.load(path) as source:
        arrays = {key: source[key].copy() for key in source.files}
    mask = arrays['face_type'] == 0
    return {key: value[mask] for key, value in arrays.items()}


def load_baseline():
    columns = ['material_code', 'gun_x_mm', 'face_type', 'event_id',
               'first_source_type']
    for prefix in ('first_cherenkov', 'first_primary_cherenkov', 'first_scint'):
        columns.extend(prefix+'_'+field for field in (
            'track_id', 't_detection_ns', 't_creation_ns', 'd_direct_mm',
            'wl_nm', 'exit_angle_deg'))
    with uproot.open(BASELINE_PAIRS) as root_file:
        arrays = root_file['event_pairs'].arrays(columns, library='np')
    mask = ((arrays['material_code'] == 1) & (arrays['gun_x_mm'] == -650)
            & (arrays['face_type'] == 0))
    require(np.count_nonzero(mask) == EVENTS, 'baseline EJ204_xm650 count mismatch')
    return {key: value[mask] for key, value in arrays.items()}


def load_rindex(path):
    wavelength, index = np.loadtxt(path, unpack=True)
    order = np.argsort(HC_EV_NM/wavelength)
    return (HC_EV_NM/wavelength)[order], index[order]


def metrics(scenario, arrays, rindex_path):
    rng = np.random.default_rng(SEED+sum(map(ord, scenario)))
    first_cherenkov = arrays['first_source_type'] == 2
    fraction = float(first_cherenkov.mean())
    result = {
        'scenario': scenario, 'events': EVENTS,
        'first_cherenkov_fraction': fraction,
        'first_cherenkov_fraction_se': math.sqrt(fraction*(1-fraction)/EVENTS),
    }
    energy, index = load_rindex(rindex_path)
    gv_energy, gv_speed = geant4_group_velocity_table(energy, index)
    result['mpt_group_velocity_408nm_mm_ns'] = float(
        np.interp(HC_EV_NM/408.0, gv_energy, gv_speed))
    for population, prefix in (
            ('cherenkov', 'first_cherenkov'),
            ('primary_cherenkov', 'first_primary_cherenkov'),
            ('scintillation', 'first_scint')):
        valid = arrays[prefix+'_track_id'] >= 0
        require(np.count_nonzero(valid) > 100, f'{scenario}/{population}: too few photons')
        propagation = (arrays[prefix+'_t_detection_ns'][valid]
                       - arrays[prefix+'_t_creation_ns'][valid])
        distance = arrays[prefix+'_d_direct_mm'][valid]
        wavelength = arrays[prefix+'_wl_nm'][valid]
        angle = folded(arrays[prefix+'_exit_angle_deg'][valid])
        result.update({
            f'{population}_n': int(valid.sum()),
            f'{population}_angle_q05_deg': float(np.quantile(angle, .05)),
            f'{population}_angle_q95_deg': float(np.quantile(angle, .95)),
            f'{population}_angle_central90_width_deg': float(
                np.quantile(angle, .95)-np.quantile(angle, .05)),
            f'{population}_angle_width_bootstrap_se_deg': bootstrap_width(angle, rng),
            f'{population}_local_direct_velocity_mm_ns': float(
                distance.mean()/propagation.mean()),
            f'{population}_local_direct_velocity_bootstrap_se_mm_ns': bootstrap_ratio(
                distance, propagation, rng),
            f'{population}_mean_detected_group_velocity_mm_ns': float(
                np.mean(np.interp(HC_EV_NM/wavelength, gv_energy, gv_speed))),
        })
    return result


def write_root(frame, path):
    with uproot.recreate(path) as root_file:
        root_file['metrics'] = {
            column: (frame[column].astype(str).to_numpy(dtype=str)
                     if frame[column].dtype == object else frame[column].to_numpy())
            for column in frame.columns
        }


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    done = json.loads((RUN/'.DONE.json').read_text())
    require(done['events_with_hits'] == EVENTS and done['exit_code'] == 0,
            'I2 DONE record invalid')
    new = load_new()
    baseline = load_baseline()
    rows = [
        metrics('constant_n_baseline', baseline,
                Path('/home/rrios/exec46_20260915/build_baseline/sslg4/data/oscnt/opsc-101/rIndex.txt')),
        metrics('bc404_3800', new, RUN/'sslg4/data/oscnt/opsc-101/rIndex.txt'),
    ]
    frame = pd.DataFrame(rows)
    baseline_row, corrected = frame.iloc[0], frame.iloc[1]
    difference = corrected.first_cherenkov_fraction-baseline_row.first_cherenkov_fraction
    difference_se = math.hypot(corrected.first_cherenkov_fraction_se,
                               baseline_row.first_cherenkov_fraction_se)
    passed = bool(difference > 0)
    frame.to_csv(OUTPUT/'i2_metrics.csv', index=False, float_format='%.12g')
    write_root(frame, OUTPUT/'i2_metrics.root')
    metadata = {
        'created_utc': datetime.now(timezone.utc).isoformat(), 'command': COMMAND,
        'input_root': done['root_path'], 'input_root_sha256': done['root_sha256'],
        'baseline_pairs': str(BASELINE_PAIRS),
        'baseline_pairs_sha256': sha256(BASELINE_PAIRS),
        'bootstrap_replicates': BOOTSTRAPS, 'bootstrap_seed': SEED,
        'caustic_width': 'q95-q05 folded exit_angle_deg, primary-like first Cherenkov',
        'local_velocity': 'mean(d_direct_mm)/mean(t_detection_ns-t_creation_ns)',
        'csv': str(OUTPUT/'i2_metrics.csv'), 'root': str(OUTPUT/'i2_metrics.root'),
    }
    (OUTPUT/'i2_metrics.meta.json').write_text(
        json.dumps(metadata, indent=2, sort_keys=True)+'\n')
    report = [
        '# EXEC_46 I2 — BC-404 validation at EJ204_xm650', '',
        'Date: 2026-09-16. Existing ROOT analysis plus one authorized 10,000-event validation cell.', '',
        'The log contains the inherited Geant4 `mat031` warning that the OPSC-101 fractional '
        'masses do not sum exactly to one. It is also present with the baseline material data and '
        'is not introduced by the RINDEX/ABSLENGTH replacement.', '',
        '| scenario | first Cherenkov [%] | primary caustic width [deg] | first-Cher local velocity [mm/ns] | first-scint local velocity [mm/ns] | MPT v_group(408 nm) [mm/ns] |',
        '|---|---:|---:|---:|---:|---:|',
    ]
    for _, value in frame.iterrows():
        report.append(
            f'| {value.scenario} | {100*value.first_cherenkov_fraction:.2f} +/- '
            f'{100*value.first_cherenkov_fraction_se:.2f} | '
            f'{value.primary_cherenkov_angle_central90_width_deg:.3f} +/- '
            f'{value.primary_cherenkov_angle_width_bootstrap_se_deg:.3f} | '
            f'{value.cherenkov_local_direct_velocity_mm_ns:.3f} +/- '
            f'{value.cherenkov_local_direct_velocity_bootstrap_se_mm_ns:.3f} | '
            f'{value.scintillation_local_direct_velocity_mm_ns:.3f} +/- '
            f'{value.scintillation_local_direct_velocity_bootstrap_se_mm_ns:.3f} | '
            f'{value.mpt_group_velocity_408nm_mm_ns:.3f} |')
    report += ['',
        f'The corrected-minus-baseline first-photon Cherenkov difference is '
        f'{100*difference:.2f} +/- {100*difference_se:.2f} percentage points. '
        f'The preregistered direction gate is **{"PASS" if passed else "FAIL"}**: '
        f'{100*corrected.first_cherenkov_fraction:.2f}% is '
        f'{"above" if passed else "below"} the measured constant-index baseline of '
        f'{100*baseline_row.first_cherenkov_fraction:.2f}% (the registered rounded value is 62.90%).', '',
        'The production campaign may proceed only when this gate is PASS and the four MPT hashes '
        'are independently verified during campaign preparation.', '',
        'Reproduce with:', '', '```bash', COMMAND, '```', ''
    ]
    report_path = OUTPUT/'REPORT_I2_BC404_VALIDATION_20260916.md'
    report_path.write_text('\n'.join(report))
    summary = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'PASS' if passed else 'FAIL',
        'first_cherenkov_difference_percentage_points': 100*difference,
        'first_cherenkov_difference_se_percentage_points': 100*difference_se,
        'report': str(report_path), 'report_sha256': sha256(report_path),
        'campaign_launch_authorized_by_I2': passed,
    }
    (OUTPUT/'analysis_summary.json').write_text(
        json.dumps(summary, indent=2, sort_keys=True)+'\n')
    print(json.dumps(summary, indent=2, sort_keys=True))
    require(passed, 'I2 direction gate failed: campaign must not launch')


if __name__ == '__main__':
    main()
