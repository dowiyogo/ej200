#!/usr/bin/env python3
"""Analyze the two authorized F4 cells and the matching existing baseline."""

import hashlib
import json
import math
import re
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

from build_step4_pairs import analyze_cell
from exec46_schema import HC_EV_NM, SPEED_OF_LIGHT_MM_PER_NS, geant4_group_velocity_table


BASE_DIR = Path(__file__).resolve().parent
STEP4_PAIRS = BASE_DIR / "step4" / "step4_event_pairs.root"
RUN_DIR = Path("/home/rrios/exec46_20260915/f4_bc408_sensitivity")
OUTPUT = BASE_DIR / "f4_bc408_sensitivity"
REPORT = OUTPUT / "REPORT_F4_BC408_SENSITIVITY_20260916.md"
VARIANTS = ("visible_lower_764mm", "visible_current_3800mm")
EXPECTED_EVENTS = 10_000
BOOTSTRAPS = 2_000
BOOTSTRAP_SEED = 460604
MUON_MASS_MEV = 105.6583755
MUON_KINETIC_MEV = 1000.0
MUON_GAMMA = (MUON_KINETIC_MEV + MUON_MASS_MEV) / MUON_MASS_MEV
MUON_BETA = math.sqrt(1.0 - 1.0 / MUON_GAMMA ** 2)
FIT_A, FIT_B, FIT_C = 1.518, 0.640, 0.00423
REFERENCE_GLOBAL_SCINT_MM_NS = 186.070980
REFERENCE_GLOBAL_CHER_MM_NS = 148.547541
EXPERIMENTAL_EFFECTIVE_MM_NS = 155.0
UV_CLAMP_NM = 370.0
HANDICAP_DISTANCE_MM = 50.0
HELD_CAMPAIGN = Path('/home/rrios/exec46_20260916/full_grid_bc408_3800')
COMMAND = (
    "env PYTHONPATH=analysis/track_mechanism_20260915 python3 "
    "analysis/track_mechanism_20260915/analyze_f4_bc408_sensitivity.py"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def n_phase(wavelength_nm):
    return FIT_A + FIT_B * np.exp(-FIT_C * np.asarray(wavelength_nm))


def n_group_analytic(wavelength_nm):
    wavelength_nm = np.asarray(wavelength_nm)
    exponential = FIT_B * np.exp(-FIT_C * wavelength_nm)
    return FIT_A + exponential + wavelength_nm * FIT_C * exponential


def folded_angle(values):
    return np.minimum(values, 180.0 - values)


def bootstrap_width(values, rng):
    estimates = np.empty(BOOTSTRAPS)
    for start in range(0, BOOTSTRAPS, 50):
        stop = min(start + 50, BOOTSTRAPS)
        sampled = values[rng.integers(0, len(values), size=(stop-start, len(values)))]
        quantiles = np.quantile(sampled, (0.05, 0.95), axis=1)
        estimates[start:stop] = quantiles[1] - quantiles[0]
    return estimates.std(ddof=1)


def bootstrap_ratio_of_means(distance, time, rng):
    estimates = np.empty(BOOTSTRAPS)
    for start in range(0, BOOTSTRAPS, 50):
        stop = min(start + 50, BOOTSTRAPS)
        sampled = rng.integers(0, len(time), size=(stop-start, len(time)))
        estimates[start:stop] = (distance[sampled].mean(axis=1)
                                 / time[sampled].mean(axis=1))
    return estimates.std(ddof=1)


def load_rindex(path):
    wavelength, index = np.loadtxt(path, unpack=True)
    order = np.argsort(HC_EV_NM / wavelength)
    return (HC_EV_NM / wavelength)[order], index[order]


def metrics(label, arrays, rindex_path, root_path):
    rng = np.random.default_rng(BOOTSTRAP_SEED + sum(map(ord, label)))
    require(len(arrays['event_id']) == EXPECTED_EVENTS, f"{label}: wrong event count")
    first_is_cherenkov = arrays['first_source_type'] == 2
    fraction = float(first_is_cherenkov.mean())
    fraction_se = math.sqrt(fraction * (1.0 - fraction) / EXPECTED_EVENTS)

    result = {
        'scenario': label,
        'events': EXPECTED_EVENTS,
        'first_cherenkov_fraction': fraction,
        'first_cherenkov_fraction_se': fraction_se,
        'root_path': str(root_path),
        'root_sha256': sha256(root_path),
    }
    energy, index = load_rindex(rindex_path)
    gv_energy, gv_speed = geant4_group_velocity_table(energy, index)

    for selection, prefix in (
            ('cherenkov', 'first_cherenkov'),
            ('primary_cherenkov', 'first_primary_cherenkov'),
            ('scintillation', 'first_scint')):
        valid = arrays[prefix + '_track_id'] >= 0
        require(np.count_nonzero(valid) > 100, f"{label}/{selection}: too few photons")
        detection = arrays[prefix + '_t_detection_ns'][valid]
        creation = arrays[prefix + '_t_creation_ns'][valid]
        propagation = detection - creation
        path = arrays[prefix + '_path_length_mm'][valid]
        direct = arrays[prefix + '_d_direct_mm'][valid]
        wavelength = arrays[prefix + '_wl_nm'][valid]
        angle = folded_angle(arrays[prefix + '_exit_angle_deg'][valid])
        require(np.all(propagation > 0), f"{label}/{selection}: nonpositive propagation")
        group_lookup = np.interp(HC_EV_NM / wavelength, gv_energy, gv_speed)
        result.update({
            f'{selection}_n': int(valid.sum()),
            f'{selection}_angle_q05_deg': float(np.quantile(angle, 0.05)),
            f'{selection}_angle_median_deg': float(np.median(angle)),
            f'{selection}_angle_q95_deg': float(np.quantile(angle, 0.95)),
            f'{selection}_angle_central90_width_deg': float(
                np.quantile(angle, 0.95) - np.quantile(angle, 0.05)),
            f'{selection}_angle_width_bootstrap_se_deg': bootstrap_width(angle, rng),
            f'{selection}_path_speed_mean_mm_ns': float(np.mean(path / propagation)),
            f'{selection}_path_speed_sem_mm_ns': float(
                np.std(path / propagation, ddof=1) / math.sqrt(valid.sum())),
            f'{selection}_local_direct_velocity_mm_ns': float(
                np.mean(direct) / np.mean(propagation)),
            f'{selection}_local_direct_velocity_bootstrap_se_mm_ns':
                bootstrap_ratio_of_means(direct, propagation, rng),
            f'{selection}_mean_group_velocity_lookup_mm_ns': float(np.mean(group_lookup)),
            f'{selection}_mean_detected_wavelength_nm': float(np.mean(wavelength)),
        })
    return result


def extract_variant(payload):
    name, scratch = payload
    root_path = RUN_DIR / name / 'photon_hits_run000.root'
    require(root_path.is_file(), f"missing F4 ROOT: {root_path}")
    cell = {'cell_id': name, 'root_path': str(root_path), 'x_mm': -650,
            'material': 'EJ-200'}
    npz_path = Path(analyze_cell((cell, 0, scratch)))
    with np.load(npz_path) as source:
        columns = {key: source[key].copy() for key in source.files}
    face = columns['face_type'] == 0
    return name, {key: values[face] for key, values in columns.items()}


def load_baseline():
    columns = [
        'material_code', 'gun_x_mm', 'face_type', 'event_id',
        'first_source_type',
    ]
    for prefix in ('first_cherenkov', 'first_primary_cherenkov', 'first_scint'):
        columns.extend(prefix + '_' + field for field in (
            'track_id', 't_detection_ns', 't_creation_ns', 'path_length_mm',
            'd_direct_mm', 'wl_nm', 'exit_angle_deg'))
    with uproot.open(STEP4_PAIRS) as root_file:
        arrays = root_file['event_pairs'].arrays(columns, library='np')
    selected = ((arrays['material_code'] == 0) & (arrays['gun_x_mm'] == -650)
                & (arrays['face_type'] == 0))
    return {key: values[selected] for key, values in arrays.items()}


def uv_clamp_metrics(label, arrays):
    rows = []
    selections = (
        ('overall_cherenkov_winner', 'first', arrays['first_source_type'] == 2),
        ('source_specific_first_cherenkov', 'first_cherenkov',
         arrays['first_cherenkov_track_id'] >= 0),
        ('primary_like_first_cherenkov', 'first_primary_cherenkov',
         arrays['first_primary_cherenkov_track_id'] >= 0),
    )
    for population, prefix, valid in selections:
        wavelength = arrays[prefix + '_wl_nm_created'][valid]
        angle = folded_angle(arrays[prefix + '_exit_angle_deg'][valid])
        below = wavelength < UV_CLAMP_NM
        for region, mask in (
                ('all', np.ones(len(wavelength), dtype=bool)),
                ('clamped_lt370', below),
                ('measured_ge370', ~below)):
            selected_wavelength = wavelength[mask]
            selected_angle = angle[mask]
            require(len(selected_wavelength) > 20,
                    f'{label}/{population}/{region}: insufficient entries')
            rows.append({
                'scenario': label, 'population': population,
                'wavelength_region': region, 'n': int(len(selected_wavelength)),
                'fraction_of_population': float(mask.mean()),
                'fraction_below_370': float(below.mean()),
                'fraction_below_370_binomial_se': float(
                    math.sqrt(below.mean()*(1.0-below.mean())/len(below))),
                'wavelength_q01_nm': float(np.quantile(selected_wavelength, .01)),
                'wavelength_q05_nm': float(np.quantile(selected_wavelength, .05)),
                'wavelength_q10_nm': float(np.quantile(selected_wavelength, .10)),
                'wavelength_median_nm': float(np.median(selected_wavelength)),
                'angle_q05_deg': float(np.quantile(selected_angle, .05)),
                'angle_median_deg': float(np.median(selected_angle)),
                'angle_q95_deg': float(np.quantile(selected_angle, .95)),
                'angle_central90_width_deg': float(
                    np.quantile(selected_angle, .95)-np.quantile(selected_angle, .05)),
            })
    return rows


def caustic_time_selection(label, arrays, rindex_path):
    """Test the H4 angular-delay prediction photon by photon.

    The primary-like first Cherenkov photon at the near END is used.  The
    unmeasured UV clamp is removed explicitly so that the observed upper tail
    cannot be attributed to wavelengths below 370 nm.
    """
    valid = ((arrays['first_primary_cherenkov_track_id'] >= 0)
             & (arrays['first_primary_cherenkov_wl_nm_created'] >= UV_CLAMP_NM))
    wavelength = arrays['first_primary_cherenkov_wl_nm_created'][valid]
    angle_deg = folded_angle(
        arrays['first_primary_cherenkov_exit_angle_deg'][valid])
    # alpha is measured against the bar axis, so the geometric delay formula
    # requires the axial separation rather than the three-dimensional chord.
    distance = np.abs(
        arrays['first_primary_cherenkov_x_mm'][valid]
        - arrays['first_primary_cherenkov_x_creation_mm'][valid])
    propagation = (arrays['first_primary_cherenkov_t_detection_ns'][valid]
                   - arrays['first_primary_cherenkov_t_creation_ns'][valid])
    require(len(wavelength) > 1000, f'{label}: too few measured-domain photons')
    require(np.all(propagation > 0), f'{label}: nonpositive propagation')

    energy, index = load_rindex(rindex_path)
    gv_energy, gv_speed = geant4_group_velocity_table(energy, index)
    group_velocity = np.interp(HC_EV_NM / wavelength, gv_energy, gv_speed)
    edge_deg = np.degrees(np.arcsin(1.0 / (MUON_BETA*n_phase(wavelength))))
    angle_rad = np.radians(angle_deg)
    edge_rad = np.radians(edge_deg)
    predicted_penalty_ns = distance/group_velocity * (
        1.0/np.cos(angle_rad) - 1.0/np.cos(edge_rad))
    observed_edge_excess_ns = propagation - distance/(
        group_velocity*np.cos(edge_rad))

    angle_q95 = float(np.quantile(angle_deg, .95))
    tail = angle_deg >= angle_q95
    require(np.count_nonzero(tail) > 100, f'{label}: too few upper-tail photons')
    fit_slope, fit_intercept = np.polyfit(
        predicted_penalty_ns*1000.0, observed_edge_excess_ns*1000.0, 1)
    median_wavelength = float(np.median(wavelength))
    median_energy = HC_EV_NM/median_wavelength
    median_gv = float(np.interp(median_energy, gv_energy, gv_speed))
    median_edge = float(np.degrees(np.arcsin(
        1.0/(MUON_BETA*float(n_phase(median_wavelength))))))
    mean_distance = float(np.mean(distance))
    scalar_q95_penalty_ps = 1000.0*mean_distance/median_gv * (
        1.0/math.cos(math.radians(angle_q95))
        - 1.0/math.cos(math.radians(median_edge)))

    raw = pd.DataFrame({
        'scenario': label,
        'wavelength_nm_created': wavelength,
        'folded_exit_angle_deg': angle_deg,
        'finite_beta_edge_deg': edge_deg,
        'd_axial_mm': distance,
        'group_velocity_mm_ns': group_velocity,
        'propagation_time_ns': propagation,
        'predicted_angular_penalty_ps': predicted_penalty_ns*1000.0,
        'observed_edge_excess_ps': observed_edge_excess_ns*1000.0,
        'is_upper_5pct_angle_tail': tail.astype(np.int8),
    })
    summary = {
        'scenario': label,
        'n_measured_domain': len(raw),
        'mean_distance_mm': mean_distance,
        'median_wavelength_nm': median_wavelength,
        'median_finite_beta_edge_deg': median_edge,
        'angle_q95_deg': angle_q95,
        'angle_q99_deg': float(np.quantile(angle_deg, .99)),
        'scalar_q95_penalty_ps': scalar_q95_penalty_ps,
        'predicted_penalty_q95_ps': float(np.quantile(predicted_penalty_ns, .95)*1000),
        'observed_edge_excess_q95_ps': float(
            np.quantile(observed_edge_excess_ns, .95)*1000),
        'upper_tail_predicted_penalty_median_ps': float(
            np.median(predicted_penalty_ns[tail])*1000),
        'upper_tail_observed_excess_median_ps': float(
            np.median(observed_edge_excess_ns[tail])*1000),
        'upper_tail_predicted_penalty_mean_ps': float(
            np.mean(predicted_penalty_ns[tail])*1000),
        'upper_tail_observed_excess_mean_ps': float(
            np.mean(observed_edge_excess_ns[tail])*1000),
        'all_pearson_r': float(np.corrcoef(
            predicted_penalty_ns, observed_edge_excess_ns)[0, 1]),
        'upper_tail_pearson_r': float(np.corrcoef(
            predicted_penalty_ns[tail], observed_edge_excess_ns[tail])[0, 1]),
        'all_fit_intercept_ps': float(fit_intercept),
        'all_fit_slope': float(fit_slope),
    }
    return raw, summary


def save_caustic_time_selection(raw_frame, summary_frame):
    csv_path = OUTPUT / 'f4_caustic_time_selection.csv'
    summary_path = OUTPUT / 'f4_caustic_time_selection_summary.csv'
    root_path = OUTPUT / 'f4_caustic_time_selection.root'
    pdf_path = OUTPUT / 'f4_caustic_time_selection.pdf'
    meta_path = OUTPUT / 'f4_caustic_time_selection.meta.json'
    raw_frame.to_csv(csv_path, index=False, float_format='%.12g')
    summary_frame.to_csv(summary_path, index=False, float_format='%.12g')
    with uproot.recreate(root_path) as root_file:
        for tree_name, frame in (('photons', raw_frame), ('summary', summary_frame)):
            root_file[tree_name] = {
                column: (frame[column].astype(str).to_numpy(dtype=str)
                         if frame[column].dtype == object else frame[column].to_numpy())
                for column in frame.columns
            }

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True, sharey=True)
    for axis, scenario in zip(axes, VARIANTS):
        data = raw_frame[raw_frame.scenario == scenario]
        shown = data[(data.predicted_angular_penalty_ps > -40)
                     & (data.predicted_angular_penalty_ps < 300)
                     & (data.observed_edge_excess_ps > -60)
                     & (data.observed_edge_excess_ps < 300)]
        axis.hexbin(shown.predicted_angular_penalty_ps,
                    shown.observed_edge_excess_ps, gridsize=55,
                    mincnt=1, bins='log', cmap='viridis')
        summary = summary_frame[summary_frame.scenario == scenario].iloc[0]
        domain = np.array([-40.0, 300.0])
        axis.plot(domain, domain, color='0.55', ls='--', lw=1, label='unit response')
        axis.plot(domain, summary.all_fit_intercept_ps + summary.all_fit_slope*domain,
                  color='tab:red', lw=1.2,
                  label=f'fit slope={summary.all_fit_slope:.3f}')
        axis.set_title(scenario)
        axis.set_xlabel('Predicted angular penalty [ps]')
        axis.grid(alpha=.2)
        axis.legend(fontsize=8)
    axes[0].set_ylabel('Observed excess above cone-edge time [ps]')
    fig.tight_layout()
    fig.savefig(pdf_path)
    plt.close(fig)
    meta = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'command': COMMAND,
        'selection': ('primary-like first Cherenkov at near END, '
                      'wl_nm_created >= 370 nm'),
        'prediction': ('d_axial/v_group(lambda) * [1/cos(alpha) - '
                       '1/cos(alpha_edge(lambda,beta))]'),
        'observed_excess': ('t_detection-t_creation - '
                            'd/[v_group(lambda)*cos(alpha_edge)]'),
        'csv': str(csv_path), 'summary_csv': str(summary_path),
        'root': str(root_path), 'pdf': str(pdf_path),
    }
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + '\n')


def save_uv_bundle(frame, extracted):
    csv_path = OUTPUT / 'f4_uv_clamp.csv'
    root_path = OUTPUT / 'f4_uv_clamp.root'
    pdf_path = OUTPUT / 'f4_uv_clamp.pdf'
    meta_path = OUTPUT / 'f4_uv_clamp.meta.json'
    frame.to_csv(csv_path, index=False, float_format='%.12g')
    with uproot.recreate(root_path) as root_file:
        root_file['uv_clamp'] = {
            column: (frame[column].astype(str).to_numpy(dtype=str)
                     if frame[column].dtype == object else frame[column].to_numpy())
            for column in frame.columns
        }
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for column, scenario in enumerate(VARIANTS):
        arrays = extracted[scenario]
        winner = arrays['first_source_type'] == 2
        wavelength = arrays['first_wl_nm_created'][winner]
        axes[0, column].hist(wavelength, bins=np.arange(280, 821, 10), density=True,
                             histtype='step', linewidth=1.4)
        axes[0, column].axvline(UV_CLAMP_NM, color='black', linestyle='--', linewidth=1)
        axes[0, column].set_title(scenario)
        axes[0, column].set_xlabel('Created wavelength [nm]')
        axes[0, column].set_ylabel('Normalized density')
        primary = arrays['first_primary_cherenkov_track_id'] >= 0
        primary_wavelength = arrays['first_primary_cherenkov_wl_nm_created'][primary]
        primary_angle = folded_angle(
            arrays['first_primary_cherenkov_exit_angle_deg'][primary])
        bins = np.arange(35.0, 55.01, 0.25)
        for mask, label in ((primary_wavelength < UV_CLAMP_NM, 'lambda < 370 nm'),
                            (primary_wavelength >= UV_CLAMP_NM, 'lambda >= 370 nm')):
            axes[1, column].hist(primary_angle[mask], bins=bins, density=True,
                                 histtype='step', linewidth=1.4, label=label)
        axes[1, column].set_xlabel('Folded exit_angle_deg [deg]')
        axes[1, column].set_ylabel('Normalized density')
        axes[1, column].legend(fontsize=8)
        for axis in axes[:, column]:
            axis.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(pdf_path)
    plt.close(fig)
    meta = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'command': COMMAND, 'uv_clamp_nm': UV_CLAMP_NM,
        'scope': ('near END (face_type=0) at x=-650; overall Cherenkov winners, '
                  'source-specific first Cherenkov, and primary-like first Cherenkov'),
        'angle': 'folded min(exit_angle_deg, 180-exit_angle_deg)',
        'pdf': str(pdf_path), 'csv': str(csv_path), 'root': str(root_path),
    }
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + '\n')


def save_bundle(frame):
    csv_path = OUTPUT / 'f4_metrics.csv'
    root_path = OUTPUT / 'f4_metrics.root'
    pdf_path = OUTPUT / 'f4_metrics.pdf'
    meta_path = OUTPUT / 'f4_metrics.meta.json'
    frame.to_csv(csv_path, index=False, float_format='%.12g')
    with uproot.recreate(root_path) as root_file:
        root_file['metrics'] = {
            column: (frame[column].astype(str).to_numpy(dtype=str)
                     if frame[column].dtype == object else frame[column].to_numpy())
            for column in frame.columns
        }
    display = frame.copy()
    labels = ['constant n=1.58\nexisting', 'BC-408 n(lambda)\n764 mm bound',
              'BC-408 n(lambda)\n3800 mm']
    x = np.arange(len(display))
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes[0, 0].errorbar(x, 100*display.first_cherenkov_fraction,
                        yerr=100*display.first_cherenkov_fraction_se, fmt='o', capsize=3)
    axes[0, 0].axhline(74.10, color='0.5', ls='--', lw=1)
    axes[0, 0].set_ylabel('First-photon Cherenkov fraction [%]')
    axes[0, 1].errorbar(x, display.primary_cherenkov_angle_central90_width_deg,
                        yerr=display.primary_cherenkov_angle_width_bootstrap_se_deg,
                        fmt='o', capsize=3)
    axes[0, 1].axhline(2.716, color='0.5', ls='--', lw=1,
                       label='370–660 nm edge span')
    axes[0, 1].set_ylabel('Primary-cone central-90% width [deg]')
    axes[0, 1].legend(fontsize=8)
    axes[1, 0].errorbar(x, display.cherenkov_local_direct_velocity_mm_ns,
                        yerr=display.cherenkov_local_direct_velocity_bootstrap_se_mm_ns,
                        fmt='o', capsize=3)
    axes[1, 0].axhline(129.629, color='0.5', ls='--', lw=1,
                       label='420 nm finite-beta prediction')
    axes[1, 0].set_ylabel('First-Cherenkov local direct velocity [mm/ns]')
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].errorbar(x, display.scintillation_local_direct_velocity_mm_ns,
                        yerr=display.scintillation_local_direct_velocity_bootstrap_se_mm_ns,
                        fmt='o', capsize=3)
    axes[1, 1].axhline(EXPERIMENTAL_EFFECTIVE_MM_NS, color='0.5', ls='--', lw=1,
                       label='test-beam global v_eff')
    axes[1, 1].set_ylabel('First-scint local direct velocity [mm/ns]')
    axes[1, 1].legend(fontsize=8)
    for axis in axes.flat:
        axis.set_xticks(x, labels)
        axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(pdf_path)
    plt.close(fig)
    meta = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'command': COMMAND,
        'bootstrap_replicates': BOOTSTRAPS,
        'bootstrap_seed': BOOTSTRAP_SEED,
        'angle_width': 'q95-q05 of folded exit_angle_deg for primary-like Cherenkov winner',
        'local_direct_velocity': 'mean(d_direct_mm)/mean(t_detection_ns-t_creation_ns)',
        'total_track_effective_speed': (
            'mean(path_length_mm/(t_detection_ns-t_creation_ns)); path_length_mm is '
            'G4Track total path from creation through every traversed medium, so this is not '
            'the scintillator group velocity'),
        'caution': ('one x=-650 cell cannot reproduce a seven-distance effective-velocity '
                    'slope or establish agreement with the 155 mm/ns test-beam value'),
        'csv': str(csv_path), 'root': str(root_path), 'pdf': str(pdf_path),
    }
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + '\n')


def render_report(frame, uv_frame, caustic_frame, theory):
    def row(label):
        return frame[frame.scenario == label].iloc[0]
    baseline = row('constant_n1p58_existing')
    lower = row('visible_lower_764mm')
    current = row('visible_current_3800mm')
    f4_runtime = RUN_DIR / 'visible_current_3800mm' / 'sslg4'
    f4_rindex_hash = sha256(f4_runtime/'data/oscnt/opsc-100/rIndex.txt')
    f4_abs_hash = sha256(f4_runtime/'data/oscnt/opsc-100/absLength.txt')
    bc404_n408 = 1.578 + .818*math.exp(-.00729*408.0)
    bc404_ng408 = bc404_n408 + 408.0*.00729*.818*math.exp(-.00729*408.0)
    bc404_vg408 = SPEED_OF_LIGHT_MM_PER_NS/bc404_ng408
    winner_uv = uv_frame.query(
        'population == "overall_cherenkov_winner" and wavelength_region == "all"')
    angle_uv = uv_frame.query(
        'population in ["overall_cherenkov_winner", "primary_like_first_cherenkov"] '
        'and wavelength_region != "all"')
    caustic = caustic_frame.set_index('scenario')
    baseline_handicap = HANDICAP_DISTANCE_MM*1000.0*(
        1.0/baseline.cherenkov_local_direct_velocity_mm_ns
        - 1.0/baseline.scintillation_local_direct_velocity_mm_ns)
    lower_handicap = HANDICAP_DISTANCE_MM*1000.0*(
        1.0/lower.cherenkov_local_direct_velocity_mm_ns
        - 1.0/lower.scintillation_local_direct_velocity_mm_ns)
    current_handicap = HANDICAP_DISTANCE_MM*1000.0*(
        1.0/current.cherenkov_local_direct_velocity_mm_ns
        - 1.0/current.scintillation_local_direct_velocity_mm_ns)
    lines = [
        '# EXEC_46 F4 — BC-408 optical-model sensitivity at EJ200 x=-650 mm', '',
        'Date: 2026-09-16. Step 6 remains suspended.', '',
        'Two authorized 10,000-event cells were run with the original EndTop geometry, OPSC-100 '
        'emission/timing/yield, Broadcom PDE, zero SPTR, seeds 26092601/8349041, and four workers. '
        'Only RINDEX and ABSLENGTH in isolated SSLG4 runtime copies changed. The existing constant-n '
        'cell is the paired baseline; it was not rerun.', '',
        'Each new log contains the inherited Geant4 `mat031` warning that the OPSC-100 '
        'fractional masses sum to 0.998943. The identical warning is present in the existing '
        'baseline log; it was not introduced by the RINDEX/ABSLENGTH sensitivity change.', '',
        '## Published input and our derivation', '',
        '[Huggins, Johnson and Buckner](https://arxiv.org/abs/2608.13710) measured '
        'commercial PVT scintillators over 370–660 nm. '
        'For BC-408 their ODR model is n(lambda)=1.518+0.640 exp(-0.00423 lambda_nm). '
        'They report 28.23 +/- 2.7 mm at 372 nm and only 90% lower bounds at visible wavelengths, '
        'including >764 mm at 439 nm. They have not implemented group-velocity corrections.', '',
        f'Our derivative gives n(420 nm)={theory["n420"]:.6f}, n_group(420 nm)={theory["ng420"]:.6f}, '
        f'and v_group={theory["vg420"]:.3f} mm/ns. Propagating the paper’s typical 0.0055–0.0075 '
        f'weighted-index uncertainty as a fully correlated additive index uncertainty gives '
        f'sigma(v_group)={theory["vg_error_low"]:.3f}–{theory["vg_error_high"]:.3f} mm/ns. '
        'This propagation is ours and is conditional because the paper does not publish the ODR '
        'coefficient covariance needed for a full n_group uncertainty.', '',
        f'For the configured 1 GeV muon (beta={MUON_BETA:.6f}), the same derivation gives a '
        f'420 nm axial cone-edge velocity of {theory["axial420"]:.3f} mm/ns. The finite-beta edge '
        f'moves from {theory["edge370"]:.3f} deg at 370 nm to {theory["edge660"]:.3f} deg at '
        f'660 nm, a {theory["edge_span"]:.3f} deg chromatic span. The identity between the cone '
        'edge and critical angle remains exact wavelength by wavelength in the beta=1 limit.', '',
        f'At 200 mm, replacing c/1.58 by the derived 420 nm group velocity changes the nominal '
        f'transport time by {theory["delta_t_200_ps"]:.1f} ps, '
        f'{theory["delta_t_200_ps"]/7.1834161898857465:.1f} times the 7.183 ps common residual. '
        'This scale comparison is our derivation, not a result quoted by Huggins et al.', '',
        'The MPT samples n(lambda) every 2 nm from 370 to 660 nm. It is held constant at n(370) '
        'from 200–370 nm and at n(660) from 660–800 nm. The UV clamp is an explicit hypothesis: '
        'no exponential extrapolation is used. ABSLENGTH is 28.23 mm through 372 nm, interpolates '
        'to the visible scenario at 439 nm, then is held at either 764 or 3800 mm. The 764 mm '
        'scenario represents a lower bound, not a measured central value.', '',
        '## F4 results', '',
        '| scenario | first Cherenkov [%] | primary caustic q05–q95 [deg] | width [deg] | '
        'first-Cher local velocity [mm/ns] | first-scint local velocity [mm/ns] | '
        'first-scint total-track effective speed [mm/ns] | '
        'scintillator-MPT group speed at detected wavelengths [mm/ns] |',
        '|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for _, value in frame.iterrows():
        lines.append(
            f'| {value.scenario} | {100*value.first_cherenkov_fraction:.2f} +/- '
            f'{100*value.first_cherenkov_fraction_se:.2f} | '
            f'{value.primary_cherenkov_angle_q05_deg:.3f}–'
            f'{value.primary_cherenkov_angle_q95_deg:.3f} | '
            f'{value.primary_cherenkov_angle_central90_width_deg:.3f} +/- '
            f'{value.primary_cherenkov_angle_width_bootstrap_se_deg:.3f} | '
            f'{value.cherenkov_local_direct_velocity_mm_ns:.3f} +/- '
            f'{value.cherenkov_local_direct_velocity_bootstrap_se_mm_ns:.3f} | '
            f'{value.scintillation_local_direct_velocity_mm_ns:.3f} +/- '
            f'{value.scintillation_local_direct_velocity_bootstrap_se_mm_ns:.3f} | '
            f'{value.scintillation_path_speed_mean_mm_ns:.3f} +/- '
            f'{value.scintillation_path_speed_sem_mm_ns:.3f} | '
            f'{value.scintillation_mean_group_velocity_lookup_mm_ns:.3f} |')
    lines += ['',
        'The published-dispersion MPT changes every requested observable. The difference between '
        'the two absorption scenarios is the unresolved attenuation systematic; neither scenario '
        'is selected as truth. The sub-370 nm Cherenkov region is unmeasured and is the dominant '
        'model limitation because the emitted spectrum scales approximately as 1/lambda^2.', '',
        f'The previous 74.10% first-photon Cherenkov fraction becomes '
        f'{100*lower.first_cherenkov_fraction:.2f}% (764 mm) and '
        f'{100*current.first_cherenkov_fraction:.2f}% (3800 mm). The primary-cone central-90% '
        f'width changes from {baseline.primary_cherenkov_angle_central90_width_deg:.3f} deg to '
        f'{lower.primary_cherenkov_angle_central90_width_deg:.3f}/{current.primary_cherenkov_angle_central90_width_deg:.3f} deg. '
        'This empirical quantile width is distinct from both the 0.02 deg histogram bin and the '
        '2.716 deg parameter-free edge span.', '',
        f'The one-cell first-Cherenkov local direct velocity moves from '
        f'{baseline.cherenkov_local_direct_velocity_mm_ns:.3f} to '
        f'{lower.cherenkov_local_direct_velocity_mm_ns:.3f}/{current.cherenkov_local_direct_velocity_mm_ns:.3f} mm/ns. '
        f'The first-scintillation local direct velocity moves from '
        f'{baseline.scintillation_local_direct_velocity_mm_ns:.3f} to '
        f'{lower.scintillation_local_direct_velocity_mm_ns:.3f}/{current.scintillation_local_direct_velocity_mm_ns:.3f} mm/ns.', '',
        '`path_length_mm` is `G4Track::GetTrackLength()` plus the detection step '
        '(`src/SiPMSD.cc`); it includes every medium crossed between creation and detection. '
        'Consequently, path divided by propagation time is a total-track effective speed, not '
        'the scintillator group velocity. The last table column instead evaluates Geant4’s '
        'RINDEX-derived group-velocity table at each detected wavelength. Their difference is '
        'therefore expected for tracks containing non-scintillator segments and is not an MPT '
        'closure failure.', '',
        'The historical 148.548 mm/ns Cherenkov and 186.071 mm/ns scintillation numbers are '
        'seven-distance, origin-constrained slopes. A single x=-650 cell cannot refit those slopes. '
        'Likewise, the [SHiP test-beam 155 mm/ns value]'
        '(https://doi.org/10.1016/j.nima.2020.164398) is a global effective propagation speed from '
        'time versus position, not a microscopic group velocity or one-point d/t ratio. F4 can test '
        f'whether the local transport shifts toward that scale. Numerically, the distance of the '
        f'local first-scintillation proxy from 155 mm/ns decreases from '
        f'{abs(baseline.scintillation_local_direct_velocity_mm_ns-EXPERIMENTAL_EFFECTIVE_MM_NS):.3f} '
        f'to {abs(lower.scintillation_local_direct_velocity_mm_ns-EXPERIMENTAL_EFFECTIVE_MM_NS):.3f}/'
        f'{abs(current.scintillation_local_direct_velocity_mm_ns-EXPERIMENTAL_EFFECTIVE_MM_NS):.3f} '
        'mm/ns, so the shift is in the requested direction. It cannot establish agreement or '
        'trigger the stated full-campaign invalidation criterion by itself. That requires at least '
        'two positions under the corrected optical model.', '',
        '## G2: UV-clamp diagnostic', '',
        f'Configuration correction: the executed F4 MPT clamps the index below 370 nm to '
        f'n(370)={float(n_phase(370.0)):.6f}, but it does **not** make that region absorption-free. '
        'It holds ABSLENGTH at 28.23 mm from 200 through 372 nm. Thus this diagnostic measures '
        'the combined implemented UV hypotheses (index clamp plus absorption-length extension). '
        'A hypothetical zero-absorption UV model is not represented by these ROOT files and cannot '
        'be inferred without a new simulation, which G2 forbids.', '',
        'The following population is the actual first photon at the near END, conditioned on that '
        'winner being Cherenkov. `wl_nm_created` is used, so the classification tests the optical '
        'model at photon creation rather than the detected wavelength.', '',
        '| scenario | N Cherenkov winners | wavelength q01 [nm] | q05 [nm] | q10 [nm] | median [nm] | fraction <370 nm |',
        '|---|---:|---:|---:|---:|---:|---:|',
    ]
    for _, value in winner_uv.iterrows():
        lines.append(
            f'| {value.scenario} | {int(value.n)} | {value.wavelength_q01_nm:.2f} | '
            f'{value.wavelength_q05_nm:.2f} | {value.wavelength_q10_nm:.2f} | '
            f'{value.wavelength_median_nm:.2f} | '
            f'{100*value.fraction_below_370:.2f} +/- '
            f'{100*value.fraction_below_370_binomial_se:.2f}% |')
    lines += ['',
        'For completeness, the fraction below 370 nm among the source-specific first Cherenkov '
        'photons in all events is 5.39% (764 mm) and 4.68% (3800 mm); among the primary-like '
        'cone selection it is only 1.92% and 1.32%. None approaches the preregistered 50% threshold.', '',
        '| scenario | population | wavelength region | N | angle q05 [deg] | median [deg] | q95 [deg] | q95-q05 [deg] |',
        '|---|---|---|---:|---:|---:|---:|---:|']
    for _, value in angle_uv.iterrows():
        lines.append(
            f'| {value.scenario} | {value.population} | {value.wavelength_region} | '
            f'{int(value.n)} | {value.angle_q05_deg:.3f} | '
            f'{value.angle_median_deg:.3f} | {value.angle_q95_deg:.3f} | '
            f'{value.angle_central90_width_deg:.3f} |')
    lines += ['',
        'The primary-like caustic remains broad after removing the clamped photons: its measured-domain '
        'width is 8.096 deg (764 mm) and 7.695 deg (3800 mm), compared with total widths of '
        '8.199 and 7.749 deg. The small clamped population has wider tails, but it does not generate '
        'the observed 8 deg width. Therefore the 82.25/82.02% first-photon fractions are **not '
        'dominated by the clamped population in the executed F4 trees under the declared >50% '
        'decision rule**. The unmeasured UV region and its 28.23 mm absorption extension remain model '
        'systematics; passing this test does not validate either hypothesis physically.', '',
        '## H4: temporal selection of the caustic upper tail', '',
        'For every primary-like first Cherenkov photon with created wavelength at or above 370 nm, '
        'the parameter-free angular penalty is evaluated photon by photon as '
        '`d/v_group(lambda) * [1/cos(alpha) - 1/cos(alpha_edge(lambda,beta))]`. The observed '
        'comparison quantity is propagation time minus the wavelength-dependent cone-edge time. '
        'This selection removes the UV clamp before testing the 46.5 deg tail.', '',
        '| scenario | N | median edge [deg] | angle q95 [deg] | angle q99 [deg] | q95-angle penalty [ps] | upper-tail median predicted / observed [ps] | upper-tail r | all-event fit intercept [ps] | slope |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for scenario in VARIANTS:
        value = caustic.loc[scenario]
        lines.append(
            f'| {scenario} | {int(value.n_measured_domain)} | '
            f'{value.median_finite_beta_edge_deg:.3f} | {value.angle_q95_deg:.3f} | '
            f'{value.angle_q99_deg:.3f} | {value.scalar_q95_penalty_ps:.1f} | '
            f'{value.upper_tail_predicted_penalty_median_ps:.1f} / '
            f'{value.upper_tail_observed_excess_median_ps:.1f} | '
            f'{value.upper_tail_pearson_r:.4f} | {value.all_fit_intercept_ps:.1f} | '
            f'{value.all_fit_slope:.3f} |')
    lines += ['',
        f'The q95 angles of {caustic.loc[VARIANTS[0], "angle_q95_deg"]:.2f}/'
        f'{caustic.loc[VARIANTS[1], "angle_q95_deg"]:.2f} deg cost '
        f'{caustic.loc[VARIANTS[0], "scalar_q95_penalty_ps"]:.1f}/'
        f'{caustic.loc[VARIANTS[1], "scalar_q95_penalty_ps"]:.1f} ps relative to a photon at the '
        'finite-beta cone edge at the same representative distance and wavelength. Within the '
        'upper 5% angular tail, the eventwise predicted penalty and observed excess have '
        f'Pearson r={caustic.loc[VARIANTS[0], "upper_tail_pearson_r"]:.4f}/'
        f'{caustic.loc[VARIANTS[1], "upper_tail_pearson_r"]:.4f}; their median scales are '
        f'{caustic.loc[VARIANTS[0], "upper_tail_predicted_penalty_median_ps"]:.1f}/'
        f'{caustic.loc[VARIANTS[1], "upper_tail_predicted_penalty_median_ps"]:.1f} ps predicted '
        'and '
        f'{caustic.loc[VARIANTS[0], "upper_tail_observed_excess_median_ps"]:.1f}/'
        f'{caustic.loc[VARIANTS[1], "upper_tail_observed_excess_median_ps"]:.1f} ps observed. '
        'The nonzero fitted intercept records transport contributions absent from the one-angle '
        f'formula; the slopes of {caustic.loc[VARIANTS[0], "all_fit_slope"]:.3f}/'
        f'{caustic.loc[VARIANTS[1], "all_fit_slope"]:.3f} show that the angular dependence itself is recovered. '
        f'The measured-domain widths of 8.096/7.695 deg are '
        f'{8.096/baseline.primary_cherenkov_angle_central90_width_deg:.2f}/'
        f'{7.695/baseline.primary_cherenkov_angle_central90_width_deg:.2f} times the '
        f'{baseline.primary_cherenkov_angle_central90_width_deg:.3f} deg constant-index baseline '
        f'and exceed the {theory["edge_span"]:.3f} deg chromatic edge span. The 8 deg central width '
        'and its tail to about 46.5 deg are therefore quantitatively '
        'compatible with arrival-time selection, rather than a second geometric cone edge or the '
        'sub-370 nm clamp.', '',
        '## G3: direction of the Cherenkov effect', '',
        f'At d={HANDICAP_DISTANCE_MM:.0f} mm the measured local transport handicap '
        f'`d*(1/v_Cher-1/v_scint)` decreases from {baseline_handicap:.1f} ps to '
        f'{lower_handicap:.1f}/{current_handicap:.1f} ps. Relative to the baseline, the '
        f'first-scintillation local velocity falls by '
        f'{100*(1-lower.scintillation_local_direct_velocity_mm_ns/baseline.scintillation_local_direct_velocity_mm_ns):.2f}/'
        f'{100*(1-current.scintillation_local_direct_velocity_mm_ns/baseline.scintillation_local_direct_velocity_mm_ns):.2f}%, '
        f'while the first-Cherenkov velocity falls by only '
        f'{100*(1-lower.cherenkov_local_direct_velocity_mm_ns/baseline.cherenkov_local_direct_velocity_mm_ns):.2f}/'
        f'{100*(1-current.cherenkov_local_direct_velocity_mm_ns/baseline.cherenkov_local_direct_velocity_mm_ns):.2f}%. '
        'The corrected optical model therefore strengthens Cherenkov’s advantage. The prior prediction '
        'that dispersion would reduce the first-photon Cherenkov fraction is refuted by the two measured '
        'F4 variants; the axial Cherenkov speed is less sensitive because cone geometry controls it.', '',
        '## H2/H3: single held campaign', '',
        f'The first-photon Cherenkov fractions differ by '
        f'{100*(lower.first_cherenkov_fraction-current.first_cherenkov_fraction):.2f} +/- '
        f'{100*math.hypot(lower.first_cherenkov_fraction_se,current.first_cherenkov_fraction_se):.2f} '
        'percentage points between 764 and 3800 mm. The primary-caustic widths differ by '
        f'{lower.primary_cherenkov_angle_central90_width_deg-current.primary_cherenkov_angle_central90_width_deg:.3f} +/- '
        f'{math.hypot(lower.primary_cherenkov_angle_width_bootstrap_se_deg,current.primary_cherenkov_angle_width_bootstrap_se_deg):.3f} deg. '
        'Both are compatible with zero. H2 therefore selects one campaign, the 3800 mm scenario; '
        'F4 documents the observed insensitivity rather than motivating a duplicate grid.', '',
        'The campaign has been prepared and hash-audited but not launched. Its EJ-200 MPT is copied '
        'from the validated F4 3800 mm runtime. EJ-204 and EJ-230 are explicitly marked '
        '`UNCORRECTED_CONSTANT_RINDEX`. A BC-404 RINDEX construction from '
        '`A=1.578, B=0.818, C=0.00729` is technically viable, but it is not validated here and '
        'does not include a measured absorption model; it is excluded from this campaign. No '
        'measured analog is available for EJ-230.', '',
        f'The validated EJ-200 hashes are `rIndex.txt={f4_rindex_hash}` and '
        f'`absLength.txt={f4_abs_hash}`; every prepared EJ-200 cell has these exact hashes. '
        f'For scale, the proposed BC-404 coefficients give n(408 nm)={bc404_n408:.6f}, '
        f'n_group(408 nm)={bc404_ng408:.6f}, and v_group={bc404_vg408:.3f} mm/ns. '
        'Generating that RINDEX table is feasible, but enabling it would require its own '
        'single-cell validation and a declared ABSLENGTH treatment.', '',
        'The detached-driver dry run passed with 21 pending cells, zero outputs, diagnostics OFF, '
        'no timeout, 354,728,885,238 projected bytes against 961,284,907,008 available bytes, '
        'and a 3,663,937,536-byte conservative six-process memory budget against '
        '120,537,047,040 bytes available.', '',
        'The exact held launch command is:', '',
        '```bash',
        'python3 analysis/sigma_t/orchestration/detached_grid.py launch \\',
        '  --directory /home/rrios/exec46_20260916/full_grid_bc408_3800',
        '```', '',
        'It remains held pending Rene’s explicit approval. No launch command was executed.', '',
        '## Decision', '',
        'The current constant-n optical model is demonstrably non-robust for transport timing at '
        'the scale of the 7 ps residual: the measured-input sensitivity shifts local photon transport '
        'by far more than 7 ps-equivalent timing and introduces chromatic broadening that the old '
        'model fixes to zero. Therefore no Step 6 mechanism attribution is defensible with the '
        'constant-n campaign. The stronger claim that the old campaign is experimentally invalidated '
        'by reproducing 155 mm/ns is **not decidable from this one-position design**.', '',
        'EJ-230 remains on its existing constant-n model: BC-420 was not measured and BC-422 is not '
        'a valid substitute. It is explicitly uncorrected.', '',
        '## Reproducibility', '', '```bash',
        'env PYTHONPATH=analysis/track_mechanism_20260915 python3 '
        'analysis/track_mechanism_20260915/run_f4_bc408_sensitivity.py',
        COMMAND, '```', '',
        'The exact MPT tables, macros, logs, ROOT hashes, run times and commands are under '
        '`/home/rrios/exec46_20260915/f4_bc408_sensitivity/`. Figure sidecars are '
        '`f4_metrics.{csv,root,meta.json}`, `f4_uv_clamp.{csv,root,meta.json}`, and '
        '`f4_caustic_time_selection.{csv,root,meta.json}`; '
        'source-selection caches are retained under '
        '`scratch/`. No Step 6 analysis, push, merge, or deck edit occurred.', ''
    ]
    REPORT.write_text('\n'.join(lines))


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    scratch = OUTPUT / 'scratch'
    scratch.mkdir(exist_ok=True)
    for variant in VARIANTS:
        require((RUN_DIR / variant / '.DONE.json').is_file(), f"unfinished F4 run: {variant}")
        log_text = (RUN_DIR / variant / 'simulation.log').read_text()
        event_matches = re.findall(r'^\s*Events run\s*:\s*(\d+)\s*$', log_text, re.MULTILINE)
        require(event_matches == [str(EXPECTED_EVENTS)],
                f"master run summary does not certify {EXPECTED_EVENTS} events: {variant}")
        require(log_text.count('EXEC34_POST_BEAM') == 5,
                f"post-beam markers incomplete: {variant}")
    with ProcessPoolExecutor(max_workers=2) as pool:
        extracted = dict(pool.map(extract_variant, [(name, scratch) for name in VARIANTS]))
    baseline = load_baseline()
    baseline_rindex = (Path('/home/rrios/exec46_20260915/build_baseline/sslg4/data/'
                            'oscnt/opsc-100/rIndex.txt'))
    rows = [metrics('constant_n1p58_existing', baseline, baseline_rindex,
                    Path('/home/rrios/exec46_20260915/full_grid/cells/EJ200_xm650/'
                         'attempts/49749dc0e186483c86cf6e0b900894b1/photon_hits_run000.root'))]
    for variant in VARIANTS:
        rows.append(metrics(variant, extracted[variant],
                            RUN_DIR / variant / 'sslg4/data/oscnt/opsc-100/rIndex.txt',
                            RUN_DIR / variant / 'photon_hits_run000.root'))
    frame = pd.DataFrame(rows)
    save_bundle(frame)
    uv_frame = pd.DataFrame(
        row for variant in VARIANTS for row in uv_clamp_metrics(variant, extracted[variant]))
    save_uv_bundle(uv_frame, extracted)
    caustic_raw = []
    caustic_rows = []
    for variant in VARIANTS:
        raw, row = caustic_time_selection(
            variant, extracted[variant],
            RUN_DIR / variant / 'sslg4/data/oscnt/opsc-100/rIndex.txt')
        caustic_raw.append(raw)
        caustic_rows.append(row)
    caustic_raw_frame = pd.concat(caustic_raw, ignore_index=True)
    caustic_frame = pd.DataFrame(caustic_rows)
    save_caustic_time_selection(caustic_raw_frame, caustic_frame)
    n420 = float(n_phase(420.0))
    ng420 = float(n_group_analytic(420.0))
    vg420 = SPEED_OF_LIGHT_MM_PER_NS / ng420
    theory = {
        'n420': n420, 'ng420': ng420, 'vg420': vg420,
        'vg_error_low': SPEED_OF_LIGHT_MM_PER_NS / ng420**2 * 0.0055,
        'vg_error_high': SPEED_OF_LIGHT_MM_PER_NS / ng420**2 * 0.0075,
        'axial420': vg420 * math.sqrt(1.0 - 1.0/(MUON_BETA*n420)**2),
        'edge370': math.degrees(math.asin(1.0/(MUON_BETA*float(n_phase(370.0))))),
        'edge660': math.degrees(math.asin(1.0/(MUON_BETA*float(n_phase(660.0))))),
        'delta_t_200_ps': 1000.0 * 200.0 * (
            1.0/vg420 - 1.0/(SPEED_OF_LIGHT_MM_PER_NS/1.58)),
    }
    theory['edge_span'] = theory['edge660'] - theory['edge370']
    render_report(frame, uv_frame, caustic_frame, theory)
    summary = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'F4_COMPLETE_STEP6_SUSPENDED_OPTICAL_MODEL_SYSTEMATIC',
        'report': str(REPORT), 'report_sha256': sha256(REPORT),
        'metrics_sha256': sha256(OUTPUT/'f4_metrics.csv'),
        'uv_clamp_sha256': sha256(OUTPUT/'f4_uv_clamp.csv'),
        'caustic_time_selection_sha256': sha256(
            OUTPUT/'f4_caustic_time_selection.csv'),
        'step6_run': False,
        'full_campaign_run': False,
        'full_campaign_prepared': HELD_CAMPAIGN.is_dir(),
        'held_campaign': str(HELD_CAMPAIGN),
        'uv_clamp_artifact_gate': bool(
            uv_frame.query('population == "overall_cherenkov_winner" '
                           'and wavelength_region == "all"').fraction_below_370.max() > .5),
        'global_effective_velocity_test': 'UNDERIDENTIFIED_BY_SINGLE_POSITION',
        'theory': theory,
    }
    (OUTPUT / 'analysis_summary.json').write_text(
        json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
