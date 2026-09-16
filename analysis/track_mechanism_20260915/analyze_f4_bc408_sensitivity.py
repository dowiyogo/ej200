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


def render_report(frame, theory):
    def row(label):
        return frame[frame.scenario == label].iloc[0]
    baseline = row('constant_n1p58_existing')
    lower = row('visible_lower_764mm')
    current = row('visible_current_3800mm')
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
        '`f4_metrics.{csv,root,meta.json}`; source-selection caches are retained under '
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
    render_report(frame, theory)
    summary = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'F4_COMPLETE_STEP6_SUSPENDED_OPTICAL_MODEL_SYSTEMATIC',
        'report': str(REPORT), 'report_sha256': sha256(REPORT),
        'metrics_sha256': sha256(OUTPUT/'f4_metrics.csv'),
        'step6_run': False,
        'global_effective_velocity_test': 'UNDERIDENTIFIED_BY_SINGLE_POSITION',
        'theory': theory,
    }
    (OUTPUT / 'analysis_summary.json').write_text(
        json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
