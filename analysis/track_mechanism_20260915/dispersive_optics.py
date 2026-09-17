"""Cell-resolved optical predictions; no fitted ODR coefficients or scalar index.

RINDEX is interpolated in photon energy. GROUPVEL uses the existing numerical
Geant4 11.4 CalculateGROUPVEL transcription, including its midpoint mesh and
normal-dispersion guard. c/GROUPVEL is the effective numerical group index.
"""
import hashlib
import json
import os
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
from exec46_schema import HC_EV_NM, SPEED_OF_LIGHT_MM_PER_NS as C, geant4_group_velocity_table

CAMPAIGN = Path(os.environ.get('EXEC46_CAMPAIGN_DIR',
    '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2'))
MATERIALS = ('EJ-200', 'EJ-204', 'EJ-230')
MEASURED_MIN_NM, MEASURED_MAX_NM = 370., 660.
MUON_MASS_MEV = 105.6583755
QUANTILES = (.05, .5, .95)


class OpticalTable:
    def __init__(self, path, beta):
        self.path = Path(path).resolve(strict=True)
        table = np.loadtxt(self.path)
        energy = HC_EV_NM / table[:, 0]
        order = np.argsort(energy)
        self.energy, self.index = energy[order], table[order, 1]
        if (not np.all(np.isfinite(table)) or np.any(table[:, 0] <= 0)
                or np.any(np.diff(self.energy) <= 0) or np.any(self.index <= 1)):
            raise ValueError(f'Invalid RINDEX: {path}')
        self.beta = beta
        self.group_energy, self.group_speed = geant4_group_velocity_table(self.energy, self.index)
        self.constant = bool(np.all(self.index == self.index[0]))
        self.sha256 = hashlib.sha256(self.path.read_bytes()).hexdigest()
        self.evaluate(table[:, 0], table[:, 0])  # identity checked on every node

    def evaluate(self, created_nm, detected_nm):
        created, detected = np.broadcast_arrays(np.asarray(created_nm, float),
                                               np.asarray(detected_nm, float))
        if (not np.all(np.isfinite(created)) or not np.all(np.isfinite(detected))
                or np.any(created <= 0) or np.any(detected <= 0)):
            raise ValueError('Nonpositive/nonfinite photon wavelength')
        n = np.interp(HC_EV_NM/created, self.energy, self.index)
        nd = np.interp(HC_EV_NM/detected, self.energy, self.index)
        vg_created = np.interp(HC_EV_NM/created, self.group_energy, self.group_speed)
        vg = np.interp(HC_EV_NM/detected, self.group_energy, self.group_speed)
        if np.any(self.beta*n <= 1):
            raise ValueError('Configured primary is below Cherenkov threshold')
        critical = np.arcsin(1/n)
        cone_beta1 = np.arccos(1/n)
        cone = np.arccos(1/(self.beta*n))
        identity = np.arccos(np.sin(cone_beta1))
        if not np.allclose(identity, critical, rtol=0, atol=2e-14):
            raise AssertionError('beta=1 cone/critical identity failed')
        outside = lambda w: (w < MEASURED_MIN_NM) | (w > MEASURED_MAX_NM)
        return {
            'n_created': n, 'n_detected': nd,
            'theta_cherenkov_beta1_deg': np.degrees(cone_beta1),
            'theta_cherenkov_deg': np.degrees(cone),
            'theta_critical_created_deg': np.degrees(critical),
            'theta_critical_detected_deg': np.degrees(np.arcsin(1/nd)),
            'edge_angle_deg': np.degrees(np.arcsin(1/(self.beta*n))),
            'phase_speed_created_mm_ns': C/n,
            'group_index_created': C/vg_created, 'group_index_detected': C/vg,
            'group_speed_created_mm_ns': vg_created, 'group_speed_mm_ns': vg,
            # Requested phase-edge expression; not a dispersive transport speed.
            'phase_edge_beta1_mm_ns': C*np.sqrt(1-1/n**2)/n,
            'transport_edge_beta1_mm_ns': vg*np.sqrt(1-1/n**2),
            'transport_edge_mm_ns': vg*np.sin(cone),
            'spectral_delay_ns_per_mm': 1/vg - 1/vg_created,
            'group_minus_phase_delay_ns_per_mm': 1/vg - nd/C,
            'created_outside_measured': outside(created).astype(float),
            'detected_outside_measured': outside(detected).astype(float),
            'created_below_measured': (created < MEASURED_MIN_NM).astype(float),
            'created_above_measured': (created > MEASURED_MAX_NM).astype(float),
            'detected_below_measured': (detected < MEASURED_MIN_NM).astype(float),
            'detected_above_measured': (detected > MEASURED_MAX_NM).astype(float),
            # Constant tables (e.g. EJ-230) have no measured dispersion interior.
            'created_clamped': (outside(created) & ~np.asarray(self.constant)).astype(float),
            'detected_clamped': (outside(detected) & ~np.asarray(self.constant)).astype(float),
        }


@lru_cache(maxsize=4)
def campaign_tables(campaign=str(CAMPAIGN)):
    """Resolve each actual cells/<id>/sslg4, never an OPSC-wide default."""
    directory = Path(campaign)
    cells = json.loads((directory/'campaign.json').read_text())['cells']
    tables, records = {}, []
    for cell in cells:
        local = directory/'cells'/cell['cell_id']
        runtime = (local/'sslg4').resolve(strict=True)
        macro = (local/'run.mac').read_text()
        angle = re.findall(r'^/muon/angle\s+([\d.eE+-]+)\s*$', macro, re.M)
        kinetic = re.findall(r'^/gun/energy\s+([\d.eE+-]+)\s+(MeV|GeV)\s*$', macro, re.M)
        if len(angle) != 1 or float(angle[0]) != 0 or len(kinetic) != 1:
            raise ValueError(f'Primary geometry/energy not supported: {cell["cell_id"]}')
        energy, unit = kinetic[0]
        gamma = 1+float(energy)*(1000 if unit == 'GeV' else 1)/MUON_MASS_MEV
        beta = np.sqrt(1-gamma**-2)
        path = runtime/'data/oscnt'/cell['opsc'].lower()/'rIndex.txt'
        optical = OpticalTable(path, beta)
        key = (MATERIALS.index(cell['material']), int(cell['x_mm']))
        if key in tables:
            raise ValueError(f'Duplicate cell {key}')
        tables[key] = optical
        absorption = runtime/'data/oscnt'/cell['opsc'].lower()/'absLength.txt'
        abs_table = np.loadtxt(absorption)
        records.append({'cell_id': cell['cell_id'], 'material': cell['material'],
            'x_mm': cell['x_mm'], 'runtime': str(runtime), 'rindex_path': str(path),
            'rindex_sha256': optical.sha256, 'rindex_constant': optical.constant,
            'rindex_min': float(optical.index.min()), 'rindex_max': float(optical.index.max()),
            'abs_length_min_mm': float(abs_table[:, 1].min()*10),
            'abs_length_max_mm': float(abs_table[:, 1].max()*10),
            'abs_length_constant': bool(np.ptp(abs_table[:, 1]) == 0),
            'beta': float(beta), 'identity_beta1': 'PASS',
            'outside_370_660': ('unmeasured constant-index model' if optical.constant
                                else 'constant clamp; not measured')})
    return tables, records


def photon_optics(frame, prefix='', position='gun_x_mm'):
    """Array dictionary aligned to frame; accepts dataframes or array mappings."""
    tables, _ = campaign_tables()
    mat = np.asarray(frame['material_code'])
    pos = np.asarray(frame[position])
    created = np.asarray(frame[prefix+'wl_nm_created'])
    detected = np.asarray(frame[prefix+'wl_nm'])
    result = {}
    for code, x in np.unique(np.column_stack([mat, pos]), axis=0):
        mask = (mat == code) & (pos == x)
        for name, value in tables[(int(code), int(x))].evaluate(created[mask], detected[mask]).items():
            if name not in result:
                result[name] = np.empty(len(mat), float)
            result[name][mask] = value
    return result


def attach_optics(frame, prefix=''):
    """New columns only; caller retains the original photon selection."""
    result = frame.copy()
    for key, value in photon_optics(frame, prefix).items():
        result['opt_'+key] = value
    return result


def distribution_summary(values):
    result = {}
    for name, value in values.items():
        value = np.asarray(value)
        if name.endswith(('clamped', '_measured')):
            result[name+'_fraction'] = float(np.mean(value))
        else:
            for label, q in zip(('q05', 'median', 'q95'), np.quantile(value, QUANTILES)):
                result[name+'_'+label] = float(q)
    return result


def optical_summary(frame):
    return distribution_summary({name[4:]: frame[name] for name in frame if name.startswith('opt_')})


def quantile_text(row, name, digits=3):
    return (f"{row[name+'_median']:.{digits}f} "
            f"[{row[name+'_q05']:.{digits}f}, {row[name+'_q95']:.{digits}f}]")


def optical_markdown(rows):
    """Population-labelled distributions, never a single material constant."""
    lines = ['Median [q05, q95]; wavelengths weighted by the selected photons. '
             'Cone: created wavelength; transport: detected wavelength. '
             'Spectral correction is (1/vg_det - 1/vg_created), in ps/mm. '
             'Constant-index models have no measured dispersion domain.', '',
        '| material | population | n(created) | theta_C(beta) [deg] | critical(created) [deg] | '
        'edge(beta) [deg] | vg(det) [mm/ns] | transport edge [mm/ns] | '
        'spectral delay [ps/mm] | group-phase delay [ps/mm] | created clamp / outside fraction | detected clamp / outside fraction |',
        '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for row in rows:
        delays = []
        for key in ('spectral_delay_ns_per_mm', 'group_minus_phase_delay_ns_per_mm'):
            delays.append(' / '.join(f'{1000*row[key+"_"+q]:.6g}'
                                    for q in ('median', 'q05', 'q95')))
        lines.append('| '+ ' | '.join([row['material'], row['population'],
            *(quantile_text(row, name) for name in ('n_created', 'theta_cherenkov_deg',
              'theta_critical_created_deg', 'edge_angle_deg', 'group_speed_mm_ns',
              'transport_edge_mm_ns')), *delays,
            *(f"{row[end+'_clamped_fraction']:.4%} / {row[end+'_outside_measured_fraction']:.4%}"
              for end in ('created', 'detected'))])+ ' |')
    lines += ['', '| material | population | phase speed(created) [mm/ns] | group index(det) | phase edge(beta=1) [mm/ns] |',
              '|---|---|---|---|---|']
    for row in rows:
        lines.append('| ' + ' | '.join([row['material'], row['population'],
            *(quantile_text(row, name) for name in ('phase_speed_created_mm_ns',
              'group_index_detected', 'phase_edge_beta1_mm_ns'))]) + ' |')
    lines += ['', 'The beta=1 identity arccos(sin(theta_C(lambda))) = arcsin(1/n(lambda)) '
              'is checked at every table node and evaluated photon wavelength; dispersion '
              'replicates the identity wavelength by wavelength. Arrival-time edge speeds use '
              'vg(lambda_det)*sin(theta_C(lambda_created)); the phase-edge formula '
              'c*sqrt(1-1/n^2)/n is retained separately. Multi-medium/reflected paths are '
              'not reconstructed by this homogeneous-scintillator reference.']
    return '\n'.join(lines)
