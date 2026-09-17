#!/usr/bin/env python3
"""Render report tables by copying registered CSV columns without recalculation."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / 'generated'
OUT.mkdir(exist_ok=True)

TABLES = [
    ('budget', '../top_npe_diag/top_npe_diag.csv', ['x', 'npe_end_per_face', 'npe_top_mean', 'npe_end_hybrid', 'npe_end_endonly', 'delta_pct'], None),
    ('optical_properties_index', 'step6_v2/step2/material_optical_properties.csv', ['material', 'opsc_code', 'rindex_wavelength_min_nm', 'rindex_wavelength_max_nm', 'rindex_min', 'rindex_max', 'rindex_constant'], None),
    ('optical_properties_absorption', 'step6_v2/step2/material_optical_properties.csv', ['material', 'opsc_code', 'abs_wavelength_min_nm', 'abs_wavelength_max_nm', 'abs_length_min_mm', 'abs_length_max_mm', 'abs_length_constant'], None),
    ('inventory', 'step6_v2/exec46_inventory.csv', ['cell_id', 'material', 'x_mm', 'events', 'photons', 'mean_npe_end', 'all_gates_pass'], None),
    ('tau_diagnostics', 'step6_v2/exec46_tau_diagnostics.csv', None, None),
    ('source_census', 'step6_v2/exec46_source_census.csv', ['scope', 'face_type', 'source_label', 'count', 'denominator', 'fraction'], None),
    ('transport_fits', 'step6_v2/step3/fit_summary.csv', None, None),
    ('transport_gd', 'step6_v2/step3/g_d.csv', ['material', 'source', 'nominal_d_mm', 'mean_tprop_ns', 'var_tprop_ns2', 'n_total'], None),
    ('transport_variance', 'step6_v2/step3/tprop_variance.csv', ['material', 'source', 'face', 'd_center_mm', 'n_photons', 'mean_tprop_ns', 'var_tprop_ns2'], None),
    ('cherenkov_counts', 'step6_v2/step4/cherenkov_counts.csv', None, None),
    ('cherenkov_angles', 'step6_v2/step4/cherenkov_angle_window.csv', None, None),
    ('cherenkov_enrichment', 'step6_v2/step4/cherenkov_enrichment.csv', None, None),
    ('cherenkov_velocity', 'step6_v2/step4/cherenkov_nc_velocity.csv', None, None),
    ('order_scaling', 'step6_v2/step4/scintillation_order_scaling.csv', ['material', 'face', 'nominal_d_mm', 'mean_n_scint', 'mean_min_creation_ns', 'mean_min_corrected_creation_ns', 'fitted_effective_fraction', 'fitted_n_eff'], None),
    ('order_handicap', 'step6_v2/step4/minimum_handicap_d50.csv', ['material', 'face', 'mean_n_scint', 'mean_n_cherenkov', 'predicted_scint_min_delay_ns', 'selection_reordering_gap_ns', 'cherenkov_transport_handicap_ns', 'total_cherenkov_minus_scint_ns'], None),
    ('within_between', 'step6_v2/step5/within_between_summary.csv', ['material', 'beta_within_ns_per_pe', 'beta_between_ns_per_pe', 'beta_gap_significance', 'registered_remnant_a2_ns_per_m2', 'descriptive_residual_a2_ns_per_m2'], None),
    ('between_specification', 'step6_v2/step5/between_common_specification.csv', None, None),
    ('f_test', 'step6_v2/step5/h1_even_f_test.csv', ['material', 'seven_position_f_value', 'seven_position_f_p_value', 'seven_position_linear_rss_ps2', 'seven_position_pol2_rss_ps2'], None),
    ('loo', 'step6_v2/step5/h1_loo_summary.csv', None, None),
    ('mixture_components', 'step6_v2/step5/mixture_components.csv', ['material', 'x_mm', 'term', 'shift_ps', 'bootstrap_se_ps'], None),
    ('t6_widths', 'step6_v2/step6/t6_widths.csv', None, None),
    ('l2_guardrail', 'step6_v2/step6/t6_l2_guardrail.csv', None, None),
    ('veff_slopes', 'step6_v2/veff/u2_slope_fits.csv', None, None),
    ('veff_conventions', 'step6_v2/veff/u3_veff.csv', None, None),
    ('veff_contrast', 'step6_v2/veff/u5_contrast.csv', None, None),
    ('cfd_fits', 'step6_v2/veff_rank/w6_cfd_fits.csv', None, None),
    ('old_new_control', 'step6_v2/t2_old_vs_new.csv', None, None),
]

def esc(value):
    return str(value).replace('&', r'\&').replace('%', r'\%').replace('_', r'\_').replace('#', r'\#')

HEADER_LABELS = {
    'rindex_wavelength_min_nm': r'$n$ $\lambda$ min [nm]',
    'rindex_wavelength_max_nm': r'$n$ $\lambda$ max [nm]',
    'rindex_min': r'$n$ min',
    'rindex_max': r'$n$ max',
    'rindex_constant': r'$n$ constante',
    'abs_wavelength_min_nm': r'$L_{abs}$ $\lambda$ min [nm]',
    'abs_wavelength_max_nm': r'$L_{abs}$ $\lambda$ max [nm]',
    'abs_length_min_mm': r'$L_{abs}$ min [mm]',
    'abs_length_max_mm': r'$L_{abs}$ max [mm]',
    'abs_length_constant': r'$L_{abs}$ constante',
}

def render(name, relative, columns, predicate):
    source = ROOT / relative
    if not source.exists():
        return False
    with source.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    if predicate:
        rows = [row for row in rows if predicate(row)]
    if not rows:
        return False
    columns = columns or list(rows[0])
    columns = [column for column in columns if column in rows[0]]
    rows = rows[:80]
    path = OUT / f'{name}.tex'
    with path.open('w') as stream:
        stream.write('% Generated from ' + relative + '; do not edit.\n')
        stream.write('\\scriptsize\n\\begin{longtable}{' + ' '.join(['p{0.13\\textwidth}'] * len(columns)) + '}\n')
        headers = [HEADER_LABELS.get(column, esc(column)) for column in columns]
        stream.write('\\toprule\n' + ' & '.join(headers) + r' \\' + '\n\\midrule\\endfirsthead\n')
        stream.write('\\toprule\n' + ' & '.join(headers) + r' \\' + '\n\\midrule\\endhead\n')
        for row in rows:
            stream.write(' & '.join(esc(row.get(column, '')) for column in columns) + r' \\' + '\n')
        stream.write('\\bottomrule\n\\end{longtable}\n\\normalsize\n')
    return True

count = sum(render(*table) for table in TABLES)
print(f'Generated {count} tables in {OUT}')