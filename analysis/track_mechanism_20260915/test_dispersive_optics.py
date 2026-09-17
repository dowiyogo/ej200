"""Focused optical-prediction tests. Never execute a campaign analysis or simulation.

Run: PYTHONPATH=analysis/track_mechanism_20260915 python3 -m unittest discover \
    -s analysis/track_mechanism_20260915 -p test_dispersive_optics.py -v
The V1/V2 integration tests read only the actual campaign macros and MPT tables.
EXEC46_TEST_STEP3_SOURCE may name an old Step 3 script to demonstrate V3 failure.
"""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import uproot

import dispersive_optics as optics
import analyze_step3 as step3
import analyze_step4 as step4


def selected_fixture():
    rows = []
    for code in range(3):
        table = optics.campaign_tables()[0][code, -650]
        prediction = table.evaluate([370., 660.], [370., 660.])
        angle = prediction['theta_critical_created_deg'].mean()
        for i, wavelength in enumerate((370., 660.)):
            rows.append(dict(material_code=code, gun_x_mm=-650, face_type=0,
                cell_code=code, source_type=2, event_id=i, track_id=10+i,
                x_creation_mm=-650., y_creation_mm=0., t_creation_ns=.01*i,
                t_detection_ns=.4+.01*i, exit_angle_deg=angle,
                wl_nm_created=wavelength, wl_nm=wavelength))
    return pd.DataFrame(rows)


class OpticalPredictions(unittest.TestCase):
    def test_v1_constant_table_exact_display_regression(self):
        values = optics.campaign_tables()[0][2, 0].evaluate([420.], [420.])
        expected = {'theta_cherenkov_beta1_deg': ('50.73475', 5),
            'theta_critical_created_deg': ('39.26525', 5),
            'group_speed_mm_ns': ('189.742062', 6),
            'phase_edge_beta1_mm_ns': ('146.902903', 6),
            'transport_edge_mm_ns': ('146.449825', 6)}
        for name, (value, digits) in expected.items():
            self.assertEqual(f'{values[name][0]:.{digits}f}', value, name)
        self.assertEqual(values['group_minus_phase_delay_ns_per_mm'][0], 0)

    def test_v2_actual_mesh_matches_f4(self):
        tables, _ = optics.campaign_tables()
        # F4 used rounded/analytic references; no alteration of the numerical mesh.
        checks = ((0, 420., {'n_created': (1.626296, 1e-6),
                    'group_index_created': (1.818695, 4e-6),
                    'group_speed_mm_ns': (164.839, .001),
                    'transport_edge_mm_ns': (129.629, .001)}),
                  (1, 408., {'n_created': (1.619785, 1e-6),
                    'group_index_created': (1.744068, 4e-6),
                    'group_speed_mm_ns': (171.893, .001)}))
        for code, wavelength, expected in checks:
            values = tables[code, 0].evaluate([wavelength], [wavelength])
            for name, (value, tolerance) in expected.items():
                self.assertAlmostEqual(values[name][0], value, delta=tolerance)
        edges = tables[0, 0].evaluate([370., 660.], [370., 660.])['edge_angle_deg']
        np.testing.assert_allclose(edges, [37.458, 40.174], rtol=0, atol=.001)

    def test_identity_energy_order_and_clamps(self):
        tables, records = optics.campaign_tables()
        self.assertEqual(len(records), 21)
        for table in tables.values():
            self.assertTrue(np.all(np.diff(table.energy) > 0))
            values = table.evaluate([300., 370., 420., 660., 700.],
                                    [300., 370., 420., 660., 700.])
            identity = np.degrees(np.arccos(np.sin(
                np.radians(values['theta_cherenkov_beta1_deg']))))
            np.testing.assert_allclose(identity, values['theta_critical_created_deg'],
                                       rtol=0, atol=1e-12)
            self.assertEqual(np.mean(values['created_clamped']), 0 if table.constant else .4)
            self.assertEqual(np.mean(values['created_below_measured']), .2)
            self.assertEqual(np.mean(values['created_above_measured']), .2)

    def test_real_runtime_per_cell_not_per_material(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            cells = []
            for x, n in ((0, 1.5), (200, 1.7)):
                name = f'cell{x}'
                runtime = base / f'runtime{x}'
                data = runtime / 'data/oscnt/opsc-100'
                data.mkdir(parents=True)
                (data/'rIndex.txt').write_text(f'200 {n}\n800 {n}\n')
                (data/'absLength.txt').write_text('200 380\n800 380\n')
                local = base/'cells'/name
                local.mkdir(parents=True)
                (local/'sslg4').symlink_to(runtime, target_is_directory=True)
                (local/'run.mac').write_text('/muon/angle 0\n/gun/energy 1 GeV\n')
                cells.append(dict(cell_id=name, material='EJ-200', opsc='OPSC-100', x_mm=x))
            (base/'campaign.json').write_text(json.dumps(dict(cells=cells)))
            tables, records = optics.campaign_tables(str(base))
            self.assertNotEqual(tables[0, 0].sha256, tables[0, 200].sha256)
            self.assertEqual(tables[0, 0].evaluate(420., 420.)['n_created'], 1.5)
            self.assertEqual(tables[0, 200].evaluate(420., 420.)['n_created'], 1.7)

    def test_created_and_detected_wavelength_have_distinct_roles(self):
        table = optics.campaign_tables()[0][0, -650]
        values = table.evaluate([370., 660.], [660., 370.])
        same = table.evaluate([370., 660.], [370., 660.])
        np.testing.assert_array_equal(values['edge_angle_deg'], same['edge_angle_deg'])
        np.testing.assert_array_equal(values['group_speed_mm_ns'], same['group_speed_mm_ns'][::-1])
        self.assertNotEqual(values['spectral_delay_ns_per_mm'][0], 0.)

    def test_v3_step3_uses_photon_wavelength_threshold(self):
        module = step3
        if os.environ.get('EXEC46_TEST_STEP3_SOURCE'):
            spec = importlib.util.spec_from_file_location('old_step3',
                os.environ['EXEC46_TEST_STEP3_SOURCE'])
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        _, summary = module.cherenkov_diagnostics(selected_fixture())
        row = summary[(summary.material == 'EJ-200')
                      & (summary.selection == 'all_source_type_2')].iloc[0]
        # Same measured angle, two wavelengths straddling their own critical angles.
        # A scalar critical angle necessarily gives 0 or 1, never 1/2.
        self.assertEqual(row.fraction_below_critical, .5)
        self.assertGreater(row.theta_critical_created_deg_q95,
                           row.theta_critical_created_deg_q05)

    def test_step4_photonwise_group_penalty(self):
        table = optics.campaign_tables()[0][0, -650]
        wl = np.tile([370., 420., 660.], 10)
        predictions = table.evaluate(wl, wl)
        angle = predictions['edge_angle_deg'] + np.linspace(.1, 3., len(wl))
        tprop = 50 / predictions['group_speed_mm_ns'] / np.cos(np.radians(angle))
        data = dict(material_code=np.zeros(len(wl), int), gun_x_mm=-650,
                    nominal_d_mm=50, npe_cherenkov=np.arange(len(wl))+1,
                    npe_primary_cherenkov=np.arange(len(wl))+1)
        for prefix in ('first_cherenkov_', 'first_primary_cherenkov_'):
            data.update({prefix+'track_id': np.arange(len(wl)),
                prefix+'t_creation_ns': np.zeros(len(wl)), prefix+'t_detection_ns': tprop,
                prefix+'x_mm': -700., prefix+'x_creation_mm': -650.,
                prefix+'exit_angle_deg': angle, prefix+'wl_nm': wl,
                prefix+'wl_nm_created': wl})
        frame = pd.DataFrame(data)
        summary, _ = step4.cherenkov_angle_window(frame)
        np.testing.assert_allclose(summary.corr_predicted_observed_penalty, 1., atol=1e-12)
        scan = step4.cherenkov_nc_scan(frame)
        self.assertEqual(len(scan), 10)
        self.assertTrue((scan.transport_edge_mm_ns_q95 > scan.transport_edge_mm_ns_q05).all())
        counts = step4.cherenkov_angle_by_count(frame)
        for _, row in counts.iterrows():
            selected = step4.quantile_bins(frame.npe_primary_cherenkov) == row.count_quantile
            self.assertAlmostEqual(row.q95_minus_edge_deg,
                np.quantile(angle[selected]-predictions['edge_angle_deg'][selected], .95))

    def test_step3_wavelength_join_preserves_selection_and_rejects_corruption(self):
        selected = selected_fixture().query('material_code == 0').copy()
        source = selected.copy()
        original = selected.drop(columns=['wl_nm', 'wl_nm_created']).iloc[::-1]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/'hits.root'
            cells = [dict(material='EJ-200', x_mm=-650, cell_id='fixture', root_path=str(root))]
            def write(frame):
                with uproot.recreate(root) as output:
                    values = {key: value.to_numpy() for key, value in frame.items()}
                    output.mktree('sipm_hits', {key: value.dtype for key, value in values.items()})
                    output['sipm_hits'].extend(values)
            write(source)
            result = step3.recover_selected_wavelengths(original, cells, step_size=1)
            pd.testing.assert_frame_equal(result[original.columns], original)
            np.testing.assert_array_equal(result.wl_nm, [660., 370.])
            write(source.iloc[:1])
            with self.assertRaisesRegex(RuntimeError, 'missing/invalid'):
                step3.recover_selected_wavelengths(original, cells, step_size=1)
            write(pd.concat([source, source.iloc[:1]]))
            with self.assertRaisesRegex(RuntimeError, 'duplicate production'):
                step3.recover_selected_wavelengths(original, cells, step_size=1)
            source.loc[0, 't_detection_ns'] += 1.
            write(source)
            with self.assertRaisesRegex(RuntimeError, 'provenance mismatch'):
                step3.recover_selected_wavelengths(original, cells, step_size=1)

    def test_step2_critical_angle_uses_detected_wavelength(self):
        import analyze_step2 as step2
        rows = []
        for code in range(3):
            for x in (-650, 650):
                critical = optics.campaign_tables()[0][code, x].evaluate(
                    [370., 660.], [370., 660.])['theta_critical_detected_deg']
                for source in (1, 2):
                    for wl in (370., 660.):
                        row = dict(material_code=code, x_mm=x)
                        for selection in ('first', 'random'):
                            for face in ('left', 'right'):
                                prefix = f'{selection}_{face}_'
                                row.update({prefix+'source_type': source,
                                    prefix+'wl_nm_created': 420., prefix+'wl_nm': wl,
                                    prefix+'exit_angle_deg': critical.mean(),
                                    prefix+'n_boundary_encounters': 3})
                        rows.append(row)
        arrays = {key: value.to_numpy() for key, value in pd.DataFrame(rows).items()}
        with tempfile.TemporaryDirectory() as temp:
            stub = Path(temp)/'fixture.root'
            stub.write_bytes(b'fixture source identity')
            with patch.object(step2, 'STEP2_DIR', Path(temp)), patch.object(step2, 'DERIVED_PATH', stub):
                result = step2.make_guiding_diagnostics(arrays, [])
            for row in result:
                if row['material'] == 'EJ-200':
                    self.assertEqual(row['fraction_exit_angle_above_TIR_critical'], .5)
            for suffix in ('.root', '.csv', '.meta.json', '.pdf'):
                self.assertTrue((Path(temp)/('cherenkov_guiding_diagnostics'+suffix)).is_file())


if __name__ == '__main__':
    unittest.main()
