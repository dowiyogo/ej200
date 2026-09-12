"""EXEC34R tests: synthetic ROOTs and mocked processes; never invoke Geant4."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, Mock

import numpy as np
import uproot
import detached_grid as d


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        # Resource gates are deterministic in unit tests (tmp may be on a
        # smaller filesystem). The real dry-run checks the actual /home disk.
        patcher = patch.object(d.shutil, 'disk_usage', return_value=Mock(free=2*1024**4))
        patcher.start()
        self.addCleanup(patcher.stop)
        patcher = patch.object(d, 'memory', return_value={'MemAvailable': 100*1024**3})
        patcher.start()
        self.addCleanup(patcher.stop)

    def fixture(self, code='OPSC-100'):
        original = d.read(d.SOURCE)
        cells = []
        for c in original['cells']:
            p = self.path/'source'/c['cell_id']
            p.mkdir(parents=True)
            (p/'run.mac').write_bytes((Path(c['output'])/'run.mac').read_bytes())
            cells.append(dict(c, output=str(p)))
        source = self.path/'source.json'
        source.write_text(json.dumps(dict(original, cells=cells)))
        h = d.read(d.HANDOFF)
        h['EJ200_OPSC_CODE'] = code
        handoff = self.path/'handoff.json'
        handoff.write_text(json.dumps(h))
        dest = self.path/'campaign'
        with contextlib.redirect_stdout(io.StringIO()):
            d.prepare(dest, source, handoff)
        return d.read(dest/'campaign.json')

    def root(self, n=10000, bad_ids=False):
        path = self.path/'test.root'
        with uproot.recreate(path) as f:
            tree = f.mktree('sipm_hits', dict(event_id='int32', face_type='int32', global_id='int32', time_ns='float64'))
            ids = np.arange(n, dtype=np.int32)
            if bad_ids:
                ids[-1] = 10000
            tree.extend(dict(event_id=ids, face_type=np.zeros(n, dtype=np.int32),
                             global_id=np.zeros(n, dtype=np.int32), time_ns=np.ones(n)))
        log = self.path/'stdout.log'
        log.write_text(f'Events run : 10000\nEvents with ≥1 hit : {n}\nEnd-left photons : {n}\nEnd-right photons : 0\nTop SiPM photons : 0\n')
        return path, log

    def test_macro_changes_only_two_numeric_tokens(self):
        raw = b'/run/numberOfThreads 1 # keep\r\n/run/eventModulo 7\r\n/random/setSeeds 1 2\r\n'
        self.assertEqual(d.adjusted_macro(raw), raw.replace(b'Threads 1', b'Threads 4').replace(b'Modulo 7', b'Modulo 1'))
        self.assertEqual(d.adjusted_macro(d.adjusted_macro(raw)), d.adjusted_macro(raw))
        with self.assertRaises(ValueError):
            d.adjusted_macro(raw + b'/run/eventModulo 1\n')

    def test_prepare_and_dry_run_never_spawn(self):
        with patch.object(d.subprocess, 'Popen', side_effect=AssertionError('MUST NOT SPAWN')):
            cfg = self.fixture()
            self.assertEqual(len(cfg['cells']), 21)
            before = {p: d.sha(p) for p in Path(cfg['output_directory']).rglob('*') if p.is_file()}
            output = io.StringIO()
            with patch('sys.argv', ['detached_grid.py', 'dry-run', '--directory', cfg['output_directory']]), contextlib.redirect_stdout(output):
                self.assertEqual(d.main(), 0)
            self.assertIn('JOBLIST 21', output.getvalue())
            self.assertIn('EXECUTED=0', output.getvalue())
            self.assertEqual(before, {p: d.sha(p) for p in Path(cfg['output_directory']).rglob('*') if p.is_file()})
            self.assertEqual(d.status(cfg)['counts'], dict(complete=0, running=0, pending=21, failed=0))

    def test_not_found_excludes_exactly_seven(self):
        cfg = self.fixture('NOT_FOUND')
        self.assertEqual(len(cfg['cells']), 14)
        self.assertEqual(len(cfg['excluded']), 7)
        self.assertFalse(any(c['material'] == 'EJ-200' for c in cfg['cells']))
        d.preflight(cfg)

    def test_preflight_aborts_disk_ram_binary_and_macro(self):
        cfg = self.fixture()
        with patch.object(d.shutil, 'disk_usage', return_value=Mock(free=0)):
            with self.assertRaisesRegex(ValueError, 'Disk'):
                d.preflight(cfg)
        with patch.object(d, 'memory', return_value={'MemAvailable': 1}):
            with self.assertRaisesRegex(ValueError, 'RAM'):
                d.preflight(cfg)
        changed = dict(cfg, binary_sha256='wrong')
        with self.assertRaisesRegex(ValueError, 'Binary'):
            d.preflight(changed)
        macro = Path(cfg['cells'][0]['output'])/'run.mac'
        macro.write_bytes(macro.read_bytes().replace(b'beamOn 10000', b'beamOn 10001'))
        with self.assertRaisesRegex(ValueError, 'Macro differs'):
            d.preflight(cfg)

    def test_readable_root_and_zero_hit_events_supported(self):
        root, log = self.root(n=9999)
        result = d.verify_output(root, log)
        self.assertEqual(result['events_run'], 10000)
        self.assertEqual(result['events_with_hits'], 9999)
        self.assertEqual(result['root_sha256'], d.sha(root))

    def test_preflight_rejects_diagnostics_on_and_unreadable_macro(self):
        cfg = self.fixture()
        original = Path.read_text
        def changed(path, *args, **kwargs):
            value = original(path, *args, **kwargs)
            return value.replace('EJ200_ENABLE_DIAGNOSTICS:BOOL=OFF',
                                 'EJ200_ENABLE_DIAGNOSTICS:BOOL=ON') if path.name == 'CMakeCache.txt' else value
        with patch.object(Path, 'read_text', changed):
            with self.assertRaisesRegex(ValueError, 'Diagnostics OFF'):
                d.preflight(cfg)
        original_bytes = Path.read_bytes
        blocked = Path(cfg['cells'][0]['output'])/'run.mac'
        def unreadable(path):
            if path == blocked:
                raise PermissionError('Unreadable macro')
            return original_bytes(path)
        with patch.object(Path, 'read_bytes', unreadable):
            with self.assertRaises(PermissionError):
                d.preflight(cfg)

    def test_reject_short_summary_wrong_ids_and_corrupt_root(self):
        root, log = self.root(bad_ids=True)
        with self.assertRaisesRegex(ValueError, 'event ID'):
            d.verify_output(root, log)
        log.write_text(log.read_text().replace('Events run : 10000', 'Events run : 9999'))
        with self.assertRaisesRegex(ValueError, '10000'):
            d.verify_output(root, log)
        root.write_bytes(b'not a ROOT')
        log.write_text(log.read_text().replace('Events run : 9999', 'Events run : 10000'))
        with self.assertRaises(Exception):
            d.verify_output(root, log)

    def test_done_integrity_and_restart_states(self):
        root, _ = self.root()
        base = self.path/'cell'
        base.mkdir()
        cell = dict(output=str(base), cell_id='C', macro_sha256='macro')
        self.assertEqual(d.state(cell), 'pending')
        d.write(base/'state.json', dict(status='RUNNING', process=d.identity(os.getpid())))
        self.assertEqual(d.state(cell), 'running')
        d.write(base/'state.json', dict(status='RUNNING', process={'pid': 99999999}))
        self.assertEqual(d.state(cell), 'failed')
        d.write(base/'.DONE', dict(macro_sha256='macro', events_run=10000, exit_code=0,
               root_path=str(root), root_size_bytes=root.stat().st_size, root_sha256=d.sha(root)))
        self.assertTrue(d.done_record(cell))
        self.assertEqual(d.state(cell), 'complete')
        # Same-length corruption must be caught on relaunch.
        with root.open('r+b') as f:
            f.write(b'xxxx')
        with self.assertRaisesRegex(ValueError, 'hash changed'):
            d.done_record(cell)

    def test_duplicate_driver_lock(self):
        with d.acquire(self.path):
            with self.assertRaisesRegex(ValueError, 'already hold'):
                d.acquire(self.path)
        with d.acquire(self.path):
            pass

    def test_launch_detaches_and_redirects_mock_only(self):
        cfg = self.fixture()
        fake = Mock(pid=os.getpid())
        with patch.object(d.subprocess, 'Popen', return_value=fake) as spawn, patch('sys.argv',
             ['detached_grid.py', 'launch', '--directory', cfg['output_directory']]), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(d.main(), 0)
        self.assertTrue(spawn.call_args.kwargs['start_new_session'])
        self.assertEqual(spawn.call_args.kwargs['stdin'], d.subprocess.DEVNULL)
        self.assertEqual(len(spawn.call_args.kwargs['pass_fds']), 1)
        self.assertIn('_driver', spawn.call_args.args[0])
        self.assertNotIn(cfg['binary'], spawn.call_args.args[0])

    def test_cell_success_done_then_skipped_and_failure_preserved(self):
        cfg = self.fixture()
        cell = cfg['cells'][0]
        root, log = self.root()
        fake = Mock(pid=os.getpid())
        fake.poll.return_value = 0
        fake.wait.return_value = 0
        def simulate(argv, **kwargs):
            attempt = Path(kwargs['cwd'])
            d.shutil.copyfile(root, attempt/'photon_hits_run000.root')
            kwargs['stdout'].write(log.read_text())
            kwargs['stdout'].write('SiPM PDE file : ' + cfg['PDE_path'] + '\n')
            kwargs['stdout'].flush()
            self.assertEqual(kwargs['env']['EJ200_DATA_DIR'], str(Path(cfg['PDE_path']).parent.parent))
            self.assertNotIn('G4FORCENUMBEROFTHREADS', kwargs['env'])
            return fake
        with d.acquire(Path(cfg['output_directory'])) as lock, patch.object(d.subprocess, 'Popen', side_effect=simulate), contextlib.redirect_stdout(io.StringIO()):
            self.assertTrue(d.run_cell(cfg, cell, lock.fileno()))
        self.assertTrue(d.done_record(cell))
        rows = [json.loads(s) for s in Path(cfg['manifest']).read_text().splitlines()]
        self.assertEqual(rows[-1]['status'], 'SIMULATION_COMPLETE')
        self.assertEqual(rows[-1]['workers'], 4)
        self.assertEqual(rows[-1]['parallelism_provenance']['exact_delta_npe_end'], 0)
        one = dict(cfg, cells=[cell])
        with patch.object(d, 'preflight', return_value={}), patch.object(d, 'run_cell', side_effect=AssertionError('Must skip DONE')), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(d.run_driver(one, 0), 0)
        # Nonzero exit never marks DONE and each retry retains old attempts.
        failed = cfg['cells'][1]
        fake.poll.return_value = 7
        fake.wait.return_value = 7
        with d.acquire(Path(cfg['output_directory'])) as lock, patch.object(d.subprocess, 'Popen', return_value=fake), contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(d.run_cell(cfg, failed, lock.fileno()))
            self.assertFalse(d.run_cell(cfg, failed, lock.fileno()))
        self.assertFalse((Path(failed['output'])/'.DONE').exists())
        self.assertEqual(len(list((Path(failed['output'])/'attempts').iterdir())), 2)
        self.assertEqual(d.state(failed), 'failed')

    def test_driver_continues_after_cell_failure_six_slots(self):
        cfg = self.fixture()
        seen = []
        def simulated(c, cell, fd):
            seen.append(cell['cell_id'])
            return len(seen) != 1
        with patch.object(d, 'preflight', return_value={}), patch.object(d, 'run_cell', side_effect=simulated), patch.object(d, 'ThreadPoolExecutor', wraps=d.ThreadPoolExecutor) as pool, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(d.run_driver(cfg, 0), 34)
        self.assertEqual(len(seen), 21)
        pool.assert_called_once_with(max_workers=6)


if __name__ == '__main__':
    unittest.main()
