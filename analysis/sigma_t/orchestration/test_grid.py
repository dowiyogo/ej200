import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import subprocess
import grid


class GridTests(unittest.TestCase):
    def test_append_only_manifest_and_memory_cap(self):
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory);gate=d/'gate.json';gate.write_text('{"gate_status":"PASS"}')
            binary=d/'binary';binary.write_text('test')
            manifest=d/'manifest.jsonl'
            grid.append_manifest(manifest,dict(cell_id='x',status='PENDING'))
            before=manifest.read_bytes()
            cfg=dict(canonical_train_parity='even',workers_per_process=1,
                     pilot_gate=str(gate),pilot_gate_sha256=grid.sha(gate),script_sha256={},
                     binary=str(binary),binary_sha256=grid.sha(binary),manifest=str(manifest),
                     cells=[dict(cell_id='x',output=str(d/'x'))],memory_reserve_fraction=.2,
                     conservative_process_budget_bytes=100,planned_concurrency=24)
            with patch.object(grid,'memory',return_value={'MemAvailable':1000}):
                todo,skipped,n=grid.plan(cfg,False)
            self.assertEqual(n,1);self.assertEqual(len(todo),1);self.assertFalse(skipped)
            self.assertEqual(manifest.read_bytes(),before)  # planning is read-only
            grid.append_manifest(manifest,dict(cell_id='x',status='RUNNING'))
            self.assertTrue(manifest.read_bytes().startswith(before))
            self.assertEqual(grid.read_manifest(manifest)['x']['status'],'RUNNING')
            with self.assertRaises(AssertionError): grid.plan(cfg,False)

    def test_nonpass_pilot_prevents_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            gate=Path(directory)/'gate.json';gate.write_text('{"gate_status":"FAIL"}')
            with self.assertRaises(AssertionError):
                grid.plan(dict(canonical_train_parity='even',workers_per_process=1,pilot_gate=str(gate)),False)

    def test_timeout_status_and_later_cells_continue(self):
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory); (d/'logs').mkdir()
            config=dict(output_directory=str(d),cell_timeout_s=.01,
                        manifest=str(d/'manifest.jsonl'),
                        python='python',run_cell='cell.py',binary='binary')
            cells=[dict(cell_id='slow',output=str(d/'slow'),material='EJ-204',opsc='OPSC-101',x_mm=0),
                   dict(cell_id='later',output=str(d/'later'),material='EJ-204',opsc='OPSC-101',x_mm=200)]
            calls=[]
            def fake_run(command,*args):
                calls.append(command[-1])
                return (124,True) if command[-1]=='0' else (1,False)
            args=type('Args',(),dict(config='unused',resume=False,execute=True))()
            with patch.object(grid.json,'loads',return_value=config), \
                 patch.object(grid.Path,'read_text',return_value='{}'), \
                 patch.object(grid,'plan',return_value=(cells,[],1)), \
                 patch.object(grid,'run_with_timeout',side_effect=fake_run), \
                 patch.object(grid,'append_manifest') as append:
                rc=grid.main(args)
            self.assertEqual(rc,34)
            self.assertEqual(calls,['0','200'])
            finished=[x.args[1] for x in append.call_args_list if x.args[1].get('exit_code') is not None]
            self.assertEqual(finished[0]['status'],'FAILED — TIMEOUT')
            self.assertEqual(finished[1]['status'],'FAILED')


if __name__=='__main__': unittest.main()
