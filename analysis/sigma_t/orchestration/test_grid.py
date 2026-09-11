import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
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


if __name__=='__main__': unittest.main()
