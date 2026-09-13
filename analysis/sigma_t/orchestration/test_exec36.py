import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
import uproot
from exec36_bridge import load,SOURCE,sensitivity_source
from run_exec36 import sanity,CONFIG


class PreservedPrimitives(unittest.TestCase):
    def test_parameterized_baseline_matches_original_exactly(self):
        root=load();rng=np.random.default_rng(36)
        self.assertIn(sensitivity_source(),SOURCE.read_text())
        for count in [0,1,3,4,10,100,1000]:
            values=root.std.vector('double')(rng.uniform(0,10,count).tolist())
            a=root.exec36.original(values);b=root.exec36.varied(values,.5,5.,4.)
            self.assertTrue(a==b or (np.isnan(a) and np.isnan(b)))
        for id in range(16):self.assertEqual(root.exec36.group(id),id//4)

    def test_cache_group_union_and_generated_denominator(self):
        root=load()
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'input.root';cache=Path(tmp)/'cache.root'
            ids=np.array([0,0,2,2,2,8,8,8,8,8,17],np.int32)
            ev=np.array([0]*10+[2],np.int32);t=np.arange(11,dtype=float)*.03
            with uproot.recreate(source) as f:
                f['sipm_hits']=dict(event_id=ev,global_id=ids,face_type=np.where(ids<8,0,np.where(ids<16,1,2)).astype(np.int32),time_ns=t)
            data=root.exec36.Data(str(source),str(cache))
            times=list(data.times(root.std.vector('int')([0,2])))
            self.assertEqual(len(times),10000)
            self.assertTrue(np.isfinite(times[0]));self.assertTrue(np.all(np.isnan(times[1:])))
            original=root.exec36.original(root.std.vector('double')(t[:5].tolist()))
            self.assertEqual(times[0],original)
            copied=root.exec36.Data(str(cache))
            self.assertEqual(copied.left,5);self.assertEqual(copied.right,5)
            self.assertTrue(np.array_equal(times,list(copied.times(root.std.vector('int')([0,2]))),equal_nan=True))

    def test_mapping_and_sanity_boundaries(self):
        self.assertEqual(CONFIG['groups']['V1']['left'],[[0,2],[4,6]])
        self.assertEqual(CONFIG['groups']['V2']['right'],[[8,9,10,11],[12,13,14,15]])
        for value,status in [(99,'FAIL_EXTREME'),(100,'INDETERMINATE'),(150,'PASS'),(350,'PASS'),(500,'INDETERMINATE'),(501,'FAIL_EXTREME')]:
            self.assertEqual(sanity(dict(sigma_core=value,primary_validity=dict(status='VALID'))),status)


if __name__=='__main__':unittest.main()
