"""Check fixed END map, simultaneous crossings, missing ends and zero events."""
import tempfile
from pathlib import Path
import unittest
import json
import numpy as np
import uproot
import ROOT


class EndTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        header=Path(__file__).with_name('end_bridge.h').resolve()
        assert ROOT.gInterpreter.Declare('#include '+json.dumps(str(header)))

    def test_native_map_reduction_and_explicit_zero_events(self):
        # event 0 invisible; event 1 both ends; event 2 left only.
        # On left the second cluster wins, testing first-crossing reduction.
        events=np.array([1]*15+[2]*5,dtype=np.int32)
        gids=np.array([0]*5+[4]*5+[12]*5+[3]*5,dtype=np.int32)
        times=np.array([4.0]*5+[2.0]*5+[2.3]*5+[1.0]*5)
        face=np.where(gids<8,0,1).astype(np.int32)
        with tempfile.TemporaryDirectory() as directory:
            raw=Path(directory)/'photon_hits_run000.root'
            with uproot.recreate(raw) as f:
                f.mktree('sipm_hits',{'event_id':'int32','global_id':'int32','face_type':'int32','time_ns':'float64'})
                f['sipm_hits'].extend(dict(event_id=events,global_id=gids,face_type=face,time_ns=times))
            output=Path(directory)/'out.root'
            root=ROOT.TFile(str(output),'RECREATE')
            result=list(ROOT.exec34.analyze_end(str(raw),3,root))
            root.Close()
            self.assertEqual(result[3],1)
            self.assertFalse(result[4])  # too few events must not look fitted
            with uproot.open(output) as f:
                rows=f['end_events'].arrays(library='np')
                self.assertEqual(len(rows['event_id']),3)
                np.testing.assert_array_equal(rows['accepted'],[False,True,False])
                np.testing.assert_array_equal(rows['npe_left'],[0,10,5])
                self.assertAlmostEqual(rows['delta_ns'][1],-.3,places=12)


if __name__=='__main__': unittest.main()
