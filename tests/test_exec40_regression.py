"""No transport: test that the exact regression cannot pass on a mean alone."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis/validation'))
from check_exec40 import physical_regression

class ExactRegression(unittest.TestCase):
    def setUp(self):
        self.reference=np.full((2000,3),397,dtype=np.int64)
        self.reference[:729,0]+=1
        self.reference[:,1]=396
        self.reference[:1881,1]+=1

    def test_identical_events_pass(self):
        self.assertEqual(physical_regression(self.reference,self.reference,self.reference)['status'],'PASS')

    def test_compensating_end_changes_abort_despite_identical_mean(self):
        changed=self.reference.copy()
        changed[0,0]+=1
        changed[0,1]-=1
        result=physical_regression(changed,self.reference,changed)
        self.assertEqual(result['difference'],0.)
        self.assertEqual(result['events_with_different_counts'],1)
        self.assertEqual(result['status'],'FAIL_ABORT')

    def test_top_change_aborts(self):
        changed=self.reference.copy()
        changed[7,2]+=1
        self.assertEqual(physical_regression(changed,self.reference,changed)['status'],'FAIL_ABORT')

    def test_wrong_event_ledger_aborts(self):
        ledger=self.reference.copy()
        ledger[3,0]+=1
        self.assertEqual(physical_regression(self.reference,self.reference,ledger)['status'],'FAIL_ABORT')

if __name__=='__main__':
    unittest.main()
