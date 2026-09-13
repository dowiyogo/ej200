import unittest
import numpy as np
from run_exec37 import crossings,c1
from exec36_bridge import load


class ExplicitPopulationTests(unittest.TestCase):
    def test_2000_generated_with_empty_events(self):
        event=np.array([0]*6+[1999]*6,dtype=np.int32)
        gid=np.array([0,1,2,3,0,1]*2,dtype=np.int32)
        ts=np.r_[np.arange(6)*.02,10+np.arange(6)*.02]
        values=crossings(event,gid,ts,[0,1,2,3],2000)
        self.assertEqual(len(values),2000)
        self.assertEqual(np.isfinite(values).sum(),2)
        self.assertTrue(np.all(np.isnan(values[1:1999])))
        root=load()
        self.assertEqual(values[0],root.exec36.original(root.std.vector('double')(ts[:6].tolist())))
        self.assertEqual(values[-1],root.exec36.original(root.std.vector('double')(ts[6:].tolist())))
        with self.assertRaises(AssertionError):crossings(np.array([2000]),np.array([0]),np.array([0.]),[0],2000)

    def test_bound_and_invalid_are_distinct(self):
        self.assertEqual(c1(88.4,True),'PASS')
        self.assertEqual(c1(88.400001,True),'FAIL')
        self.assertEqual(c1(70.,False),'INVALID')
        self.assertEqual(c1(float('nan'),True),'INVALID')


if __name__=='__main__':unittest.main()
