import unittest
import numpy as np
from run_exec35 import end_callback, annotate
from robust_widths import measure


class EndBootstrapTest(unittest.TestCase):
    def test_unchanged_end_primitives_bootstrap_units(self):
        values = np.random.default_rng(719).normal(2000,120,1000)
        callback, version = end_callback()
        row, boot, edges, extra = measure(values,replicas=5,callback=callback)
        self.assertEqual(np.asarray(extra).shape,(5,2))
        self.assertTrue(np.all((np.asarray(extra)>80)&(np.asarray(extra)<160)))
        self.assertTrue(100 < row['sigma_core'] < 140)
        self.assertTrue(version)

    def test_bad_gaussian_does_not_veto_valid_robust(self):
        values = np.random.default_rng(100).exponential(5,1000)
        row, *_ = measure(values,replicas=20)
        old = dict(sigma_ps=90.,uncertainty_ps=4.,bootstrap_error_ps=3.,fit_status=4,chi2_ndf=190.)
        r = annotate(row,old,'TOP','c','EVAL',3)
        self.assertEqual(r['primary_validity']['status'],'VALID')
        self.assertEqual(r['gaussian_validity']['status'],'INVALID')


if __name__=='__main__':
    unittest.main()
