import unittest
from validity_exec35 import verdict


class ValidityTests(unittest.TestCase):
    def sample(self):
        return dict(provenance_valid=True, n_eff=1000, efficiency=1.,
                    sigma_core=10., sigma_core_se=.2, sigma_core_bootstrap_fraction=1.,
                    sigma_gauss=10., sigma_gauss_se=.2,
                    sigma_gauss_bootstrap_fraction=1., gauss_fit_status=0, gauss_chi2_ndf=1.)

    def test_failed_gaussian_does_not_veto_distribution_free_width(self):
        r = self.sample(); r.update(gauss_fit_status=4, gauss_chi2_ndf=100.)
        self.assertEqual(verdict(r, True)['status'], 'INVALID')
        self.assertEqual(verdict(r)['status'], 'VALID')

    def test_primary_rejects_insufficient_sample_and_unstable_bootstrap(self):
        r = self.sample(); r['n_eff'] = 199
        self.assertIn('n_eff_below_200', verdict(r)['reasons'])
        r = self.sample(); r['sigma_core_bootstrap_fraction'] = .94
        self.assertEqual(verdict(r)['status'], 'INVALID')
        r = self.sample(); r['sigma_core_se'] = 3.
        self.assertEqual(verdict(r)['status'], 'INVALID')

    def test_good_fit_and_exact_boundaries(self):
        r = self.sample(); r.update(n_eff=200, efficiency=.05, gauss_chi2_ndf=2.)
        self.assertEqual(verdict(r, True)['status'], 'VALID')
        r['sigma_core_se'] = 0.
        self.assertEqual(verdict(r)['status'], 'INVALID')


if __name__ == '__main__':
    unittest.main()
