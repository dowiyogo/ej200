import unittest
import numpy as np
from scipy.stats import norm
from robust_widths import (exact_robust, legacy_robust, rsigma, measure,
                           interpolated_fwhm, fixed_timestamps, function_only, PACKAGE)


class WidthTests(unittest.TestCase):
    def test_exact_levels_and_legacy_difference(self):
        v = np.linspace(0, 1, 10001)
        self.assertAlmostEqual(legacy_robust(v)['sigma_core'], .34)
        self.assertAlmostEqual(exact_robust(v)['sigma_core'], .34135)
        self.assertAlmostEqual(rsigma(v), .5/1.349)

    def test_normal_population_quantiles(self):
        v = norm.ppf((np.arange(10000)+.5)/10000)
        self.assertAlmostEqual(exact_robust(v)['sigma_core'], 1., places=3)
        self.assertAlmostEqual(rsigma(v), 1., places=3)

    def test_interpolation_and_boundary_convention(self):
        self.assertEqual(interpolated_fwhm([0, 2, 4, 2, 0], np.arange(6)), 2.)
        self.assertEqual(interpolated_fwhm([4, 0], np.arange(3)), 1.)

    def test_bootstrap_is_reproducible_and_keeps_generated_denominator(self):
        events = np.r_[np.arange(100), np.full(20, np.nan)]
        a, samples, _, _ = measure(events, replicas=30)
        b, again, _, _ = measure(events, replicas=30)
        np.testing.assert_array_equal(samples, again)
        self.assertEqual(a['N_generated_partition'], 120)
        self.assertEqual(a['efficiency'], 100/120)
        self.assertGreater(a['sigma_core_se'], 0.)

    def test_fixed_timestamps_match_original_primitive_with_missing_events(self):
        original = function_only(PACKAGE/'upstream/analysis_core/timing_fit_pipeline.py',
                                 'compute_tN', {'np': np})
        rng = np.random.default_rng(57)
        ev = rng.integers(0, 30, 500); gids = rng.integers(16, 25, 500); ts = rng.random(500)
        channels = [16, 19, 20, 21]
        result = fixed_timestamps(ev, gids, ts, channels, 40)
        for n in range(1, 21):
            expected = original(ev, ts, np.isin(gids, channels), n, np.arange(40))
            np.testing.assert_array_equal(result[np.isfinite(result[:, n-1]), n-1], expected)


if __name__ == '__main__':
    unittest.main()
