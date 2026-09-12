"""Apply the published EXEC35 rule without changing any fit or selection."""
import json
import math
from pathlib import Path

RULE_PATH = Path(__file__).resolve().parents[1] / 'validity_exec35.json'
RULE = json.loads(RULE_PATH.read_text())


def positive(x):
    return x is not None and math.isfinite(x) and x > 0


def verdict(row, gaussian=False):
    limits = RULE['validity']
    reasons = []
    if not row.get('provenance_valid', False):
        reasons.append('provenance_or_frozen_selection')
    if row['n_eff'] < limits['n_eff_min']:
        reasons.append('n_eff_below_200')
    if row['efficiency'] < limits['efficiency_min']:
        reasons.append('efficiency_below_0.05')
    key = 'sigma_gauss' if gaussian else 'sigma_core'
    width, error = row.get(key), row.get(key + '_se')
    if not positive(width):
        reasons.append('nonpositive_or_nonfinite_width')
    if not positive(error) or (positive(width) and error / width > limits['bootstrap_relative_se_max']):
        reasons.append('invalid_or_large_bootstrap_SE')
    fraction = row.get(key + '_bootstrap_fraction')
    # Old TOP bootstrap success fractions were not archived. Keep that fact
    # explicit, while evaluating the archived finite bootstrap SE as declared.
    if fraction is not None and fraction < limits['bootstrap_finite_fraction_min']:
        reasons.append('bootstrap_finite_fraction_below_0.95')
    if not gaussian and fraction is None:
        reasons.append('missing_bootstrap_fraction')
    if gaussian:
        if row.get('gauss_fit_status') != 0:
            reasons.append('fit_status_nonzero_or_unsuccessful')
        chi = row.get('gauss_chi2_ndf')
        if chi is None or not math.isfinite(chi) or not 0 <= chi <= limits['gaussian_chi2_ndf_max']:
            reasons.append('Gaussian_model_inadequate')
    return {'status': 'INVALID' if reasons else 'VALID', 'reasons': reasons,
            'bootstrap_fraction_verified': fraction is not None}
