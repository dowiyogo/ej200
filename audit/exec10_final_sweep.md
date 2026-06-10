# EXEC_10 final consistency sweep

Date: 2026-06-11

| Item | Status | Evidence |
|---|---|---|
| Explicit units on figure axes | OK | EXEC_10 Landau, timing-density, and velocity figures label Npe/event, ns, cm, cm/ns, or dimensionless F explicitly. Existing channel-ID axes are dimensionless identifiers. |
| Attenuation-fit legend precision | CORRECTED | `analysis/exec07_photon_budget.py` now prints lambda and its uncertainty with two decimal places. |
| x=300 mm statistical tie wording | CORRECTED | Position table labels it `statistical tie between two nearest candidates`; backup gate slide retains the measured 0.08 sigma statement. |
| FPT uncertainties rounded to 0.000 ns | CORRECTED | `analysis/exec07/make_beamer.py` uses adaptive 3--5 decimal precision. |
| Configuration commits | CORRECTED | EXEC_10 deck cites real analysis commits `7b4bc80`, `ad26d78`, and `bcefa24`, plus the prior published checkpoint `0201413`. |
| English consistency | OK | New captions, glossary, and slides are in English. |

## Additional findings documented without correction

1. The historical unweighted mean-time fit in `analysis/exec07_photon_budget.py`
   gives about 27.67 cm/ns, while the EXEC_10 statistically weighted fit gives
   25.616 +/- 0.005 cm/ns. Both global linear models are rejected by the
   structured residuals. Proposed correction: replace global straight-line
   velocity language with a future dedicated GROUPVEL/path-length audit.
2. In `analysis/exec07/exec10_landau_mpv.csv`, 35 Moyal-Gaussian fits place
   `sigma_gauss` at its lower bound with large uncertainty, and 43 high-light
   channels have chi2/ndf above 5. Proposed correction: future test-beam
   comparison should rely first on MPV and Landau width, and revisit the
   Gaussian-resolution component only with a constrained detector-response
   model.
3. `analysis/exec07_photon_budget.py` still generates the historical
   `P3_fano_vs_x.png` and Poisson/Gaussian `P4` diagnostics. They remain useful
   low-level diagnostics but are no longer used as the main statistical
   conclusion in the EXEC_10 deck.
