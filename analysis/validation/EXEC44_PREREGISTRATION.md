# EXEC44 preregistration — qualitative H3' acceptance replacement

Published before evaluating H3' or changing the golden acceptance result.
EXEC44 uses only already-published EXEC42/43 results. It does not recalculate
an observable, reopen a simulation ROOT, or authorize transport.

The earlier records remain unchanged:

- `EXEC40_PREREGISTRATION.md`: SHA256
  `cd13d28db1811e5c74d9b17d60d1be7aa0971f02ff8768d48528126602fb5531`
- `EXEC41_PREREGISTRATION.md`: SHA256
  `933433a2882c4753301c675913df1f37db614e4afc8316ae87df6c122021d418`
- `EXEC42_H3_PREREGISTRATION.md`: SHA256
  `37e15cf4a26271075babc90a1efc2f1c0d0a8ce1a0be6ea64e48e1d78053d6f4`
- `EXEC43_PREREGISTRATION.md`: SHA256
  `159ca702a34836bf9e94dac01e5a8e1af02e294d76bbd15087a552b78a15f96e`

The evidence source is `/home/rrios/exec43_20260914/results.json`, SHA256
`c9544351da187a692c34705a6669753e3b5195caa352006b1864625c73f8c6c9`.

## Retraction of H3

The EXEC42 material-invariance H3 is removed from the active acceptance
criteria. Its file, result, and refuting campaign remain in the record. H3 is
classified as a retracted, physically inappropriate hypothesis because its
first-encounter population excludes scintillation photons absorbed in the bulk
before reaching a surface, and that selection depends on material.

## H3' — qualitative absorption-filtering hypothesis

First-encounter escape depends on material only through absorption filtering
along the path before the first surface. H3' passes if and only if all four
conditions pass:

1. The emission-spectrum-weighted effective refractive indices of EJ-200,
   EJ-204 and EJ-230 are mutually compatible within `1e-4`, excluding the
   dispersion mechanism implemented in the model.
2. At every one of the seven common positions, the fraction of absent photons
   (`produced_scint - registered first encounters`) increases monotonically in
   the order EJ-200, EJ-204, EJ-230, the order of `1/lambda_att`.
3. At every common position, first-encounter escape increases monotonically in
   that same order: EJ-200 < EJ-204 < EJ-230.
4. For every available directional comparison, the measured median
   `sec(theta)=1/cos(theta)` is larger for non-escaping photons than for
   escaping photons.

Each condition is evaluated by reading the published EXEC43 result and its
named CSV/JSON evidence only. No quantitative WorldPV-envelope agreement is
required by H3'. A condition with missing published evidence is
`NOT_EVALUABLE`, which prevents H3' from passing.

## H3'' — future quantitative hypothesis, outside acceptance

H3'' would predict each material-pair escape difference from attenuation
lengths and the measured pre-encounter path distribution within a tolerance
declared before evaluation. It is future work and is not an active acceptance
criterion. The current ROOT omits accumulated track length, creation
coordinates, and first-encounter coordinates. EXEC43 found that its central
WorldPV assignment overpredicted all 14 target differences by 3.73–7.63 paired
SE; the broad ambiguity envelope is not used by H3'.

## Acceptance decision

After H3 is retired, `ready_for_acceptance` becomes true if V1 primary,
corrected V2-H1, and matched V5 pass all 21 cells, all four H3' conditions pass,
and no active criterion is `FAIL` or `NOT_EVALUABLE`. H3'' and known model
limitations are registered outside this Boolean rule.
