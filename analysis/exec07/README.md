# EXEC_09 / completed EXEC_08 analysis

This directory contains the reproducible EXEC_08 products for the completed
EXEC_07 EndTop scan:

- 31 beam positions, 2,000 events per position.
- SSLG4 OPSC-101 / EJ-204.
- Broadcom AFBR-S4N66P024M.
- 16 End plus 70 Top channels.
- `/sipm/jitterSigma 0`: intrinsic timing only.

The physical Top hardware has 32 SiPMs. The 70-channel Top geometry analyzed
here is deliberately a simulation coverage study and does not replicate it.

## Regenerate everything

From the repository root:

```bash
python analysis/exec07_photon_budget.py \
  --data-dir /home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000 \
  --output-dir analysis/exec07

python analysis/exec07/make_beamer.py \
  --output-dir analysis/exec07 \
  --positions both
```

The EXEC_09 timing-mechanism comparison uses existing EndTop and End-only
data only:

```bash
python analysis/exec07/exec09_timing_mechanism.py \
  --endtop-dir /home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000 \
  --end-only-dir /home/reriosto/SHiP/t0minidaq/scan_end_wrapped_2026-06-09
```

It writes a block-jackknife metric table and the normalized log-y comparison
used to establish preferential late-photon drainage through the Top windows.
The full trace and anti-artifact check are in
`audit/exec09_timing_mechanism.md`.

The first command begins with a blocking validation of all 31 ROOT files. It
requires exactly 2,000 unique event IDs, matching `gun_x_mm`, and aggregate
channel IDs 0 through 85 in every file.

For a quicker development iteration that omits the two 31-page per-channel
PDFs, add `--skip-channel-pdfs`. Do not use that option for the final
deliverable.

## Single geometry source

`common.py` is the only source for scan positions, channel groups, Top
positions, nearest-channel selection and EJ-204 constants. Both the analysis
and the 31 per-position geometry figures import it.

At an exact geometric tie, such as x=0 between IDs 50 and 51, the analysis
selects the higher-N_pe member of the tied nearest set. This keeps the
highlighted channel geometrically valid while making the required
nearest-channel/maximum-profile cross-check unambiguous.

The localization gate is strict outside the known window-track pattern. At
positions `x % 20 == 10`, the maximum may be the nearest or second-nearest
window when the nearest deficit is at most 15%; see
`audit/exec08b_window_dip.md`. Outside that pattern, a non-nearest maximum is
accepted only when statistically compatible with the nearest at one sigma.

## SUM4 leading edge

The Python analysis ports the existing implementation from
`analysis/congruent_sum4_timing.C`:

- clusters: `{0..3}`, `{4..7}`, `{8..11}`, `{12..15}`;
- normalized double-exponential single-PE pulse, 0.5/5 ns;
- absolute leading-edge threshold of 4 PE;
- no time-walk correction.

The reported End `sigma_group = sigma(DeltaT)/sqrt(2)` is intrinsic and must
not be interpreted as including sensor SPTR or electronics jitter.

## Main products

- `summary_exec07.csv`
- `exec07_photon_budget_report.pdf`
- `exec09_report_full.pdf` (51 slides, all 31 positions)
- `exec09_report_key.pdf` (27 slides, 7 key positions)
- `figs/P1_...png` through `figs/P7_...png`
- `figs/exec09_tail_comparison.png`
- `exec09_tail_metrics.csv`, `exec09_tail_comparison.csv`
- `top_localization_gate.csv`
- 31 `figs/muon_{position}mm_geometry.png`
