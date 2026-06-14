# End-only + Mylar EJ-204 — Beamer presentation

Simulation analysis and Beamer deck for the SHiP Timing Detector
"End-only + Mylar" scan (EJ-204, Geant4 11.4.0 / SSLG4-OPSim).

## Quick regeneration

```bash
python build_deck.py \
    --data-dir /home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614 \
    --out-dir  /path/to/output
```

This runs the full pipeline: analysis → figures/tables/`deck_values.json` → `main.tex` → `main.pdf`.

To skip re-analysis (use existing outputs):

```bash
python build_deck.py --skip-analysis
```

To generate `main.tex` without compiling:

```bash
python build_deck.py --skip-analysis --no-compile
```

## Traceability check

```bash
python verify_deck.py --out-dir /path/to/output
```

Verifies that every number in `main.tex` (via LaTeX macros) traces back to
`deck_values.json`, and that `deck_values.json` matches the source CSVs.
Exit code 0 = all checks pass.

## Files

| File | Description |
|------|-------------|
| `analysis_endonly_mylar.py` | Reads ROOT + CSV data, writes figures, tables, `deck_values.json` |
| `build_deck.py` | Reads `deck_values.json`, generates `main.tex`, compiles PDF |
| `verify_deck.py` | Numerical traceability checker |
| `deck_values.json` | Single source of truth for all numbers in the deck |
| `main.tex` | Auto-generated LaTeX; do not edit by hand |
| `main.pdf` | Compiled 16-page Beamer deck |
| `figures/` | 7 PNG figures (DPI 150) |
| `tables/` | 3 LaTeX/CSV tables |

## Data directory (read-only)

Default: `/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614/`

Required inputs:
- `run_metadata.txt` — simulation parameters
- `analysis/attenuation_curve.csv` — Npe vs position
- `analysis/attenuation_fits.csv` — λ_eff fit results
- `analysis/sigma_t_sum4.csv` — σ_t(x) from SUM4 trigger
- `*.root` — raw Geant4 hits (tree `sipm_hits`)

## Key physics caveats

- **σ_t is intrinsic only**: optical propagation + SUM4 estimator.
  SPTR (~106 ps) and FastIC (~10 ps) are NOT folded in.
- **No test-beam counterpart**: April-2026 TB used top readout + Tyvek.
- **SiPM naming**: AFBR-S4N66P014M is the physical device;
  AFBR-S4N66P024M is used only for the PDE curve (63 % @ 420 nm).
- **Longitudinal axis = X** (bar half-length 700 mm in ±X).
