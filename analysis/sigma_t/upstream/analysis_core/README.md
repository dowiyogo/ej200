# analysis_core — EXEC_16 Timing-Fit Pipeline

> **Campaign:** EXEC_16 — seeded Gaussian fit pipeline for temporal resolution  
> **Author:** René Ríos Torres (ULS, SHiP T0 group)  
> **Rollback tag:** `pre-exec16-fit-pipeline` on the `ej200` repo

---

## What this package is

`analysis_core` contains a **new, self-contained timing-resolution pipeline** for the
SHiP T0 scintillator-bar characterization. It processes the ROOT hit-trees produced by the
Geant4 simulations in `validation_skinfix_msi_300ev_31pos/` and computes the temporal
resolution σ_t(x) for each detector-group configuration (TOP SUM4/SUM8, END SUM8, etc.)
using a **Gaussian fit seeded at the histogram peak with a local window**.

This package does **not** modify any existing scripts. It is a clean alternative to
`analyze_skinfix_300ev_31pos.py`, retained in its original location as reference.

---

## Why this is different from the previous script

| Aspect | `analyze_skinfix_300ev_31pos.py` | This pipeline |
|--------|----------------------------------|---------------|
| Resolution metric | `np.std(v, ddof=1)` — global RMS | σ from Gaussian fit on core only |
| Treatment of tails | Tails inflate σ via quadratic weight | Fit window excludes tails |
| Seed for fit | — | MAD×1.4826 (robust, asymptotically Gaussian) |
| Error on σ | Bootstrap on std | Fit error + independent bootstrap |
| Tail diagnostic | `robust_sigma` computed but never plotted | χ²/NDF vs position figure |
| Visual verification | No overlay figures | Histogram + fit overlay for 3 representative positions |

**The key insight** (from Constanza Valdivieso's methodology): the crude RMS includes
photon-arrival-time tails (secondary scintillation components, reflections, delayed
photons visible as bumps at ~1.4–1.6 ns in the TOP histograms). These tails are real but
do not represent the core Gaussian that governs the timing precision of a fast trigger.
The Gaussian fit on `[peak ± 2·σ_seed]` measures the **physically relevant core**.

---

## Directory layout

```
analysis_core/
├── README.md                     # this file
├── timing_fit_pipeline.py        # main entry point (orchestrates everything)
├── lib/
│   ├── fit_engine.py             # PyROOT Gaussian fit, seeded-at-peak
│   ├── robust_seeds.py           # MAD×1.4826, median, histogram peak finder
│   ├── sidecar.py                # writes .root / .csv / .meta.json sidecars
│   └── gates.py                  # QA gates (abort-on-anomaly)
├── config/
│   └── exec16_config.yaml        # ALL tunable constants and hooks
└── outputs/                      # NOT versioned; contains sidecars and figures
    └── .gitkeep
```

---

## How to run

```bash
# From any working directory:
python3 /home/reriosto/SHiP/analysis_core/timing_fit_pipeline.py \
        --config /home/reriosto/SHiP/analysis_core/config/exec16_config.yaml

# Run on a single material for testing (Checkpoint 3):
python3 timing_fit_pipeline.py --config config/exec16_config.yaml --materials EJ-230
```

Requirements: Python ≥ 3.9, ROOT (with PyROOT), uproot, numpy, scipy, matplotlib, pyyaml.

---

## What it produces

For each material × group × N_value, three sidecar files under `outputs/`:

| Extension | Contents |
|-----------|----------|
| `.root` | TH1 (time distribution), TF1 (fitted Gaussian), TGraphErrors (σ_fit vs x) |
| `.csv`  | x_gun_mm, group, N, mu_fit_ns, sigma_fit_ns, sigma_fit_err_ps, bootstrap_err_ps, sigma_mad_ps, chi2_ndf, fit_status, efficiency, n_events |
| `.meta.json` | SHA-256 of input ROOT, fit window, branch+SHA (runtime-verified), timestamp |

Plus figures:
- **σ_fit(x) per group** with bootstrap error bars
- **χ²/NDF vs x per group** (diagnostic: should be ~1 if core is Gaussian)
- **Overlay histograms** at x = −690, 0, +690 mm: histogram + fitted Gaussian + fit window shading

---

## QA gates

| Gate | Trigger | Behaviour |
|------|---------|-----------|
| QA-0 | Runtime SHA ≠ declared SHA | **ABORT** |
| QA-1 | TTree missing expected branch | **ABORT** |
| QA-2 | μ_fit or σ_fit outside physical bounds | **ABORT** per-point |
| QA-3 | χ²/NDF > 3.0 at >30% of positions in a group | Flag group `core_not_gaussian` (warn, continue) |
| QA-3d | SHA-256 manifest mismatch | **ABORT** |

---

## Checkpoints (implementation sequence)

- **CP1** ← *you are here*: structure + config + README
- **CP2**: `fit_engine.py` + `robust_seeds.py` + unit test on one distribution
- **CP3**: full run on EJ-230 only; inspect sidecars and figures
- **CP4**: EJ-204 + comparisons + final report + SHA-256 manifest

---

## Methodology reference

Gaussian-fit-at-peak approach: Constanza Valdivieso, `congruent_sum4_timing.C` (SHiP T0).  
MAD×1.4826 as robust σ estimator: Rousseeuw & Croux (1993).  
Rollback tag `pre-exec16-fit-pipeline` on repo `ej200` at SHA `a0368c4`.
