# BUILD.md — How to compile talk_v6

## Quick build (host AlmaLinux → WSL2 → retrieve PDF)

```bash
# 1. [Optional] Regenerate figures on host (if data changed)
source /opt/root/bin/thisroot.sh
cd /home/rrios/ej200
python3 analysis/presentation_v5/scripts/analysis_v5.py  # v5 figs (same data)

# 2. Copy .tex and figs to WSL2
scp -P 9022 analysis/presentation_v6/talk_v6.tex \
    reriosto@127.0.0.1:/home/reriosto/SHiP/ej200/presentations/presentation_v6/

scp -P 9022 analysis/presentation_v6/figs/*.pdf \
    reriosto@127.0.0.1:/home/reriosto/SHiP/ej200/presentations/presentation_v6/figs/

# 3. Compile on WSL2 (two passes for cross-references)
ssh -p 9022 reriosto@127.0.0.1 \
    "cd /home/reriosto/SHiP/ej200/presentations/presentation_v6 && \
     pdflatex -interaction=nonstopmode talk_v6.tex && \
     pdflatex -interaction=nonstopmode talk_v6.tex"

# 4. Retrieve PDF
scp -P 9022 \
    reriosto@127.0.0.1:/home/reriosto/SHiP/ej200/presentations/presentation_v6/talk_v6.pdf \
    /home/rrios/ej200/analysis/presentation_v6/talk_v6.pdf
```

## Visual QA

```bash
mkdir -p /tmp/qa_v6
pdftoppm -png -r 120 analysis/presentation_v6/talk_v6.pdf /tmp/qa_v6/slide
ls /tmp/qa_v6/  # should have 36 files
```

## TeX packages required (WSL2 TeX Live 2020)

| Package | Used for |
|---------|---------|
| beamer | presentation class |
| Warsaw theme | color/layout |
| amsmath, amssymb | equations |
| booktabs, tabularx, colortbl | tables |
| graphicx | figures |
| tikz | geometry + physics diagrams |
| tikz libraries: arrows.meta, decorations.pathreplacing, patterns, calc, shapes.geometric | TikZ features |

**NOT used** (not available or not needed):
- multirow — removed
- pgfplots — removed (replaced with hand-drawn TikZ)
- XeLaTeX — not needed, pdflatex only

## Figure dependencies

All figures in `figs/`:

| File | Origin | Description |
|------|--------|-------------|
| `v5_veff_fit.pdf` | analysis_v5.py §1 | Δt vs x linear fit |
| `v5_veff_residual.pdf` | analysis_v5.py §1 | Nonlinearity residual |
| `v5_sigma_vs_x.pdf` | analysis_v5.py §2-3 | Global+oracle σ(x) |
| `v5_top_position_loo.pdf` | analysis_v5.py §4 | TOP LOO position |
| `v5_pareto.pdf` | analysis_v5.py §6 | N_TOP scan |
| `fig1_kscan.pdf` | root_best_est | k-scan at x=0 |
| `fig_end_mscan.pdf` | root_best_est | m-scan at x=0 |
| `fig_ntop_scan_full.pdf` | root_best_est | k-scan all positions |
| `fig_npe_x.pdf` | napkin figs | N_pe vs x |
| `fig_sigma_t_x.pdf` | napkin figs | σ_t vs x |
| `fig_bulk_survival.pdf` | napkin figs | Bulk att. survival |
| `fig_survival_RN.pdf` | napkin figs | Reflector bounce survival |
| `figM1_mat_sigma_end.pdf` | materials analysis | σ_END material comparison |
| `figM2_mat_npe.pdf` | materials analysis | N_pe material comparison |

## Input data

```
results/scan_end_vik_sparse_top_v2/photon_hits_ej230_ntop{4,8,14,20}_x{-600,...,+600}.root
```

13 positions × 4 N_TOP values × 5000 events each. EJ-230, GEN-3 geometry.

## Known compile characteristics

- Zero Overfull \hbox warnings (verified 2026-08-21)
- One suppressed Overfull \vbox on BLUE slide (handled via `[shrink=3]`)
- 36 pages total (2 title/outline + 22 main + 1 appendix divider + 13 appendix backup)
- Output: ~480 KB PDF
