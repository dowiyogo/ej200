# SIM_TRUST_VERDICT — EXEC_17 — 2026-06-20 15:20

## Veredicto global: **CONFIABLE**

| Test | Resultado | Evidencia |
|------|----------|-----------|
| V1 | ✓ **PASS** | n_found=86/86; END_L x_mean=-699.5 mm (pred=-699.8); END_R x_mean=699.5 mm (pred=699.8); TOP y_mean=29.6 mm (pred=29.8); TOP x-map rms_err=0.27 mm |
| V2 | ✓ **PASS** | ~MPV sum_eV=300 eV (~99 photons); skewness=3.73 (Landau>0 expected); frac_zero=0.0000 |
| V3 | ✓ **PASS** | Observed peak=407.8 nm; Predicted=408.8 nm; Δ=1.0 nm |
| V4 | ✓ **PASS** | END_L near/far asymmetry=0.999 (>0.1 expected); TOP ridge median \|Δx\|=2.0 mm (<40 mm expected) |
| V5 | ✓ **PASS** | τ_d fit = 1.835 ± 0.013 ns (config=1.8 ns); predicted far-end bump at 93333 ps |
| V6 | ✓ **PASS** | Mean fractional asymmetry=0.0299 (3.0%) over 15 mirror pairs |
| V7 | ✓ **PASS** | n=133 (channel,position) pairs; fit: a=774.3 ps, b=0.0 ps (expect b≈0), R²=0.953 |
| V8 | ✓ **PASS** | max \|bias\| across σ_true=[30,200] ps: 1.4% (< 5% → PASS, < 15% → CONCERN) |
| V9 | ✓ **PASS** | Dead events=0.000% (<1% expected); anomalous (>5σ)=0.060%; t_mean=611 ps, t_std=153 ps |

## Totales: 9 PASS / 0 CONCERN / 0 FAIL

## Constantes verificadas (leídas de fuente esta sesión)
- `kNEndSiPMs = 8` (DetectorConstruction.hh:37)
- `kNTopSiPMs = 70` (DetectorConstruction.hh:38)
- `kBarHalfX = 700.0` mm, `kBarHalfY = 30.0` mm, `kBarHalfZ = 5.0` mm
- `kEndPitch = 7.5` mm; END_L at x=-699.75 mm, END_R at +699.75 mm
- TOP: 35 SiPMs each side, step 20 mm, cx = -692..+692 mm
- `SCINTILLATIONTIMECONSTANT1 = 1.8 ns` (opsc-101.mac:11)
- `SCINTILLATIONYIELD = 10400/MeV` (opsc-101.mac:9)
- `ABSLENGTH = 160.0 cm` (DetConstr.cc override for OPSC-101)
- Emission peak = 408.8 nm (scntComp1.txt maximum)

## `VAL_STOP = True` — DETENTE aquí, espera OK de René antes de Fase 2.

*(Números generados programáticamente desde los sidecars CSV — ninguno fue tecleado.)*