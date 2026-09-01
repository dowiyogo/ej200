# presentations/ — Índice de presentaciones

| Deck | Fecha | Estado | Fase óptica | Una línea |
|------|-------|--------|-------------|-----------|
| [v6/](v6/) | 2026-09-01 | **ACTUAL** | Fase 7 (`5576687`) | Presentación principal: Fase 7, GEN-2/GEN-3, best estimate BLUE |
| [v5/](v5/) | agosto 2026 | OBSOLETO | Fase 7 (`5576687`) | Estimador BLUE END+TOP; EJ-230, R=0.95; supersedido por v6 |
| [v4/](v4/) | agosto 2026 | OBSOLETO | Fase 7 (`5576687`) | Trayectoria napkin → Geant4; IQR+BLUE; supersedido por v5/v6 |
| [napkin_first_principles/](napkin_first_principles/) | junio–agosto 2026 | OBSOLETO | Analítico R=0.98 / simulación Fase 7 | Derivación analítica de primer principios (survival R^N, estimador); absorbido por v6 |
| [best_est_2026-08-17/](best_est_2026-08-17/) | 2026-08-17 | OBSOLETO | Fase 7 (`5576687`) | Best estimate con ROOT verificado; σ_END=53.68 ps, σ_BLUE=15.21 ps |
| [optim_2026-08-17/](optim_2026-08-17/) | 2026-08-17 | OBSOLETO | Fase 7 (`5576687`) | Optimización Phase II: END+TOP sparse combinado |
| [optim_2026-08-16/](optim_2026-08-16/) | 2026-08-16 | OBSOLETO | Fase 7 (`5576687`) | Primera iteración optimización α/β (END-only) |
| [summary_2026-08-16/](summary_2026-08-16/) | 2026-08-16 | OBSOLETO | Fase 7 (`5576687`) | Síntesis multimaterial multimodo del día del scan |
| [end_vikuiti_2026-08-16/](end_vikuiti_2026-08-16/) | 2026-08-16 | OBSOLETO | Fase 7 (`5576687`) | Template End-only+Vikuiti+air gap; placeholders; absorbido por v6 |
| [exec14/](exec14/) | 2026-06-15/16 | OBSOLETO | Fuera de cronología principal | Comparación EJ-204 vs EJ-230; feat/endonly-mylar @ `3ae135f` |
| [sim_status_hi_2026-06-08/](sim_status_hi_2026-06-08/) | 2026-06-08 | OBSOLETO | INDETERMINADO | Estado simulación alta estadística; predates Fase 3 (2026-06-10) |
| [sim_status_2026-06-08/](sim_status_2026-06-08/) | 2026-06-08 | OBSOLETO | INDETERMINADO | Primer informe de estado; predates Fase 3 (2026-06-10) |

## Cronología óptica de referencia

| Fase | Commit | Fecha | Descripción |
|------|--------|-------|-------------|
| Fase 3 | `1f6aca1` | 2026-06-10 | Near-GEN-1: `dielectric_metal`, sin air gap — DESCARTADA |
| Fase 4 | `f39b84c` | 2026-06-18 | GEN-1: skin surface `dielectric_metal`, elimina TIR — INVÁLIDA |
| Fase 5 | `fc5d8b0` | 2026-06-21 | SK-DD: skin surface `dielectric_dielectric` R=0.95 |
| Fase 6 | `4cb2c23` | 2026-08-12 | 2T-MYLAR: air gap explícito, bug de velocidad de grupo — DATOS INCORRECTOS |
| Fase 7 | `5576687` | 2026-08-14 | 2T-MYLAR corregido: `GROUPVEL_air = c/n_bar = 18.97 cm/ns` — VÁLIDA |

**Defectos por fase (usar al declarar OBSOLETO):**
- **GEN-1 (Fase 4):** skin `dielectric_metal` elimina TIR; ~0.37 pe/evento en END. Todo resultado de esa fase es inválido.
- **Ventana Fase 6 (2026-08-12/14):** velocidad de grupo incorrecta en air gap; tiempos de propagación erróneos.
- **Fase 3:** `dielectric_metal` sin air gap; internamente consistente pero óptica descartada.

## Notas

- Los PDF compilados se conservan si existen (registro de lo presentado).
- Los archivos auxiliares LaTeX (`*.aux`, `*.nav`, `*.toc`, `*.snm`, `*.fls`, `*.fdb_latexmk`, `*.synctex.gz`) están excluidos por `.gitignore`.
- El framework de análisis de exec14 está en `analysis/exec14/`.
- La documentación QA de exec14 está en `docs/branch_diagnosis/exec14_qa/`.
