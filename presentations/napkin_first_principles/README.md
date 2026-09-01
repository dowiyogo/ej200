# napkin_first_principles — junio–agosto 2026

**Estado:** OBSOLETO (absorbido por v6, sección "Napkin First, Geant4 Second")
**Fase óptica:**
- Modelo analítico (napkin): R=0.98 (especificación del fabricante Vikuiti 3M ESR) — intencional;
  el napkin usa la spec, no el valor simulado.
- Figuras de simulación (fig_gen.py → fig_sigma_t_x.pdf, fig_npe_x.pdf): Fase 7 —
  datos de `results/scan_end_vikuiti/` vía `analysis/optim/phase_ab_optimal.csv` (2026-08-16).
**Dataset:**
- Analítico: cálculo cerrado (R^N survival, N_pe promedio)
- Simulación: `analysis/optim/phase_ab_optimal.csv` ← `results/scan_end_vikuiti/` (Fase 7)

## Defectos conocidos

**Reflectividad R = 0.95 en figuras de simulación (corregida en `c7acb7a`):** El modelo
analítico (napkin) usó R=0.98 (spec Vikuiti 3M ESR) — intencional y correcto. Las figuras
de simulación (`fig_gen.py`, `fig_gen.C`) usaron datos de `scan_end_vikuiti` producidos
con `REFLECTIVITY=0.95` — error de configuración. La discrepancia fue declarada en
`talk_v6.tex:A1/A2`.

Corregido en `Materials.cc:354` (commit `c7acb7a`, 2026-09-01). Nuevas simulaciones con
R=0.98 serán directamente comparables con el napkin sin factor de escala correctivo.
Las figuras actuales del deck corresponden a R=0.95 y son históricas.
Ver `docs/branch_diagnosis/REFLECTIVITY_CHANGE.md §5`.

## Por qué se conserva

Código fuente del razonamiento "napkin first": derivación analítica de la supervivencia
de fotones y el estimador de timing, independiente de la simulación Geant4. Sirve como
referencia de primer principios y base pedagógica del deck v6.

## Contenido

| Archivo | Función |
|---------|---------|
| `tex/talk.tex` | Versión original del deck napkin |
| `tex/talk_v3.tex` | Versión 3 del deck (más completa) |
| `fig_gen.py` | Genera fig_sigma_t_x.pdf, fig_npe_x.pdf desde phase_ab_optimal.csv |
| `fig_gen.C` | Equivalente en ROOT C++ |
| `napkin.py` | Cálculo analítico (survival, estimador) |
| `values/` | Constantes y macros LaTeX compartidas |
