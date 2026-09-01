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

Ningún defecto de fase óptica en las figuras de simulación (Fase 7 correcta).
La discrepancia R=0.98 (napkin) vs R=0.95 (código) es declarada explícitamente
en la diapositiva A1/A2 de talk_v6.tex como parte del análisis de consistencia.

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
