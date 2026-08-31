# REFLECTOR_LABEL.md — Propuesta de etiquetas para el reflector en talk_v6

Contexto del problema:
- Volumen simulado: `VikuitiXXXLV`, material `G4_MYLAR`, superficie `dielectric_dielectric`, R=0.95
- Nombre de producto de referencia: Vikuiti 3M ESR (spec R≈0.98)
- El deck usa "Vikuiti ESR" en todo el documento, lo que puede inducir a pensar que R≈0.98

El texto sustituto propuesto es:
  **"specular reflector R = 0.95 (Vikuiti ESR as physical candidate)"**

---

## Inventario completo de ocurrencias en talk_v6.tex

### Grupo A — Ya consistentes (no requieren cambio)

Estas líneas ya muestran R=0.95 explícitamente o están en contexto de auditoría.
No es necesario cambiar el texto.

| Línea | Slide | Texto actual | Estado |
|-------|-------|-------------|--------|
| 33 | (macro) | `\newcommand{\Rvikuiti}{0.95}` | Correcto — define R=0.95 |
| 741 | Design Decision | `Reflector: \textbf{Vikuiti ESR} ($R=0.95$ constant in code)` | Consistente |
| 790 | Summary | `Vikuiti (R=0.95) recovers non-TIR photons (+23%)` | Consistente |
| 859 | A1 — Geant4 Surfaces | `Manufacturer Vikuiti ESR spec: $R \approx 0.98$` | Correcto — describe la discrepancia |
| 860 | A1 — Geant4 Surfaces | `The code comment says "R=0.98 Vikuiti ESR" but the actual constant set is 0.95.` | Correcto — nota (comentario corregido por C1) |
| 887 | A3 — Air-Gap Geometry | `Vikuiti reflectivity at air–Mylar interface ($R=0.95$)` | Consistente |
| 936 | B2 — 1-PE Anomaly | `absorbed by the "Vikuiti" metal surface` | Contexto histórico (GEN-1); entrecomillado ✓ |
| 937 | B2 — 1-PE Anomaly | `Vikuiti/TIR ratio of $\sim 10^3$` | Histórico GEN-1, no afirmación del modelo actual |
| 1192 | J1 — Provenance | `Reflector & Vikuiti ESR ($R=0.95$) \\` | Consistente |

### Grupo B — Usos narrativos cortos / etiquetas TikZ (prioridad baja)

"Vikuiti" funciona como nombre propio de la familia de reflectores sin implicar R=0.98.
El riesgo de confusión es bajo porque no acompañan valores numéricos.
Se puede dejar como está o sustituir puntualmente si el presentador prefiere uniformidad.

| Línea | Slide | Texto actual | Propuesta opcional |
|-------|-------|-------------|-------------------|
| 97, 156, 159 | The Bar… / TIR | Comentarios TikZ: `% Vikuiti top/bottom` | Sin cambio necesario |
| 121 | The Bar… | `\node{Vikuiti}` (etiqueta de diagrama) | `\node{Reflector}` |
| 171 | TIR Channels | `non-TIR $\to$ Vikuiti` | `non-TIR $\to$ reflector` |
| 191, 199 | Two Populations | `never hit Vikuiti` / `hit Vikuiti` | `never hit reflector` / `hit reflector` |
| 252 | Bounces | `did not bounce off the $H$-face Vikuiti` | `did not bounce off the $H$-face reflector` |
| 369 | GEN-3 Model | `Vikuiti only acts on photons that escape the bar.` | Sin cambio o `reflector only acts on...` |

### Grupo C — Ocurrencias que justifican la sustitución propuesta

Estas líneas nombran "Vikuiti ESR" en contexto cuantitativo o de especificación,
sin aclarar inmediatamente R=0.95. Son las candidatas prioritarias al texto sustituto.

| Línea | Slide | Texto actual | Propuesta de sustitución |
|-------|-------|-------------|-------------------------|
| **53** | Subtítulo | `EJ-230 + Vikuiti ESR \| Napkin First, Geant4 Second` | `EJ-230 + specular reflector R = 0.95 (Vikuiti ESR as physical candidate) \| Napkin First, Geant4 Second` **— o más corto:** `EJ-230 + Reflector R=0.95 (Vikuiti ESR) \| Napkin First, Geant4 Second` |
| **82** | The Bar… | `\textbf{Reflector:} Vikuiti ESR on all non-SiPM surfaces ($R = \Rvikuiti$).` | `\textbf{Reflector:} specular wrap, $R = \Rvikuiti$ (Vikuiti ESR as physical candidate).` |
| **250** | Bounces | `Napkin: $R=0.98$ (Vikuiti spec).` | Sin cambio — describe correctamente la especificación del fabricante |
| **148** | TIR | `recovered by the Vikuiti reflector (Population II)` | `recovered by the specular reflector (Population II)` _(o dejar como está)_ |

### Prioridad de cambio recomendada

1. **Alta** (líneas 53, 82): son las primeras apariciones de "Vikuiti ESR" en el deck y fijan la expectativa del lector sobre R. Sustituir aquí tiene el mayor impacto en claridad.
2. **Media** (líneas 148, 171, 191, 199, 252, 369): usos narrativos; no urgentes pero contribuyen a coherencia terminológica.
3. **Baja** (etiquetas TikZ 97, 121, 156, 159): en diagramas la etiqueta "Vikuiti" o "Reflector" no cambia la física. Decisión estética.

---

## Nota sobre la línea 860 tras el commit C1

El commit `71eeb20` (C1) corrigió el comentario en `DC.cc:310` de `dielectric_metal R=0.98` a
`dielectric_dielectric R=0.95`. La línea 860 del deck dice:

> `The code comment says "R=0.98 Vikuiti ESR" but the actual constant set is 0.95.`

Con el comentario ya corregido, esta frase describe un estado que ya no existe.
Propuesta de corrección al deck (requiere aprobación):

```latex
% Antes (línea 860):
The code comment says ``R=0.98 Vikuiti ESR'' but the actual constant set is 0.95.

% Después:
The code simulates $R=0.95$ (constant); Vikuiti 3M ESR manufacturer spec is $R\approx 0.98$.
```

---

*Informe de propuesta — ninguna línea del deck fue editada. Editar requiere aprobación.*
