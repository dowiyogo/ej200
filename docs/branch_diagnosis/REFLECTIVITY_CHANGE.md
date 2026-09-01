# REFLECTIVITY_CHANGE.md — Corrección R=0.95→0.98 en CreateBarSkinReflector

Fecha: 2026-09-01  
Commit §2: `c7acb7a`

---

## 1. Cambio realizado

`src/Materials.cc:354` — `CreateBarSkinReflector()`:

```diff
- const std::vector<G4double> refl = {0.95, 0.95};
+ const std::vector<G4double> refl = {0.98, 0.98};  // Vikuiti 3M ESR spec (R≈0.98, 3M product sheet)
```

Comentario de cabecera actualizado (`Materials.cc:332`):

```diff
- //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate.
+ //   (2) angle < theta_c → non-TIR; REFLECTIVITY=0.98 (Vikuiti 3M ESR spec) models the reflective wrap.
```

**Conservado sin cambios:** `G4_MYLAR` como material volumétrico, `VikuitiXXXLV` como nombre
del volumen. La investigación §1 confirmó que con `REFLECTIVITY` explícita en superficie
`dielectric_dielectric`, las ecuaciones de Fresnel no se calculan — el RINDEX del volumen Mylar
(1.65) no entra en ninguna decisión óptica. Solo el valor de la surface MPT importa.

---

## 2. Justificación

`REFLECTIVITY=0.95` era un error de configuración. La hoja de datos del fabricante (Vikuiti 3M ESR)
especifica R≈0.98. El valor 0.95 nunca fue una decisión de diseño documentada; los comentarios
originales del código citaban 0.98 mientras el valor codificado era 0.95.

En `G4OpBoundaryProcess`, para `dielectric_dielectric` con `REFLECTIVITY` explícita en el MPT
de la superficie, Geant4 usa esa probabilidad directamente — sin calcular Fresnel. El RINDEX
del volumen Mylar no se consulta para la decisión de reflexión/transmisión.

---

## 3. Datasets invalidados

Todos los datasets producidos con la Fase 7 fueron simulados con `REFLECTIVITY=0.95`.
Necesitan re-simulación para ser válidos con R=0.98.

| Dataset | Fecha | Estado |
|---------|-------|--------|
| `scan_end_vikuiti` (End-only GEN-3) | 2026-08-16 | Invalidado — R=0.95 |
| `scan_end_vik_sparse_top_v2` (END+TOP) | 2026-08-16/17 | Invalidado — R=0.95 |

Los números publicados en `presentations/v6/` (σ_END=53.68 ps, σ_BLUE=15.21 ps) corresponden
a R=0.95 y son provisionalmente válidos como referencia interna hasta la re-simulación.

---

## 4. Dirección del efecto (cualitativo — sin estimación de magnitud)

El cambio afecta únicamente a los fotones non-TIR (ángulo θ < 39.3°); los fotones TIR
(θ > 39.3°) se reflejan al 100% independientemente de `REFLECTIVITY`.

| Variable | Dirección esperada |
|----------|-------------------|
| N_pe | Mayor (más fotones non-TIR recuperados) |
| σ_t | Menor (más fotones → menor incertidumbre estadística en el estimador) |

No se estima la magnitud del cambio. La cuantificación requiere re-simulación con R=0.98.

---

## 5. Alineación napkin / simulación tras la corrección

El modelo analítico (`presentations/napkin_first_principles/`) siempre usó R=0.98
(especificación del fabricante). La simulación Fase 7 usaba R=0.95, generando la
discrepancia documentada en `TALKV6_CONSISTENCY.md` y declarada en la diapositiva
A1 de `talk_v6.tex`.

Tras el commit `c7acb7a`, nuevas simulaciones usarán R=0.98. La comparación
napkin-vs-Geant4 que el deck describe como diagnóstico de discrepancia será,
con los nuevos datos, una validación directa sin factor de corrección.

---

## 6. REFLECTOR_LABEL.md — superado (§3a)

`docs/branch_diagnosis/REFLECTOR_LABEL.md` proponía sustituir "Vikuiti ESR" en el deck
por "specular reflector R=0.95 (Vikuiti ESR as physical candidate)" para señalar la
discrepancia R=0.95 (código) vs R≈0.98 (spec).

Esa discrepancia queda eliminada por el commit `c7acb7a`. El deck puede mantener
"Vikuiti ESR (R=0.98)" sin discrepancia respecto al código.
`REFLECTOR_LABEL.md` queda **superado en su recomendación principal**.

Quedan pendientes para cuando el deck se regenere con datos R=0.98:
- `\newcommand{\Rvikuiti}{0.95}` (línea 33 de `talk_v6.tex`) → actualizar a `0.98`
- Ocurrencias hardcoded `R=0.95` en líneas 741, 790, 887, 1192 → actualizar
- Líneas 859–860 (descripción de la discrepancia como defecto) → obsoletas; reescribir

El documento `REFLECTOR_LABEL.md` permanece como registro histórico con el inventario
de ocurrencias; marcado como superado en su cabecera.

---

## 7. Refinamiento óptico pendiente: PolishedESR_LUT (§3b)

`src/external/OPSimTool/` implementa el finish `PolishedESR_LUT` (y variantes
`RoughESR_LUT`, `PolishedESRGrease_LUT`, `RoughESRGrease_LUT`), que encapsulan la
BRDF medida del Vikuiti 3M ESR con dependencia angular completa (modelo LUT de Geant4).
Contrasta con el modelo actual: `REFLECTIVITY=0.98` escalar uniforme, que aplica la
misma probabilidad de reflexión para todos los ángulos de incidencia.

**Relevancia física:** los fotones non-TIR que llegan al reflector cubren un rango
angular amplio; los fotones con múltiples rebotes TIR llegan en ángulos muy rasantes
(θ ≈ θ_c). La BRDF medida del ESR tiene dependencia angular documentada; un modelo
escalar promedia sobre esa variación y puede introducir sesgo diferente por población
de fotones.

**Implicaciones de adoptar PolishedESR_LUT:**
1. `surf->SetFinish(polished)` → `surf->SetFinish(PolishedESR_LUT)` en
   `CreateBarSkinReflector()` (`Materials.cc:349`).
2. Retirar la propiedad `REFLECTIVITY` del MPT de la superficie: con finish LUT,
   Geant4 consulta la tabla medida directamente; `REFLECTIVITY` explícita sobrescribiría el LUT.
3. Verificar que la geometría lateral del bar (superficie plana pulida, air gap 0.10 mm)
   corresponde a la configuración de la muestra con que se construyó el LUT de OPSimTool.
4. Re-simulación completa.

**Por qué queda fuera de alcance ahora:** cambiar el finish LUT altera el modelo
físico más allá de una corrección de parámetro escalar y requiere validación
independiente de la correspondencia LUT–geometría. La corrección R=0.95→0.98 alinea
la simulación con la especificación escalar del fabricante y es el paso mínimo
justificado por datos disponibles.

Se documenta como **candidato a evaluar** en una fase futura de validación óptica avanzada;
no es una tarea pendiente de este ciclo.
