# TALKV6_CONSISTENCY.md — Auditoría de consistencia interna de talk_v6

Auditoría realizada: 2026-08-31  
Rama auditada: `feat/bar-end-vikuiti` (HEAD `b007750`)  
Fichero auditado: `analysis/presentation_v6/talk_v6.tex` (1225 líneas)  
Alcance: solo lectura — git, código fuente, scripts de análisis  

---

## A1 — Tabla de proveniencia de figuras

> Para cada figura objetivo: qué script la generó, de qué datos, y en qué generación óptica cae.

| Figura | Script generador | CSV/datos intermedios | Directorio de datos ROOT | Generación óptica | Commit hash |
|--------|-----------------|----------------------|--------------------------|-------------------|-------------|
| `fig_sigma_t_x.pdf` | `talks/napkin_first_principles/fig_gen.py` | `analysis/optim/phase_ab_optimal.csv` | `results/scan_end_vikuiti/` | GEN-3 física¹ (End-only) | **sin rastreo git** |
| `fig_npe_x.pdf` | `talks/napkin_first_principles/fig_gen.py` | `analysis/optim/phase_ab_optimal.csv` | `results/scan_end_vikuiti/` | GEN-3 física¹ (End-only) | **sin rastreo git** |
| `figM1_mat_sigma_end.pdf` | `analysis/presentation_v4/scripts/fig_materials.py` | `analysis/optim/phase_ab_optimal.csv` | `results/scan_end_vikuiti/` | GEN-3 física¹ (End-only) | **sin rastreo git** |
| `figM2_mat_npe.pdf` | `analysis/presentation_v4/scripts/fig_materials.py` | `analysis/optim/phase_ab_optimal.csv` | `results/scan_end_vikuiti/` | GEN-3 física¹ (End-only) | **sin rastreo git** |
| `fig_end_mscan.pdf` | `analysis/optim/root_best_est/new_analysis_plots.C` | datos embebidos (de `phase_ab_kscan_full.csv`) | `results/scan_end_vikuiti/` | GEN-3 física¹ (End-only) | **sin rastreo git** |

**¹ Nota sobre la clasificación GEN-2 vs GEN-3:**  
El `CONFIGURATION_AUDIT.md` etiqueta `phase_ab_optimal.csv` (y por tanto figM1/M2/M3 y las figuras de napkin) como datos "GEN-2" (sin cámara de aire). Sin embargo:

- `phase_ab.py` lee explícitamente de `results/scan_end_vikuiti/` (variable `ENDVIK`)
- El binario `build_end_vikuiti/ej200_bar_sim` contiene la cadena `AirGapToVikuiti_`, confirming Fase 7 con cámara de aire
- El scan `scan_end_vikuiti` se ejecutó el 2026-08-16 (verificado en `results/scan_end_vikuiti/scan_run.log`), dos días después del commit de Fase 7
- Ergo: todos los datos de `scan_end_vikuiti` son **Fase 7 = física GEN-3** (con cámara de aire)

El `CONFIGURATION_AUDIT.md` confunde "geometría End-only (sin SiPMs TOP)" con "GEN-2 (sin cámara de aire)". Son dos distinciones ortogonales.

**Comandos de verificación:**
```bash
# Confirmar que phase_ab.py usa scan_end_vikuiti
grep -n "ENDVIK\|scan_end_vik" analysis/optim/phase_ab.py

# Confirmar que el binario tiene cámara de aire (Fase 7)
strings build_end_vikuiti/ej200_bar_sim | grep AirGapToVikuiti_

# Confirmar fecha del scan
head -5 results/scan_end_vikuiti/scan_run.log
```

**Por qué no hay commits:** Los directorios `analysis/`, `talks/`, `results/` y `presentations/` no están rastreados en git (aparecen como `??` en `git status`). No existe hash de commit para ninguno de los scripts o figuras auditados.

---

## A2 — Afirmación B1: atribución del delta 50.5 → 53.7 ps

**Afirmación en el deck** (talk_v6.tex, línea 924):  
> "GEN-2 END-only gives σ_END=50.5 ps (smaller than GEN-3 53.7 ps) **because TOP SiPMs in GEN-3 intercept ~3 ps worth of early photons** that would otherwise have reached the END face."

**Tabla B1** (líneas 906–921):

| Parámetro | GEN-2 | GEN-3 |
|-----------|-------|-------|
| Cámara de aire | **No** | **Sí** (0.10 mm) |
| Reflector | skin surface, R=0.95 | border surface, R=0.95 |
| SiPMs TOP | No | Sí (N=4,8,14,20) |
| σ_END | 50.5 ps | 53.7 ps |

**Problema identificado:**  
La comparación GEN-2 → GEN-3 cambia **dos variables simultáneamente**:
1. Cámara de aire: 0 → 0.10 mm (+ volumen Mylar 0.05 mm)
2. SiPMs TOP: apagados → encendidos

La atribución en el pie de la diapositiva asigna el delta de +3.2 ps íntegramente a la intercepción de fotones por los SiPMs TOP. Esto no puede verificarse con la comparación mostrada porque los cambios de geometría (cámara de aire) y de configuración (SiPMs TOP) son simultáneos.

**Evidencia parcialmente mitigante:**  
El `CONFIGURATION_AUDIT.md` registra que GEN-3 End-only (con cámara de aire, sin SiPMs TOP) también da σ_END=50.5 ps, lo que implicaría que la cámara de aire sola no modifica σ_END. Si esto es correcto, el delta 50.5→53.7 ps sería efectivamente atribuible a los SiPMs TOP. Sin embargo, esta afirmación no tiene un run de referencia explícito en el deck que lo demuestre directamente; es una aserción del CONFIGURATION_AUDIT, no un dato mostrado al lector.

**Nota adicional (C1, línea 973):**  
El pie de C1 dice "σ_t values above are from END-only geometry (GEN-2/GEN-3 no-TOP)". La equivalencia entre GEN-2 y GEN-3 sin TOP se asume (50.5 ps para ambos) pero no se muestra como comparación explícita en ninguna diapositiva.

**Veredicto para A2:** La afirmación **es plausible** según el CONFIGURATION_AUDIT, pero **no está demostrada internamente** en el deck. Para blindarla, hace falta un run GEN-3 (con cámara de aire) en geometría End-only y mostrarlo junto al run GEN-3 con TOP.

---

## A3 — R=0.95 (código) y R=0.98 (napkin): ubicaciones exactas

| Valor | Ubicación | Línea | Contexto |
|-------|-----------|-------|---------|
| R=0.95 | `src/Materials.cc` | 354 | `const std::vector<G4double> refl = {0.95, 0.95};` en `CreateBarSkinReflector()` |
| R=0.95 | `src/Materials.cc` | 357 | `mpt->AddProperty("REFLECTIVITY", energy, refl);` — surface MPT efectiva |
| R=0.95 | `src/Materials.cc` | 332 | Comentario: "(2) angle < theta_c → non-TIR; REFLECTIVITY=0.95 models Mylar/ESR substrate" |
| R=0.98 | `talks/napkin_first_principles/fig_gen.py` | 8 | `fig_survival_RN.pdf — R^N vs N (analytic, R=0.98)` |
| R=0.98 | `analysis/presentation_v6/FINAL_NUMBERS.md` | tabla | "Λ_refl (H=10mm face, R=0.98): 404.6 mm — napkin, R=0.98" |
| R=0.98 (erróneo) | `src/DetectorConstruction.cc` | 310 | `//   air→Vik : dielectric_metal R=0.98 (Vikuiti 3M ESR)` — comentario desactualizado |

**El deck reconoce explícitamente la discrepancia** en el frame A1/A2 (líneas 830–867):
- "The code comment says 'R=0.98 Vikuiti ESR' but the actual constant set is 0.95" (paráfrasis del deck)
- El frame confirma que R=0.95 es el valor simulado y R≈0.98 es la spec del fabricante

**Defecto adicional en DC.cc línea 310:**  
El comentario dice `dielectric_metal R=0.98` pero el código en la línea siguiente (313) llama a `CreateBarSkinReflector()` que retorna `dielectric_dielectric`. Son dos errores en un solo comentario: tipo de superficie equivocado **y** valor de R equivocado.

---

## A4 — "Vikuiti ESR" vs material simulado

**Afirmación del deck** (talk_v6.tex, línea 82):  
> "Reflector: Vikuiti ESR on all non-SiPM surfaces (R=0.95)"

**Jerarquía de volúmenes en código** (verificado en `src/DetectorConstruction.cc`):

| Nivel | Nombre en código | Material | Superficie óptica |
|-------|-----------------|----------|------------------|
| Volumen bar | `BarLV` | scint (`CreateScintillator()`) | — |
| Cámara de aire | `AirGapXXXLV` | `G4_AIR` (worldMat, DC.cc:256) | border `CreateBarSurface()` → `dielectric_dielectric` |
| Reflector | `VikuitiXXXLV` | **`Materials::CreateMylar()`** → G4_MYLAR (DC.cc:257,372) | border `CreateBarSkinReflector()` → `dielectric_dielectric`, R=0.95 |

**Contradicciones identificadas:**

1. **Nombre vs material**: El volumen se llama `VikuitiXXXLV` (DC.cc:372) pero el material es `G4_MYLAR` (no existe material "Vikuiti ESR" en Geant4 ni en el código).

2. **Comentario DC.cc:310 doblemente incorrecto**:
   - Dice `dielectric_metal` → código usa `dielectric_dielectric` (DC.cc:313)
   - Dice `R=0.98` → código usa `R=0.95` (Materials.cc:354)

3. **Comentario Materials.cc:205 incorrecto**:
   > "NOTE: CreateBarSurface() and CreateMylar() are not used in the current geometry"  
   Pero DC.cc:312 llama `CreateBarSurface()` y DC.cc:257 llama `CreateMylar()`. Ambas funciones **sí se usan**.

4. **Spec del fabricante vs simulación**:
   - Vikuiti 3M ESR especificación: R≈0.98 (típico a λ>400 nm)
   - Código simula: R=0.95 constante en todas las longitudes de onda

**Qué se simula realmente:**  
Una superficie `dielectric_dielectric` con REFLECTIVITY=0.95 constante entre el volumen de cámara de aire (G4_AIR) y un volumen de 0.05 mm de G4_MYLAR. El G4_MYLAR tiene RINDEX=1.65 y ABSLENGTH=1 µm (Materials.cc:265,267), lo que lo hace efectivamente opaco: los fotones que atraviesan la superficie entran en el Mylar y se absorben. La reflexión ocurre por la propiedad de superficie (R=0.95), no por las propiedades del bulk del material.

El deck identifica este material como "Vikuiti ESR" a lo largo de la presentación. La identificación es aspiracional (Vikuiti ESR es el producto de referencia que se pretende modelar) pero la simulación usa R=0.95, no R=0.98.

---

## Resumen: afirmaciones del deck que no se sostienen tal como están escritas

| # | Slide | Afirmación | Problema |
|---|-------|-----------|---------|
| **I** | B1 (línea 924) | "difference because TOP SiPMs in GEN-3 intercept ~3 ps worth of early photons" | La comparación GEN-2→GEN-3 cambia dos variables (cámara de aire + SiPMs TOP). La atribución es plausible pero no demostrada en el deck. |
| **II** | C1 (línea 973) | "GEN-2/GEN-3 no-TOP" como equivalentes | No se muestra un run de comparación directa GEN-2 vs GEN-3 sin TOP que valide la equivalencia. |
| **III** | DC.cc comentario (línea 310) | "air→Vik : dielectric_metal R=0.98" | Doble error: superficie es `dielectric_dielectric` (Materials.cc:313,329) y R=0.95 (Materials.cc:354). |
| **IV** | Materials.cc comentario (línea 205) | "CreateBarSurface() y CreateMylar() not used in current geometry" | Ambas funciones SÍ se usan (DC.cc:257,312). Comentario obsoleto. |
| **V** | Deck general | "Vikuiti ESR" como reflector | El material simulado es G4_MYLAR. La superficie usa R=0.95, no R≈0.98 (spec del fabricante). |
| **VI** | A1 (implícito) | figM1/M2/M3 descritas como "GEN-2/GEN-3 no-TOP" | Los datos provienen de `scan_end_vikuiti` con `build_end_vikuiti` Fase 7 (física GEN-3 con cámara de aire), no de un build sin cámara de aire. El CONFIGURATION_AUDIT confunde geometría End-only con generación GEN-2. |

---

## Apéndice: cadena de datos para las cinco figuras

```
results/scan_end_vikuiti/   (5000 ev × 13 pos × 3 mats, Fase 7, build_end_vikuiti)
    └─► analysis/optim/phase_ab.py
            └─► analysis/optim/phase_ab_optimal.csv
                    ├─► talks/napkin_first_principles/fig_gen.py
                    │       └─► fig_sigma_t_x.pdf, fig_npe_x.pdf
                    └─► analysis/presentation_v4/scripts/fig_materials.py
                            └─► figM1_mat_sigma_end.pdf, figM2_mat_npe.pdf

    └─► analysis/optim/phase_ab.py
            └─► analysis/optim/phase_ab_kscan_full.csv
                    └─► analysis/optim/root_best_est/new_analysis_plots.C
                            └─► fig_end_mscan.pdf  (datos embebidos al compilar)
```

Todos los directorios (`analysis/`, `talks/`, `results/`) están sin rastrear en git. No hay hashes de commit disponibles para ninguno de estos artefactos.

---

*Auditoría read-only. Ningún archivo fue modificado. Ninguna simulación fue ejecutada.*
