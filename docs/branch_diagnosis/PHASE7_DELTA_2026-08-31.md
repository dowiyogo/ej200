# PHASE7_DELTA_2026-08-31

Diagnóstico acotado: cuánto cambia el resultado al pasar de Fase 4 (GEN-1) a Fase 7 (2T-Vikuiti).
Corrida pequeña, en scratch y **no autoritativa**. Propósito: decidir alcance de regeneración.

---

## 1. Hash compilado y verificación de condiciones ópticas (§2 del PROMPT)

**Binario Phase 7**: `build_end_vikuiti/ej200_bar_sim`
compilado desde `/home/rrios/ej200` en rama `feat/bar-end-vikuiti` HEAD `b00775081063b7d37cf786b2df45b8503fc4dabd`.

Fase 7 NO es ancestro de `feat/bar-end-vikuiti` (divergieron en `4cb2c23`).
Las 4 condiciones se verifican directamente en el código de la rama:

| Condición | Archivo:Línea | Valor | Estado |
|-----------|--------------|-------|--------|
| bar→air: `dielectric_dielectric`, `polished` | [src/DetectorConstruction.cc:312](../../src/DetectorConstruction.cc#L312), [src/Materials.cc:214](../../src/Materials.cc#L214) | `surf->SetType(dielectric_dielectric)` en `CreateBarSurface()`, colocada como `scintAirSurface` en `ScintToAirGap_*` border surfaces | ✓ |
| air→reflector: `dielectric_dielectric` (no `dielectric_metal`) | [src/Materials.cc:347](../../src/Materials.cc#L347) | `surf->SetType(dielectric_dielectric)` en `CreateBarSkinReflector()`, usada como `airReflectorSurface` en `AirGapToVikuiti_*` border surfaces | ✓ |
| air gap explícito con `RINDEX` declarado | [src/DetectorConstruction.cc:256](../../src/DetectorConstruction.cc#L256), [src/DetectorConstruction.cc:266](../../src/DetectorConstruction.cc#L266) | `airGapMat = worldMat (G4_AIR)`, `mpt->AddProperty("RINDEX", {2.0*eV,4.0*eV}, {1.0,1.0})` | ✓ |
| `GROUPVEL_air = c/n_bar` | [src/DetectorConstruction.cc:264](../../src/DetectorConstruction.cc#L264) | `const G4double vg = CLHEP::c_light / 1.58` → 18.97 cm/ns | ✓ |

Las 4 condiciones se cumplen. Se procede.

> **Nota**: El comentario en DC.cc:310 dice `"air→Vik : dielectric_metal R=0.98"` — es un comentario desactualizado.
> El código real usa `CreateBarSkinReflector()` que retorna `dielectric_dielectric`. No afecta la verificación.

---

## 2. Paridad de parámetros con Dataset F

Dataset F: `runs/t0minidaq_endtop_scan_20260618_204959`, commit `a0368c4` en
`feat/endtop-sslg4`, Fase 4 (GEN-1). Los parámetros se extraen de
[macros/t0minidaq_scan_5000/endtop_5000_31_x+690mm.mac](../../macros/t0minidaq_scan_5000/endtop_5000_31_x+690mm.mac)
(macro original del scan).

| Parámetro | Valor en Dataset F | Origen | Valor usado aquí |
|-----------|-------------------|--------|-----------------|
| Readout | `EndTop` | `endtop_5000_31_x+690mm.mac:7` | `EndTop` ✓ |
| Scintillator | `OPSC-101` (EJ-204) | `endtop_5000_31_x+690mm.mac:8` | `OPSC-101` ✓ |
| SiPM model | `AFBR-S4N66P024M` | `endtop_5000_31_x+690mm.mac:9` | `AFBR-S4N66P024M` ✓ |
| jitterSigma | `0 ns` | `endtop_5000_31_x+690mm.mac:11` | `0 ns` ✓ |
| Partícula | `mu-` | `endtop_5000_31_x+690mm.mac:12` | `mu-` ✓ |
| Energía | `1 GeV` | `endtop_5000_31_x+690mm.mac:13` | `1 GeV` ✓ |
| Ángulo | `0` (vertical, -Z) | `endtop_5000_31_x+690mm.mac:14` | `0` ✓ |
| Posición gun (por defecto) | `{0, 0, 60 mm}` dirección `-z` | [src/PrimaryGeneratorAction.cc:21-22](../../src/PrimaryGeneratorAction.cc#L21) | misma default ✓ |
| Cruzamiento bar | 10 mm (kBarHalfZ=5 mm, gun entra por +Z sale por -Z) | [src/DetectorConstruction.cc:34](../../src/DetectorConstruction.cc#L34) | idéntico ✓ |
| Posiciones | x=0, x=400, x=690 mm | Dataset F cubre −690 a +690 en 31 pos | 3 posiciones de ese rango ✓ |

**Binario Fase 4**: `build_t0minidaq/ej200_bar_sim`. Verificación de que es GEN-1:
`strings build_t0minidaq/ej200_bar_sim | grep G4LogicalSkinSurface` → `_ZN20G4LogicalSkinSurfaceC1E...`
(constructor llamado → crea skin surface en barLV, comportamiento Fase 4).
`strings build_end_vikuiti/ej200_bar_sim | grep G4LogicalSkinSurface` → solo `GetSurface` (no constructor).

Diferencia respecto al Dataset F original:
- Semillas: Dataset F usó semillas específicas por posición; aquí sin `random/setSeeds` (semilla del sistema).
- Threads: Dataset F usó 24 threads; aquí ejecución serial. Sin efecto en física.
- Eventos: Dataset F tiene 5000 ev/pos; aquí 300 ev/pos (propósito declarado del diagnóstico).

---

## 3. Tabla comparativa por posición y grupo

**Metodología**: ajuste gaussiano con ventana ±2σ_MAD alrededor del pico
(`σ_MAD = 1.4826·median(|x−median(x)|)`), bin width adaptativo (objetivo: ~20 bins en ±2σ_MAD).
**Todos los valores son intrínsecos**: sin SPTR, sin jitter electrónico
(`/sipm/jitterSigma 0 ns` en macro). N=300 eventos por posición y fase.

Script de análisis: `scratch_nonauthoritative/analyze_delta.py` (metodología idéntica
a `analyze_t0minidaq_endtop_corefit.py` — misma función `fit_gaussian_core`, mismo estimador MAD).

### 3.1 N_pe por grupo (hits detectados por evento)

| Grupo | x=0 P4 | x=0 P7 | ratio | x=400 P4 | x=400 P7 | ratio | x=690 P4 | x=690 P7 | ratio |
|-------|--------|--------|-------|---------|---------|-------|---------|---------|-------|
| END_L | 0.3 | 83.9 | — | 0.1 | 39.4 | — | 0.1 | 25.8 | — |
| END_R | 0.4 | 83.3 | — | 2.7 | 251.3 | — | 897.7 | 2278.9 | 2.5× |
| TOP_SUM4† | 114.2 | 1679.5 | **14.7×** | 125.2 | 1647.2 | **13.2×** | 85.1 | 1013.2 | **11.9×** |
| TOP_SUM8† | 114.2 | 1679.5 | **14.7×** | 125.2 | 1647.2 | **13.2×** | 85.1 | 1013.2 | **11.9×** |

† TOP_SUM4 y TOP_SUM8 usan el mismo grupo de canales (face_type=2, Top SiPMs), difieren solo en N.
La razón END_L y END_R Fase4/P7 no se calcula porque Fase 4 tiene Npe≈0 (GEN-1 elimina TIR → fotones no llegan a extremos).

### 3.2 Eficiencia por grupo (fracción de eventos con ≥N hits)

| Grupo | N | x=0 P4 | x=0 P7 | x=400 P4 | x=400 P7 | x=690 P4 | x=690 P7 |
|-------|---|--------|--------|---------|---------|---------|---------|
| END_L | 4 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 1.000 |
| END_R | 4 | 0.007 | 1.000 | 0.303 | 1.000 | 1.000 | 1.000 |
| TOP_SUM4 | 4 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| TOP_SUM8 | 8 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

La firma LOW_EFF del GEN-1 en END_L (reportada en ANALYSIS_CORE_FIT_REPORT.md) queda confirmada:
eff=0.000 en las 3 posiciones para Fase 4. Fase 7 recupera eff=1.000 en todas.

### 3.3 σ_t por grupo con error (ajuste gaussiano al core)

Unidades: picosegundos. Error = σ del parámetro en el ajuste (1σ).

#### x = 0 mm

| Grupo | N | P4 σ (ps) | P4 err | P7 σ (ps) | P7 err | Δσ (ps) | σ_comb | Δ/σ_comb |
|-------|---|-----------|--------|-----------|--------|---------|--------|---------|
| END_L | 4 | — (eff=0) | — | 161.3 | ±11.5 | — | — | — |
| END_R | 4 | — (eff≈0) | — | 141.1 | ±11.3 | — | — | — |
| TOP_SUM4 | 4 | 89.0 | ±5.4 | 3.2 | ±0.2 | **−85.8** | 5.4 | **15.9σ** |
| TOP_SUM8 | 8 | 98.4 | ±6.0 | 4.6 | ±0.3 | **−93.9** | 6.0 | **15.5σ** |

#### x = 400 mm

| Grupo | N | P4 σ (ps) | P4 err | P7 σ (ps) | P7 err | Δσ (ps) | σ_comb | Δ/σ_comb |
|-------|---|-----------|--------|-----------|--------|---------|--------|---------|
| END_L | 4 | — (eff=0) | — | 271.5 | ±25.9 | — | — | — |
| END_R | 4 | 1336.9 | ±751 | 70.3 | ±4.7 | **−1266.6** | 751 | 1.7σ† |
| TOP_SUM4 | 4 | 72.2 | ±5.0 | 11.8 | ±1.0 | **−60.4** | 5.1 | **12.0σ** |
| TOP_SUM8 | 8 | 81.9 | ±5.8 | 5.5 | ±0.3 | **−76.3** | 5.8 | **13.1σ** |

† END_R x=400 Fase4 tiene eff=0.303 y el ajuste solo tiene 91 eventos válidos → σ_err muy grande. El delta de 1266 ps es real pero el error es grande. Veredicto basado en eficiencia cualitativa, no en Δ/σ.

#### x = 690 mm

| Grupo | N | P4 σ (ps) | P4 err | P7 σ (ps) | P7 err | Δσ (ps) | σ_comb | Δ/σ_comb |
|-------|---|-----------|--------|-----------|--------|---------|--------|---------|
| END_L | 4 | — (eff=0) | — | 413.5 | ±44.1 | — | — | — |
| END_R | 4 | 25.7 | ±1.4 | 1.5 | ±0.1 | **−24.2** | 1.4 | **16.8σ** |
| TOP_SUM4 | 4 | 86.7 | ±5.6 | 28.6 | ±7.0 | **−58.0** | 8.9 | **6.5σ** |
| TOP_SUM8 | 8 | 108.2 | ±7.2 | 9.5 | ±0.6 | **−98.7** | 7.2 | **13.6σ** |

---

## 4. Veredicto por readout

### END_L

```
REGENERAR
```

Fase 4 tiene eff=0.000 en las 3 posiciones — ningún evento con ≥4 hits End-Left.
La comparación cuantitativa de σ_t no es posible porque Fase 4 no tiene señal.
Fase 7 recupera eff=1.000 en todas las posiciones.
Es un cambio cualitativo, no un delta pequeño.

### END_R

```
REGENERAR
```

- x=0: eff 0.007→1.000 (cualitativo)
- x=400: eff 0.303→1.000 con σ que colapsa de 1336 ps a 70 ps
- x=690: σ=25.7 ps → 1.5 ps a **16.8σ** (con nota: en Fase 4, END_R a x=690 funciona porque el cañón está a 10 mm del extremo; en las demás posiciones, END_R también está muerto)

La recuperación del End readout en las posiciones centrales (x=0, x=400) es la diferencia física crítica.

### TOP_SUM4

```
REGENERAR
```

| Posición | Δσ (ps) | Significancia |
|----------|---------|--------------|
| x=0 | −85.8 | **15.9σ** |
| x=400 | −60.4 | **12.0σ** |
| x=690 | −58.0 | **6.5σ** |

En todas las posiciones el delta supera el error combinado de forma clara.

### TOP_SUM8

```
REGENERAR
```

| Posición | Δσ (ps) | Significancia |
|----------|---------|--------------|
| x=0 | −93.9 | **15.5σ** |
| x=400 | −76.3 | **13.1σ** |
| x=690 | −98.7 | **13.6σ** |

---

## 5. Recomendación sobre alcance de la regeneración

**Regenerar el escaneo completo de 31 posiciones (End y Top).**

Justificación:
1. El efecto es sistemático en todas las posiciones probadas. No hay ninguna posición donde
   el delta sea pequeño o compatible con cero. El rango probado (x=0 = centro, x=400 = media
   longitud, x=690 = extremo) cubre los tres regímenes físicos del bar.
2. End readout pasa de estar muerto (eff≈0) a eff=1.000 en todas las posiciones centrales.
   Las conclusiones sobre End de EXEC_16-22 no tienen ninguna validez física bajo Fase 4.
3. Top readout mejora σ_t entre 58 y 99 ps (6-16σ). Los valores absolutos de EXEC_16-22
   para σ_t(TOP) son incorrectos por este margen.
4. El ratio EJ-230/EJ-204 no se puede asumir conservado: la dependencia de σ en N_pe
   no es lineal, y N_pe cambia ×12-15× entre fases. Se necesita una corrida EJ-230
   equivalente antes de afirmar si la razón se preserva o no (hipótesis a comprobar,
   no conclusión de este diagnóstico).

**Esta recomendación es solo de texto. No se ejecuta nada aquí.**

---

## 6. Ruta del scratch y declaración de no autoridad

```
Scratch: /tmp/rrios/claude-1002/-home-rrios-ej200/b814b8c0-d020-491e-907a-c940267e5fe6/
         scratchpad/phase7_delta_scratch_nonauthoritative/
```

Contenido del scratch:
- `run_diagnostic.sh` — script de simulación (300 ev × 3 pos × 2 fases)
- `phase4/photon_hits_ej204_x{0,400,690}.root` — 3 archivos ROOT Fase 4
- `phase7/photon_hits_ej204_x{0,400,690}.root` — 3 archivos ROOT Fase 7
- `analyze_delta.py` — script de análisis (metodología MAD-Gaussian, misma que corefit)
- `delta_results.json` — resultados completos en JSON

**Declaración**: Este diagnóstico es no autoritativo y no está registrado en
`orchestrator/datasets.py`. No reemplaza al escaneo completo. Los valores numéricos
son orientativos para estimar el tamaño del efecto a 300 eventos; no son los valores
de publicación. El scratch es temporal y puede ser eliminado sin consecuencia.

---

## Notas de auditoría

- **Sesión**: 2026-08-31 (generado en esta sesión de diagnóstico)
- **Ramas involucradas**: `docs/branch-diagnosis-2026-08-31` (documento),
  `feat/bar-end-vikuiti` HEAD `b007750` (código Phase 7)
- **Commit hash compilado**: verificado con `strings build_end_vikuiti/ej200_bar_sim`
- **Dataset F commit (`a0368c4`)** confirmado GEN-1 por `git merge-base --is-ancestor f39b84c a0368c4 → True`
  (ver DATA_AUDIT_2026-08-31.md §Dataset F)
- **Commit `5576687`** (Fase 7 original en `feat/endtop-sslg4`) no es ancestro de `feat/bar-end-vikuiti`.
  Las 4 condiciones físicas de Fase 7 se reimplementaron independientemente en `feat/bar-end-vikuiti`
  (verificadas en código, §1 de este documento).
