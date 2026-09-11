# MORNING_REPORT — EXEC_16 — 2026-06-20 (autonomous run)

## 1. Veredicto

**COMPLETADO-CON-FLAGS**

Todos los niveles 0–4 ejecutados. Beamer compilado (16 páginas, sin errores ni overfull boxes).
Flags son físicamente esperados (ver §4); ninguno indica corrupción de datos.

---

## 2. Fix de identidad aplicado

**Causa raíz (§A):** `face_type` ternario {0,1,2}; `TOP_IDS=range(16,36)` erróneo.

**Fix aplicado — leído de `DetectorConstruction.hh:37-38`:**
- `K_END = 8` (kNEndSiPMs)
- `N_TOP = 70` (kNTopSiPMs)
- `BOUNDARY = 16`
- classify_gid(): g<8→END_L, 8≤g<16→END_R, g≥16→TOP
- Gate de partición: **PASSED** (86 canales únicos, sin colisiones)
- global_id observado: 0..85 — conforme

---

## 3. Decisiones autónomas (política §B)

- **TOP_IDS** corregido a range(16,86) — leído del header
- **WRAPPING** label: "reflector_skin_surface (Mylar-equivalent)" — wrap volume removido (DetConstr.cc:89)
- **τ_d = 1.8 ns** — de opsc-101.mac:11, no heredado de EJ-230
- **Dataset no registrado en datasets.py** — gap 2×2 conocido; DATA_DIR leído directamente
- **ORCH_DIR**: datasets.py está en {ORCH_DIR}/analysis/ (no en raíz)
- **31 posiciones × 5000 evt** — conteo directo; paso medio 46 mm
- **Canal nearest TOP**: máx ⟨Npe⟩ programático → gid=50 @ x=0 (cx=−12 mm ≈ beam; coherente)
- **NPE(x) double-exp**: convergencia falló (maxfev=5000); datos graficados sin ajuste
- **Beamer overfull**: 3 frames corregidos con height/keepaspectratio; compilación limpia

---

## 4. Flags emitidos

| Flag | Detalle | Acción |
|------|---------|--------|
| CURVE-ABORT (>30% sigma_outlier) | END_L_SUM4, END_R_SUM4, END_L_SUM8, END_R_SUM8 | Excluidas de figuras; conservadas en CSV; no fabricados resultados aguas abajo |
| double-exp fit failed | TOP NPE(x) atenuación | Datos graficados sin ajuste; anotado |

**Física de las curvas END abortadas:** σ(END) = 30 ps @ x=−690 mm (near-end) → 909 ps @ x=0 → >1000 ps @ x=+690 mm. Es físicamente correcto — distribuciones anchas con pocos PE lejos del extremo. No es fallo de ajuste. Ajustar SIGMA_BAND_PS=(5,2000) para visualizar curva completa.

---

## 5. Números destacados (desde results.json + CSVs — no tecleados)

### Resolución intrínseca @ x=0 mm

| Estimador | σ_fit [ps] | bootstrap ± [ps] |
|-----------|-----------|-----------------|
| **TOP_SUM4** | **72.6** | ±1.6 |
| **TOP_SUM8** | **74.4** | ±1.7 |
| END_COMBINED_SUM4 | 866 | ±376 |
| END_COMBINED_SUM8 | 883 | ±101 |

**TOP_SUM4 y TOP_SUM8 < 100 ps en TODAS las 31 posiciones del scan.**

### TOP_SUM4 σ(x) (posiciones representativas)

| x [mm] | σ [ps] |
|--------|--------|
| −690 | 81.9 |
| −350 | 70.8 |
| 0 | **72.6** |
| +350 | 67.0 |
| +690 | 72.6 |

### Beneficio EndTop

| x [mm] | σ_END-only | σ_EndTop | ratio |
|--------|-----------|----------|-------|
| −690 | 609 ps | 295 ps | 0.484 |
| 0 | 883 ps | 520 ps | **0.589** |
| +690 | 660 ps | 264 ps | **0.400** |

**Ratio < 1 en TODAS las posiciones → TOP siempre ayuda en EJ-204.**

### END near-end performance

- END_L_SUM4 @ x=−690 mm: **σ = 30.0 ps** (intrínseco)
- TOP_SUM4 @ x=−690 mm: 81.9 ps

### Resolución espacial

- **v_eff = 15.02 cm/ns** (slope de ⟨Δt⟩(x); referencia ≈ 15.5 cm/ns — consistente)

### Estimado con SPTR + FastIC (solo cuadratura informativa)

√(72.6² + 106² + 10²) ≈ **128.8 ps**

### Contexto SHiP

Objetivo: 100 ps (preferible 50 ps). Intrínseco: ✓ (TOP_SUM4/SUM8).
Con SPTR+FastIC: ≈ 129 ps estimado — sobre el objetivo. EXEC_02b cierra este gap.

---

## 6. Rutas de salida

| Tipo | Ruta |
|------|------|
| Figuras L0–L4 (PDF+PNG) | `analysis_core/out/EXEC_16/L0..L4/` — 11 figuras |
| CSVs sidecars | `L0/inventory.csv`, `L1/npe_per_face.csv`, `L2/sigma_t_all_estimators.csv`, `L3/endonly_vs_endtop.csv + spatial_resolution.csv` |
| results.json | `out/EXEC_16/results.json` |
| Beamer PDF | `beamer/EXEC_16/EXEC_16_endtop_ej204.pdf` (16 páginas, 8.3 MB) |
| Beamer .tex | `beamer/EXEC_16/EXEC_16_endtop_ej204.tex` (generado por build_beamer.py) |

---

## 7. Hashes de commit y tags (analysis_core master)

| CP | Hash | Tag de rollback |
|----|------|----------------|
| init | 892fb8b | EXEC_16-pre-cp0 |
| CP1 (L0+L1) | ad34864 | EXEC_16-pre-cp1, EXEC_16-pre-cp2 |
| CP2 (L2) | 8d5e398 | EXEC_16-pre-cp3 |
| CP3 (L3+L4) | 63c2481 | EXEC_16-pre-cp4 |
| CP4 (Beamer+report) | *(este commit)* | EXEC_16-cp4 |

**git push (NO ejecutado):** `git push origin master`

---

## 8. Para revisar René

1. **Curvas END abortadas**: Ajustar `SIGMA_BAND_PS = (5, 2000)` en exec16_endtop_ej204.py y re-correr para ver la curva completa END_L/R vs x (30 ps→1000 ps).

2. **TOP_SUM4 vs SUM8 anomalía**: σ(SUM4)=72.6 < σ(SUM8)=74.4 — pequeña pero contraintuitiva. SUM8 puede agregar canales lejanos con fotones tardíos. Confirmar con análisis de gid vs ⟨t_1⟩.

3. **double-exp NPE(x)**: Ajuste no convergió. Re-intentar con p0 más fino.

4. **Clustering**: id//4 fijo vs ventana móvil → **decisión de Gerardo**.

5. **HOOK_WALK**: Corrección time-walk THR_k — stub reservado en código (línea ~203 de exec16_endtop_ej204.py).

6. **EXEC_02b**: Reintroducir SPTR≈106 ps + FastIC≈10 ps en cuadratura.

7. **Cierre 2×2**: EJ-230 EndTop pendiente (scan alta estadística).

---

**CAVEAT: Resolución INTRÍNSECA (sin SPTR≈106 ps ni jitter FastIC≈10 ps; sin gate). Wrapping: reflector skin surface on BarLV (Mylar-equivalent). τ_d=1.8 ns, τ_r=0.7 ns, λ_abs=160 cm (OPSC-101/EJ-204).**
