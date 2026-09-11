# analysis/timing/ — Estimadores de tiempo de llegada

Scripts de reconstrucción de timestamp y análisis de resolución temporal σ_t.
Todos leen el TTree `sipm_hits` de los archivos ROOT de Fase 7 (ver §2.1 de
`docs/branch_diagnosis/ANALYSIS_CONSOLIDATION_AUDIT.md`).

---

## Tabla de estimadores

| Estimador | Script(s) | Modelo de pulso τ_r / τ_f | σ SPTR incluido | Dispositivo de referencia | Fuente |
|-----------|-----------|--------------------------|-----------------|--------------------------|--------|
| Waveform + dCFD 14% | `sipm_waveform_dcfd.py`, `sipm_waveform_dcfd.cpp` | Elegido por `--pulse-model` (ver tabla de modelos abajo; ningún par mezcla τ de fuentes distintas) | Sí — σ = `--transit-sigma-ps` por fotón (sin default; ver Lee 2025) | Varía según modelo | `pulse_models.py` (única fuente de verdad) |
| dCFD directo sobre hits | `analyze_dCFD.py`, `analyze_dCFD_fraction.C` | Sin modelo de pulso (opera sobre timestamps de fotones) | No (aplica jitter de electrónica fijo 30 ps; **ver nota**) | — | FastIC+ ASIC (sin cita formal; 30 ps no es jitter del ASIC) |
| FPT (primer fotón) | `resolution_vs_x_FPT.py`, `analyze_basic.py` | Sin modelo de pulso | No | — | Estadístico de orden puro (mínimo de hit times) |
| k-ésimo fotón | `analyze_dCFD_5thPhoton.C` | Sin modelo de pulso | No | — | 5.º hit time ordenado |
| SUM4 leading-edge | `analysis/exec14/engine.py` (en main, fuera de este directorio) | 0.5 ns / 5.0 ns (SPR de ej200_endonly) | No | SiPM de ej200_endonly (no identificado); **incomparable con sipm_waveform_dcfd** | `common.py` vía `/home/reriosto/SHiP/ej200_endonly/analysis/exec07` |
| FPT por SiPM individual | `fpt_vs_n_profile.C`, `fpt_vs_n_profile_batch.C` | Sin modelo de pulso | No | — | Perfil ⟨t_n⟩ vs n promediado sobre eventos |
| Múltiples estimadores (barrido) | `TimeMarkScan.C` | Sin modelo de pulso | No | — | FPT, media de k primeros, media pesada |
| Ranking de SiPMs | `SiPMRankingScan_v2.C`, `SiPMRankingScan_RMS.C`, `SiPMRankingScan_coreSigma.C` | Sin modelo de pulso | No | — | FPT por canal individual |

### Modelos de pulso disponibles en `sipm_waveform_dcfd.py`

`--pulse-model` es requerido; sin él el script aborta. Fuente de verdad: `pulse_models.py`.

| `--pulse-model` | τ_r | τ_f | Sistema que describe | Fuente | Estado |
|----------------|-----|-----|----------------------|--------|--------|
| `penarodriguez_shortened` | 2.0 ns | 3.0 ns | Montaje con circuito de acortamiento (cancelación polo-cero) | arXiv:2411.16710 §4 | PUBLICADO (otro montaje) |
| `broadcom_intrinsic` | NO MEDIDO | 55.0 ns | SiPM intrínseco AFBR-S4N66P024M | DS105 + arXiv:2411.16710 §2 | DATASHEET (τ_r NO MEDIDO) |
| `fastic_measured` | NO MEDIDO | NO MEDIDO | Montaje real FastIC+ de este banco | PENDIENTE | **NO MEDIDO** |

El par `(engine.py: 0.5/5.0 ns)` y el par `(sipm_waveform_dcfd.py: elegido por modelo)` **no
son comparables entre sí**: describen montajes distintos y posiblemente dispositivos distintos.
Los resultados de σ_t de estos dos estimadores no deben combinarse ni compararse directamente.

**σ SPTR (para `--transit-sigma-ps`):**  
Sin default desde 2026-09-01. Lee et al. (IEEE TRPMS 2025) mide σ_intrínseco ≈ 58.2 ps
(137 FWHM / 2.355) y σ_detector ≈ 73.0 ps (172 FWHM / 2.355) para AFBR-S4N66P014M a
OV ≈ 15.5 V. Este banco opera a OV = 10 V: **valores de Lee son cotas optimistas**.
El script imprime esta advertencia en cada ejecución.

**Nota sobre 30 ps en `analyze_dCFD.py`:**  
El valor de 30 ps era la resolución temporal esperada de un detector completo diferente
(teja EJ-228 con SiPM FBK de 2 mm²), no el jitter del FastIC+. Ver
`docs/branch_diagnosis/ELECTRONICS_PARAMETERS.md`.

---

## Descripción de scripts

### `sipm_waveform_dcfd.py` / `sipm_waveform_dcfd.cpp`

Reconstruye la waveform del SiPM sumando respuestas SPE (Lv et al. Eq. 1) con pulso
Peña-Rodríguez `v(t) = A(1−e^{−t/τ_r})e^{−t/τ_f}`, aplica dCFD al 14%.
**No migrar** el binario `sipm_waveform_dcfd` (83 KB, compilado en W2).

Uso:
```bash
python analysis/timing/sipm_waveform_dcfd.py photon_hits_run*.root --face 0
```

### `analyze_dCFD.py`

dCFD al 14% aplicado directamente a los timestamps de fotones (sin reconstrucción de
waveform). Produce timestamp por evento y cara; añade jitter gaussiano de electrónica
(30 ps FastIC+). Requiere solo `event_id`, `face_type`, `time_ns`, `gun_x_mm`.

### `analyze_dCFD_fraction.C`

Barrido de fracción dCFD (parámetro configurable, default 14%). Análisis de cara 0 (end
left) por defecto.

### `analyze_dCFD_5thPhoton.C`

Trigger por el 5.º fotón ordenado. Cara 0. Sin dCFD; sirve para comparar con el umbral
de 5 pe.

### `analyze_basic.py`

FPT básico (mínimo de `time_ns` por evento por cara). Sin reconstrucción de waveform.
Compatibilidad completa con Fase 7. Punto de partida para análisis nuevos.

### `resolution_vs_x_dCFD.py`

σ_t vs posición x usando dCFD al 14%. Jitter electrónica 30 ps FastIC+. Produce figura
PDF de σ_t vs x.

### `resolution_vs_x_FPT.py`

σ_t vs posición x usando FPT (primer fotón). Sin electrónica. Útil como referencia para
cuantificar el coste del estimador.

### `fpt_vs_n_profile.C` / `fpt_vs_n_profile_batch.C`

Perfil ⟨t_n⟩ vs n (número de fotón llegado) por SiPM individual. Visualización del
timing prompt photon budget por canal.

### `fpt_vs_n_profile_batch_slides.C` + `fpt_manifest_to_beamer.py` + `png_to_beamer.py`

Pipeline de generación de presentación Beamer desde un scan batch. Ejecutar desde
`analysis/timing/`:
```bash
root -l -b -q fpt_vs_n_profile_batch_slides.C
```
Genera un manifiesto TSV y lo convierte a `.tex` via `fpt_manifest_to_beamer.py`.

### `TimeMarkScan.C`

Barrido de modos de time mark (FPT, media de k primeros, media pesada) vs posición x.

### `SiPMRankingScan_v2.C` / `SiPMRankingScan_RMS.C` / `SiPMRankingScan_coreSigma.C`

Ranking de SiPMs individuales por resolución temporal (FPT). `v2`: ajuste gaussiano
al núcleo; `RMS`: RMS puro; `coreSigma`: σ del núcleo (sin colas).

---

## Datos disponibles vs geometría

Ver también `docs/branch_diagnosis/ANALYSIS_CONSOLIDATION.md §5`.

Los scripts que hacen referencia a `N_TOP_SIPMS` o IDs de top SiPMs asumen la geometría
con la que fueron desarrollados. Main tiene `kNTopSiPMs = 70` (IDs 16–85). Los scripts
en este directorio fueron escritos en una rama con `kNTopSiPMs = 20` (IDs 16–35).

**No modificar** los scripts al migrarlos. Documentado como geometría divergente; la
propuesta de derivación en runtime está en `ANALYSIS_CONSOLIDATION_AUDIT.md §2.5b`.
