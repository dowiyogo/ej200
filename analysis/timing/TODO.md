# analysis/timing/ — Trabajo pendiente

Lista en orden de valor para el análisis de resolución temporal.  
No es un plan de ejecución; es un registro de qué falta y qué datos son necesarios.

---

## 1. Medir el SPTR con los datos de EOS (alta prioridad)

**Qué resuelve:** convertir el rango derivado 85–92 ps RMS (barrido de intensidad láser,
midpoints de rangos) en una medida propia con incertidumbre formal.

**Datos disponibles:**
- ROOT del barrido de intensidad: `/eos/experiment/ship/DSTD/`
- Decodificador: `https://gitlab.cern.ch/cvaldivi/test-beam-ship-april-2026/`
- Los archivos tienen ToA y ToT evento a evento (a diferencia del barrido actual con
  valores brutos del osciloscopio sin resolución evento a evento)

**Procedimiento:**
1. Para cada nivel de intensidad: histograma de ToA → σ_ToA; espectro de carga → N_pe
   (valor central, no rango)
2. Ajuste OLS de σ²(N) = SPTR²/N + σ_elec² con tres o más puntos, con propagación de
   incertidumbre (errores en σ_ToA desde el ajuste gaussiano; errores en N_pe desde el
   espectro de carga)
3. Resultado: SPTR y σ_elec con errores formales; comparación con Lee 2025 y con el rango
   derivado aquí

**Limitaciones actuales del barrido de osciloscopio:**
- N_pe reportados como rangos (no valores centrales)
- σ_ToA brutos (sin ajuste evento a evento)
- Límite instrumental del osciloscopio de ~21 ps que impide resolver el suelo electrónico

---

## 2. Verificar la fracción óptima del dCFD (prioridad media)

**Qué resuelve:** determinar si 14% es el óptimo para el AFBR-S4N66P024M en este montaje,
o si el óptimo real es 3–6% (Cattaneo, SiPMs de 3×3 mm²) o mayor.

**Datos disponibles:**
- `analyze_dCFD_fraction.C` implementa el barrido de fracción (ya existe en `analysis/timing/`)
- Datos de Fase 7 END-only: `results/scan_end_vikuiti/` o `results/scan_end_vik_sparse_top_v2/`

**Procedimiento:**
- Correr `analyze_dCFD_fraction.C` sobre datos de Fase 7 con el modelo `fastic_measured`
  (τ_r=2 ns, τ_f=3 ns), SPTR del rango propio (85–92 ps)
- Registrar la salida como output versionado en `results/` o en `docs/`

**Referencias:**
- Cattaneo et al. (arXiv:1402.1404 §II-F): óptimo 3–6% para SiPMs de 3×3 mm²
- Derenzo et al. (PMC nihms596188): óptimo sube con mayor jitter de tránsito (plausible
  para 6×6 mm², no verificado)

---

## 3. Confirmar τ_r/τ_f con el FastIC+ en la cadena completa (prioridad media)

**Qué resuelve:** confirmar el par (2, 3) ns con medida propia con el FastIC+ en cadena,
en lugar de tomarlo de las capturas del osciloscopio sin el ASIC en la señal.

**Datos disponibles:**
- Capturas de osciloscopio del banco (1 GHz, 6.4 GS/s) caracterizan el pulso libre y
  el pulso tras la red de cancelación polo-cero → dan (τ_fall=55 ns, τ_rise=2 ns / τ_fall=3 ns)
- El FastIC+ (TIA ~20 Ω + discriminador de corriente) no debería alterar el flanco de
  subida de 2 ns, pero no se ha verificado directamente

**Procedimiento:**
- Captura de pulso del SiPM con el FastIC+ en la cadena y el osciloscopio en paralelo
- Ajuste exponencial al flanco de subida y al de caída del pulso capturado
- Si τ_r medido difiere de 2 ns: actualizar `pulse_models.py` → `fastic_measured.tau_r_ns`

---

## 4. Resolver la ambigüedad FWHM/sigma del 106 ps (prioridad baja, requiere persona)

**Qué resuelve:** determinar si el valor de 106 ps inyectado en las simulaciones Geant4
de la colaboración es σ o FWHM.

**Por qué importa:** si es FWHM (110 ps FWHM ≡ ~47 ps σ), el jitter está sobreestimado
en un factor ~2.355 respecto a la convención estándar. Si es σ, es coherente con el rango
85–92 ps derivado del banco.

**Las notas del grupo describen** "100–120 ps FWHM" y "~106 ps RMS" como equivalentes.
Son incompatibles; alguien mezcló convenciones al fijar el parámetro.

**Cómo resolverlo:** consultar con Gerardo la fuente original del valor; identificar si
vino de una medida propia (en qué unidades) o de una referencia publicada.

**Sin resolver:** no usar 106 ps en análisis propios sin aclarar la convención.
Ver `docs/branch_diagnosis/ELECTRONICS_PARAMETERS.md §Referencia de la colaboración`.

---

*Registrado: 2026-09-01. Actualizar a medida que las tareas se completen.*
