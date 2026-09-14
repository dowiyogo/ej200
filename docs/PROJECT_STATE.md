## 1. Establecido

- El guard destruía reflexiones internas totales: 3 750 479 de 4 930 224 encuentros barra→World, N=500; conteos, no fotones únicos. [EXEC_27: semilla/comandos](/home/rrios/REPORT_EXEC27_20260910.md:65), [sidecar](/home/rrios/ej200_exec26_20260909/build_exec27_20260910/run500/tables/bar_world_states.meta.json).
- Borde `dielectric_dielectric`+`polished`, entrada 0,98: retorno 13,371145 %, absorción 2,008837 % frente a 2 %; error NO REGISTRADO. Supervivencia→Fresnel: hipótesis de código; «1,2σ» no calibrado. [EXEC_26](/home/rrios/REPORT_EXEC26_PHASE1B_20260909.md:47), [EXEC_29](/home/rrios/REPORT_EXEC29_20260910.md:332), [sidecar](/home/rrios/ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/tables/prediction_comparison.meta.json).
- Borde `dielectric_metal`: retorno 98,000104 % para entrada 0,98; error NO REGISTRADO; muestra V1 END-only, no comparación emparejada con la anterior. [EXEC_30](/home/rrios/REPORT_EXEC30_20260910.md:219), [sidecar](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/cells/V1/reflection_panels.meta.json).
- Ganancia del reflector: END-only 1,09895552 ± 0,00852071; EndTop 1,01154746 ± 0,00830277; errores aproximados sin covarianza entre corridas. [EXEC_30](/home/rrios/REPORT_EXEC30_20260910.md:268), [sidecar](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/factorial_gains.meta.json).
- Encuentros aire→wrap por óptico iniciado: dieléctrico 0,48896774; metálico 10,66476630; error NO REGISTRADO; no son rebotes por centelleo. [EXEC_31, fuente EXEC_30](/home/rrios/REPORT_EXEC31_20260911.md:98), [D3](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/cells/D3/encounter_multiplicity.meta.json), [V1](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/cells/V1/encounter_multiplicity.meta.json).
- Balance terminal: residual 0,00000000 %, identidad contable en END-only D0/V1; atribución a EXEC_34C NO REGISTRADA. [EXEC_30](/home/rrios/REPORT_EXEC30_20260910.md:54), [rectificación EXEC_38](/home/rrios/REPORT_EXEC38_20260913.md:155), [sidecar](/home/rrios/ej200_exec26_20260909/build_exec30_20260910/g0_components.meta.json).
- Luz central: EndTop 397,152500 ± 2,409927 frente a END-only 1173,337500 ± 6,639454 pe/end; factor literal 2,95 NO REGISTRADO. [EXEC_29](/home/rrios/REPORT_EXEC29_20260910.md:48), [D1](/home/rrios/ej200_exec26_20260909/build_exec29_20260910/cells/D1/cell.meta.json), [D0](/home/rrios/ej200_exec26_20260909/build_exec29_20260910/cells/D0/cell.meta.json).
- Correlación del mismo extremo positiva en 88/88 casos; D0 izquierda ρ=0,262 ± 0,023; σ directa/σ normalizada Q=1,170 ± 0,022. [EXEC_38 V7](/home/rrios/REPORT_EXEC38_20260913.md:176), [sidecar](/home/rrios/exec38_20260913/V7_correlations.meta.json).
- Escalado temporal EndTop/END-only: extremos 1,235557 ± 0,030355 y 1,336580 ± 0,029581 frente a 1,719; errores aproximados, sin PASS/FAIL. [EXEC_37 C2](/home/rrios/REPORT_EXEC37_20260913.md:53), [sidecar](/home/rrios/exec37_20260913/same_end_results.meta.json).
- Cotas descriptivas: λ=422,25–601,67 mm y v=146,655–154,916 mm/ns, inferiores al nominal y a 189,742; errores por extremo publicados, modelos rechazados. [EXEC_38](/home/rrios/REPORT_EXEC38_20260913.md:79), [λ](/home/rrios/exec38_20260913/V3_attenuation.meta.json), [v](/home/rrios/exec38_20260913/V4_velocity.meta.json).
- Workers 1/4/12/24: diferencia exacta cero en conteos por evento y pe/end; EJ-204 central, N=2000, sin generalización a otras configuraciones. [EXEC_33](/home/rrios/REPORT_EXEC33_20260911.md:21), [sidecar](/home/rrios/exec33_20260911/S1/invocation.meta.json).
- Simetrías universales 0,5 %/0,3 %: NO REGISTRADAS como resultados establecidos; EXEC_36 no fijó umbral y reportó máximo especular 2,664 %. [EXEC_38 V6](/home/rrios/REPORT_EXEC38_20260913.md:159), [sidecar](/home/rrios/exec36_20260913/symmetry.meta.json).
- V1 pasa en 21/21 celdas, con razón producción/(yield·energía depositada) 0,999862–1,000178. [EXEC_42](../analysis/reports/exec42/REPORT_EXEC42_20260914.md), [sidecar](../analysis/reports/exec42/campaign/analysis/contract_results.json).
- V2-H1 pasa en 21/21 celdas; la fracción de escape al primer encuentro abarca 0,256163–0,263241. [EXEC_42](../analysis/reports/exec42/REPORT_EXEC42_20260914.md), [sidecar](../analysis/reports/exec42/campaign/analysis/contract_results.json).
- V5 emparejado pasa en 21/21 celdas comparando detecciones de frontera con incidentes observados y PDE superficial. [EXEC_42](../analysis/reports/exec42/REPORT_EXEC42_20260914.md), [sidecar](../analysis/reports/exec42/campaign/analysis/contract_results.json).
- M1 queda confirmado y M2 descartado: n_eff=1,58 en los tres materiales, con dispersión 4,44×10⁻¹⁶. [EXEC_43](../analysis/reports/exec43/REPORT_EXEC43_20260914.md), [sidecar](../analysis/reports/exec43/campaign/results.json).
- H3' cualitativa pasa sus cuatro condiciones preregistradas. [preregistro](../analysis/validation/EXEC44_PREREGISTRATION.md), [evaluación](../analysis/validation/EXEC44_H3_PRIME_EVALUATION.json).

## 2. Pendiente

- σ_t con electrónica validada: ausente; al menos cuatro versiones vivas, ambigüedad FWHM/σ y calibración sin resolver. [EXEC_34C](/home/rrios/REPORT_EXEC34C_20260912.md:186), [registro de electrónica](../analysis/sigma_t/ELECTRONICS_PROVENANCE.md).
- V1/V2/V5: NOT EVALUABLE en EXEC_38; estado histórico superado por la instrumentación de EXEC_40 y la grilla de EXEC_42. [EXEC_38](/home/rrios/REPORT_EXEC38_20260913.md:1), [EXEC_42](../analysis/reports/exec42/REPORT_EXEC42_20260914.md).
- ρ experimental y cuantificación de un modo común electrónico adicional: NO REGISTRADOS; EXEC_38 establece únicamente correlaciones de simulación. [EXEC_38 V7](/home/rrios/REPORT_EXEC38_20260913.md:167).
- Cuadratura: residuo indirecto 28,409–41,948 ps frente a referencias SPTR mayores; faltan equivalencia de observables y calibración, no constituye medición electrónica. [EXEC_37 C3](/home/rrios/REPORT_EXEC37_20260913.md:66).
- Sensor/PDE: simulación AFBR-S4N66P024M a 12 V frente a FBK NUV-MT 14M a OV=10 V; transferencia de respuesta sin establecer. [EXEC_36](/home/rrios/REPORT_EXEC36_20260913.md:484).
- Contrato dorado `ready_for_acceptance=false` en EXEC_38: estado histórico superado por EXEC_44. [EXEC_38](/home/rrios/REPORT_EXEC38_20260913.md:17), [contrato vigente](../analysis/validation/GOLDEN_REFERENCE_20260913.json).
- H3'' cuantitativa: faltan longitud de pista acumulada y coordenadas xyz de creación y primer encuentro; coste bruto estimado 56 bytes/fila, sin overhead ni compresión ROOT. [contrato](../analysis/validation/GOLDEN_REFERENCE_20260913.json).
- El índice SSLG4 es constante, n=1,58, mientras el PVT real varía aproximadamente 1,57–1,61 entre 380–450 nm; falta validar la traducción a otro motor. [EXEC_43](../analysis/reports/exec43/REPORT_EXEC43_20260914.md), [contrato](../analysis/validation/GOLDEN_REFERENCE_20260913.json).
- La fusión de la rama de validación a `main` queda pendiente de decisión explícita por el aumento de salida de aproximadamente 700 MB a 7 GB por celda N=10000. [archivo](../analysis/reports/EXEC40_43_ARCHIVE.json).

## 3. Retractado

- `CONTENT_AUDIT.md` atribuyó todo el exceso a recuperación por reflector; autor individual NO REGISTRADO; retirada documentada en EXEC_29, no en EXEC_34C. [EXEC_29](/home/rrios/REPORT_EXEC29_20260910.md:334).
- Informe EXEC_34C publicó σ gaussianos TOP; informe EXEC_35 declara los 21 inválidos bajo su regla preregistrada. [EXEC_35](/home/rrios/REPORT_EXEC35_20260912.md:7).
- Informe EXEC_34C: 17 signos negativos/4 positivos; p≈0,007 NO REGISTRADO. EXEC_35: 14/7, p nominal=0,18924713; no prueba ausencia de sesgo. [EXEC_34C](/home/rrios/REPORT_EXEC34C_20260912.md:124), [EXEC_35](/home/rrios/REPORT_EXEC35_20260912.md:342).
- «N=1 gana robustamente en 21 celdas» y retractación de todos los ganadores interiores: NO REGISTRADO; EXEC_35 conservó los índices y no hizo argmin robusto. [EXEC_35](/home/rrios/REPORT_EXEC35_20260912.md:344).
- Comparación EndTop del informe EXEC_36: violación de cota ausente en END-only central según EXEC_37; no validación global ni reescritura del guard anterior. [EXEC_37](/home/rrios/REPORT_EXEC37_20260913.md:49).
- Banda de cordura para SiPM aislado independiente de posición: EXEC_36 registra fallos a distintas distancias; retractación formal «mal formulada» NO REGISTRADA. [EXEC_36](/home/rrios/REPORT_EXEC36_20260913.md:51).
- «TOP destruyen guiado: 25–30 % frente a 2,24 %, refutado en EXEC_34C»: afirmación, autor y refutación NO REGISTRADOS; la razón TOP/centelleo documentada mezcla poblaciones. [EXEC_29](/home/rrios/REPORT_EXEC29_20260910.md:332).
- Informe EXEC_38 rechaza sus perfiles mono/biexponenciales, factor común entre materiales y recta temporal global; las cotas pasan solo descriptivamente. [EXEC_38](/home/rrios/REPORT_EXEC38_20260913.md:5).
- Genealogía inicial corregida por informe EXEC_31: `bd78211` introdujo metal; `f3a8062` seleccionó dieléctrico; autor individual del relato inicial NO REGISTRADO. [EXEC_31](/home/rrios/REPORT_EXEC31_20260911.md:96).
- H3 de invariancia material era físicamente inapropiada: la población de primeros encuentros excluye fotones absorbidos en bulk. Afirmada y fallida en EXEC_42; mecanismo refutado en EXEC_43 y sustitución formalizada en EXEC_44. [H3 original](../analysis/validation/EXEC42_H3_PREREGISTRATION.md), [resolución](../analysis/validation/GOLDEN_REFERENCE_20260913.json).

## 4. No reproducible

- `napkin.py:71`: `sim_npe={1310.8,941.0,701.3}`, diccionario literal sin procedencia de corrida; se reproduce la cita, no el resultado físico. [EXEC_32](/home/rrios/REPORT_EXEC32_20260911.md:237).
- `napkin_values.csv`: 52,1 ps rotulado x=0; macros históricos centro 53,36 y mínimo 52,07; discrepancia de etiqueta documentada, sin sustituir valores. [EXEC_32](/home/rrios/REPORT_EXEC32_20260911.md:237).
- `CLAIMS_AUDIT_20260910.md`: 778 de 967 filas `unverified`; son ocurrencias documentales, no mediciones independientes. [EXEC_29](/home/rrios/REPORT_EXEC29_20260910.md:322).
- `talk_v6.tex` en MSI: edición por tubería SSH sin commit; rama y PDF no comparten una fuente versionada allí. [EXEC_24_REPORT.md:131](../presentations/v6/EXEC_24_REPORT.md:131).

### 5. Estado del contrato de aceptación

- Pruebas activas: V1, V2-H1, V5 emparejado y H3' cualitativa; todas usan observables físicos portables.
- V1, V2-H1 y V5 pasan 21/21 celdas; H3' pasa 4/4 condiciones.
- No queda ninguna prueba activa en FAIL ni NOT EVALUABLE.
- `ready_for_acceptance=true`: las cuatro pruebas activas pasan con su cobertura requerida.
- H3 original se conserva retractada; no participa en el estado actual.
- H3'' queda como trabajo futuro, fuera del contrato de aceptación.
- Cambiar el estado requiere nueva evidencia con sidecar verificable o una modificación preregistrada del contrato.
