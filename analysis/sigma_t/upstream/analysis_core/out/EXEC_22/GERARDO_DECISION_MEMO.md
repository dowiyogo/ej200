# GERARDO_DECISION_MEMO — EXEC_22 (actualizado)
Fecha: 2026-06-21

## Estado del fix optico

Se aplico el fix TIR (dielectric_dielectric + REFLECTIVITY=0.95) en:
- ej200_endonly (branch exec21-optfix)
- ej200 EndTop build (branch exec22-endtop-optfix)

**Auditoria de velocidad (T2): PASA** — cero fotones superluminicos.

## Que funciona y que no

| | Antes del fix | Despues del fix (parcial) |
|--|--------------|--------------------------|
| npe_L(x=0) | 0.37 PE | ~0.4 PE (apenas cambia) |
| npe_L(x=-600) | ~0.01 PE | 26.1 PE (gran mejora!) |
| npe_L(x=-690) | 903 PE | 896 PE (sin cambio) |
| Estructura bimodal | ausente | 48% prompt @ x=-690 |
| σ_END(centro) | 882 ps | ~985 ps (aun no fisico) |

La TIR SI encendio (bimodal visible, perfil mas suave). El problema es que
los fotones no-TIR (~23% por rebote) siguen escapando sin recuperacion.
El centro sigue fallando.

## Por que el "aporte del TOP 85%" sigue invalido

La comparacion EndTop vs END-only de EXEC_20 usaba σ_END = 882 ps (no fisico).
Hasta que el END sea fisico (~40 PE, σ ~ decenas de ps), la metrica de aporte
no es significativa.

## Proximos pasos tecnicos (para decidir con Rene)

1. **Fix completo (geometria)**: anadir volumen de aire + panel Mylar alrededor de la barra.
   ETA: ~1-2 dias de implementacion en branch exec21-optfix.

2. **Re-simulacion EndTop**: una vez el fix completo este validado, re-simular la campana
   de 5000 ev/pos con el build corregido. El comando seria (NO ejecutar sin aprobacion):
   ```
   [ver EXEC_22_REPORT.md §T4 para el comando completo]
   ```

3. **Comparacion con Gerardo**: solo posible con datos re-simulados y fix completo.

## Decisiones pendientes (sin cambio respecto a EXEC_21)

1. Estimador END-only de Gerardo: t_avg vs single-end (confirmar)
2. Topologia SUM TOP: dynamic-4 vs id//4 (decision de Gerardo)
3. Electronica (SPTR/FastIC): diferida por decision de Rene
