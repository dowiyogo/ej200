# GERARDO_DECISION_MEMO — EXEC_20 (reencuadrado, intrínseco/PDE only)
Fecha: 2026-06-20
Nota de René: la electrónica (SPTR, FastIC) queda diferida. Solo física intrínseca.

## 1. Reconciliación del baseline END-only

### Por qué nuestro END-only (882 ps @ centro) es peor que el de Gerardo
- Nuestro estimador: t_avg = (t_L + t_R)/2, requiere que AMBOS extremos disparen
- En x=0: END_L ve 0.37 PE/ev, END_R ve 0.38 PE/ev → casi ningún evento dispara ambos
- σ_END(t_avg, x=0) = 882 ps — timing puro de fluctuación Poisson

- Alternativa (un solo extremo, cerca de ese extremo):
  END_L_SUM4 @ x=-690 mm: σ ≈ 30.5 ps (near-end, EXEC_19 T4 canónico)
  Pero esto no es position-independent

### ¿Qué usa Gerardo exactamente? (confirmar)
  □ ¿Suma analógica de todos los SiPMs de un extremo?
  □ ¿Promedio (t_L+t_R)/2 o solo el mejor extremo?
  □ ¿Posición evaluada (centro vs extremo vs promedio del scan)?
  → Una vez confirmado, podemos hacer la comparación apples-to-apples.

## 2. Aporte cuantificado del TOP (intrínseco)

- σ_END-only (t_avg) @ x=0: 882 ps
- σ_TOP_nearest @ x=0: 105.7 ps
- σ_EndTop (GLS) @ x=0: 122.2 ps
- **Mejora @ x=0: 86.2%**
- **Mejora media (31 posiciones): 84.5%**

  Física: el TOP ayuda MUCHO al centro (END casi ciego) y POCO en los extremos
  (donde END ya da σ~30 ps). El aporte del TOP es posición-dependiente.

## 3. Decisión de topología SUM [DECISIÓN GERARDO]

Con sqrt_n binning + walk correction (EXEC_20 canónico):
  - dynamic nearest-4: σ = 69.4 ps (selecciona 4 con más luz)
  - id//4 fijo cluster:  σ = 76.4 ps (4 adyacentes)

  Nota: el 'co-localizado' aplica a END (misma cara, misma luz).
  El TOP es SIEMPRE distribuido (fila a 20 mm).
  Ambas opciones TOP son distribuidas; la diferencia es si el cluster es fijo o adaptativo.

## 4. Electrónica diferida (por decisión de René)
  SPTR ≈ 106 ps, FastIC ≈ 10 ps → NO en headline de esta campaña.
  Ver EXEC_19 para los números completos de σ_tot bajo Modelo A/B.
  Resumen: solo bajo Modelo B (FastIC+ digital por canal) se alcanzaría < 100 ps.