# GERARDO_DECISION_MEMO — EXEC_21 actualizado
Fecha: 2026-06-21
**Para**: Gerardo Vásquez | **De**: EXEC_21 autonomous diagnosis

## URGENTE: el baseline END-only de la sim NO es físico

Los resultados de EXEC_20 ("aporte del TOP 85%") comparaban:
- σ_TOP ≈ 106 ps (físico)
- σ_END ≈ 882 ps (NO físico — déficit de luz ×300)

**No presentar el "aporte del TOP 85%" hasta que END sea físico.**

---

## Qué encontramos: culpable identificado

El código usa `dielectric_metal` para la superficie del skin de la barra.
Esto elimina la Reflexión Total Interna (TIR). Consecuencias:

| Modelo | <Npe>_END @ centro | Predicción teórica | Real detector |
|--------|-------------------|-------------------|---------------|
| dielectric_metal, R=0.98 (actual) | 0.37 PE | 0.04 PE | × |
| dielectric_metal, R=1.0 (perfecto) | ~1.2 PE | ~1.2 PE | × |
| dielectric_dielectric + Mylar (real) | pendiente | ~749 PE | ~40 PE |

Incluso con R=1.0 la absorción en trayecto largo (164 rebotes × 20mm = 328 cm >> λ_abs=160 cm) limita a ~1 PE.

**La única solución es TIR (dielectric_dielectric + envolvente Mylar explícita).**

---

## Estado del fix (branch exec21-optfix)

**Fix parte 1 (DONE, commit af5ddb7):** `dielectric_metal` → `dielectric_dielectric`
- Resultado @ x=-690: σ_END: 609 → 259 ps; ¡estructura bimodal visible! (47.5% fotones rápidos)
- Resultado @ x=0: σ_END: 882 → 862 ps (sin cambio — fotones no-TIR escapan sin Mylar)

**Fix parte 2 (PENDIENTE):** Añadir volumen Mylar explícito
- Los fotones a ángulo < θ_c = 39.3° no hacen TIR → actualmente escapan → perdidos
- Con envolvente Mylar (R=0.97-0.99): esos fotones se recuperan
- Implementación: geometría explícita (barra + gap de aire + paneles Mylar)
- ETA: ~1-2 días de implementación + validación

---

## Consecuencia para T3 (aporte del TOP)

El "aporte del TOP 85%" del EXEC_20 no es válido porque el baseline END es incorrecto.
Una vez que el fix Mylar esté completo:
- σ_END(centro) debería bajar de 882 ps a O(100 ps) si el modelo es correcto
- La comparación con Gerardo (85-90 ps de dos extremos) se podrá hacer de verdad
- El "aporte del TOP" se recalculará con un baseline END físico

---

## Decisiones pendientes (sin cambio)

1. Estimador END-only de Gerardo: t_avg o single-end? — confirmar
2. Topología SUM TOP: dynamic-4 vs id//4 — decisión de Gerardo
3. Electrónica (SPTR/FastIC): diferida por decisión de René
