# Comparación simulación vs test beam reportado

La posición longitudinal exacta de la medición de Gerardo sigue
**PENDIENTE de confirmar**. Por ello, el punto central de simulación se muestra
como referencia, no como una comparación posición-a-posición definitiva.

| observable | simulación Geant4 End | test beam (reportado) | nota |
|---|---:|---:|---|
| `sigma_LR`, centro (`x=0 mm`) | 138.48 ± 1.04 ps | 125 ps | Metodología temporal espejo; posición TB por confirmar. |
| `sigma_single`, centro (`x=0 mm`) | 97.92 ± 0.73 ps | ≈88 ps | `sigma_LR/sqrt(2)`; sin SPTR/jitter y sin corrección de walk. |
| `sigma_single`, punto más cercano a la medición | PENDIENTE | ≈88 ps | Falta confirmar la posición longitudinal de Gerardo. |
| `sigma_single`, centro, corte `NPE >= p25` | 94.64 ps | no reportado | Proxy del corte de pulse width; no equivalente sin `ToT(NPE)`. |
| `sigma_single`, centro, corte `NPE >= p50` | 93.05 ps | no reportado | Proxy del corte de pulse width; no equivalente sin `ToT(NPE)`. |
| MPV de NPE, centro | 360.58 NPE | no comparable | Fit Landau de `NPE_L+NPE_R` simulado. |
| MPV de ToT / pulse width | no disponible | PENDIENTE/TODO | Debe reportarse separado; requiere calibración `ToT(NPE)`. |

## Alcance

- Test beam: número reportado, no re-análisis de datos.
- Simulación: EJ-204, End readout, `+Y` wrapped, muón de 1 GeV, 10 000
  eventos por posición.
- La cadena SUM4 usa leading-edge con walk inherente, pero no aplica una
  corrección de time-walk. No incluye SPTR ni jitter electrónico.
