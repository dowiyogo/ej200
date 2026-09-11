# SPTR_PROVENANCE.md — Procedencia del valor de SPTR para AFBR-S4N66P024M

Fecha: 2026-09-01

---

## 1. Hallazgo: el datasheet no especifica SPTR

El datasheet **AFBR-S4N66P024M-DS105** (Broadcom, julio 2025) lista "Excellent SPTR
and CRT" en la sección de Features pero **no incluye ninguna fila de SPTR en la
tabla de Optical and Electrical Features**. El valor "~100 ps" que apareció en la
diapositiva I1 de `talk_v6.tex` no tenía procedencia documental y fue eliminado.

---

## 2. Referencia completa — medida publicada

**Artículo:**
> S. Lee, W.-S. Choong, R. Heller, J. W. Cates,
> "Timing Performance with Broadcom Metal Trench Silicon Photomultipliers",
> *IEEE Transactions on Radiation and Plasma Medical Sciences* **9**(4), 406–411 (2025).
> DOI: 10.1109/TRPMS.2024.3518479

**Dispositivo medido:**  
AFBR-S4N66P014M — $6\times6$\,mm² NUV-MT (Metal Trench SiPM).  
Ópticamente idéntico al **AFBR-S4N66P024M** (mismo active area, misma familia NUV-MT);
la diferencia es el encapsulado/conector, no el semiconductor.

**Condiciones de medida:**
- $V_\text{bias} = 48\,\text{V}$ ($\approx15.5\,\text{V}$ de sobretensión; $V_\text{bd}\approx32.5\,\text{V}$)
- Láser 408 nm, pulso ~15 ps, single-photon events
- Front-end LNHF $\approx2\,\text{GHz}$
- Temperatura $\approx20\,^\circ\text{C}$

**Resultados para el dispositivo de $6\times6$\,mm²:**

| Magnitud | Valor | Unidad |
|----------|-------|--------|
| SPTR intrínseco | $137\pm4$ | ps FWHM |
| SPTR de detector (incl. electrónica) | 172 | ps FWHM |

---

## 3. Conversiones FWHM → σ

Para una distribución gaussiana: $\sigma = \text{FWHM} / 2.355$.

| SPTR FWHM | σ equivalente |
|-----------|---------------|
| $137\pm4\,\text{ps}$ (intrínseco) | $\sigma\approx58\,\text{ps}$ |
| $172\,\text{ps}$ (detector) | $\sigma\approx73\,\text{ps}$ |

**Por qué esta conversión es obligatoria:** $\sigma_\text{TOP}=15.2\,\text{ps}$ es una
desviación estándar. Combinar en cuadratura una cantidad FWHM con una σ introduce un
error de factor $2.355$. FWHM debe convertirse a σ antes de cualquier suma en cuadratura.

---

## 4. Escalado con área — misma familia, mismas condiciones

Datos del mismo artículo (Lee et al. 2025), todos a $V_\text{bias}=48\,\text{V}$:

| Área activa | SPTR intrínseco | Capacitancia típica | Causa según autores |
|-------------|----------------|---------------------|---------------------|
| $2\times2\,\text{mm}^2$ | $45\pm1\,\text{ps FWHM}$ | $\approx160\,\text{pF}$ | — |
| $\approx4\times4\,\text{mm}^2$ | $55\pm1\,\text{ps FWHM}$ | $\approx580\,\text{pF}$ | capacitancia + transit-time skew |
| $6\times6\,\text{mm}^2$ | $137\pm4\,\text{ps FWHM}$ | $\approx1550\,\text{pF}$ | idem — escala con área |

La degradación del SPTR con el área es atribuida a: (1) capacitancia de drenaje que
alarga la respuesta del front-end, y (2) dispersión del tiempo de tránsito de la avalancha
a través del área del píxel.

**Nota de operación:** A sobretensión menor de $\approx40\,\text{V}$ el SPTR del dispositivo
de $6\,\text{mm}$ no era medible con el setup del paper. $V_\text{bias}=48\,\text{V}$
($\approx15.5\,\text{V}$ de sobretensión) está próximo al máximo absoluto del datasheet
($16\,\text{V}$ de sobretensión). El SPTR citado aplica a ese punto de operación;
a menor sobretensión el SPTR empeora de forma no cuantificada en este paper.

---

## 5. Punto abierto — propagación de SPTR en estadísticos de orden

La diapositiva I1 usaba la fórmula:

$$\sigma_\text{SPTR} \approx \frac{100\,\text{ps}}{\sqrt{k \cdot N_\text{active}}}$$

Esta fórmula es válida para el **promedio** de $k \cdot N$ observaciones independientes
(cada una con dispersión $\sigma_\text{SPTR}$). El estimador TOP es el **$k$-ésimo
estadístico de orden** de un conjunto de tiempos de llegada de fotones, **no un promedio**.

Para un estadístico de orden $X_{(k)}$, la propagación de la incertidumbre de timing de
cada SiPM ($\sigma_\text{SPTR}$) al estadístico final no sigue la raíz de $k \cdot N$;
la forma exacta depende de la distribución de los tiempos de llegada y requiere derivación
estadística formal (función de distribución del $k$-ésimo mínimo de una muestra).

**Estado:** pendiente de derivación. La fórmula de promedio fue eliminada de la diapositiva
y reemplazada por un enunciado cualitativo con nota explícita de que la propagación correcta
está por establecer. No se fabricó ninguna fórmula alternativa ni se recalculó $\sigma_t$.

---

## 6. PDE — consistencia con el datasheet

El archivo `data/sipm/AFBR-S4N66P024M_pde.txt` tiene en la cabecera:

```
# Source: AFBR-S4N66P024M-DS105, Figure 6 "PDE vs. Wavelength",
# typical operation at 12 V above breakdown, digitized from the datasheet plot.
# Datasheet cross-check: spectral range 250-900 nm; peak PDE 0.63 at 420 nm.
```

El valor implementado en el punto de emisión máxima de EJ-230 (~420 nm):
`0.630` — **consistente con el datasheet (63% típico a 420 nm, 12 V de sobretensión)**.
No se modifica nada.

Nota adicional del header: "PDE excludes optical crosstalk (~23%) and afterpulsing (<1%)."
Esto confirma que el crosstalk del 23% (declarado en el datasheet a 12 V de sobretensión)
no está incluido en la PDE usada por la simulación.
