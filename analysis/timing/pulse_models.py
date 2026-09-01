"""
pulse_models.py — Única fuente de verdad para configuraciones de pulso SPE del SiPM.

Cada entrada es un par (tau_r, tau_f) medido en UNA MISMA configuración física.
Está PROHIBIDO combinar tau_r de una configuración con tau_f de otra: los dos
parámetros están acoplados por C_J (capacitancia de la microcelda) y describen
sistemas distintos si se toman de montajes distintos.

Referencia de diseño del pulso (Peña-Rodríguez, arXiv:2411.16710 §2):
    v(t) = A * (1 - exp(-t/tau_r)) * exp(-t/tau_f)
    tau_r = R_S * C_J   (resistencia de microcelda, subida rápida)
    tau_f = R_Q * C_J   (resistencia de quenching, caída lenta)
Como R_S << R_Q, el tau_r intrínseco debe ser << tau_f — décimas de ns, no 2 ns.

Origen probable del error en el código anterior:
    El Apéndice A de arXiv:2411.16710 trae un ejemplo de código con Rt=2e-9 y Ft=50e-9.
    Alguien tomó el Rt del ejemplo (que corresponde al montaje acortado de §4)
    y lo combinó con el Ft=55 ns del datasheet DS105 (valor intrínseco),
    obteniendo un par (2 ns, 55 ns) que no describe ninguna configuración real.
"""

FWHM_TO_SIGMA = 1.0 / 2.355  # conversión FWHM → sigma gaussiana; explícita, no precalculada

PULSE_MODELS = {
    "penarodriguez_shortened": {
        "tau_r_ns": 2.0,
        "tau_f_ns": 3.0,
        "source": "arXiv:2411.16710 §4 — circuito de acortamiento con cancelación polo-cero",
        "note": (
            "NO es la respuesta intrínseca del SiPM. "
            "Es la respuesta del montaje de lectura de ese trabajo (amplificador RF clase A "
            "con cancelación polo-cero para evitar apilamiento). "
            "tau_f intrínseco del AFBR-S4N66P024M = 55 ns (DS105); "
            "el acortamiento lo reduce a 3 ns en ese montaje."
        ),
    },
    "broadcom_intrinsic": {
        "tau_r_ns": None,
        "tau_f_ns": 55.0,
        "source": (
            "DS105 tabla 'Optical and Electrical Features' (tau_fall = R_Q·C_J). "
            "tau_r = R_S·C_J no publicado."
        ),
        "note": (
            "tau_r no medido ni publicado por Broadcom. "
            "El script abortará a menos que se pase --tau-r-ns con un valor medido."
        ),
    },
    "fastic_measured": {
        "tau_r_ns": 2.0,
        "tau_f_ns": 3.0,
        "source": (
            "Capturas de osciloscopio del banco (1 GHz, 6.4 GS/s): "
            "ajuste exponencial da tau_fall = 55 ns (pulso libre) y par (2 ns, 3 ns) "
            "tras la red de cancelación polo-cero. "
            "La rama de tiempo del FastIC+ (TIA de ~20 Ohm + discriminador de corriente) "
            "no altera el flanco de subida (2 ns)."
        ),
        "note": (
            "Configuración real de este detector. "
            "El par (2, 3) ns es el pulso que ve la electrónica real de este banco. "
            "Requiere confirmación formal midiendo tau_r/tau_f directamente con el FastIC+ "
            "en la cadena completa (ver analysis/timing/TODO.md)."
        ),
    },
}

# ── SPTR — Lee et al., IEEE TRPMS 9(4) 406–411 (2025) ──────────────────────────
# DOI 10.1109/TRPMS.2024.3518479
# Dispositivo: AFBR-S4N66P014M (6×6 mm² NUV-MT), V_bias = 48 V, V_OV ≈ 15.5 V
#
# ADVERTENCIA DE PUNTO DE OPERACIÓN:
#   Este detector (W4) opera a V_OV = 10 V (W4_Weekly_meeting_Timinig_Detector_Status.pdf).
#   Lee midió a V_OV ≈ 15.5 V. El SPTR empeora al bajar el sobrevoltaje.
#   Los valores de Lee son cotas OPTIMISTAS para el punto de operación de este banco.
#
SPTR_LEE_INTRINSIC_FWHM_PS = 137.0   # ±4 ps; SiPM solo, sin electrónica
SPTR_LEE_DETECTOR_FWHM_PS  = 172.0   # FWHM incluyendo electrónica del banco de Lee

# Conversión explícita FWHM → sigma usando la constante definida arriba
SPTR_LEE_INTRINSIC_SIGMA_PS = SPTR_LEE_INTRINSIC_FWHM_PS * FWHM_TO_SIGMA   # ≈ 58.2 ps
SPTR_LEE_DETECTOR_SIGMA_PS  = SPTR_LEE_DETECTOR_FWHM_PS  * FWHM_TO_SIGMA   # ≈ 73.0 ps

SPTR_OV_WARNING = (
    f"ADVERTENCIA SPTR: Lee et al. (IEEE TRPMS 2025) midió a V_OV ≈ 15.5 V. "
    f"Este banco opera a V_OV = 10 V (W4 Weekly meeting). "
    f"El SPTR empeora al bajar el sobrevoltaje: los valores de Lee "
    f"(sigma_intrinsico ≈ {SPTR_LEE_INTRINSIC_SIGMA_PS:.1f} ps, "
    f"sigma_detector ≈ {SPTR_LEE_DETECTOR_SIGMA_PS:.1f} ps) "
    f"son cotas optimistas para el punto de operación real."
)

# ── Catálogo completo de SPTR ────────────────────────────────────────────────
# Cada entrada: dispositivo, sobrevoltaje, sigma_ps, fuente, estado.
# Estados: PUBLICADO · DERIVADO (banco propio) · PUBLICADO (otro OV) · PUBLICADO (otro dispositivo)
#
# Ajuste σ²(N) = SPTR²/N + σ_elec² sobre el barrido de intensidad láser (2026):
#   N=2,  σ=63 ps → 1/N=0.500, σ²=3969 ps²
#   N=4,  σ=43 ps → 1/N=0.250, σ²=1849 ps²   (midpoint rango 3–5)
#   N=7.5,σ=33.5ps→ 1/N=0.133, σ²=1122 ps²   (midpoints de rangos)
#   Pendiente (OLS) = 7877 ps²  →  SPTR = 88.8 ≈ 89 ps RMS
#   Intercepto = −6 ps² (→ σ_elec ≈ 0; cota superior ≤ 18 ps)
#
SPTR_CATALOG = {
    "lee_2025_intrinsic": {
        "device": "AFBR-S4N66P014M (6×6 mm² NUV-MT)",
        "ov_v": 15.5,
        "sigma_ps": SPTR_LEE_INTRINSIC_SIGMA_PS,   # 137 FWHM / 2.355 ≈ 58.2 ps
        "source": "Lee et al., IEEE TRPMS 9(4) 406–411 (2025), DOI 10.1109/TRPMS.2024.3518479",
        "status": "PUBLICADO (otro punto de operacion)",
        "note": "SiPM solo, sin electrónica del banco de Lee.",
    },
    "lee_2025_detector": {
        "device": "AFBR-S4N66P014M (6×6 mm² NUV-MT) + electrónica del banco de Lee",
        "ov_v": 15.5,
        "sigma_ps": SPTR_LEE_DETECTOR_SIGMA_PS,    # 172 FWHM / 2.355 ≈ 73.0 ps
        "source": "Lee et al., IEEE TRPMS 9(4) 406–411 (2025), DOI 10.1109/TRPMS.2024.3518479",
        "status": "PUBLICADO (otro punto de operacion)",
        "note": "Incluye electrónica del banco de Lee (no FastIC+).",
    },
    "banco_propio_low": {
        "device": "AFBR-S4N66P014M (6×6 mm² NUV-MT) + FastIC+",
        "ov_v": 10.0,
        "sigma_ps": 85.0,
        "source": (
            "Banco propio SHiP, barrido de intensidad láser, V_bias=38 V, umbral FastIC+=35. "
            "Ajuste OLS sigma^2=SPTR^2/N+sigma_elec^2; pendiente=7877 ps^2, SPTR=88.8 ps; "
            "rango declarado 85–92 ps (extremos del rango de N y sigma medidos)."
        ),
        "status": "DERIVADO (banco propio)",
        "note": (
            "Cota inferior del rango derivado. Requiere medida formal con datos de EOS "
            "(ver analysis/timing/TODO.md). "
            "Coherente con Lee 2025 (89 > 58 ps a OV menor, dirección esperada)."
        ),
    },
    "banco_propio_high": {
        "device": "AFBR-S4N66P014M (6×6 mm² NUV-MT) + FastIC+",
        "ov_v": 10.0,
        "sigma_ps": 92.0,
        "source": (
            "Banco propio SHiP, barrido de intensidad láser, V_bias=38 V, umbral FastIC+=35. "
            "Mismo ajuste que banco_propio_low; cota superior del rango."
        ),
        "status": "DERIVADO (banco propio)",
        "note": "Cota superior del rango derivado. Ver banco_propio_low.",
    },
    "colaboracion_geant4_ref": {
        "device": "Referencia Geant4 colaboración SHiP (origen: Hamamatsu S13360-3050CS, 3–4.5 V OV)",
        "ov_v": None,
        "sigma_ps": 106.0,
        "source": "Estándar colaboración SHiP; heredado de Hamamatsu S13360-3050CS.",
        "status": "PUBLICADO (otro dispositivo)",
        "note": (
            "ADVERTENCIA DE UNIDADES (ambigüedad sin resolver): "
            "las notas del grupo describen '100–120 ps FWHM' y '~106 ps RMS' como equivalentes. "
            "Son incompatibles: 110 ps FWHM = ~47 ps sigma. "
            "Si el valor inyectado en Geant4 (106 ps) se interpreta como sigma cuando "
            "la medida original era FWHM, el jitter está sobreestimado en factor ~2.355. "
            "Requiere consulta con Gerardo para resolver la convención. "
            "No usar este valor sin aclarar la ambigüedad."
        ),
    },
}
