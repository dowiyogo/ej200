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
        "tau_r_ns": None,
        "tau_f_ns": None,
        "source": "PENDIENTE — medida en banco de láser con FastIC+",
        "note": (
            "Configuración real de este detector. Aún no medida. "
            "Es la que debe usarse al finalizar la caracterización. "
            "Requiere --tau-r-ns y --tau-f-ns explícitos."
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
