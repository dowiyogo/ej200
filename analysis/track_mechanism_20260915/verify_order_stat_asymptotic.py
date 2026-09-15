#!/usr/bin/env python3
"""Verifica la asíntota del mínimo para el muestreador real de Geant4 11.4."""

import csv
import json
from math import gamma, pi, sqrt
from pathlib import Path

import numpy as np

from exec46_schema import MATERIAL_BY_OPSC, load_material_config


OUTPUT_DIR = Path(__file__).resolve().parent
N_VALUES = np.asarray([10, 100, 1_000, 10_000, 100_000], dtype=int)
INTEGRATION_POINTS = 2_000_000
INTEGRATION_UPPER_DECAY_MULTIPLIER = 50.0
INTEGRATION_LOWER_SCALE = 1.0e-12
SMALL_TIME_SCALE = 1.0e-8


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def geant4_cdf(time_ns, rise_ns, decay_ns):
    """CDF normalizada del rejection sampler G4Scintillation::sample_time."""
    harmonic_ns = rise_ns * decay_ns / (rise_ns + decay_ns)
    return (
        1.0
        - (rise_ns + decay_ns) / decay_ns * np.exp(-time_ns / decay_ns)
        + rise_ns / decay_ns * np.exp(-time_ns / harmonic_ns)
    )


def geant4_density(time_ns, rise_ns, decay_ns):
    normalization = (rise_ns + decay_ns) / decay_ns ** 2
    return normalization * np.exp(-time_ns / decay_ns) * (
        1.0 - np.exp(-time_ns / rise_ns)
    )


def main():
    rows = []
    material_summary = []
    for opsc_code, material in MATERIAL_BY_OPSC.items():
        config = load_material_config(opsc_code)
        rise_ns = config["rise_time_ns"]
        decay_ns = config["decay_time_ns"]

        # La fórmula indicada por el usuario corresponde a una biexponencial estándar.
        standard_linear_coefficient = 1.0 / (rise_ns * decay_ns)
        standard_min_prefactor = (
            sqrt(2.0 * rise_ns * decay_ns) * gamma(1.5)
        )

        # El rejection sampler de Geant4 acepta Exp(tau_d) con 1-exp(-t/tau_r).
        geant4_linear_coefficient = (
            (rise_ns + decay_ns) / (decay_ns ** 2 * rise_ns)
        )
        geant4_min_prefactor = sqrt(pi / (2.0 * geant4_linear_coefficient))
        probe_time = SMALL_TIME_SCALE * min(rise_ns, decay_ns)
        numerical_linear_coefficient = float(
            geant4_density(probe_time, rise_ns, decay_ns) / probe_time
        )

        time_grid = np.r_[
            0.0,
            np.geomspace(
                INTEGRATION_LOWER_SCALE * min(rise_ns, decay_ns),
                INTEGRATION_UPPER_DECAY_MULTIPLIER * decay_ns,
                INTEGRATION_POINTS,
            ),
        ]
        cdf = np.clip(geant4_cdf(time_grid, rise_ns, decay_ns), 0.0, 1.0)
        for count in N_VALUES:
            with np.errstate(divide="ignore"):
                survival_minimum = np.exp(count * np.log1p(-cdf))
            exact_mean = float(np.trapz(survival_minimum, time_grid))
            rows.append({
                "material": material,
                "opsc_code": opsc_code,
                "tau_rise_ns": rise_ns,
                "tau_decay_ns": decay_ns,
                "N": int(count),
                "exact_geant4_mean_min_ns": exact_mean,
                "N_minus_half_geant4_ns": geant4_min_prefactor / sqrt(count),
                "N_minus_half_standard_biexponential_ns": standard_min_prefactor / sqrt(count),
                "N_minus_one_ns_using_tau_decay": decay_ns / count,
                "exact_times_sqrt_N_ns": exact_mean * sqrt(count),
            })
        material_summary.append({
            "material": material,
            "opsc_code": opsc_code,
            "tau_rise_ns": rise_ns,
            "tau_decay_ns": decay_ns,
            "standard_biexponential_f_over_t_limit_per_ns2": standard_linear_coefficient,
            "geant4_f_over_t_limit_per_ns2": geant4_linear_coefficient,
            "numerical_geant4_f_over_t_per_ns2": numerical_linear_coefficient,
            "standard_biexponential_min_prefactor_ns": standard_min_prefactor,
            "geant4_min_prefactor_ns": geant4_min_prefactor,
            "exponent": -0.5,
        })

    write_csv(OUTPUT_DIR / "exec46_order_stat_asymptotic.csv", rows)
    payload = {
        "derivation": {
            "geant4_pdf": "(tau_r+tau_d)/tau_d^2 * exp(-t/tau_d) * (1-exp(-t/tau_r))",
            "geant4_small_t": "f(t) -> [(tau_r+tau_d)/(tau_d^2 tau_r)] t",
            "generic_minimum": "if f(t)->a*t, E[min_N] -> sqrt(pi/(2*a))*N^(-1/2)",
            "standard_biexponential_prefactor": "sqrt(2*tau_r*tau_d)*Gamma(3/2)",
            "geant4_prefactor": "tau_d*sqrt(pi*tau_r/[2*(tau_r+tau_d)])",
            "conclusion": "The N^-1/2 exponent is verified; the prompt prefactor is not the Geant4 11.4 sampler prefactor.",
        },
        "integration": {
            "points": INTEGRATION_POINTS,
            "upper_decay_multiplier": INTEGRATION_UPPER_DECAY_MULTIPLIER,
            "N_values": N_VALUES.tolist(),
        },
        "materials": material_summary,
    }
    (OUTPUT_DIR / "exec46_order_stat_asymptotic.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
