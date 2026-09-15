#!/usr/bin/env python3
"""Esquema y configuración compartidos por el análisis fotón-a-fotón EXEC_46."""

from collections import OrderedDict
from pathlib import Path
import re

import numpy as np


TREE_NAME = "sipm_hits"
CAMPAIGN_DIR = Path("/home/rrios/exec46_20260915/full_grid")
SSLG4_DIR = Path("/home/rrios/exec46_20260915/build_baseline/sslg4")
MATERIAL_BY_OPSC = {
    "OPSC-100": "EJ-200",
    "OPSC-101": "EJ-204",
    "OPSC-106": "EJ-230",
}
SOURCE_LABELS = {0: "primary", 1: "scintillation", 2: "Cherenkov", 3: "other"}

# ROOT expone Int_t como int32_t y Double_t como double a través de uproot.
BRANCH_TYPES = OrderedDict([
    ("event_id", "int32_t"),
    ("face_type", "int32_t"),
    ("global_id", "int32_t"),
    ("local_id", "int32_t"),
    ("time_ns", "double"),
    ("energy_eV", "double"),
    ("wl_nm", "double"),
    ("pde", "double"),
    ("x_mm", "double"),
    ("y_mm", "double"),
    ("z_mm", "double"),
    ("gun_x_mm", "double"),
    ("track_id", "int32_t"),
    ("t_detection_ns", "double"),
    ("t_creation_ns", "double"),
    ("x_creation_mm", "double"),
    ("y_creation_mm", "double"),
    ("z_creation_mm", "double"),
    ("wl_nm_created", "double"),
    ("path_length_mm", "double"),
    ("exit_angle_deg", "double"),
    ("n_boundary_encounters", "int32_t"),
    ("source_type", "int32_t"),
])

# Geometría copiada de DetectorConstruction.cc en el commit de simulación 4967ec8.
BAR_HALF_X_MM = 700.0
BAR_HALF_Y_MM = 30.0
BAR_HALF_Z_MM = 5.0
END_HALF_X_MM = 0.25
END_HALF_Y_MM = 3.0
END_HALF_Z_MM = 3.0
END_PITCH_MM = 7.5
TOP_HALF_X_MM = 3.0
TOP_HALF_Y_MM = 0.25
TOP_HALF_Z_MM = 3.0
N_END_PER_FACE = 8
N_TOP = 70
N_SENSORS = 2 * N_END_PER_FACE + N_TOP
LEFT_FACE = 0
RIGHT_FACE = 1
TOP_FACE = 2

SPEED_OF_LIGHT_MM_PER_NS = 299.792458
HC_EV_NM = 1239.8419843320026


def sensor_geometry(global_id):
    """Devuelve (face_type, local_id, x_mm, y_mm, z_mm)."""
    if not 0 <= global_id < N_SENSORS:
        raise ValueError(f"global_id fuera de rango: {global_id}")
    if global_id < N_END_PER_FACE:
        face_type = LEFT_FACE
        local_id = global_id
        x_mm = -(BAR_HALF_X_MM - END_HALF_X_MM)
        y_mm = (local_id - 3.5) * END_PITCH_MM
        return face_type, local_id, x_mm, y_mm, 0.0
    if global_id < 2 * N_END_PER_FACE:
        face_type = RIGHT_FACE
        local_id = global_id - N_END_PER_FACE
        x_mm = BAR_HALF_X_MM - END_HALF_X_MM
        y_mm = (local_id - 3.5) * END_PITCH_MM
        return face_type, local_id, x_mm, y_mm, 0.0
    local_id = global_id - 2 * N_END_PER_FACE
    x_mm = (-692.0 + 20.0 * local_id) if local_id < 35 else (12.0 + 20.0 * (local_id - 35))
    return TOP_FACE, local_id, x_mm, BAR_HALF_Y_MM - TOP_HALF_Y_MM, 0.0


SENSOR_MAP = {global_id: sensor_geometry(global_id) for global_id in range(N_SENSORS)}


def material_paths(opsc_code):
    """Resuelve los archivos de propiedades desde la instalación usada por la campaña."""
    stem = opsc_code.lower()
    return {
        "macro": SSLG4_DIR / "macros" / "oscnt" / f"{stem}.mac",
        "rindex": SSLG4_DIR / "data" / "oscnt" / stem / "rIndex.txt",
        "emission": SSLG4_DIR / "data" / "oscnt" / stem / "scntComp1.txt",
        "absorption": SSLG4_DIR / "data" / "oscnt" / stem / "absLength.txt",
    }


def _const_property_ns(macro_path, property_name):
    text = Path(macro_path).read_text()
    pattern = re.compile(
        r"^\s*/mpt/\{scnt\}/addConstProperty\s+" + re.escape(property_name)
        + r"\s+([0-9.eE+-]+)\s+(ps|ns)\b",
        re.MULTILINE,
    )
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise RuntimeError(f"propiedad {property_name} ausente o duplicada en {macro_path}")
    value, unit = matches[0]
    return float(value) * (1.0e-3 if unit == "ps" else 1.0)


def load_material_config(opsc_code):
    """Lee tiempos, RINDEX y ABSLENGTH sin heredar valores entre materiales."""
    paths = material_paths(opsc_code)
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    wavelength_nm, rindex = np.loadtxt(paths["rindex"], unpack=True)
    absorption_wavelength_nm, absorption_cm = np.loadtxt(
        paths["absorption"], unpack=True)
    order = np.argsort(HC_EV_NM / wavelength_nm)
    energy_ev = (HC_EV_NM / wavelength_nm)[order]
    rindex = rindex[order]
    return {
        "opsc_code": opsc_code,
        "material": MATERIAL_BY_OPSC[opsc_code],
        "decay_time_ns": _const_property_ns(paths["macro"], "SCINTILLATIONTIMECONSTANT1"),
        "rise_time_ns": _const_property_ns(paths["macro"], "SCINTILLATIONRISETIME1"),
        "energy_ev": energy_ev,
        "rindex": rindex,
        "rindex_wavelength_nm": wavelength_nm,
        "absorption_wavelength_nm": absorption_wavelength_nm,
        "absorption_length_mm": absorption_cm * 10.0,
        "paths": paths,
    }


def geant4_group_velocity_table(energy_ev, rindex):
    """Reproduce CalculateGROUPVEL de G4MaterialPropertiesTable 11.4."""
    if len(energy_ev) != len(rindex) or len(energy_ev) == 0:
        raise ValueError("tabla RINDEX vacía o inconsistente")
    if len(energy_ev) == 1:
        return np.asarray(energy_ev), np.asarray([SPEED_OF_LIGHT_MM_PER_NS / rindex[0]])
    e0, e1 = energy_ev[0], energy_ev[1]
    n0, n1 = rindex[0], rindex[1]
    speed = SPEED_OF_LIGHT_MM_PER_NS / (n0 + (n1 - n0) / np.log(e1 / e0))
    if speed < 0.0 or speed > SPEED_OF_LIGHT_MM_PER_NS / n0:
        speed = SPEED_OF_LIGHT_MM_PER_NS / n0
    group_energy = [e0]
    group_speed = [speed]
    for index in range(2, len(energy_ev)):
        mean_n = 0.5 * (n0 + n1)
        speed = SPEED_OF_LIGHT_MM_PER_NS / (mean_n + (n1 - n0) / np.log(e1 / e0))
        if speed < 0.0 or speed > SPEED_OF_LIGHT_MM_PER_NS / mean_n:
            speed = SPEED_OF_LIGHT_MM_PER_NS / mean_n
        group_energy.append(0.5 * (e0 + e1))
        group_speed.append(speed)
        e0, n0 = e1, n1
        e1, n1 = energy_ev[index], rindex[index]
    speed = SPEED_OF_LIGHT_MM_PER_NS / (n1 + (n1 - n0) / np.log(e1 / e0))
    if speed < 0.0 or speed > SPEED_OF_LIGHT_MM_PER_NS / n1:
        speed = SPEED_OF_LIGHT_MM_PER_NS / n1
    group_energy.append(e1)
    group_speed.append(speed)
    return np.asarray(group_energy), np.asarray(group_speed)


def group_velocity_mm_per_ns(wavelength_nm, material_config):
    energy_ev = HC_EV_NM / wavelength_nm
    table_energy, table_speed = geant4_group_velocity_table(
        material_config["energy_ev"], material_config["rindex"])
    return np.interp(energy_ev, table_energy, table_speed)
