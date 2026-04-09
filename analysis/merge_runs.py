#!/usr/bin/env python3
"""
merge_runs.py — fusiona archivos photon_hits_run*.root generados por scan.mac

El problema con hadd es que G4AnalysisManager escribe metadatos internos de ROOT
con ownership compartido entre runs, lo que corrompe el archivo al hacer hadd.
Este script evita ese problema leyendo con uproot y reescribiendo con awkward/uproot.

Uso:
    python merge_runs.py                          # busca photon_hits_run*.root en .
    python merge_runs.py --pattern "run_*.root"   # patron personalizado
    python merge_runs.py --dir /ruta/a/runs       # directorio de entrada
    python merge_runs.py --out merged.root        # nombre de salida
    python merge_runs.py --check                  # imprime resumen y sale sin escribir
"""

import argparse
import glob
import pathlib
import sys

import numpy as np
import pandas as pd
import uproot


TREE_NAME   = "sipm_hits"
BASE_COLS   = [
    "event_id", "face_type", "global_id", "local_id",
    "time_ns", "energy_eV", "wl_nm", "pde",
    "x_mm", "y_mm", "z_mm",
]
OPTIONAL_COLS = ["gun_x_mm"]

FACE_NAMES = {0: "End left  (−X)", 1: "End right (+X)", 2: "Top  (+Y)"}


# ── Lectura ───────────────────────────────────────────────────────────────────

def read_one_file(path: str, event_offset: int) -> tuple[pd.DataFrame, int]:
    """
    Lee el TTree sipm_hits de un archivo ROOT.
    Suma event_offset al event_id para que los IDs sean únicos entre runs.
    Retorna (dataframe, nuevo_offset).
    """
    with uproot.open(path) as rf:
        if TREE_NAME not in rf:
            raise KeyError(
                f"'{path}' no contiene el TTree '{TREE_NAME}'.\n"
                f"  Llaves disponibles: {list(rf.keys())}"
            )
        tree = rf[TREE_NAME]
        available = set(tree.keys())
        cols = BASE_COLS + [c for c in OPTIONAL_COLS if c in available]
        df = tree.arrays(cols, library="pd")

    if "gun_x_mm" not in df.columns:
        df["gun_x_mm"] = 0.0

    # Desplazar event_id para que no colisione con runs anteriores.
    # Importante: calculamos el max ANTES de sumar el offset.
    n_events_this_run = int(df["event_id"].max()) + 1
    df["event_id"] = df["event_id"] + event_offset

    new_offset = event_offset + n_events_this_run
    return df, new_offset


def load_all(files: list[str]) -> pd.DataFrame:
    dfs = []
    event_offset = 0

    for i, fpath in enumerate(files):
        try:
            df_run, event_offset = read_one_file(fpath, event_offset)
            n_ph   = len(df_run)
            n_evts = df_run["event_id"].nunique()
            x_vals = sorted(df_run["gun_x_mm"].unique())
            print(f"  [{i+1:>3}/{len(files)}] {pathlib.Path(fpath).name:<35} "
                  f"{n_ph:>8,} fotones  {n_evts:>5} eventos  "
                  f"x ∈ {{{', '.join(f'{x:.0f}' for x in x_vals[:3])}"
                  f"{'...' if len(x_vals) > 3 else ''}}} mm")
            dfs.append(df_run)
        except Exception as exc:
            print(f"  [!] ERROR en '{fpath}': {exc}", file=sys.stderr)
            print(f"      Archivo saltado.", file=sys.stderr)

    if not dfs:
        raise RuntimeError("Ningún archivo se pudo leer. Abortando.")

    df_all = pd.concat(dfs, ignore_index=True)
    return df_all


# ── Escritura ─────────────────────────────────────────────────────────────────

def write_root(df: pd.DataFrame, out_path: str):
    """Escribe el DataFrame fusionado como TTree en un archivo ROOT nuevo."""

    # Separar columnas enteras y flotantes para uproot
    int_cols   = ["event_id", "face_type", "global_id", "local_id"]
    float_cols = ["time_ns", "energy_eV", "wl_nm", "pde",
                  "x_mm", "y_mm", "z_mm", "gun_x_mm"]

    data = {}
    for c in int_cols:
        data[c] = df[c].to_numpy(dtype=np.int32)
    for c in float_cols:
        data[c] = df[c].to_numpy(dtype=np.float64)

    with uproot.recreate(out_path) as fout:
        fout[TREE_NAME] = data

    print(f"\n  ROOT escrito → {out_path}")
    print(f"  TTree '{TREE_NAME}' con {len(df):,} entradas")


# ── Resumen ───────────────────────────────────────────────────────────────────

def print_summary(df: pd.DataFrame):
    total  = len(df)
    n_evts = df["event_id"].nunique()
    n_xpos = df["gun_x_mm"].nunique()

    print(f"\n{'='*58}")
    print(f"  Total fotones detectados : {total:>10,}")
    print(f"  Eventos únicos           : {n_evts:>10,}")
    print(f"  Fotones / evento (media) : {total/n_evts:>10.1f}")
    print(f"  Posiciones x del gun     : {n_xpos:>10,}")
    print("-"*58)
    for face, grp in df.groupby("face_type"):
        name = FACE_NAMES.get(int(face), f"face {face}")
        print(f"  {name:<22} : {len(grp):>8,}  ({100*len(grp)/total:.1f}%)")
    print("="*58)

    if n_xpos > 1:
        x_sorted = sorted(df["gun_x_mm"].unique())
        print(f"\n  Posiciones gun_x_mm [mm]:")
        print(f"    {x_sorted}")
    print()


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Fusiona archivos photon_hits_run*.root (evita bugs de hadd con G4AnalysisManager)"
    )
    parser.add_argument(
        "--pattern", default="photon_hits_run*.root",
        help="Glob pattern de archivos de entrada (default: photon_hits_run*.root)"
    )
    parser.add_argument(
        "--dir", default=".",
        help="Directorio donde buscar los archivos (default: directorio actual)"
    )
    parser.add_argument(
        "--out", default="photon_hits_merged.root",
        help="Archivo ROOT de salida (default: photon_hits_merged.root)"
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Solo imprime resumen, no escribe el archivo de salida"
    )
    args = parser.parse_args()

    # Buscar archivos
    search_pattern = str(pathlib.Path(args.dir) / args.pattern)
    files = sorted(glob.glob(search_pattern))

    if not files:
        print(f"ERROR: no se encontraron archivos con patrón '{search_pattern}'")
        print("  Verifica --dir y --pattern")
        sys.exit(1)

    print(f"Encontrados {len(files)} archivos con patrón '{search_pattern}':")
    print()

    # Leer y fusionar
    df = load_all(files)

    # Resumen
    print_summary(df)

    # Escribir (salvo que solo se pida el check)
    if args.check:
        print("  [--check] No se escribe archivo de salida.")
    else:
        out_path = str(pathlib.Path(args.dir) / args.out) \
                   if not pathlib.Path(args.out).is_absolute() else args.out
        write_root(df, out_path)
        print(f"\n  Listo. Puedes leerlo con:")
        print(f"    python analyze.py {out_path}")
        print(f"    python resolution_vs_x.py {out_path}")


if __name__ == "__main__":
    main()
