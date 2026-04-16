#!/usr/bin/env python3
"""
sipm_waveform_dcfd.py

Reconstruye la waveform de un SiPM a partir de los tiempos de fotoelectrones
simulados en Geant4, usando:

1. La superposicion de respuestas SPE de Lv et al. (2026), Eq. (1):
       V(t) = sum_i v_i(t - t_transit_i - t0_i)
2. El pulso base de SiPM de Peña-Rodriguez (2024), Eq. (2.1):
       v(t) = A * (1 - exp(-t/tau_r)) * exp(-t/tau_f)
3. Un discriminador dCFD al 14% sobre el flanco ascendente.

El script opera sobre el TTree `sipm_hits` producido por esta simulacion y
estima la resolucion temporal ajustando una Gaussiana a la distribucion de
tiempos reconstruidos evento a evento.
"""

from __future__ import annotations

import argparse
import math
import pathlib
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit


def gauss(x: np.ndarray, mu: float, sigma: float, amplitude: float) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


@dataclass
class WaveformResult:
    event_uid: str
    source_file: str
    event_id: int
    face_type: int
    n_pe: int
    t_dcfd_ns: float
    amplitude_pe: float
    t_peak_ns: float
    gun_x_mm: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruye waveforms de SiPM desde sipm_hits ROOT y aplica "
            "dCFD al 14% para estimar resolucion temporal."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        default=["photon_hits_run*.root"],
        help="Archivos ROOT o globs a procesar.",
    )
    parser.add_argument("--tree", default="sipm_hits", help="Nombre del TTree.")
    parser.add_argument(
        "--face",
        default="0",
        help="Cara a analizar: 0, 1, 2 o 'all'. Default: 0.",
    )
    parser.add_argument(
        "--time-branch",
        default="time_raw_ns",
        help="Branch base para t0 de Geant4. Default: time_raw_ns.",
    )
    parser.add_argument(
        "--fraction",
        type=float,
        default=0.14,
        help="Fraccion del dCFD sobre la amplitud maxima. Default: 0.14.",
    )
    parser.add_argument(
        "--dt-ps",
        type=float,
        default=10.0,
        help="Paso temporal del waveform en ps. Default: 10 ps.",
    )
    parser.add_argument(
        "--transit-sigma-ps",
        type=float,
        default=200.0,
        help="Sigma del tiempo de transito de Lv et al. en ps. Default: 200 ps.",
    )
    parser.add_argument(
        "--tau-r-ns",
        type=float,
        default=2.0,
        help="Constante de subida Broadcom en ns. Default: 2 ns.",
    )
    parser.add_argument(
        "--tau-f-ns",
        type=float,
        default=55.0,
        help="Constante de bajada Broadcom en ns. Default: 55 ns.",
    )
    parser.add_argument(
        "--pulse-sigma-pe",
        type=float,
        default=0.1,
        help="Sigma de la amplitud SPE en pe. Default: 0.1 pe.",
    )
    parser.add_argument(
        "--electronics-sigma-ps",
        type=float,
        default=0.0,
        help="Jitter extra gaussiano de electronica en ps. Default: 0 ps.",
    )
    parser.add_argument(
        "--tail-window-factor",
        type=float,
        default=8.0,
        help="Longitud del kernel en multiplos de tau_f. Default: 8.",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Limita el numero de eventos procesados.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Semilla RNG. Default: 12345.",
    )
    parser.add_argument(
        "--plot-waveform",
        action="store_true",
        help="Guarda una waveform ejemplo con el cruce dCFD.",
    )
    parser.add_argument(
        "--waveform-out",
        default="sipm_waveform_example.png",
        help="Salida de la waveform ejemplo.",
    )
    parser.add_argument(
        "--hist-out",
        default="sipm_dcfd_resolution.png",
        help="Salida del histograma de tiempos reconstruidos.",
    )
    parser.add_argument(
        "--csv-out",
        default="sipm_dcfd_times.csv",
        help="CSV con tiempos reconstruidos por evento.",
    )
    return parser.parse_args()


def resolve_inputs(patterns: Iterable[str]) -> list[pathlib.Path]:
    files: list[pathlib.Path] = []
    for pattern in patterns:
        candidate = pathlib.Path(pattern)
        if candidate.is_absolute():
            if any(ch in pattern for ch in "*?[]"):
                matches = sorted(candidate.parent.glob(candidate.name))
                files.extend(matches)
            elif candidate.exists():
                files.append(candidate)
            continue

        matches = sorted(pathlib.Path(".").glob(pattern))
        if matches:
            files.extend(matches)
            continue
        if candidate.exists():
            files.append(candidate)
    unique = []
    seen = set()
    for path in files:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(resolved)
    return unique


def read_schema(root_file: pathlib.Path, tree_name: str, time_branch: str) -> tuple[str, list[str]]:
    with uproot.open(root_file) as handle:
        tree = handle[tree_name]
        available = set(tree.keys())
    base_branch = time_branch if time_branch in available else "time_ns"
    columns = ["event_id", "face_type", base_branch]
    if "gun_x_mm" in available:
        columns.append("gun_x_mm")
    return base_branch, columns


def prepare_chunk(df: pd.DataFrame, base_branch: str, face: str) -> pd.DataFrame:
    if "gun_x_mm" not in df.columns:
        df["gun_x_mm"] = 0.0
    df = df.rename(columns={base_branch: "t0_ns"})
    return select_faces(df, face)


def build_spe_kernel(
    dt_ns: float,
    tau_r_ns: float,
    tau_f_ns: float,
    tail_window_factor: float,
) -> tuple[np.ndarray, np.ndarray]:
    kernel_span_ns = max(tail_window_factor * tau_f_ns, 5.0 * dt_ns)
    time_axis_ns = np.arange(0.0, kernel_span_ns + dt_ns, dt_ns)
    kernel = (1.0 - np.exp(-time_axis_ns / tau_r_ns)) * np.exp(-time_axis_ns / tau_f_ns)
    return time_axis_ns, kernel


def reconstruct_waveform(
    hit_times_ns: np.ndarray,
    dt_ns: float,
    kernel_t_ns: np.ndarray,
    kernel: np.ndarray,
    transit_sigma_ns: float,
    pulse_sigma_pe: float,
    electronics_sigma_ns: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    if hit_times_ns.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    transit_jitter = rng.normal(0.0, transit_sigma_ns, size=hit_times_ns.size)
    amplitudes = rng.normal(1.0, pulse_sigma_pe, size=hit_times_ns.size)
    amplitudes = np.clip(amplitudes, 0.0, None)
    shifted_times = np.asarray(hit_times_ns, dtype=float) + transit_jitter

    pre_ns = max(5.0 * transit_sigma_ns, dt_ns)
    start_ns = shifted_times.min() - pre_ns
    stop_ns = shifted_times.max() + kernel_t_ns[-1] + dt_ns
    n_bins = int(np.ceil((stop_ns - start_ns) / dt_ns)) + 1
    waveform_t_ns = start_ns + np.arange(n_bins, dtype=float) * dt_ns

    hit_bins = np.rint((shifted_times - start_ns) / dt_ns).astype(int)
    hit_bins = np.clip(hit_bins, 0, n_bins - 1)

    delta_train = np.zeros(n_bins, dtype=float)
    np.add.at(delta_train, hit_bins, amplitudes)

    waveform = np.convolve(delta_train, kernel, mode="full")[:n_bins]

    if electronics_sigma_ns > 0.0:
        waveform += rng.normal(0.0, electronics_sigma_ns, size=waveform.size)

    return waveform_t_ns, waveform


def dcfd_time_from_waveform(
    time_ns: np.ndarray,
    waveform: np.ndarray,
    fraction: float,
) -> float:
    if waveform.size == 0:
        return float("nan")

    peak_idx = int(np.argmax(waveform))
    peak_amp = float(waveform[peak_idx])
    if not np.isfinite(peak_amp) or peak_amp <= 0.0:
        return float("nan")

    threshold = fraction * peak_amp
    rising = waveform[: peak_idx + 1]
    crossing = np.where(rising >= threshold)[0]
    if crossing.size == 0:
        return float("nan")

    idx = int(crossing[0])
    if idx == 0:
        return float(time_ns[0])

    x0, x1 = time_ns[idx - 1], time_ns[idx]
    y0, y1 = waveform[idx - 1], waveform[idx]
    if np.isclose(y1, y0):
        return float(x1)
    return float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))


def fit_resolution(times_ns: np.ndarray) -> tuple[float, float, float]:
    if times_ns.size < 5:
        return float("nan"), float("nan"), float("nan")

    lo, hi = np.percentile(times_ns, [1, 99])
    width = max(hi - lo, 0.2)
    lo -= 0.25 * width
    hi += 0.25 * width
    bins = min(80, max(25, int(np.sqrt(times_ns.size) * 2)))

    counts, edges = np.histogram(times_ns, bins=bins, range=(lo, hi))
    centers = 0.5 * (edges[:-1] + edges[1:])

    mu0 = float(np.mean(times_ns))
    sigma0 = float(np.std(times_ns))
    sigma0 = max(sigma0, 1e-3)
    amp0 = float(max(counts.max(), 1.0))

    try:
        popt, pcov = curve_fit(
            gauss,
            centers,
            counts,
            p0=[mu0, sigma0, amp0],
            bounds=([lo, 0.0, 0.0], [hi, hi - lo, 10.0 * amp0]),
            maxfev=10000,
        )
        perr = np.sqrt(np.diag(pcov))
        return float(popt[0]), abs(float(popt[1])), float(perr[1])
    except Exception:
        sigma = float(np.std(times_ns))
        sigma_err = sigma / math.sqrt(max(times_ns.size - 1, 1))
        return float(np.mean(times_ns)), sigma, sigma_err


def select_faces(df: pd.DataFrame, face: str) -> pd.DataFrame:
    if face.lower() == "all":
        return df.copy()
    face_int = int(face)
    return df[df["face_type"] == face_int].copy()


def process_file(
    root_file: pathlib.Path,
    args: argparse.Namespace,
    rng: np.random.Generator,
    kernel_t_ns: np.ndarray,
    kernel: np.ndarray,
    max_events_left: int | None,
) -> list[WaveformResult]:
    base_branch, columns = read_schema(root_file, args.tree, args.time_branch)
    dt_ns = args.dt_ps * 1e-3
    transit_sigma_ns = args.transit_sigma_ps * 1e-3
    electronics_sigma_ns = args.electronics_sigma_ps * 1e-3

    results: list[WaveformResult] = []
    carry = pd.DataFrame()

    def consume_groups(frame: pd.DataFrame) -> bool:
        nonlocal results
        if frame.empty:
            return False
        grouped = frame.groupby(["event_id", "face_type"], sort=True)
        for (event_id, face_type), grp in grouped:
            waveform_t_ns, waveform = reconstruct_waveform(
                grp["t0_ns"].to_numpy(),
                dt_ns=dt_ns,
                kernel_t_ns=kernel_t_ns,
                kernel=kernel,
                transit_sigma_ns=transit_sigma_ns,
                pulse_sigma_pe=args.pulse_sigma_pe,
                electronics_sigma_ns=electronics_sigma_ns,
                rng=rng,
            )
            t_dcfd_ns = dcfd_time_from_waveform(waveform_t_ns, waveform, args.fraction)
            if not np.isfinite(t_dcfd_ns):
                continue

            peak_idx = int(np.argmax(waveform))
            results.append(
                WaveformResult(
                    event_uid=f"{root_file.name}:{int(event_id)}:{int(face_type)}",
                    source_file=root_file.name,
                    event_id=int(event_id),
                    face_type=int(face_type),
                    n_pe=int(len(grp)),
                    t_dcfd_ns=float(t_dcfd_ns),
                    amplitude_pe=float(waveform[peak_idx]),
                    t_peak_ns=float(waveform_t_ns[peak_idx]),
                    gun_x_mm=float(grp["gun_x_mm"].iloc[0]),
                )
            )
            if max_events_left is not None and len(results) >= max_events_left:
                return True
        return False

    for chunk in uproot.iterate(
        f"{root_file}:{args.tree}",
        expressions=columns,
        library="pd",
        step_size="25 MB",
    ):
        chunk = prepare_chunk(chunk, base_branch, args.face)
        if chunk.empty:
            continue
        if not carry.empty:
            chunk = pd.concat([carry, chunk], ignore_index=True)

        last_event = chunk["event_id"].iloc[-1]
        last_face = chunk["face_type"].iloc[-1]
        is_last_group = (chunk["event_id"] == last_event) & (chunk["face_type"] == last_face)
        carry = chunk[is_last_group].copy()
        completed = chunk[~is_last_group].copy()

        if consume_groups(completed):
            return results

    if not carry.empty:
        consume_groups(carry)

    return results


def plot_example_waveform(
    root_file: pathlib.Path,
    args: argparse.Namespace,
    rng: np.random.Generator,
    kernel_t_ns: np.ndarray,
    kernel: np.ndarray,
) -> None:
    import matplotlib.pyplot as plt

    base_branch, columns = read_schema(root_file, args.tree, args.time_branch)
    df = None
    for chunk in uproot.iterate(
        f"{root_file}:{args.tree}",
        expressions=columns,
        library="pd",
        step_size="5 MB",
    ):
        chunk = prepare_chunk(chunk, base_branch, args.face)
        if not chunk.empty:
            df = chunk
            break
    if df is None:
        return
    if df.empty:
        return

    event_key, grp = next(iter(df.groupby(["event_id", "face_type"], sort=True)))
    waveform_t_ns, waveform = reconstruct_waveform(
        grp["t0_ns"].to_numpy(),
        dt_ns=args.dt_ps * 1e-3,
        kernel_t_ns=kernel_t_ns,
        kernel=kernel,
        transit_sigma_ns=args.transit_sigma_ps * 1e-3,
        pulse_sigma_pe=args.pulse_sigma_pe,
        electronics_sigma_ns=args.electronics_sigma_ps * 1e-3,
        rng=rng,
    )
    t_dcfd_ns = dcfd_time_from_waveform(waveform_t_ns, waveform, args.fraction)
    peak_amp = float(np.max(waveform))
    threshold = args.fraction * peak_amp

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(waveform_t_ns, waveform, color="#1f4e79", lw=1.8)
    ax.axhline(threshold, color="#c0392b", ls="--", lw=1.2, label=f"{args.fraction*100:.0f}% CFD")
    if np.isfinite(t_dcfd_ns):
        ax.axvline(t_dcfd_ns, color="#c0392b", ls=":", lw=1.2, label=f"t_dCFD = {t_dcfd_ns:.3f} ns")
    ax.set_title(
        f"Waveform SiPM sintetizada ({root_file.name}, evento {event_key[0]}, cara {event_key[1]})"
    )
    ax.set_xlabel("Tiempo [ns]")
    ax.set_ylabel("Amplitud [pe arb.]")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.waveform_out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_time_histogram(df: pd.DataFrame, args: argparse.Namespace, mu_ns: float, sigma_ns: float) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.5))
    times_ns = df["t_dcfd_ns"].to_numpy()
    bins = min(80, max(25, int(np.sqrt(times_ns.size) * 2)))
    counts, edges, _ = ax.hist(
        times_ns,
        bins=bins,
        color="#5dade2",
        alpha=0.8,
        edgecolor="#1b4f72",
    )
    if np.isfinite(mu_ns) and np.isfinite(sigma_ns) and sigma_ns > 0.0:
        x = np.linspace(edges[0], edges[-1], 500)
        area = np.sum(counts) * (edges[1] - edges[0])
        ax.plot(x, gauss(x, mu_ns, sigma_ns, area / (sigma_ns * np.sqrt(2.0 * np.pi))), color="#922b21", lw=2.0)
    ax.set_xlabel("Tiempo dCFD [ns]")
    ax.set_ylabel("Eventos / bin")
    ax.set_title("Distribucion temporal reconstruida con waveform SiPM + dCFD 14%")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.hist_out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    files = resolve_inputs(args.inputs)
    if not files:
        raise SystemExit("No se encontraron archivos ROOT de entrada.")

    rng = np.random.default_rng(args.seed)
    kernel_t_ns, kernel = build_spe_kernel(
        dt_ns=args.dt_ps * 1e-3,
        tau_r_ns=args.tau_r_ns,
        tau_f_ns=args.tau_f_ns,
        tail_window_factor=args.tail_window_factor,
    )

    results: list[WaveformResult] = []
    max_events_left = args.max_events
    for root_file in files:
        chunk_results = process_file(
            root_file,
            args=args,
            rng=rng,
            kernel_t_ns=kernel_t_ns,
            kernel=kernel,
            max_events_left=max_events_left,
        )
        results.extend(chunk_results)
        if max_events_left is not None:
            max_events_left -= len(chunk_results)
            if max_events_left <= 0:
                break

    if not results:
        raise SystemExit("No se pudieron reconstruir eventos para la seleccion pedida.")

    df = pd.DataFrame([vars(item) for item in results])
    df.to_csv(args.csv_out, index=False)

    mu_ns, sigma_ns, sigma_err_ns = fit_resolution(df["t_dcfd_ns"].to_numpy())
    plot_time_histogram(df, args, mu_ns, sigma_ns)

    if args.plot_waveform:
        plot_example_waveform(files[0], args, rng, kernel_t_ns, kernel)

    face_label = args.face
    print(f"Archivos procesados: {len(files)}")
    print(f"Seleccion de cara: {face_label}")
    print(f"Eventos reconstruidos: {len(df)}")
    print(f"Pulso Broadcom: tau_r = {args.tau_r_ns:.3f} ns, tau_f = {args.tau_f_ns:.3f} ns")
    print(f"dCFD: {args.fraction * 100.0:.1f}%")
    print(f"Paso temporal: {args.dt_ps:.1f} ps")
    print(f"Sigma transit-time: {args.transit_sigma_ps:.1f} ps")
    print(f"Jitter electronica: {args.electronics_sigma_ps:.1f} ps")
    print(f"Tiempo medio reconstruido: {mu_ns:.6f} ns")
    print(f"Resolucion temporal sigma_t: {sigma_ns * 1e3:.2f} ps")
    print(f"Error estadistico del sigma: {sigma_err_ns * 1e3:.2f} ps")
    print(f"CSV: {args.csv_out}")
    print(f"Histograma: {args.hist_out}")
    if args.plot_waveform:
        print(f"Waveform ejemplo: {args.waveform_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
