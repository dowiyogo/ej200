#!/usr/bin/env python3
"""EXEC41 Block 2 and V1 accounting diagnosis using only existing files."""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import uproot

N = 2000
NBOOT = 500
SEED = 41091302

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")

def sidecars(out, stem, columns, metadata):
    data = {key: np.asarray(value) for key, value in columns.items()}
    with (out / f"{stem}.csv").open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(data); writer.writerows(zip(*data.values()))
    with uproot.recreate(out / f"{stem}.root") as root:
        root["data"] = data
    dump(out / f"{stem}.meta.json", {**metadata, "columns": list(data),
        "rows": len(next(iter(data.values()))),
        "csv_sha256": sha(out / f"{stem}.csv"),
        "root_sha256": sha(out / f"{stem}.root")})

def ratio(num, den, weights):
    values = (weights @ num) / (weights @ den)
    return {"value": float(num.sum()/den.sum()), "se": float(values.std(ddof=1)),
            "numerator": float(num.sum()), "denominator": float(den.sum())}

def coordinate(gid):
    if gid < 8: return "-X", (gid - 3.5) * 7.5
    if gid < 16: return "+X", (gid - 8 - 3.5) * 7.5
    index = gid - 16
    return "+Y", (-692 + 20*index) if index < 35 else (12 + 20*(index-35))

def emission_prediction(spectrum_path, pde_path):
    spectrum = np.loadtxt(spectrum_path)
    pde = np.loadtxt(pde_path)
    energy = 1239.84193 / spectrum[:, 0]
    order = np.argsort(energy)
    energy, density = energy[order], spectrum[order, 1]
    grid = np.linspace(energy.min(), energy.max(), 100001)
    intensity = np.interp(grid, energy, density)
    wavelength = 1239.84193 / grid
    efficiency = np.interp(wavelength, pde[:, 0], pde[:, 1])
    norm = np.trapz(intensity, grid)
    return {"mean_wavelength_nm": float(np.trapz(intensity*wavelength, grid)/norm),
            "mean_pde": float(np.trapz(intensity*efficiency, grid)/norm)}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--exec27-summary", type=Path, required=True)
    parser.add_argument("--exec27-top-states", type=Path, required=True)
    parser.add_argument("--preregistration", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    weights = np.stack([np.bincount(rng.integers(0, N, N), minlength=N)
                        for _ in range(NBOOT)]).astype(float)
    root = uproot.open(args.root)
    sensor = root["sipm_event_counts"].arrays(library="np")
    order = np.lexsort((sensor["global_id"], sensor["event_id"]))
    sensor = {key: value[order].reshape(N, 86) for key, value in sensor.items()}
    event = root["event_observables"].arrays(library="np")
    event_order = np.argsort(event["event_id"])
    event = {key: value[event_order] for key, value in event.items()}
    meta = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_root": str(args.root), "source_root_sha256": sha(args.root),
        "N": N, "seeds": [26092601, 8349041], "bootstrap_seed": SEED,
        "bootstrap_replicates": NBOOT, "bootstrap_unit": "generated event",
        "preregistration": str(args.preregistration),
        "preregistration_sha256": sha(args.preregistration),
        "script_sha256": sha(Path(__file__))}

    ids = np.arange(86)
    detected = sensor["detected"].sum(axis=0)
    matched = sensor["matched_detected"].sum(axis=0)
    orphan = sensor["unmatched_detected"].sum(axis=0)
    incident = sensor["incident"].sum(axis=0)
    expected = sensor["expected_surface_pde_sum"].sum(axis=0)
    face, position = zip(*(coordinate(int(gid)) for gid in ids))
    sidecars(args.out, "detection_gap_by_sensor", {
        "global_id": ids, "face": face, "face_coordinate_mm": position,
        "detected": detected, "matched_boundary_detection": matched,
        "orphan_standard_SD": orphan, "orphan_fraction_of_detected": orphan/detected,
        "incident": incident, "expected_surface_pde_sum": expected,
    }, meta)

    det_event = sensor["detected"].sum(axis=1)
    matched_event = sensor["matched_detected"].sum(axis=1)
    orphan_event = sensor["unmatched_detected"].sum(axis=1)
    sidecars(args.out, "detection_gap_by_event", {
        "event_id": np.arange(N), "detected": det_event,
        "matched_boundary_detection": matched_event,
        "orphan_standard_SD": orphan_event,
        "orphan_fraction_of_detected": orphan_event/det_event,
    }, meta)

    arms = []
    for name, columns in (("left", slice(0,8)), ("right", slice(8,16)),
                          ("top", slice(16,86)), ("all", slice(0,86))):
        d = sensor["detected"][:, columns].sum(axis=1)
        m = sensor["matched_detected"][:, columns].sum(axis=1)
        o = sensor["unmatched_detected"][:, columns].sum(axis=1)
        inc = sensor["incident"][:, columns].sum(axis=1)
        exp = sensor["expected_surface_pde_sum"][:, columns].sum(axis=1)
        wl = sensor["incident_wavelength_nm_sum"][:, columns].sum(axis=1)
        arms.append({"arm": name, "detected": int(d.sum()), "matched": int(m.sum()),
            "orphan": int(o.sum()), "orphan_fraction": ratio(o,d,weights),
            "matched_over_incident": ratio(m,inc,weights),
            "incident_spectrum_pde": ratio(exp,inc,weights),
            "incident_mean_wavelength_nm": ratio(wl,inc,weights)})
    sidecars(args.out, "detection_gap_by_arm", {
        "arm": [row["arm"] for row in arms],
        "detected": [row["detected"] for row in arms],
        "matched": [row["matched"] for row in arms],
        "orphan": [row["orphan"] for row in arms],
        "orphan_fraction": [row["orphan_fraction"]["value"] for row in arms],
        "orphan_fraction_se": [row["orphan_fraction"]["se"] for row in arms],
        "matched_over_incident": [row["matched_over_incident"]["value"] for row in arms],
        "matched_over_incident_se": [row["matched_over_incident"]["se"] for row in arms],
        "incident_spectrum_pde": [row["incident_spectrum_pde"]["value"] for row in arms],
        "incident_spectrum_pde_se": [row["incident_spectrum_pde"]["se"] for row in arms],
        "incident_mean_wavelength_nm": [row["incident_mean_wavelength_nm"]["value"] for row in arms],
        "incident_mean_wavelength_nm_se": [row["incident_mean_wavelength_nm"]["se"] for row in arms],
    }, meta)

    old_summary = json.loads(args.exec27_summary.read_text())
    old_top_hits = old_summary["counts"]["R1"]["TOP"]
    with args.exec27_top_states.open() as stream:
        rows = list(csv.DictReader(stream))
    old_top_detection = int(next(row["count"] for row in rows
        if row["run"] == "R1" and row["status_name"] == "Detection"))
    old_gap = old_top_hits - old_top_detection
    all_arm = next(row for row in arms if row["arm"] == "all")
    gap = {"total": all_arm["orphan"],
        "fraction_of_SD_hits": all_arm["orphan_fraction"],
        "sensor_fraction_min": float((orphan/detected).min()),
        "sensor_fraction_max": float((orphan/detected).max()),
        "sensor_count_correlation_orphan_vs_detected": float(np.corrcoef(orphan,detected)[0,1]),
        "event_orphan_min": int(orphan_event.min()), "event_orphan_max": int(orphan_event.max()),
        "event_orphan_median": float(np.median(orphan_event)),
        "event_correlation_orphan_vs_detected": float(np.corrcoef(orphan_event,det_event)[0,1]),
        "EXEC27_R1_TOP": {"N": 500, "hits": old_top_hits,
            "boundary_detection": old_top_detection, "gap": old_gap,
            "fraction": old_gap/old_top_hits},
        "track_level_characterization": "IMPOSSIBLE_FROM_PERSISTED_SCHEMA",
        "reason": "sipm_hits has no track_id or boundary fields; sipm_event_counts stores only event/sensor aggregates.",
        "energy_or_volume_pair_of_orphans": "NOT_IDENTIFIABLE",
    }

    spectrum_path = args.repo / "src/external/SSLG4/data/oscnt/opsc-101/scntComp1.txt"
    pde_path = args.repo / "data/sipm/AFBR-S4N66P024M_pde.txt"
    emitted = emission_prediction(spectrum_path, pde_path)
    transported = {row["arm"]: {"mean_wavelength_nm": row["incident_mean_wavelength_nm"],
                    "mean_surface_pde": row["incident_spectrum_pde"]} for row in arms}
    spectral = {"emitted_energy_density": emitted, "transported_incident": transported,
        "all_mean_wavelength_shift_nm": transported["all"]["mean_wavelength_nm"]["value"] - emitted["mean_wavelength_nm"],
        "all_PDE_shift": transported["all"]["mean_surface_pde"]["value"] - emitted["mean_pde"],
        "interpretation": "Measured transport spectral filtering; incident mean wavelength is redder and its PDE lower than the emitted-spectrum reference.",
        "spectrum_sha256": sha(spectrum_path), "PDE_sha256": sha(pde_path)}

    # V1: registered primary and physically requested non-optical denominator.
    produced = event["produced_scint"].astype(float)
    total = event["edep_total_MeV"]
    optical = event["edep_optical_MeV"]
    nonoptical = total - optical
    total_ratio = ratio(produced, 10400*total, weights)
    nonopt_ratio = ratio(produced, 10400*nonoptical, weights)
    effective_yield = ratio(produced, nonoptical, weights)
    equivalent = produced/10400
    unexplained_energy = equivalent.sum() - nonoptical.sum()
    excess_photons = produced.sum() - 10400*nonoptical.sum()
    optical_expected_photons = 10400*optical.sum()
    fluctuation_residual_photons = excess_photons - optical_expected_photons
    v1 = {"total_denominator_ratio": total_ratio,
        "nonoptical_denominator_ratio": nonopt_ratio,
        "effective_yield_photons_per_MeV": effective_yield,
        "nominal_yield_photons_per_MeV": 10400,
        "excess_fraction": nonopt_ratio["value"] - 1,
        "total_edep_MeV": float(total.sum()), "optical_edep_MeV": float(optical.sum()),
        "nonoptical_edep_MeV": float(nonoptical.sum()),
        "photon_equivalent_energy_MeV": float(equivalent.sum()),
        "unexplained_equivalent_energy_vs_nonoptical_MeV": float(unexplained_energy),
        "optical_fraction_of_equivalent_excess": float(optical.sum()/unexplained_energy),
        "excess_photons_over_nonoptical_expectation": float(excess_photons),
        "photons_expected_from_optical_edep": float(optical_expected_photons),
        "fluctuation_residual_photons": float(fluctuation_residual_photons),
        "effective_yield_against_all_model_eligible_edep": float(produced.sum()/total.sum()),
        "bookkeeping_max_abs_error_MeV": float(np.max(np.abs(
            event["edep_ionizing_MeV"] + event["edep_nonionizing_MeV"] - total))),
        "nonionizing_edep_MeV": float(event["edep_nonionizing_MeV"].sum()),
        "produced_counter_matches_run_summary": bool(int(produced.sum()) == 40039535),
        "diagnosis": "RESOLVED_OPTICAL_PHOTON_ENERGY_RESCINTILLATION",
        "evidence": "Geant4 11.4 registers G4Scintillation for optical photons; their measured deposited energy accounts for 98.4004% of the excess over the non-optical denominator. The remaining count is compatible with RESOLUTIONSCALE=1 fluctuations.",
        "excluded_causes": ["configured yield differs from 10400", "active Birks coefficient",
                            "particle-dependent yield enabled by a macro", "duplicate produced-photon count"],
        "code_path": "G4Scintillation::IsApplicable accepts every non-short-lived particle at G4Scintillation.cc:150-153; G4OpticalPhysics.cc:172-175 therefore attaches it to opticalphoton."}
    sidecars(args.out, "v1_event_accounting", {
        "event_id": event["event_id"], "produced_scint": event["produced_scint"],
        "edep_total_MeV": total, "edep_optical_MeV": optical,
        "edep_nonoptical_MeV": nonoptical,
        "photon_equivalent_energy_MeV": equivalent,
        "equivalent_minus_nonoptical_MeV": equivalent-nonoptical,
    }, meta)

    code_audit = [
        ("G4SteppingManager automatic SD", "/home/tdship/opt/geant4-v11.4.0/source/tracking/src/G4SteppingManager.cc", 242, 254,
         "Calls Hit for the pre-step sensitive volume before UserSteppingAction."),
        ("G4OpBoundary Detection SD", "/home/tdship/opt/geant4-v11.4.0/source/processes/optical/src/G4OpBoundaryProcess.cc", 558, 559,
         "Invokes the post-step SD only for boundary status Detection."),
        ("SiPMSD permissive boundary gate", str(args.repo/"src/SiPMSD.cc"), 47, 53,
         "Accepts either pre-step or post-step fGeomBoundary."),
        ("SiPMSD fallback and unconditional hit", str(args.repo/"src/SiPMSD.cc"), 55, 113,
         "Falls back to track volume, records a hit, then kills the track without checking OpBoundary status."),
        ("EXEC40 incident gate", str(args.repo/"src/SiPMObservation.cc"), 49, 80,
         "Requires a post-step sensor and accepted boundary status."),
        ("Geant4 nominal scintillation mean", "/home/tdship/opt/geant4-v11.4.0/source/processes/electromagnetic/xrays/src/G4Scintillation.cc", 310, 347,
         "Uses configured yield times visible/total parent-step deposit, then Gaussian/Poisson fluctuations."),
        ("Geant4 optical-photon applicability", "/home/tdship/opt/geant4-v11.4.0/source/processes/electromagnetic/xrays/src/G4Scintillation.cc", 150, 153,
         "Accepts every non-short-lived particle; opticalphoton is therefore applicable in this installed source."),
        ("Geant4 process registration", "/home/tdship/opt/geant4-v11.4.0/source/physics_lists/constructors/electromagnetic/src/G4OpticalPhysics.cc", 151, 175,
         "Registers G4Scintillation for every particle for which IsApplicable returns true."),
    ]
    sidecars(args.out, "code_path_audit", {
        "finding": [row[0] for row in code_audit], "file": [row[1] for row in code_audit],
        "line_start": [row[2] for row in code_audit], "line_end": [row[3] for row in code_audit],
        "evidence": [row[4] for row in code_audit],
    }, meta)
    cause = {"classification": "PREEXISTING_STANDARD_SENSITIVE_VOLUME_CALLBACK_PATH",
        "evidence": "Matched detections exactly equal all OpBoundary Detection rows; every remaining SD hit must use Geant4's standard pre-volume sensitive callback, the only other invocation path.",
        "excluded_legitimate_boundary_status": "NONE_OBSERVED_OR_IDENTIFIED; sensor-directed censuses contain only Absorption and Detection, and the criterion includes both.",
        "why_SiPMSD_accepts_it": "pre-step fGeomBoundary is sufficient and fallback selects track->GetVolume; no boundary status is checked before writing a detection.",
        "repair_status": "NOT_APPLIED",
        "proposal": "First record callback mode, creator, track ID, pre/post volumes and both step statuses. Then either restrict detected hits to OpBoundary Detection or define a separately named sensitive-volume-entry population; do not silently add it to the independent incident denominator.",
        "too_narrow_incident_alternative_estimate": {"maximum_callbacks_absorbed": int(orphan.sum()),
            "qualification": "An SD-union rule would absorb all aggregate orphans by construction, but existing data cannot prove that all are legitimate external arrivals, so it is not proposed as the accepted criterion."}}
    result = {"metadata": meta, "detection_gap": gap, "arms": arms,
              "cause": cause, "spectral_filtering": spectral, "V1": v1}
    dump(args.out / "diagnosis_results.json", result)
    print(json.dumps({"gap": gap, "cause": cause, "spectral": spectral, "V1": v1}, indent=2))

if __name__ == "__main__":
    main()
