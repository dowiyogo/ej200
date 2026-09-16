#!/usr/bin/env python3
"""Prepare and run the authorized EXEC_46 I2 EJ204_xm650 validation cell."""

import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from datetime import datetime, timezone

import numpy as np
import uproot


BASE = Path('/home/rrios/exec46_20260915')
OUTPUT = Path('/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650')
BINARY = BASE/'build_baseline/ej200_bar_sim'
# Start from the validated 3800 mm BC-408 runtime so the result is also the
# exact combined runtime used by the subsequent production campaign.
SOURCE_SSLG4 = BASE/'f4_bc408_sensitivity/visible_current_3800mm/sslg4'
SOURCE_MACRO = BASE/'full_grid/cells/EJ204_xm650/run.mac'
EXPECTED_EVENTS = 10_000
EXPECTED_SEEDS = '26092601 8349041'
A, B, C = 1.578, 0.818, 0.00729
MEASURED_MIN_NM, MEASURED_MAX_NM, STEP_NM = 370.0, 660.0, 2.0
UV_ANCHOR_NM, UV_ABS_MM = 372.0, 26.58
VISIBLE_ANCHOR_NM, VISIBLE_ABS_MM = 439.0, 3800.0


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def n_bc404(wavelength_nm):
    return A+B*math.exp(-C*wavelength_nm)


def rindex_rows():
    measured = np.arange(MEASURED_MIN_NM, MEASURED_MAX_NM+.5*STEP_NM, STEP_NM)
    return ([(200.0, n_bc404(MEASURED_MIN_NM))]
            + [(float(wavelength), n_bc404(float(wavelength))) for wavelength in measured]
            + [(800.0, n_bc404(MEASURED_MAX_NM))])


def absorption_rows_mm():
    return [(200.0, UV_ABS_MM), (UV_ANCHOR_NM, UV_ABS_MM),
            (VISIBLE_ANCHOR_NM, VISIBLE_ABS_MM), (800.0, VISIBLE_ABS_MM)]


def write_table(path, rows, value_scale=1.0):
    with path.open('w') as stream:
        for wavelength, value in rows:
            stream.write(f'{wavelength:.2f} {value/value_scale:.12g}\n')


def prepare():
    require(not OUTPUT.exists(), f'refusing to overwrite {OUTPUT}')
    require(BINARY.is_file() and os.access(BINARY, os.X_OK), 'binary unavailable')
    require('EJ200_ENABLE_DIAGNOSTICS:BOOL=OFF' in
            (BINARY.parent/'CMakeCache.txt').read_text(), 'diagnostics must be OFF')
    OUTPUT.mkdir(parents=True)
    shutil.copyfile(SOURCE_MACRO, OUTPUT/'run.mac')
    shutil.copytree(SOURCE_SSLG4, OUTPUT/'sslg4')
    macro = (OUTPUT/'run.mac').read_text()
    for token in ('/run/numberOfThreads 4', '/run/eventModulo 1',
                  f'/random/setSeeds {EXPECTED_SEEDS}', '/run/beamOn 10000',
                  '/muon/gunX -650 mm', '/det/scintillator OPSC-101'):
        require(macro.count(token) == 1, f'macro contract failed: {token}')
    mpt = OUTPUT/'sslg4/data/oscnt/opsc-101'
    write_table(mpt/'rIndex.txt', rindex_rows())
    # opsc-101.mac declares ABSLENGTH values in cm.
    write_table(mpt/'absLength.txt', absorption_rows_mm(), value_scale=10.0)
    require(np.allclose(np.loadtxt(mpt/'rIndex.txt'), np.asarray(rindex_rows()),
                        rtol=0, atol=5e-12), 'RINDEX round-trip mismatch')
    expected_abs_cm = np.asarray([(wavelength, value/10.0)
                                  for wavelength, value in absorption_rows_mm()])
    require(np.allclose(np.loadtxt(mpt/'absLength.txt'), expected_abs_cm,
                        rtol=0, atol=5e-12), 'ABSLENGTH round-trip mismatch')
    metadata = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'cell_id': 'EJ204_xm650', 'opsc': 'OPSC-101', 'x_mm': -650,
        'events': EXPECTED_EVENTS, 'workers': 4, 'eventModulo': 1,
        'seeds': [26092601, 8349041],
        'binary': str(BINARY), 'binary_sha256': sha256(BINARY),
        'source_macro': str(SOURCE_MACRO),
        'source_macro_sha256': sha256(SOURCE_MACRO),
        'macro_sha256': sha256(OUTPUT/'run.mac'),
        'source_sslg4': str(SOURCE_SSLG4),
        'rindex': {
            'model': 'n(lambda_nm)=1.578+0.818*exp(-0.00729*lambda_nm)',
            'sample_step_nm': STEP_NM, 'measured_domain_nm': [370.0, 660.0],
            'below_domain': 'constant n(370 nm) from 200 to 370 nm',
            'above_domain': 'constant n(660 nm) from 660 to 800 nm',
            'sha256': sha256(mpt/'rIndex.txt'),
        },
        'absorption': {
            'uv_measurement': '26.58 +/- 3.0 mm at 372 nm',
            'below_372_nm': 'constant 26.58 mm',
            'between_372_and_439_nm': 'Geant4 material-table interpolation',
            'at_and_above_439_nm': 'constant 3800 mm through 800 nm',
            'sha256': sha256(mpt/'absLength.txt'),
        },
        'combined_runtime_status': {
            'EJ-200': 'CORRECTED_BC408_3800_VALIDATED_F4',
            'EJ-204': 'CORRECTED_BC404_3800_PENDING_I2',
            'EJ-230': 'UNCORRECTED_CONSTANT_RINDEX',
        },
    }
    (OUTPUT/'configuration.meta.json').write_text(
        json.dumps(metadata, indent=2, sort_keys=True)+'\n')


def run():
    command = [str(BINARY), '-m', str(OUTPUT/'run.mac')]
    environment = os.environ.copy()
    environment.update({'EJ200_DATA_DIR': '/home/rrios/ej200/data',
                        'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
                        'MKL_NUM_THREADS': '1'})
    started = datetime.now(timezone.utc)
    start = time.monotonic()
    with (OUTPUT/'simulation.log').open('xb') as log:
        process = subprocess.run(command, cwd=OUTPUT, env=environment,
                                 stdin=subprocess.DEVNULL, stdout=log,
                                 stderr=subprocess.STDOUT, check=False)
    wall = time.monotonic()-start
    require(process.returncode == 0, f'simulation exit code {process.returncode}')
    root_path = OUTPUT/'photon_hits_run000.root'
    require(root_path.is_file(), 'missing ROOT output')
    seen = np.zeros(EXPECTED_EVENTS, dtype=bool)
    with uproot.open(root_path) as root_file:
        tree = root_file['sipm_hits']
        entries = int(tree.num_entries)
        for batch in tree.iterate(['event_id'], step_size='256 MB', library='np'):
            events = batch['event_id'].astype(np.int64, copy=False)
            require(np.all((events >= 0) & (events < EXPECTED_EVENTS)), 'bad event ID')
            seen[events] = True
    require(np.all(seen), 'not all generated events represented')
    done = {
        'status': 'COMPLETE', 'start_utc': started.isoformat(),
        'end_utc': datetime.now(timezone.utc).isoformat(),
        'command': command, 'cwd': str(OUTPUT), 'exit_code': process.returncode,
        'wall_s': wall, 'events_with_hits': int(seen.sum()),
        'root_entries': entries, 'root_path': str(root_path),
        'root_size_bytes': root_path.stat().st_size,
        'root_sha256': sha256(root_path),
    }
    (OUTPUT/'.DONE.json').write_text(json.dumps(done, indent=2, sort_keys=True)+'\n')
    return done


def main():
    prepare()
    done = run()
    print(json.dumps(done, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
