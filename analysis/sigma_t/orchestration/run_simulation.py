#!/usr/bin/env python
"""Run one native-name simulation, preserving stage provenance and direct RSS."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import threading
import time


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def memory():
    text = Path('/proc/meminfo').read_text()
    return {key: int(re.search(r'^'+key+r':\s+(\d+)', text, re.M)[1])*1024
            for key in ('MemTotal', 'MemAvailable')}


def run(args):
    dest = Path(args.output).resolve()
    dest.mkdir(parents=True, exist_ok=False)
    binary = Path(args.binary).resolve()
    cache = (binary.parent/'CMakeCache.txt').read_text()
    assert 'EJ200_ENABLE_DIAGNOSTICS:BOOL=OFF' in cache
    (dest/'sslg4').symlink_to(binary.parent/'sslg4', target_is_directory=True)
    macro = f'''/run/numberOfThreads {args.workers}
/run/eventModulo 1
/random/setSeeds 26092601 8349041
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/det/readout EndTop
/det/scintillator {args.opsc}
/sipm/model AFBR-S4N66P024M
/run/initialize
/sipm/jitterSigma 0 ns
/gun/particle mu-
/gun/energy 1 GeV
/muon/angle 0
/muon/gunX {args.x} mm
/control/echo EXEC34_PRE_BEAM
/run/beamOn {args.events}
/control/echo EXEC34_POST_BEAM
'''
    (dest/'run.mac').write_text(macro)
    command = ['/usr/bin/time', '-v', '-o', str(dest/'resource_usage.txt'),
               str(binary), '-m', str(dest/'run.mac')]
    meta = dict(stage='simulation', start_utc=now(), command=shlex.join(command),
                cwd=str(dest), simulation_commit=args.simulation_commit,
                binary=str(binary), binary_sha256=sha(binary),
                macro_sha256=sha(dest/'run.mac'), geant4_version='11.4.0',
                seeds=[26092601, 8349041], N_generated=args.events,
                material=args.material, opsc_code=args.opsc, configuration='EndTop',
                N_TOP=70, x_mm=args.x, workers=args.workers, eventModulo=1,
                diagnostics=False, sptr_ns=0, particle='mu-', energy_GeV=1,
                direction='vertical', memory_before=memory(),
                rss_method='GNU /usr/bin/time -v wait4 maximum resident set size for simulation process',
                rng='Geant4 per-event seeds; fixed master seeds; eventModulo=1; ROOT hit order not used for split')
    (dest/'simulation.meta.json').write_text(json.dumps(meta, indent=2)+'\n')
    started = time.monotonic()
    markers = {}
    with (dest/'stdout.log').open('w') as log:
        process = subprocess.Popen(command, cwd=dest, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in process.stdout:
            log.write(line)
            log.flush()
            if line.strip() in ('EXEC34_PRE_BEAM', 'EXEC34_POST_BEAM'):
                markers[line.strip()] = time.monotonic()-started
        rc = process.wait()
    meta.update(end_utc=now(), exit_code=rc, markers_elapsed_s=markers,
                startup_s=markers.get('EXEC34_PRE_BEAM'),
                startup_definition='process launch to master marker after /run/initialize, before /run/beamOn; includes geometry/physics initialization')
    usage = (dest/'resource_usage.txt').read_text()
    elapsed = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)', usage)[1]
    meta['wall_s'] = sum(float(v)*60**i for i,v in enumerate(reversed(elapsed.split(':'))))
    meta['max_rss_bytes'] = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', usage)[1])*1024
    if rc == 0:
        log = (dest/'stdout.log').read_text()
        meta['events_run'] = int(re.search(r'Events run\s*:\s*(\d+)', log)[1])
        assert meta['events_run'] == args.events
        meta['left_total'] = int(re.findall(r'End-left\s+photons\s*:\s*(\d+)', log)[-1])
        meta['right_total'] = int(re.findall(r'End-right\s+photons\s*:\s*(\d+)', log)[-1])
        meta['npe_end_mean'] = (meta['left_total']+meta['right_total'])/(2*args.events)
        meta['PDE_path'] = re.search(r'SiPM PDE file\s*:\s*(.*)', log)[1].strip()
        meta['PDE_sha256'] = sha(meta['PDE_path'])
        meta['root_path'] = str(dest/'photon_hits_run000.root')
        meta['root_sha256'] = sha(meta['root_path'])
    (dest/'simulation.meta.json').write_text(json.dumps(meta, indent=2)+'\n')
    print(json.dumps(meta, indent=2), flush=True)
    return rc


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--simulation-commit', default='420addf0fd6029d5b2f0e235f472a8ae47f31fac')
    p.add_argument('--events', type=int, default=10000)
    p.add_argument('--workers', type=int, default=24)
    p.add_argument('--material', default='EJ-204')
    p.add_argument('--opsc', default='OPSC-101')
    p.add_argument('--x', type=int, default=0)
    raise SystemExit(run(p.parse_args()))
