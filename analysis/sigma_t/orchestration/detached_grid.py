#!/usr/bin/env python3
"""EXEC34R: simulation-only, durable 6 processes x 4 Geant4 workers.

No timeout exists. A future optional timeout must use a measurement at FOUR
workers and the same event count, never the 24-worker EXEC34A pilot wall time.
prepare/dry-run/status never execute the simulation binary. Only launch starts
an independent session (setsid); all output is redirected before detachment.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import difflib
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import uuid

from run_simulation import memory, now, sha

HERE = Path(__file__).resolve()
SOURCE = Path('/home/rrios/exec34b_20260911/campaign.json')
HANDOFF = Path('/home/rrios/EXEC34_HANDOFF_20260911.json')
SCALING = Path('/home/rrios/exec33_20260911/scaling.json')
REPRO = Path('/home/rrios/exec33_20260911/event_reproducibility.json')
DEFAULT = Path('/home/rrios/exec34r_20260912')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    """Atomic, fsynced disk state; never truncate a previous attempt's data."""
    path = Path(path)
    tmp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    with tmp.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)
    fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def append(path, row):
    with Path(path).open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        stream.write(json.dumps(dict(utc=now(), **row)) + '\n')
        stream.flush()
        os.fsync(stream.fileno())


def adjusted_macro(raw):
    """Replace only numeric thread/modulo tokens, preserving all other bytes."""
    for command, value in (('/run/numberOfThreads', '4'), ('/run/eventModulo', '1')):
        pattern = rb'(?m)^([ \t]*' + re.escape(command.encode()) + rb'[ \t]+)\d+([ \t]*(?:#[^\r\n]*)?\r?$)'
        raw, count = re.subn(pattern, lambda m: m[1] + value.encode() + m[2], raw)
        require(count == 1, 'Expected exactly one ' + command)
    return raw


def macro_value(raw, command):
    values = re.findall(r'^\s*' + re.escape(command) + r'[ \t]+([^\r\n#]+)', raw.decode(), re.M)
    require(len(values) == 1, 'Missing/duplicate macro command ' + command)
    return values[0].strip()


def prepare(dest, source=SOURCE, handoff=HANDOFF, scaling=SCALING, repro=REPRO):
    cfg, h = read(source), read(handoff)
    code = h.get('EJ200_OPSC_CODE')
    require(code == 'NOT_FOUND' or isinstance(code, str) and code.startswith('OPSC-'),
            'Missing explicit EJ200_OPSC_CODE in handoff')
    rows = read(scaling)
    four = next(r for r in rows if r['workers'] == 4)
    require(four['status'] == 'PASS', 'Four-worker benchmark failed')
    require({r['workers'] for r in rows} == {1, 4, 12, 24} and
            all(r['delta_S1'] == 0 for r in rows), 'Scaling N_pe evidence mismatch')
    require(all(r['different_events_from_S1'] == 0 for r in read(repro)),
            'Event reproducibility evidence mismatch')
    planned = []
    excluded = []
    for cell in cfg['cells']:
        if cell['material'] == 'EJ-200' and code == 'NOT_FOUND':
            excluded.append(cell['cell_id'])
            continue
        require(cell['material'] != 'EJ-200' or cell['opsc'] == code, 'EJ200 OPSC mismatch')
        src = Path(cell['output'])/'run.mac'
        raw = src.read_bytes()
        modified = adjusted_macro(raw)
        require(macro_value(raw, '/run/beamOn') == '10000', 'Wrong event count')
        require(macro_value(raw, '/det/scintillator') == cell['opsc'], 'Wrong material')
        require(re.fullmatch(r'[A-Za-z0-9_]+', cell['cell_id']), 'Unsafe cell id')
        planned.append((cell, src, raw, modified))
    require(len(planned) == (14 if code == 'NOT_FOUND' else 21), 'Unexpected job count')
    require(len({c['cell_id'] for c, *_ in planned}) == len(planned), 'Duplicate cells')
    # Preparation does not modify source macros or the EXEC34B manifest/partials.
    require(not (dest/'campaign.json').exists(), 'Already prepared; use dry-run/status/launch')
    dest.mkdir(parents=True, exist_ok=True)
    evidence = {str(p): sha(p) for p in (source, handoff, scaling, repro)}
    out = dict(schema='EXEC34R simulation launcher v1', prepared_utc=now(),
               output_directory=str(dest), manifest=str(dest/'manifest.jsonl'),
               binary=cfg['binary'], binary_sha256=cfg['binary_sha256'],
               simulation_commit=cfg['simulation_commit'],
               pilot_ROOT=h['pilot_ROOT'], pilot_ROOT_size=Path(h['pilot_ROOT']).stat().st_size,
               PDE_path=h['PDE_path'], PDE_sha256=h['PDE_sha256'],
               geant4_version=h['Geant4_version'], evidence_sha256=evidence,
               runtime_environment={k: v for k, v in os.environ.items()
                                    if k.startswith('G4') or k == 'LD_LIBRARY_PATH'},
               handoff=str(handoff), EJ200_OPSC_CODE=code, excluded=excluded,
               concurrency=6, workers=4, eventModulo=1, N_generated=10000,
               timeout_s=None, cells=[],
               expected_rss_per_process_bytes=int(four['max_rss_kib'])*1024,
               # 4-worker benchmark is N=2000. Budget 2x the larger measured
               # simulation RSS (including N=10000 pilot) plus 256 MiB for
               # streaming ROOT validation, per active cell; reserve 20% RAM.
               memory_budget_per_process_bytes=2*max(int(four['max_rss_kib'])*1024,
                   h['pilot_configuration']['max_rss_bytes']) + 256*1024**2,
               parallelism_provenance=dict(source='EXEC33 S1/S2/S3/S4',
                   scaling=str(scaling), reproducibility=str(repro), workers_tested=[1,4,12,24],
                   exact_delta_npe_end=0, different_events=0,
                   scope='Measured EXEC33 benchmark; not a new test of all grid cells',
                   previous_grid_workers=1, pilot_workers=24, selected_workers=4,
                   efficiency_four_workers=four['efficiency']))
    for cell, src, raw, modified in planned:
        target = dest/'cells'/cell['cell_id']
        target.mkdir(parents=True, exist_ok=False)
        (target/'run.mac').write_bytes(modified)
        (target/'sslg4').symlink_to(Path(cfg['binary']).parent/'sslg4', target_is_directory=True)
        out['cells'].append(dict(cell, source_macro=str(src), source_macro_sha256=sha(src),
                                 output=str(target), macro_sha256=sha(target/'run.mac')))
    (dest/'grid_logs').mkdir(exist_ok=True)
    write(dest/'campaign.json', out)
    cell, src, raw, modified = planned[0]
    print(''.join(difflib.unified_diff(raw.decode().splitlines(True), modified.decode().splitlines(True),
                  fromfile=str(src), tofile=str(dest/'cells'/cell['cell_id']/'run.mac'))))
    print(f'Prepared {len(planned)} cells; excluded: {excluded}; NO SIMULATION STARTED')


def preflight(cfg):
    require((cfg['concurrency'], cfg['workers'], cfg['eventModulo'], cfg['N_generated'], cfg['timeout_s'])
            == (6, 4, 1, 10000, None), 'Frozen EXEC34R configuration changed')
    for p, digest in cfg['evidence_sha256'].items():
        require(sha(p) == digest, 'Evidence changed: ' + p)
    code = read(cfg['handoff'])['EJ200_OPSC_CODE']
    require(code == cfg['EJ200_OPSC_CODE'], 'Handoff OPSC mismatch')
    require(len(cfg['cells']) == (14 if code == 'NOT_FOUND' else 21), 'Job count mismatch')
    binary = Path(cfg['binary'])
    require(os.access(binary, os.X_OK) and sha(binary) == cfg['binary_sha256'], 'Binary mismatch')
    require('EJ200_ENABLE_DIAGNOSTICS:BOOL=OFF' in (binary.parent/'CMakeCache.txt').read_text(),
            'Diagnostics OFF not verified')
    require((binary.parent/'sslg4').is_dir(), 'Missing sslg4 runtime data')
    require(sha(cfg['PDE_path']) == cfg['PDE_sha256'], 'PDE changed')
    for path, digest in cfg.get('material_resources_sha256', {}).items():
        require(Path(path).is_file() and sha(path) == digest,
                'Material resource missing/changed: ' + path)
    for cell in cfg['cells']:
        raw = (Path(cell['output'])/'run.mac').read_bytes()
        src = Path(cell['source_macro'])
        require(sha(src) == cell['source_macro_sha256'], 'Source macro changed')
        require(raw == adjusted_macro(src.read_bytes()) and sha(Path(cell['output'])/'run.mac') == cell['macro_sha256'],
                'Macro differs beyond authorized tokens: ' + cell['cell_id'])
        require(cell['material'] != 'EJ-200' or code != 'NOT_FOUND' and cell['opsc'] == code,
                'EJ200 not authorized by handoff')
        local_sslg4 = Path(cell['output'])/'sslg4'
        runtime_sslg4 = Path(cell.get('sslg4_runtime', binary.parent/'sslg4')).resolve()
        require(runtime_sslg4.is_dir(), 'Missing cell SSLG4 runtime: ' + cell['cell_id'])
        require(local_sslg4.is_symlink() and local_sslg4.resolve() == runtime_sslg4,
                'Missing/wrong local sslg4 link: ' + cell['cell_id'])
    import uproot  # Fail before detaching if ROOT reader dependencies are absent.
    require(Path(cfg['pilot_ROOT']).stat().st_size == cfg['pilot_ROOT_size'], 'Pilot size changed')
    pilot_events = cfg.get('pilot_N_generated', cfg['N_generated'])
    require(isinstance(pilot_events, int) and pilot_events > 0, 'Invalid pilot event count')
    projected = int(cfg['pilot_ROOT_size']*cfg['N_generated']/pilot_events*len(cfg['cells'])*1.3)
    free = shutil.disk_usage(cfg['output_directory']).free
    available = memory()['MemAvailable']
    expected = cfg['expected_rss_per_process_bytes']*6
    budget = cfg['memory_budget_per_process_bytes']*6
    require(free >= projected, f'Disk: {free} available < {projected} projected')
    require(available*.8 >= budget, f'RAM: {available} available < budget {budget} plus 20% reserve')
    fresh = sum(not (Path(cell['output'])/'attempts').exists() for cell in cfg['cells'])
    if not Path(cfg['manifest']).exists():
        require(fresh == len(cfg['cells']),
                'A first launch requires nonexistent attempt output directories')
    return dict(utc=now(), status='PASS', disk_projected_bytes=projected,
                disk_projection_rule='pilot_ROOT_size * target_N/pilot_N * job_count * 1.3',
                pilot_ROOT_size=cfg['pilot_ROOT_size'], pilot_N_generated=pilot_events,
                target_N_generated=cfg['N_generated'], job_count=len(cfg['cells']),
                disk_available_bytes=free,
                expected_six_process_rss_bytes=expected, conservative_six_cell_budget_bytes=budget,
                memory_available_bytes=available, macros_readable=len(cfg['cells']),
                fresh_attempt_output_directories=fresh,
                diagnostics=False, binary_sha256=cfg['binary_sha256'], timeout_s=None)


def identity(pid):
    try:
        stat = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
        if stat[0] == 'Z':
            return None
        return dict(pid=pid, start_ticks=stat[19], boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip())
    except (OSError, IndexError):
        return None


def alive(token):
    return bool(token) and identity(token['pid']) == token


def done_record(cell, verify_hash=True):
    marker = Path(cell['output'])/'.DONE'
    if not marker.exists():
        return None
    row = read(marker)
    require(row['macro_sha256'] == cell['macro_sha256'] and row['events_run'] == 10000 and row['exit_code'] == 0,
            'Invalid DONE marker: ' + cell['cell_id'])
    root = Path(row['root_path'])
    require(root.is_file() and root.stat().st_size == row['root_size_bytes'], 'DONE ROOT missing/size changed')
    if verify_hash:
        require(sha(root) == row['root_sha256'], 'DONE ROOT hash changed')
    return row


def state(cell):
    try:
        if done_record(cell, verify_hash=False):
            return 'complete'
        path = Path(cell['output'])/'state.json'
        if not path.exists():
            return 'pending'
        row = read(path)
        if row['status'] == 'RUNNING' and (alive(row.get('process')) or alive(row.get('driver'))):
            return 'running'
        return 'failed'
    except (OSError, ValueError, KeyError):
        return 'failed'


def status(cfg):
    result = {s: [] for s in ('complete', 'running', 'pending', 'failed')}
    for cell in cfg['cells']:
        result[state(cell)].append(cell['cell_id'])
    return dict(counts={k: len(v) for k, v in result.items()}, cells=result,
                note='Disk markers + PID/start time/boot identity; launch additionally rehashes DONE ROOTs')


def verify_output(root, log):
    """Generated denominator comes from the master run summary, NOT unique hits.

The native ROOT has only sipm_hits, so zero-hit events cannot be counted there.
Read every branch/basket; require all hit IDs in the generated 0..9999 range.
Cross-check all three face totals and unique hit-event count against the log.
"""
    import numpy as np
    import uproot
    text = Path(log).read_text(errors='replace')
    def total(pattern):
        values = re.findall(pattern, text)
        require(len(values) == 1, 'Missing/ambiguous master summary: ' + pattern)
        return int(values[0])
    events = total(r'Events run\s*:\s*(\d+)')
    require(events == 10000, 'Run did not generate 10000 events')
    expected = [total(r'End-left\s+photons\s*:\s*(\d+)'),
                total(r'End-right\s+photons\s*:\s*(\d+)'),
                total(r'Top SiPM\s+photons\s*:\s*(\d+)')]
    expected_hit_events = total(r'Events with ≥1 hit\s*:\s*(\d+)')
    require(Path(root).is_file() and Path(root).stat().st_size > 0, 'Missing ROOT')
    seen, counts, entries = np.zeros(events, dtype=bool), np.zeros(3, dtype=np.int64), 0
    with uproot.open(root) as file:
        tree = file['sipm_hits']
        require({'event_id', 'face_type', 'global_id', 'time_ns'} <= set(tree.keys()), 'Missing ROOT branches')
        for arrays in tree.iterate(step_size='32 MB', library='np'):
            ids, faces = arrays['event_id'], arrays['face_type']
            require(np.all((ids >= 0) & (ids < events)), 'ROOT event ID outside generated range')
            require(np.all(np.isin(faces, [0, 1, 2])), 'Invalid face ID')
            require(np.all(np.isfinite(arrays['time_ns'])), 'Nonfinite hit times')
            seen[ids] = True
            counts += np.bincount(faces, minlength=3)
            entries += len(ids)
        require(entries == tree.num_entries, 'Incomplete ROOT read')
    require(counts.tolist() == expected, 'ROOT/log photon totals disagree')
    require(int(seen.sum()) == expected_hit_events, 'ROOT/log hit-event totals disagree')
    return dict(events_run=events, generated_count_source='master Events run summary, ROOT hit IDs/totals cross-checked',
                root_entries=entries, events_with_hits=int(seen.sum()), left_total=expected[0],
                right_total=expected[1], top_total=expected[2], npe_end_mean=sum(expected[:2])/(2*events),
                root_path=str(root), root_sha256=sha(root), root_size_bytes=Path(root).stat().st_size)


def run_cell(cfg, cell, lock_fd):
    base = Path(cell['output'])
    # A new directory per attempt means relaunch never truncates a partial ROOT.
    attempt = base/'attempts'/uuid.uuid4().hex
    attempt.mkdir(parents=True)
    shutil.copyfile(base/'run.mac', attempt/'run.mac')
    runtime_sslg4 = Path(cell.get(
        'sslg4_runtime', Path(cfg['binary']).parent/'sslg4')).resolve()
    (attempt/'sslg4').symlink_to(runtime_sslg4, target_is_directory=True)
    logfile = Path(cfg['output_directory'])/'grid_logs'/(cell['cell_id']+'.log')
    argv = [cfg['binary'], '-m', str(attempt/'run.mac')]
    raw = (attempt/'run.mac').read_bytes()
    # Reproduce the prepared Geant4 environment even from a different SSH
    # shell. Do not let a new G4FORCENUMBEROFTHREADS or PDE override leak in.
    runtime_env = dict(cfg['runtime_environment'], EJ200_DATA_DIR=str(Path(cfg['PDE_path']).parent.parent))
    env = {k: v for k, v in os.environ.items()
           if not k.startswith('G4') and k not in ('LD_LIBRARY_PATH', 'EJ200_DATA_DIR')}
    row = dict(cell_id=cell['cell_id'], status='RUNNING', stage='simulation',
               output=str(attempt), cwd=str(attempt), cell_directory=str(base),
               command=shlex.join(argv), argv=argv, workers=4, eventModulo=1,
               exit_code=None, wall_s=None, root_sha256=None, start_utc=now(),
               driver=identity(os.getpid()), log=str(logfile), timeout_s=None,
               material=cell['material'], opsc_code=cell['opsc'], x_mm=cell['x_mm'],
               seeds=[int(v) for v in macro_value(raw, '/random/setSeeds').split()],
               N_generated=10000, configuration='EndTop', N_TOP=70, sptr_ns=0,
               binary=cfg['binary'], binary_sha256=cfg['binary_sha256'], diagnostics=False,
               sslg4_runtime=str(runtime_sslg4),
               optical_model_status=cell.get('optical_model_status', 'UNSPECIFIED'),
               macro_sha256=cell['macro_sha256'], source_macro_sha256=cell['source_macro_sha256'],
               simulation_commit=cfg['simulation_commit'], geant4_version=cfg['geant4_version'],
               PDE_path=cfg['PDE_path'], PDE_sha256=cfg['PDE_sha256'],
               parallelism_provenance=cfg['parallelism_provenance'],
               evidence_sha256=cfg['evidence_sha256'], launcher_sha256=sha(HERE),
               runtime_environment=runtime_env,
               rng='Geant4 per-event seeds; fixed master seeds; eventModulo=1',
               thread_limits=dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    started = time.monotonic()
    rc = None
    try:
        with (attempt/'stdout.log').open('x') as stream:
            # The shared driver lock is inherited by the binary. If the driver
            # dies, surviving simulations still prevent a duplicate launch.
            process = subprocess.Popen(argv, cwd=attempt, stdin=subprocess.DEVNULL,
                stdout=stream, stderr=subprocess.STDOUT,
                env=dict(env, **(row['runtime_environment'] | row['thread_limits'])), pass_fds=(lock_fd,))
            row['process'] = identity(process.pid)
            write(base/'state.json', row)
            write(attempt/'simulation.meta.json', row)
            append(cfg['manifest'], row)
            # Append-only public log is fed from a durable per-attempt stdout.
            # Copy incrementally while waiting, with NO deadline or timeout.
            with logfile.open('a') as public, (attempt/'stdout.log').open() as reader:
                public.write(f'\nEXEC34R attempt {attempt.name} {row["start_utc"]} {row["command"]}\n')
                while process.poll() is None:
                    public.write(reader.read())
                    public.flush()
                    time.sleep(1)
                public.write(reader.read())
                public.flush()
            rc = process.wait()
        row.update(exit_code=rc, wall_s=time.monotonic()-started, end_utc=now())
        require(rc == 0, f'Simulation exit {rc}')
        pdes = re.findall(r'SiPM PDE file\s*:\s*(.*)', (attempt/'stdout.log').read_text(errors='replace'))
        require(pdes and all(Path(p.strip()).resolve() == Path(cfg['PDE_path']).resolve() for p in pdes),
                'Simulator PDE path differs from pinned provenance')
        require(sha(cfg['PDE_path']) == cfg['PDE_sha256'], 'PDE changed during simulation')
        row.update(verify_output(attempt/'photon_hits_run000.root', attempt/'stdout.log'))
        # New status avoids redefining EXEC34B PASS (which included analysis).
        row['status'] = 'SIMULATION_COMPLETE'
    except Exception as error:
        row.update(status='FAILED', error=str(error), exit_code=rc,
                   wall_s=time.monotonic()-started, end_utc=now())
        partial = attempt/'photon_hits_run000.root'
        if partial.exists() and rc is not None:
            row.update(root_path=str(partial), root_sha256=sha(partial),
                       root_size_bytes=partial.stat().st_size, partial_output_preserved=True)
    write(attempt/'simulation.meta.json', row)
    write(base/'state.json', row)
    append(cfg['manifest'], row)
    if row['status'] == 'SIMULATION_COMPLETE':
        write(base/'.DONE', row)
    with logfile.open('a') as stream:
        stream.write(json.dumps(dict(utc=now(), status=row['status'], error=row.get('error'), wall_s=row['wall_s']))+'\n')
    print(f'{now()} {cell["cell_id"]} {row["status"]}', flush=True)
    return row['status'] == 'SIMULATION_COMPLETE'


def acquire(dest):
    lock = (dest/'driver.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise ValueError('Driver or surviving simulations already hold campaign lock')
    return lock


def run_driver(cfg, lock_fd):
    write(Path(cfg['output_directory'])/'preflight.json', preflight(cfg))
    todo = [c for c in cfg['cells'] if not done_record(c)]
    print(f'{now()} START 6x4 no timeout; pending={len(todo)}', flush=True)
    with ThreadPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(lambda c: run_cell(cfg, c, lock_fd), todo))
    print(f'{now()} END {json.dumps(status(cfg))}', flush=True)
    return 0 if all(results) else 34


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'dry-run', 'status', 'launch', '_driver'])
    parser.add_argument('--directory', type=Path, default=DEFAULT)
    parser.add_argument('--lock-fd', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    dest = args.directory.resolve()
    if args.action == 'prepare':
        prepare(dest)
        return 0
    cfg = read(dest/'campaign.json')
    require(cfg['output_directory'] == str(dest), 'Campaign directory mismatch')
    if args.action == 'status':
        print(json.dumps(status(cfg), indent=2))
        return 0
    if args.action == '_driver':
        require(args.lock_fd is not None, 'Internal driver requires inherited lock')
        os.fstat(args.lock_fd)
        return run_driver(cfg, args.lock_fd)
    check = preflight(cfg)
    print(json.dumps(check, indent=2), flush=True)
    if args.action == 'dry-run':
        for cell in cfg['cells']:
            print(json.dumps(dict(cell_id=cell['cell_id'], material=cell['material'], opsc=cell['opsc'],
                x_mm=cell['x_mm'], workers=4, eventModulo=1, events=10000, state=state(cell),
                macro=str(Path(cell['output'])/'run.mac'),
                command_template=shlex.join([cfg['binary'], '-m', str(Path(cell['output'])/'attempts/<unique-attempt>/run.mac')]))))
        print(f'JOBLIST {len(cfg["cells"])}; excluded={cfg["excluded"]}; EXECUTED=0', flush=True)
        return 0
    with acquire(dest) as lock:
        # Hash-check completed ROOTs before detach; fail closed on changed data.
        todo = [c for c in cfg['cells'] if not done_record(c)]
        if not todo:
            print('All cells DONE; nothing launched')
            return 0
        write(dest/'preflight.json', check)
        with (dest/'grid_logs/driver.log').open('a') as log:
            child = subprocess.Popen([sys.executable, '-u', str(HERE), '_driver', '--directory', str(dest),
                '--lock-fd', str(lock.fileno())], stdin=subprocess.DEVNULL, stdout=log,
                stderr=subprocess.STDOUT, start_new_session=True, pass_fds=(lock.fileno(),), cwd=dest)
        write(dest/'driver.json', dict(utc=now(), process=identity(child.pid), detached=True,
                                      method='setsid via Popen(start_new_session=True)', timeout_s=None))
        print(f'Detached driver PID {child.pid}; log={dest}/grid_logs/driver.log', flush=True)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f'ABORT: {error}', file=sys.stderr, flush=True)
        raise SystemExit(2)
