#!/usr/bin/env python3
"""Run the single authorized EXEC40 cell exactly once; never launch a grid."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--d1', type=Path, required=True)
    a = p.parse_args()
    resources = a.binary.parent/'sslg4'
    required = [Path('macros/oscnt/opsc-101.mac'), *[Path('data/oscnt/opsc-101')/n
        for n in ('scntComp1.txt', 'rIndex.txt', 'absLength.txt')]]
    for rel in required:
        if sha(resources/rel) != sha(a.d1/'sslg4'/rel):
            raise SystemExit(f'SSLG4 resource differs from D1: {rel}')
    a.out.mkdir(parents=True, exist_ok=False)  # fail closed against a second run
    (a.out/'sslg4').symlink_to(resources, target_is_directory=True)
    repo = Path(__file__).resolve().parents[2]
    original = (a.d1/'run.mac').read_text()
    assert original.count('/run/numberOfThreads 1') == 1
    assert '/run/eventModulo' not in original
    macro = original.replace('/run/numberOfThreads 1', '/run/numberOfThreads 4\n/run/eventModulo 1')
    (a.out/'run.mac').write_text(macro)
    command = ['/usr/bin/time', '-v', '-o', str(a.out/'resource_usage.txt'),
               str(a.binary), '-m', str(a.out/'run.mac')]
    meta = dict(N=2000, seeds=[26092601, 8349041], workers=4, eventModulo=1,
        material='EJ-204 / OPSC-101', readout='EndTop', N_TOP=70, x_mm=0,
        particle='mu-', energy_GeV=1, angle=0, jitter_ns=0, diagnostics=False,
        base='420addf', source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        command=shlex.join(command), cwd=str(a.out), binary_sha256=sha(a.binary),
        macro_sha256=sha(a.out/'run.mac'), d1_root=str(a.d1/'photon_hits_run000.root'),
        preregistration_sha256=sha(repo/'analysis/validation/EXEC40_PREREGISTRATION.md'),
        worker_invariance_evidence='EXEC33 S1/S2/S3/S4 exact zero Npe/end difference at 1/4/12/24 workers',
        start_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    meta['sslg4_resources_sha256'] = {str(rel):sha(resources/rel) for rel in required}
    diff = subprocess.check_output(['git','diff','HEAD'],cwd=repo,text=True)
    (a.out/'source.patch').write_text(diff)
    meta['source_patch_sha256'] = sha(a.out/'source.patch')
    meta['source_files_sha256'] = {str(x.relative_to(repo)):sha(x) for pattern in
        ('src/*.cc', 'include/*.hh', 'data/sipm/*.txt') for x in repo.glob(pattern)}
    target = a.out/'invocation.meta.json'
    target.write_text(json.dumps(meta,indent=2)+'\n')
    print(meta['command'], flush=True)
    start = time.monotonic()
    with (a.out/'stdout.log').open('w') as log:
        code = subprocess.call(command,cwd=a.out,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy())
    meta.update(exit_code=code,wall_s=time.monotonic()-start,end_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    root = a.out/'photon_hits_run000.root'
    if root.exists():
        meta.update(root_bytes=root.stat().st_size,root_sha256=sha(root))
    target.write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta),flush=True)
    raise SystemExit(code)

if __name__ == '__main__':
    main()
