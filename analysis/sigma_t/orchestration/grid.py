#!/usr/bin/env python
"""EXEC34B only: bounded independent cells, append-only manifest and resume.

This program is prepared by EXEC34A but not launched with --execute there.
Without --execute it validates and prints the plan, creating no cell outputs.
"""
import argparse
import fcntl
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from run_simulation import now, sha, memory


def append_manifest(path, record):
    with Path(path).open('a') as f:
        f.write(json.dumps(dict(utc=now(),**record))+'\n')
        f.flush(); os.fsync(f.fileno())


def read_manifest(path):
    latest = {}
    if Path(path).exists():
        for line in Path(path).read_text().splitlines():
            if line.strip():
                row=json.loads(line)
                if 'cell_id' in row: latest[row['cell_id']]=row
    return latest


def verified_pass(cell_output):
    output=Path(cell_output)
    p=output/'analysis/gate.json'
    if not p.exists(): return False
    gate=json.loads(p.read_text())
    return gate['gate_status']=='PASS' and all(sha(path)==value for path,value in gate['sidecar_sha256'].items())


def plan(config, resume):
    assert config['canonical_train_parity']=='even' and config['workers_per_process']==1
    assert json.loads(Path(config['pilot_gate']).read_text())['gate_status']=='PASS'
    assert sha(config['pilot_gate'])==config['pilot_gate_sha256']
    for path,digest in config['script_sha256'].items():
        assert sha(path)==digest, 'Campaign script changed: '+path
    assert sha(config['binary'])==config['binary_sha256']
    previous=read_manifest(config['manifest'])
    todo=[]; skipped=[]
    for cell in config['cells']:
        state=previous.get(cell['cell_id'],{})
        if state.get('status')=='PASS':
            assert resume, 'Existing completed campaign requires --resume'
            assert verified_pass(cell['output']), 'Completed cell sidecars no longer verify'
            skipped.append(cell['cell_id'])
        else:
            if state.get('status') not in (None,'PENDING'):
                assert resume, 'Existing started/failed campaign requires --resume'
            todo.append(cell)
    available=memory()['MemAvailable']
    live_cap=int((available*(1-config['memory_reserve_fraction']))//config['conservative_process_budget_bytes'])
    concurrency=min(config['planned_concurrency'],live_cap,os.cpu_count() or 1,max(1,len(todo)))
    assert concurrency>=1, 'Insufficient free RAM for even one conservative process budget'
    return todo,skipped,concurrency


def command(config,cell,resume):
    result=[config['python'],config['run_cell'],'--binary',config['binary'],'--output',cell['output'],
            '--workers','1','--material',cell['material'],'--opsc',cell['opsc'],'--x',str(cell['x_mm'])]
    if resume: result.append('--resume')
    return result


def main(args):
    config=json.loads(Path(args.config).read_text())
    todo,skipped,concurrency=plan(config,args.resume)
    print(json.dumps(dict(execute=args.execute,concurrency=concurrency,skipped=skipped,
                          pending=[c['cell_id'] for c in todo],commands=[shlex.join(command(config,c,args.resume)) for c in todo]),indent=2),flush=True)
    if not args.execute: return 0
    output=Path(config['output_directory']);output.mkdir(parents=True,exist_ok=True)
    lock=(output/'grid.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
    # Hold the process-wide lock until exit; a second runner cannot duplicate work.
    (output/'logs').mkdir(exist_ok=True)
    def run(cell):
        env=os.environ.copy()
        env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        logfile=output/'logs'/(cell['cell_id']+'_'+now().replace(':','')+'.log')
        with logfile.open('w') as log:
            rc=subprocess.call(command(config,cell,args.resume),stdout=log,stderr=subprocess.STDOUT,env=env)
        return dict(cell_id=cell['cell_id'],status='PASS' if rc==0 and verified_pass(cell['output']) else 'FAIL',
                    exit_code=rc,output=cell['output'],log=str(logfile),
                    gate=str(Path(cell['output'])/'analysis/gate.json'))
    failure=False; iterator=iter(todo)
    with ThreadPoolExecutor(max_workers=concurrency) as pool:
        active={}
        def launch():
            cell=next(iterator,None)
            if cell is None: return False
            append_manifest(config['manifest'],dict(cell_id=cell['cell_id'],status='RUNNING',
                command=shlex.join(command(config,cell,args.resume)),output=cell['output']))
            active[pool.submit(run,cell)]=cell
            return True
        for _ in range(concurrency): launch()
        while active:
            finished,_=wait(active,return_when=FIRST_COMPLETED)
            for future in finished:
                cell=active.pop(future)
                try: result=future.result()
                except Exception as error: result=dict(cell_id=cell['cell_id'],status='FAIL',exit_code=-1,error=str(error))
                append_manifest(config['manifest'],result)
                if result['status']!='PASS': failure=True
            if not failure:
                for _ in finished: launch()
            # After the first failure: no new cells. Existing processes finish;
            # their evidence is retained, with no kill/delete or physics retry.
    return 34 if failure else 0


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',required=True);p.add_argument('--resume',action='store_true')
    p.add_argument('--execute',action='store_true')
    raise SystemExit(main(p.parse_args()))
