#!/usr/bin/env python
"""One reproducible cell: simulation -> adapter/primitives -> sidecars -> G-P."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys
from run_simulation import now, sha

HERE = Path(__file__).resolve().parent


def main(args):
    dest = Path(args.output).resolve()
    simulation_meta = dest/'simulation.meta.json'
    analysis = Path(args.analysis_output).resolve() if args.analysis_output else dest/'analysis'
    journal = []
    def stage(name, command):
        start=now(); print(shlex.join(command),flush=True)
        rc=subprocess.call(command)
        journal.append(dict(stage=name,command=shlex.join(command),start_utc=start,end_utc=now(),exit_code=rc))
        if dest.exists():
            # Unique invocation journal; no overwrite of previous resume attempts.
            log=dest/('cell_invocation_'+start.replace(':','').replace('+','_')+'.json')
            log.write_text(json.dumps(journal,indent=2)+'\n')
        if rc: raise SystemExit(rc)
    if simulation_meta.exists():
        if not args.resume: raise RuntimeError('Existing simulation requires --resume; nothing overwritten')
        sim=json.loads(simulation_meta.read_text())
        assert sim['exit_code']==0 and sha(sim['root_path'])==sim['root_sha256']
        assert sim['N_generated']==10000 and sim['material']==args.material and sim['opsc_code']==args.opsc and sim['x_mm']==args.x
        assert sim['workers']==args.workers and sim['eventModulo']==1
        assert sim['seeds']==[26092601,8349041] and sim['N_TOP']==70 and sim['diagnostics'] is False
    else:
        stage('simulation',[sys.executable,str(HERE/'run_simulation.py'),'--binary',args.binary,
            '--output',str(dest),'--workers',str(args.workers),'--material',args.material,'--opsc',args.opsc,'--x',str(args.x)])
    if (analysis/'analysis.meta.json').exists():
        if not args.resume: raise RuntimeError('Existing analysis requires --resume')
        data=json.loads((analysis/'analysis.meta.json').read_text())
        assert all(sha(path)==digest for path,digest in data['script_sha256'].items()), 'Changed analysis scripts; choose a new --analysis-output'
        assert data['simulation']['root_sha256']==json.loads(simulation_meta.read_text())['root_sha256']
    else:
        stage('analysis',[sys.executable,str(HERE/'analyze.py'),'--simulation',str(dest),'--output',str(analysis)])
    stage('gate',[sys.executable,str(HERE/'gate.py'),'--analysis',str(analysis)])
    return 0


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',required=True);p.add_argument('--output',required=True)
    p.add_argument('--workers',type=int,default=24)
    p.add_argument('--material',default='EJ-204');p.add_argument('--opsc',default='OPSC-101')
    p.add_argument('--x',type=int,default=0);p.add_argument('--resume',action='store_true')
    p.add_argument('--analysis-output')
    raise SystemExit(main(p.parse_args()))
