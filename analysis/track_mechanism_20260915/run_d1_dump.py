#!/usr/bin/env python3
"""Compila un lector del MPT usando objetos de producción; inicialización sin eventos."""
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

REPO = Path('/home/rrios/ej200')
BUILD = Path('/home/rrios/exec46_20260915/build_baseline')
OUTPUT = BUILD.parent/'decay_diagnostics'
CELL_DIRECTORY = BUILD.parent/'full_grid/cells/EJ200_xp0'
SOURCE = Path(__file__).with_name('dump_effective_mpt.cc').resolve()
INITIALIZATION_LIMIT_S = 120

def digest(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024**2),b''): result.update(chunk)
    return result.hexdigest()

def main():
    OUTPUT.mkdir(exist_ok=True)
    done=json.loads((CELL_DIRECTORY/'.DONE').read_text())
    if digest(done['binary']) != done['binary_sha256']:
        raise RuntimeError('Production binary hash mismatch')
    flags=(BUILD/'CMakeFiles/ej200_bar_sim.dir/flags.make').read_text()
    options=[]
    for line in flags.splitlines():
        if line.startswith(('CXX_DEFINES =','CXX_INCLUDES =','CXX_FLAGS =')):
            options+=shlex.split(line.split('=',1)[1])
    compile_command=['/usr/bin/c++',*options,'-c',str(SOURCE),'-o',str(OUTPUT/'dump_effective_mpt.o')]
    subprocess.run(compile_command,check=True)
    command=shlex.split((BUILD/'CMakeFiles/ej200_bar_sim.dir/link.txt').read_text())
    command=[str(OUTPUT/'dump_effective_mpt.o') if x.endswith('/main.cc.o') else x for x in command]
    command[command.index('-o')+1]=str(OUTPUT/'dump_effective_mpt')
    subprocess.run(command,cwd=BUILD,check=True)
    macro=Path(done['cwd'])/'run.mac'
    run=[str(OUTPUT/'dump_effective_mpt'),str(macro),str(OUTPUT/'effective_mpt.json')]
    env=os.environ.copy()
    env.update(json.loads((BUILD.parent/'full_grid/campaign.json').read_text())['runtime_environment'])
    with (OUTPUT/'effective_mpt.log').open('w') as f:
        subprocess.run(run,cwd=done['cwd'],env=env,stdout=f,stderr=subprocess.STDOUT,
                       check=True,timeout=INITIALIZATION_LIMIT_S)
    inputs=[SOURCE,macro,*list((BUILD/'sslg4/data/oscnt/opsc-100').glob('*.txt')),
            BUILD/'sslg4/macros/oscnt/opsc-100.mac']
    inputs += [BUILD/x for x in command if x.endswith(('.o','.a'))]
    provenance=dict(compile_argv=compile_command,link_argv=command,link_cwd=str(BUILD),
                    run_argv=run,run_cwd=done['cwd'],simulation_commit=done['simulation_commit'],
                    binary_sha256_expected=done['binary_sha256'],binary_sha256_now=digest(done['binary']),
                    files_sha256={str(p):digest(p) for p in inputs},events_generated=0)
    (OUTPUT/'d1_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (OUTPUT/'opsc-100.mac.txt').write_bytes((BUILD/'sslg4/macros/oscnt/opsc-100.mac').read_bytes())
    print('D1 completed: full effective MPT saved; zero generated events')

if __name__=='__main__': main()
