#!/usr/bin/env python
"""EXEC37: unchanged same-end V2 timing on the two existing END-only ROOTs."""
import csv
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys
import time

import numpy as np
import uproot
from exec36_bridge import load,SOURCE
from run_exec36 import width
from robust_widths import NAMES
from validity_exec35 import RULE,RULE_PATH
from run_exec35 import sha,now,clean,write_json,csv_rows,latest

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
BASE=Path('/home/rrios/exec37_20260913')
CONFIG_PATH=HERE.parent/'exec37_config.json'
CONFIG=json.loads(CONFIG_PATH.read_text())


def c1(normalized_width,valid):
    if not valid or not math.isfinite(normalized_width) or normalized_width<=0:return 'INVALID'
    return 'PASS' if normalized_width<=CONFIG['C1']['bound_ps'] else 'FAIL'


def crossings(event,gid,time_ns,ids,n_generated):
    """Preserve exact LeadingEdgeTime; explicit population size from manifest."""
    assert np.all((event>=0)&(event<n_generated))
    mask=np.isin(gid,ids);ev=event[mask];values=time_ns[mask]
    counts=np.bincount(ev,minlength=n_generated)
    values=values[np.argsort(ev,kind='stable')]
    offsets=np.r_[0,np.cumsum(counts)]
    root=load();result=np.full(n_generated,np.nan)
    for i in range(n_generated):
        v=root.std.vector('double')(values[offsets[i]:offsets[i+1]].tolist())
        result[i]=root.exec36.original(v)
    return result


def resolve(cell):
    spec=CONFIG['cells'][cell];manifest=Path(spec['manifest'])
    meta=json.loads(manifest.read_text())
    sim=meta['cells'][cell]
    assert sim['source_commit']==spec['commit'] and sim['exit_code']==0
    assert sim['N']==2000 and sim['N_TOP']==0 and sim['readout']=='End' and sim['x_mm']==0
    assert sim['seeds']==[26092601,8349041] and sim['material']=='EJ-204 / OPSC-101'
    raw=Path(sim['cwd'])/'photon_hits_run000.root'
    if not raw.is_file():raise FileNotFoundError(f'{cell}: no archived ROOT; STOP, no simulation: {raw}')
    hitmeta=raw.with_suffix('.meta.json');hm=json.loads(hitmeta.read_text())
    assert hm['sidecars']['.root']['path']==str(raw)
    assert sha(raw)==hm['sidecars']['.root']['sha256']
    with manifest.with_name('manifest.csv').open() as f:
        row=next(r for r in csv.DictReader(f) if r['cell_id']==cell)
    assert row['commit']==sim['source_commit'] and int(row['N'])==sim['N']
    assert row['start_utc']==sim['start_utc'] and row['end_utc']==sim['end_utc'] and int(row['exit_code'])==0
    assert sha(Path(sim['cwd'])/'run.mac')==sim['macro_sha256']
    return dict(simulation=sim,manifest_row=row,raw_ROOT=str(raw),root_sha256=sha(raw),
        input_sha256={str(p):sha(p) for p in [manifest,manifest.with_name('manifest.csv'),hitmeta,Path(sim['cwd'])/'run.mac']},
        raw_sidecar_metadata=hm,reflector=spec['reflector'])


def analyze(cell,source,reference,provenance):
    start=now();clock=time.monotonic();n=source['simulation']['N']
    out=BASE/'cells'/cell;out.mkdir(parents=True,exist_ok=False)
    with uproot.open(source['raw_ROOT']) as f:
        a=f['sipm_hits'].arrays(['event_id','global_id','face_type','time_ns'],library='np')
    ev=a['event_id'];gid=a['global_id'];ts=a['time_ns'];face=a['face_type']
    assert np.all((ev>=0)&(ev<n)) and np.all((gid>=0)&(gid<16))
    assert np.all(np.isin(face,[0,1])) and np.all((gid<8)==(face==0)) and np.all(np.isfinite(ts))
    left=np.bincount(ev[face==0],minlength=n);right=np.bincount(ev[face==1],minlength=n)
    npe=(left+right)/2.;mean=float(npe.mean());sem=float(npe.std(ddof=1)/math.sqrt(n))
    assert mean==float(source['manifest_row']['npe_end_mean'])
    assert np.isclose(sem,float(source['manifest_row']['npe_end_sem']),rtol=1e-13)
    rows=[]
    with uproot.recreate(out/'results.root') as root:
        root['event_yields']=dict(event_id=np.arange(n,dtype=np.int32),npe_left=left,npe_right=right,npe_end=npe)
        for end,groups in CONFIG['groups'].items():
            t1=crossings(ev,gid,ts,groups[0],n);t2=crossings(ev,gid,ts,groups[1],n)
            accepted=np.isfinite(t1)&np.isfinite(t2);delta=np.where(accepted,t1-t2,np.nan)
            r,boot,edges,gauss=width(delta,cell+'_'+end)
            ref=next(r for r in reference['metadata']['rows'] if r['variant']=='V2' and r['end']==end)
            assert ref['group_ids']==groups and ref['rise_ns']==.5 and ref['fall_ns']==5 and ref['threshold_PE']==4
            assert r['N_generated_partition']==2000
            normalized=r['sigma_core']/math.sqrt(2);bnormalized=boot[:,0]/math.sqrt(2)
            r.update(cell_id=cell,end=end,observable='same_end_T1_minus_T2',variant='V2',N_generated=n,
                group_ids=groups,CH1='T1 low cluster',CH5='T2 high cluster',rise_ns=.5,fall_ns=5.,threshold_PE=4.,
                N_TOP=0,x_mm=0,geometry='END-only',units='ps',npe_end_mean=mean,npe_end_sem=sem,
                source_commit=source['simulation']['source_commit'],simulation_seeds=source['simulation']['seeds'],
                reflector=source['reflector'],endtop_reference=ref,
                C1=c1(normalized,r['primary_validity']['status']=='VALID'),
                C1_bootstrap_fraction_le_bound=float(np.mean(bnormalized<=88.4)),
                normalized_core_ci_low=r['sigma_core_ci_low']/math.sqrt(2),normalized_core_ci_high=r['sigma_core_ci_high']/math.sqrt(2))
            for key in ['sigma_core','sigma_gauss']:
                ratio=ref[key]/r[key]
                r['C2_'+key+'_ratio']=ratio
                r['C2_'+key+'_ratio_se_zero_covariance']=ratio*math.hypot(ref[key+'_se']/ref[key],r[key+'_se']/r[key])
                r['C2_'+key+'_minus_1p719']=ratio-1.719
                r['C2_'+key+'_departure_percent']=100*(ratio/1.719-1)
            r['C2_interpretation']=CONFIG['C2']['uncertainty']+'; no PASS/FAIL'
            if r['C1']=='PASS':
                real=bnormalized<=88.4
                residual=np.sqrt(88.4**2-bnormalized[real]**2)
                r.update(C3_residual_ps=math.sqrt(88.4**2-normalized**2),
                    C3_conditional_bootstrap_se=float(np.std(residual,ddof=1)) if len(residual)>1 else None,
                    C3_ci_low=float(np.quantile(residual,.025)) if len(residual) else None,
                    C3_ci_high=float(np.quantile(residual,.975)) if len(residual) else None,
                    C3_real_domain_fraction=float(real.mean()))
            else:
                r.update(C3_residual_ps=None,C3_conditional_bootstrap_se=None,C3_ci_low=None,C3_ci_high=None,C3_real_domain_fraction=None)
            r['C3_interpretation']=CONFIG['C3']['interpretation']
            if r['C3_residual_ps'] is not None:
                for label,value in [('Lee_intrinsic',58.2),('Lee_detector',73.),('bank_low',85.),('bank_high',92.)]:
                    r['C3_minus_'+label+'_ps']=r['C3_residual_ps']-value
            rows.append(r)
            root[end+'_events']=dict(event_id=np.arange(n,dtype=np.int32),T1_ns=t1,T2_ns=t2,delta_ns=delta,accepted=accepted)
            root[end+'_bootstrap']=dict(**{name:boot[:,j] for j,name in enumerate(NAMES)},sigma_gauss=gauss)
            root[end+'_hist']=(np.histogram(delta[accepted]*1000,edges)[0].astype(float),edges)
        root['rows_json']=json.dumps(clean(rows),allow_nan=False)
    csv_rows(out/'results.csv',rows)
    assert sha(source['raw_ROOT'])==source['root_sha256']
    meta=dict(schema='EXEC37 v1',cell_id=cell,start_utc=start,end_utc=now(),wall_s=time.monotonic()-clock,
        source=source,reference=reference,configuration=CONFIG,validity_rule=RULE,provenance=provenance,
        rows=rows,npe_end_mean=mean,npe_end_sem=sem,ROOT_version=load().gROOT.GetVersion(),numpy_version=np.__version__,
        uproot_version=uproot.__version__,python_version=sys.version,
        files={name:sha(out/name) for name in ['results.root','results.csv']})
    write_json(out/'results.meta.json',meta)
    return dict(utc=now(),cell_id=cell,status='COMPLETE',output=str(out),
        sidecar_sha256={name:sha(out/name) for name in ['results.root','results.csv','results.meta.json']})


def main():
    # Resolve BOTH archived ROOTs before calculating either result. Missing -> STOP.
    sources={cell:resolve(cell) for cell in ['D0','D3']}
    assert sha(CONFIG_PATH)==sha(BASE/'preregistration.json')
    manifest=Path('/home/rrios/exec36_20260913/groups_manifest.jsonl')
    rec=latest(manifest,'COMPLETE')['EJ204_xp0'];refdir=Path(rec['output'])
    assert all(sha(refdir/name)==digest for name,digest in rec['sidecar_sha256'].items())
    reference=dict(record=rec,metadata=json.loads((refdir/'results.meta.json').read_text()),manifest_sha256=sha(manifest))
    scripts=[Path(__file__),CONFIG_PATH,RULE_PATH,HERE/'run_exec36.py',HERE/'exec36_bridge.py',SOURCE,
        HERE/'robust_widths.py',HERE/'validity_exec35.py',HERE/'run_exec35.py']
    provenance=dict(command=shlex.join([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]]),cwd=str(Path.cwd()),
        commit=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip(),
        source_sha256={str(p.resolve()):sha(p) for p in scripts},
        preregistration_sha256=sha(BASE/'preregistration.md'),
        sampling='2000 generated IDs from each original manifest; no fixed EXEC36 10000-event adapter used',
        normalization=CONFIG['normalization'],bootstrap=CONFIG['bootstrap'])
    journal=BASE/'manifest.jsonl';done=latest(journal,'COMPLETE') if journal.exists() else {}
    for cell,source in sources.items():
        if cell in done:
            rec=done[cell];out=Path(rec['output'])
            assert all(sha(out/name)==digest for name,digest in rec['sidecar_sha256'].items())
            assert json.loads((out/'results.meta.json').read_text())['provenance']['source_sha256']==provenance['source_sha256']
            continue
        rec=analyze(cell,source,reference,provenance)
        with journal.open('a') as f:f.write(json.dumps(rec)+'\n');f.flush()
        print(cell,'COMPLETE',flush=True)


if __name__=='__main__':main()
