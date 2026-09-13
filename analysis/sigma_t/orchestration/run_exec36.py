#!/usr/bin/env python
"""Existing-ROOT analysis only: single SiPM -> same-end groups -> sensitivity."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

import numpy as np
import uproot
from robust_widths import measure,NAMES
from validity_exec35 import verdict,RULE,RULE_PATH
from exec36_bridge import load,SOURCE,sensitivity_source
from run_exec35 import sha,now,write_json,csv_rows,latest,clean

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
BASE=Path('/home/rrios/exec36_20260913')
CONFIG_PATH=HERE.parent/'exec36_config.json'
CONFIG=json.loads(CONFIG_PATH.read_text())
SIM_MANIFEST=Path('/home/rrios/exec34r_20260912/manifest.jsonl')
PRIOR_MANIFEST=Path('/home/rrios/exec35_20260912/analysis_manifest.jsonl')


def sanity(row):
    if row['primary_validity']['status']!='VALID':return 'INDETERMINATE'
    width=row['sigma_core']
    if 150<=width<=350:return 'PASS'
    if width<100 or width>500:return 'FAIL_EXTREME'
    return 'INDETERMINATE'


def fit(values_ps,name):
    root=load(); v=root.std.vector('double')((np.asarray(values_ps)/1000).tolist())
    r=list(root.exec36.fit(v,name))
    return dict(sigma_gauss=r[0],sigma_gauss_covariance_se=r[1],gauss_chi2_ndf=r[2],
                gauss_n_eff=int(r[3]),gauss_fit_status=0 if r[4] else 1,gauss_used_fit=bool(r[4]),gauss_mean_ns=r[5])


def width(values_ns,label):
    values_ps=np.asarray(values_ns)*1000
    point=fit(values_ps[np.isfinite(values_ps)],label+'_point')
    def callback(v,i):
        r=fit(v,label+'_boot_'+str(i))
        return r['sigma_gauss'] if r['gauss_used_fit'] else np.nan
    row,boot,edges,extra=measure(values_ps,callback=callback)
    extra=np.asarray(extra);valid=extra[np.isfinite(extra)]
    row.update(point,provenance_valid=True,sigma_gauss_se=float(np.std(valid,ddof=1)) if len(valid)>1 else np.nan,
        sigma_gauss_bootstrap_fraction=len(valid)/len(extra),
        sigma_gauss_ci_low=float(np.quantile(valid,.025)) if len(valid) else np.nan,
        sigma_gauss_ci_high=float(np.quantile(valid,.975)) if len(valid) else np.nan)
    row['core_over_gauss']=row['sigma_core']/row['sigma_gauss'] if row['sigma_gauss']>0 else np.nan
    ratios=boot[:,0]/extra
    row['core_over_gauss_se']=float(np.std(ratios[np.isfinite(ratios)],ddof=1))
    row['ratio_uncertainty']='paired event-bootstrap core/gauss ratios on successful Gaussian replicas'
    row['primary_validity']=verdict(row);row['gaussian_validity']=verdict(row,gaussian=True)
    for name in ['sigma_core','sigma_gauss']:
        row[name+'_over_sqrt2']=row[name]/math.sqrt(2)
        row[name+'_over_sqrt2_se']=row[name+'_se']/math.sqrt(2)
    return row,boot,edges,extra


def verify_record(rec):
    out=Path(rec['output'])
    assert all(sha(out/name)==digest for name,digest in rec['sidecar_sha256'].items())
    return json.loads((out/'results.meta.json').read_text())


def campaign_sanity():
    single=latest(BASE/'single_manifest.jsonl','COMPLETE');assert len(single)==21
    checks=[r for rec in single.values() for r in verify_record(rec)['rows']]
    failures=[dict(cell_id=r['cell_id'],end=r['end'],width_ps=r['sigma_core']) for r in checks if r['sanity']=='FAIL_EXTREME']
    return dict(all_singles_finished_before_grouping=True,rows=checks,extreme_failures=failures,
        grouped_status='NOT COMPARABLE' if failures else 'NO EXTREME CAMPAIGN FAILURE; retain cell/end sanity status and hardware caveats')


def analyze(job):
    cell_id,phase,simulation,previous,provenance,guard=job
    start=now();clock=time.monotonic();root=load()
    # Previous EXEC35 analysis is a frozen reference, never an input to a new selection.
    old_dir=Path(previous['output'])
    assert all(sha(old_dir/name)==digest for name,digest in previous['sidecar_sha256'].items())
    old=json.loads((old_dir/'widths.meta.json').read_text())
    sim=old['simulation']
    raw=Path(simulation['output'])/'photon_hits_run000.root'
    assert str(raw)==sim['root_path'] and sim['root_sha256']==simulation['root_sha256']
    assert sim['N_generated']==sim['events_run']==10000 and sim['workers']==4 and sim['eventModulo']==1
    assert sim['seeds']==[26092601,8349041] and sim['sptr_ns']==0 and sim['diagnostics'] is False
    assert sha(sim['PDE_path'])==sim['PDE_sha256']
    out=BASE/'cells'/cell_id/phase;out.mkdir(parents=True,exist_ok=False)
    resultroot=out/'results.root'
    if phase=='single':
        assert sha(raw)==sim['root_sha256']
        data=root.exec36.Data(str(raw),str(resultroot))
        input_evidence=dict(raw_ROOT=str(raw),sha256=sha(raw))
        update=uproot.update
    else:
        rec=latest(BASE/'single_manifest.jsonl','COMPLETE')[cell_id]
        cached=verify_record(rec)
        assert cached['simulation']['root_sha256']==sim['root_sha256']
        cache=Path(rec['output'])/'results.root'
        data=root.exec36.Data(str(cache))
        input_evidence=dict(exact_END_cache=str(cache),sha256=sha(cache),source_root=str(raw),source_sha256=sim['root_sha256'])
        update=uproot.recreate
    assert data.left==sim['left_total'] and data.right==sim['right_total']
    rows=[];checks={};events=np.arange(10000,dtype=np.int32)
    configs=[CONFIG['baseline']] if phase!='sensitivity' else [dict(rise_ns=p[0],fall_ns=p[1],threshold_PE=t)
        for p in CONFIG['sensitivity_pulses_ns'] for t in CONFIG['sensitivity_thresholds_PE']]
    with update(resultroot) as f:
        def save(values,identity,t1=None,t2=None):
            key=identity['distribution_id']
            row,boot,edges,gaussian=width(values,key)
            row.update(identity,cell_id=cell_id,phase=phase,units='ps',N_generated=10000)
            if identity['observable']=='single_SiPM_minus_gun':
                row['sanity']=sanity(row)
                row['normalization_note']='Primary raw single-SiPM width relative to t_gun=0; sqrt2 fields are not a physical single-SiPM observable and must not be reported as such'
            else:
                row['normalization_note']=CONFIG['normalization']
                baseline_single=next(r for r in guard['rows'] if r['cell_id']==cell_id and r['end']==row['end'])
                row['baseline_single_sanity']=baseline_single['sanity']
                row['comparability']='NOT COMPARABLE' if guard['extreme_failures'] else baseline_single['sanity']+' SANITY ONLY; hardware/discriminator calibration pending'
            rows.append(row)
            values_ps=values*1000
            payload=dict(event_id=events,time_ps=values_ps,accepted=np.isfinite(values))
            if t1 is not None:payload.update(T1_ps=t1*1000,T2_ps=t2*1000)
            f[key+'_events']=payload
            f[key+'_bootstrap']=dict(**{name:boot[:,j] for j,name in enumerate(NAMES)},sigma_gauss=gaussian)
            f[key+'_hist']=(np.histogram(values_ps[np.isfinite(values_ps)],edges)[0].astype(float),edges)
        for cfg in configs:
            tag=f"r{cfg['rise_ns']:g}_f{cfg['fall_ns']:g}_th{cfg['threshold_PE']:g}".replace('.','p')
            varied=phase=='sensitivity'
            times=lambda ids:np.asarray(list(data.times(root.std.vector('int')(ids),cfg['rise_ns'],cfg['fall_ns'],cfg['threshold_PE'],varied)))
            if phase in ['single','sensitivity']:
                for end,ids in CONFIG['single_ids'].items():
                    save(times(ids),dict(observable='single_SiPM_minus_gun',variant='single',end=end,
                        group_ids=[ids],excluded_ids=[i for i in range(0 if end=='left' else 8,8 if end=='left' else 16) if i not in ids],
                        distribution_id=f'single_{end}_{tag}',**cfg))
            if phase!='single':
                timestamps={}
                for variant,ends in CONFIG['groups'].items():
                    for end,groups in ends.items():
                        t1,t2=times(groups[0]),times(groups[1]);timestamps[(variant,end)]=(t1,t2)
                        accepted=np.isfinite(t1)&np.isfinite(t2)
                        delta=np.where(accepted,t1-t2,np.nan)
                        save(delta,dict(observable='same_end_T1_minus_T2',variant=variant,end=end,group_ids=groups,
                            excluded_ids=[i for i in range(0 if end=='left' else 8,8 if end=='left' else 16) if i not in groups[0]+groups[1]],
                            CH1='T1: low cluster',CH5='T2: high cluster',distribution_id=f'{variant}_{end}_{tag}',**cfg),t1,t2)
                if cfg==CONFIG['baseline']:
                    # Gold control: reconstruct prior opposite-end reduction exactly,
                    # without refitting or substituting its historical result.
                    reconstructed=np.fmin(*timestamps[('V2','left')])-np.fmin(*timestamps[('V2','right')])
                    with uproot.open(old_dir/'widths.root') as previous_root:
                        previous_delta=previous_root['end_events']['delta_ps'].array(library='np')
                    assert np.array_equal(reconstructed*1000,previous_delta,equal_nan=True), 'Changed END timestamps/reduction'
                    checks['exact_prior_opposite_end_event_array_match']=True
            print(cell_id,phase,tag,'complete',flush=True)
        f['rows_json']=json.dumps(clean(rows),allow_nan=False)
    csv_rows(out/'results.csv',rows)
    metadata=dict(schema='EXEC36 v1',cell_id=cell_id,phase=phase,simulation=sim,
        previous_EXEC35_record=previous,previous_END_rows=[r for r in old['rows'] if r['arm']=='END'],
        inherited_EXEC35_provenance={k:old[k] for k in ['preregistration','provenance','previous_sidecar_sha256']},
        input_evidence=input_evidence,configuration=CONFIG,validity_rule=RULE,
        config_sha256=sha(CONFIG_PATH),validity_sha256=sha(RULE_PATH),provenance=provenance,
        guard=guard,checks=checks,rows=rows,start_utc=start,end_utc=now(),wall_s=time.monotonic()-clock,
        ROOT_version=root.gROOT.GetVersion(),numpy_version=np.__version__,uproot_version=uproot.__version__,
        python_version=sys.version,bootstrap=CONFIG['bootstrap'],units='ps',
        label='Intrinsic detected-photon discriminator sensitivity; no SPTR, jitter or TDC quantization',
        files={name:sha(out/name) for name in ['results.root','results.csv']})
    write_json(out/'results.meta.json',metadata)
    del data
    return dict(utc=now(),cell_id=cell_id,status='COMPLETE',phase=phase,output=str(out),
        sidecar_sha256={name:sha(out/name) for name in ['results.root','results.csv','results.meta.json']})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phase',choices=['single','groups','sensitivity'],required=True)
    p.add_argument('--jobs',type=int,default=4);args=p.parse_args()
    simulations=latest(SIM_MANIFEST,'SIMULATION_COMPLETE');previous=latest(PRIOR_MANIFEST,'COMPLETE')
    assert len(simulations)==len(previous)==21 and simulations.keys()==previous.keys()
    assert sha(CONFIG_PATH)==sha(BASE/'preregistration.json')
    guard={} if args.phase=='single' else campaign_sanity()
    sources=[Path(__file__),HERE/'exec36_bridge.py',SOURCE,CONFIG_PATH,RULE_PATH,HERE/'validity_exec35.py',HERE/'robust_widths.py',HERE/'run_exec35.py']
    provenance=dict(command=shlex.join([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]]),
        cwd=str(Path.cwd()),commit=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip(),
        source_sha256={str(p.resolve()):sha(p) for p in sources},
        preregistration_sha256=sha(BASE/'preregistration.md'),
        manifest_sha256={str(p):sha(p) for p in [SIM_MANIFEST,PRIOR_MANIFEST]},
        sensitivity_body_identity='Verbatim source slice SprPeakTime through LeadingEdgeTime; constants supplied by named namespace, never edit upstream',
        bootstrap_pairing='Same event-index resamples in all distributions; correlations retained in stored replica vectors')
    journal=BASE/(args.phase+'_manifest.jsonl')
    done=latest(journal,'COMPLETE') if journal.exists() else {}
    ids=[CONFIG['sensitivity_cell']] if args.phase=='sensitivity' else sorted(simulations)
    jobs=[]
    for cell_id in ids:
        if cell_id in done:
            stored=verify_record(done[cell_id])
            assert stored['provenance']['source_sha256']==provenance['source_sha256']
        else:jobs.append((cell_id,args.phase,simulations[cell_id],previous[cell_id],provenance,guard))
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures=[pool.submit(analyze,job) for job in jobs]
        for future in as_completed(futures):
            rec=future.result()
            with journal.open('a') as f:f.write(json.dumps(rec)+'\n');f.flush();os.fsync(f.fileno())
            print(rec['cell_id'],args.phase,'COMPLETE',flush=True)


if __name__=='__main__':main()
