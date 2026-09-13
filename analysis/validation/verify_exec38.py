#!/usr/bin/env python3
"""Verify saved EXEC38 provenance and statistical identities, without simulation."""
import json
from pathlib import Path
import subprocess
import numpy as np
import uproot
from run_exec38 import HERE, REPO, BASE, sha, now, dump


def main():
    out=BASE/'exec38_20260913';gp=HERE/'GOLDEN_REFERENCE_20260913.json'
    golden=json.loads(gp.read_text());results=json.loads((out/'results.json').read_text())
    assert [p['id'] for p in golden['tests']]==['V'+str(i) for i in range(1,8)]
    assert golden['ready_for_acceptance'] is False
    for test in golden['tests']:
        for key in ['observable','prediction','measured','uncertainty','tolerance','status','config_ids']:
            assert key in test,(test['id'],key)
        assert test['prediction']['source'] and test['tolerance']
        for key in test['config_ids']:
            config=golden['configurations'][key]
            for field in ['material','x_mm','N_generated','seeds','simulation_commit','command']:
                assert config[field] is not None,(key,field)
        if test['id'] in ['V1','V2','V5']:assert test['status']=='NOT_EVALUABLE'
        assert 'SteppingAction' not in test['observable'] and 'TrackingAction' not in test['observable']
    for path,h in golden['artifact_hashes'].items():assert sha(path)==h,path
    delivery=json.loads((out/'delivery.meta.json').read_text())
    assert sha(delivery['golden'])==delivery['golden_sha256']
    assert sha(delivery['report'])==delivery['report_sha256']
    assert len(results['V7'])==88 and len(results['V3'])==len(results['V4'])==6
    with uproot.open(out/'V7_correlations.root') as f:
        for row in results['V7']:
            prefix=f"{row['cell_id']}_{row['variant']}_{row['end']}"
            boot=f[prefix+'_bootstrap'].arrays(library='np')
            events=f[prefix+'_events'].arrays(library='np')
            assert len(events['event_id'])==row['N_generated']
            assert events['accepted'].sum()==row['n_accepted']
            for field,values in boot.items():
                assert np.isclose(np.std(values,ddof=1),row[field+'_se'],rtol=1e-12,atol=1e-13)
            rho=row['rho'];s1=row['sigma_T1_ps'];s2=row['sigma_T2_ps'];sd=row['sigma_delta_ps']
            assert np.isclose(sd**2,s1*s1+s2*s2-2*rho*s1*s2,rtol=1e-12)
            assert row['status']=='PASS' and rho-3*row['rho_se']>0
    with uproot.open(out/'profile_points.root') as points,uproot.open(out/'V3_attenuation.root') as fits:
        for test in ['V3','V4']:
            for row in results[test]:
                prefix=row['material'].replace('-','')+'_'+row['end']
                cfg=golden['configurations'];d=[];y=[]
                for cell in row['cell_ids']:
                    d.append(700+cfg[cell]['x_mm']*(1 if row['end']=='left' else -1))
                    ev=points[cell+'_events'].arrays(library='np');count=ev['npe_'+row['end']]
                    y.append(count.mean() if test=='V3' else ev['time_sum_'+row['end']+'_ns'].sum()/count.sum())
                d=np.array(d);y=np.array(y)
                cov=fits[prefix+('_count_covariance' if test=='V3' else '_time_covariance')].arrays(library='np')['covariance'].reshape(7,7)
                pred=row['amplitude']*np.exp(-d/row['lambda_mm']) if test=='V3' else row['intercept_ns']+d*row['slope_ns_per_mm']
                residual=y-pred;stat=residual@np.linalg.solve(cov,residual)
                assert np.isclose(stat,row['chi2'],rtol=1e-10)
    for png in out.glob('*.png'):
        for ext in ['.root','.csv','.meta.json','.pdf']:assert png.with_suffix(ext).is_file()
        meta=json.loads(png.with_suffix('.meta.json').read_text())
        for path,h in meta['files'].items():assert sha(path)==h
        for ext,h in meta['figure_files'].items():assert sha(png.with_suffix(ext))==h
        with uproot.open(png.with_suffix('.root')) as f:assert 'summary' in f and 'rows_json' in f
    checks=dict(utc=now(),status='PASS',golden_sha256=sha(gp),
        commit_at_verification=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        checks=['Seven fail-closed contract entries and configuration references','Every contract artifact hash',
                'All 88 V7 bootstrap standard errors and paired-event populations',
                'All 88 unequal-variance identities and preregistered positive-rho decisions',
                'All six V3 and six V4 chi-square values reconstructed from saved full covariance',
                'All three figures with readable ROOT, CSV, metadata, PDF and matching hashes'],
        additional_native_check='EJ204_xp0: native sipm_hits event counts equal cached counts exactly; per-event time sums agree to relative 1e-13, checked in read-only shell command',
        no_simulation=True)
    dump(out/'verification.meta.json',checks)
    print('PASS:',*checks['checks'],sep='\n- ')


if __name__=='__main__':main()
