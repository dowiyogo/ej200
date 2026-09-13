#!/usr/bin/env python3
"""Read-only transport regression and preregistered one-cell physical checks.

Never runs Geant4. Outputs paired CSV/ROOT/meta sidecars and a machine verdict.
Requires numpy, scipy, uproot; analysis resamples events, not independent photons.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import uproot
from scipy.stats import chi2

N = 2000
SEED = 40091301
REPLICATES = 500

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def arrays(tree):
    return tree.arrays(library='np')

def hit_counts(tree):
    counts = np.zeros((N, 3), dtype=np.int64)
    for a in tree.iterate(['event_id','face_type'], step_size='64 MB', library='np'):
        assert np.all((a['event_id'] >= 0) & (a['event_id'] < N))
        assert np.all((a['face_type'] >= 0) & (a['face_type'] < 3))
        counts += np.bincount(3*a['event_id']+a['face_type'],minlength=3*N).reshape(N,3)
    return counts

def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')

def sidecar(out, stem, data, meta):
    columns = {key:np.asarray(value) for key,value in data.items()}
    with (out/f'{stem}.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(columns)
        w.writerows(zip(*columns.values()))
    with uproot.recreate(out/f'{stem}.root') as f:
        f['data'] = columns
    dump(out/f'{stem}.meta.json', dict(meta, columns=list(columns), rows=len(next(iter(columns.values()))),
        csv_sha256=sha(out/f'{stem}.csv'), root_sha256=sha(out/f'{stem}.root')))

def ratio(numerator, denominator, weights):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    estimate = float(numerator.sum()/denominator.sum())
    boots = (weights@numerator)/(weights@denominator)
    return dict(value=estimate, se=float(boots.std(ddof=1)), numerator=float(numerator.sum()),
                denominator=float(denominator.sum()))

def compare(value, expected):
    return 'PASS' if abs(value['value']-expected) <= 3*value['se'] else 'FAIL'

def physical_regression(measured, archived, ledger):
    diff=measured-archived
    npe=measured[:,:2].sum(axis=1)/2
    return dict(status='PASS' if np.array_equal(diff,np.zeros_like(diff)) and np.array_equal(ledger,measured)
              and float(npe.mean()) == 397.1525 else 'FAIL_ABORT',
        npe_end_mean=float(npe.mean()), npe_end_sem=float(npe.std(ddof=1)/np.sqrt(N)),
        reference=397.1525, difference=float(npe.mean()-397.1525),
        events_with_different_counts=int(np.any(diff!=0,axis=1).sum()),
        event_ledger_matches_hits=bool(np.array_equal(ledger,measured)), totals=measured.sum(axis=0).tolist())

def run(a):
    a.out.mkdir(parents=True, exist_ok=True)
    invocation = json.loads((a.cell/'invocation.meta.json').read_text())
    assert invocation['exit_code'] == 0, 'Simulation did not complete'
    root_path = a.cell/'photon_hits_run000.root'
    root = uproot.open(root_path)
    meta = dict(created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), N=N,
        seeds=invocation['seeds'], simulation_command=invocation['command'],
        source_root=str(root_path), source_root_sha256=sha(root_path),
        analysis_command=' '.join(sys.argv), analysis_sha256=sha(Path(__file__)),
        bootstrap_seed=SEED, bootstrap_replicates=REPLICATES, bootstrap_unit='generated event',
        preregistration_sha256=invocation['preregistration_sha256'])
    events = arrays(root['event_observables'])
    order = np.argsort(events['event_id'])
    events = {k:v[order] for k,v in events.items()}
    assert np.array_equal(events['event_id'],np.arange(N)), 'Missing or duplicate generated event'
    measured = hit_counts(root['sipm_hits'])
    archived = hit_counts(uproot.open(a.d1/'photon_hits_run000.root')['sipm_hits'])
    d1_csv = np.genfromtxt(a.d1/'event_yields.csv', delimiter=',', names=True)
    d1_csv = d1_csv[np.argsort(d1_csv['event_id'])]
    csv_counts = np.column_stack([d1_csv[x] for x in ('END_L','END_R','TOP')]).astype(np.int64)
    assert np.array_equal(archived,csv_counts), 'D1 ROOT and archived CSV disagree'
    ledger = np.column_stack([events[x] for x in ('detected_left','detected_right','detected_top')])
    diff = measured-archived
    g1 = physical_regression(measured, archived, ledger)
    sidecar(a.out,'event_regression',dict(event_id=np.arange(N),
        **{f'{name}_{label}':matrix[:,i] for name,matrix in [('new',measured),('D1',archived),('delta',diff)]
           for i,label in enumerate(('left','right','top'))}),meta)
    d1_meta=json.loads((a.d1/'cell_summary.json').read_text())
    d1_bytes=(a.d1/'photon_hits_run000.root').stat().st_size
    g3=dict(status='MEASURED',root_bytes=root_path.stat().st_size,D1_root_bytes=d1_bytes,
        root_ratio=root_path.stat().st_size/d1_bytes,wall_s=invocation['wall_s'],
        D1_wall_s=d1_meta['wall_s'],wall_ratio=invocation['wall_s']/d1_meta['wall_s'],
        D1_workers=1,D1_diagnostics=True,new_workers=4,new_diagnostics=False,
        matching_S2_wall_s=358.38,matching_S2_wall_ratio=invocation['wall_s']/358.38)
    g3['proposal_required']=g3['root_ratio']>2
    result=dict(G_I_1=g1,G_I_3=g3,metadata=meta,ready_for_acceptance=False)
    sizes=[(t,k,int(branch.compressed_bytes),int(branch.uncompressed_bytes))
           for t in root.keys(cycle=False) for k,branch in root[t].items()]
    sidecar(a.out,'root_branch_sizes',dict(tree=[r[0] for r in sizes],branch=[r[1] for r in sizes],
        compressed_bytes=[r[2] for r in sizes],uncompressed_bytes=[r[3] for r in sizes]),meta)
    g3['tree_compressed_bytes']={t:int(sum(r[2] for r in sizes if r[0]==t)) for t in root.keys(cycle=False)}
    print('G-I.1',json.dumps(g1),flush=True)
    if g1['status'] != 'PASS':
        result.update(G_I_2={'status':'NOT_EVALUATED_AFTER_ABORT'},
                      V1={'status':'NOT_EVALUATED_AFTER_ABORT'},V2={'status':'NOT_EVALUATED_AFTER_ABORT'},
                      V5={'status':'NOT_EVALUATED_AFTER_ABORT'})
        dump(a.out/'results.json',result)
        return 2

    rng=np.random.default_rng(SEED)
    weights=np.stack([np.bincount(rng.integers(0,N,N),minlength=N) for _ in range(REPLICATES)]).astype(float)
    sidecar(a.out,'event_observables',events,meta)
    v1=ratio(events['produced_scint'],10400*events['edep_total_MeV'],weights)
    v1.update(status=compare(v1,1.),expected=1.,yield_photons_per_MeV=10400,
        denominator='total deposited energy of all particles in BarLV',
        total_edep_MeV=float(events['edep_total_MeV'].sum()),
        optical_edep_MeV=float(events['edep_optical_MeV'].sum()),
        nonionizing_edep_MeV=float(events['edep_nonionizing_MeV'].sum()),
        produced_optical=int(events['produced_optical'].sum()),produced_scint=int(events['produced_scint'].sum()))
    v1['secondary_nonoptical_denominator']=ratio(events['produced_scint'],
        10400*(events['edep_total_MeV']-events['edep_optical_MeV']),weights)

    hist=np.zeros((N,20),dtype=np.int64)
    populations={key:np.zeros((N,2),dtype=np.int64) for key in ('primary_scint_exiting','scint_bar_air','all_exiting')}
    checks=dict(rows=0,invalid_surface=0,invalid_orientation=0,invalid_cosine=0,invalid_norm=0,
                nonfinite=0,unexpected_outcome=0)
    pair_counts={}
    keys=[]
    min_cos=1.;max_cos=-1.
    for b in root['first_bar_encounters'].iterate(step_size='64 MB',library='np'):
        ev=b['event_id'].astype(np.int64)
        assert np.all((ev>=0)&(ev<N))
        keys.append((ev.astype(np.uint64)<<32)|b['track_id'].astype(np.uint64))
        mu=b['cos_incidence']
        checks['rows']+=len(ev)
        checks['invalid_surface']+=int((b['normal_valid']!=1).sum())
        checks['invalid_orientation']+=int((b['normal_orientation_valid']!=1).sum())
        checks['invalid_cosine']+=int(((mu < -1e-9)|(mu>1+1e-9)).sum())
        checks['invalid_norm']+=int((np.abs(b['normal_norm']-1)>1e-9).sum())
        checks['nonfinite']+=int((~np.isfinite(mu)).sum())
        checks['unexpected_outcome']+=int((b['outcome']==0).sum())
        min_cos=min(min_cos,float(mu.min()));max_cos=max(max_cos,float(mu.max()))
        exiting=b['exiting_bar']==1
        scint=b['source']==1
        # Escape means leaving the bar into air/world, not transfer into a SiPM.
        air=np.array([str(s).startswith('AirGap') or str(s)=='WorldPV' for s in b['post_volume']])
        escaped=(b['outcome']==2)&air
        masks={'primary_scint_exiting':scint&exiting,'scint_bar_air':scint&exiting&air,'all_exiting':exiting}
        for name,mask in masks.items():
            populations[name][:,0]+=np.bincount(ev[mask],minlength=N)
            populations[name][:,1]+=np.bincount(ev[mask&escaped],minlength=N)
        mask=masks['primary_scint_exiting']&np.isfinite(mu)&(mu>=0)&(mu<=1)
        bins=np.minimum((mu[mask]*20).astype(int),19)
        hist+=np.bincount(20*ev[mask]+bins,minlength=N*20).reshape(N,20)
        # Compact physical-volume pair, source, and outcome audit, not a bounce census.
        labels=np.array([f'{p}|{q}|{s}|{r}' for p,q,s,r in zip(b['pre_volume'],b['post_volume'],b['source'],b['boundary_status'])])
        unique,counts=np.unique(labels,return_counts=True)
        for label,count in zip(unique,counts):pair_counts[label]=pair_counts.get(label,0)+int(count)
    joined=np.concatenate(keys)
    del keys
    joined.sort()
    checks['duplicate_photon_keys']=int(np.count_nonzero(joined[1:]==joined[:-1]))
    del joined
    checks.update(cos_min=min_cos,cos_max=max_cos)
    v2={key:ratio(values[:,1],values[:,0],weights) for key,values in populations.items()}
    primary=v2['primary_scint_exiting']
    v2['escape_status']='PASS' if primary['value']+3*primary['se']>=.36 and primary['value']-3*primary['se']<=.40 else 'FAIL'
    v2['escape_expected_interval']=[.36,.40]
    edges=np.linspace(0,1,21);expected=np.diff(edges**2)
    observed=hist.sum(axis=0)/hist.sum()
    boot_hist=weights@hist
    boot_hist/=boot_hist.sum(axis=1)[:,None]
    covariance=np.cov(boot_hist,rowvar=False,ddof=1)
    rank=int(np.linalg.matrix_rank(covariance))
    delta=observed-expected
    statistic=float(delta@np.linalg.pinv(covariance)@delta)
    pvalue=float(chi2.sf(statistic,rank))
    pearson=float(np.sum((hist.sum(axis=0)-hist.sum()*expected)**2/(hist.sum()*expected)))
    v2['cosine_gof']=dict(chi_square_event_bootstrap=statistic,rank=rank,p_value=pvalue,
        status='PASS' if pvalue>=.01 else 'FAIL',threshold_p=.01,
        pearson_chi_square_photons=pearson,pearson_df=19,
        pearson_p_value=float(chi2.sf(pearson,19)),photons=int(hist.sum()),
        mean_cosine_binned=float(np.sum(observed*(edges[:-1]+edges[1:])/2)),
        null='p(mu)=2mu',scope='first exiting scintillation encounter, includes sensor surfaces',
        applicability_caveat='First hits from a localized source in a finite bar need not follow equilibrium isotropic flux.')
    v2['status']='PASS' if v2['escape_status']=='PASS' and v2['cosine_gof']['status']=='PASS' else 'FAIL'
    v2['observation_checks']=checks
    sidecar(a.out,'first_encounter_cosine',dict(bin_low=edges[:-1],bin_high=edges[1:],
        photons=hist.sum(axis=0),observed_probability=observed,expected_probability=expected,
        probability_se=np.sqrt(np.diag(covariance))),meta)
    sidecar(a.out,'first_encounter_events',dict(event_id=np.arange(N),
        **{f'bin_{i:02}':hist[:,i] for i in range(20)},
        **{f'{name}_{kind}':values[:,i] for name,values in populations.items() for i,kind in enumerate(('encounters','escaped'))}),meta)
    sidecar(a.out,'first_encounter_pairs',dict(pair_source_status=list(pair_counts),
        photons=list(pair_counts.values())),meta)

    sensors=arrays(root['sipm_event_counts'])
    sidecar(a.out,'sipm_event_counts',sensors,meta)
    sensor_keys=86*sensors['event_id']+sensors['global_id']
    channel_shape_ok=np.array_equal(np.sort(sensor_keys),np.arange(N*86))
    channel_data={k:v[np.argsort(sensor_keys)].reshape(N,86) for k,v in sensors.items()}
    v5={}
    for name,ids in [('left',range(8)),('right',range(8,16)),('top',range(16,86)),('all',range(86))]:
        inc=channel_data['incident'][:,ids].sum(axis=1)
        det=channel_data['detected'][:,ids].sum(axis=1)
        expected_pde=channel_data['expected_surface_pde_sum'][:,ids].sum(axis=1)
        value=ratio(det,inc,weights)
        value.update(status=compare(value,.619397),emission_pde_reference=.619397,
            incident_spectrum_surface_pde=ratio(expected_pde,inc,weights),
            incident_scint_fraction=float(channel_data['incident_scint'][:,ids].sum()/inc.sum()),
            unmatched_detected=int(channel_data['unmatched_detected'][:,ids].sum()),
            unknown_surface_pde=int(channel_data['unknown_surface_pde'][:,ids].sum()),
            surface_detected=int(channel_data['surface_detection'][:,ids].sum()),
            surface_absorbed=int(channel_data['surface_absorption'][:,ids].sum()),
            transmitted=int(channel_data['transmitted'][:,ids].sum()),
            reflected_encounters_excluded=int(channel_data['reflected_encounters_excluded'][:,ids].sum()))
        value['secondary_surface_residual']=ratio(det-expected_pde,inc,weights)
        v5[name]=value
    v5['status']='PASS' if all(v5[x]['status']=='PASS' for x in ('left','right','top','all')) else 'FAIL'
    sums={k:v.sum(axis=0) for k,v in channel_data.items() if k not in ('event_id','global_id','face_type')}
    channel_ratio=sums['detected']/sums['incident']
    channel_boot=(weights@channel_data['detected'])/(weights@channel_data['incident'])
    sidecar(a.out,'sipm_channel_summary',dict(global_id=np.arange(86),**sums,
        detection_ratio=channel_ratio,ratio_se=channel_boot.std(axis=0,ddof=1)),meta)
    energy_ok=bool(np.all(np.isfinite(events['edep_total_MeV'])) and np.all(events['edep_total_MeV']>0)
        and np.all(events['produced_scint']>0) and np.all(events['produced_optical']>=events['produced_scint'])
        and np.allclose(events['edep_ionizing_MeV']+events['edep_nonionizing_MeV'],events['edep_total_MeV'],rtol=1e-14))
    first_ok=all(checks[x]==0 for x in ('invalid_surface','invalid_orientation','invalid_cosine','invalid_norm',
        'nonfinite','unexpected_outcome','duplicate_photon_keys')) and bool(np.all(hist.sum(axis=1)>0))
    channel_hits=np.column_stack([channel_data['detected'][:,:8].sum(axis=1),
        channel_data['detected'][:,8:16].sum(axis=1),channel_data['detected'][:,16:].sum(axis=1)])
    sensor_ok=channel_shape_ok and bool(np.array_equal(channel_hits,measured)) and \
        int(sensors['unmatched_detected'].sum())==0 and bool(np.all(sensors['detected']==sensors['detected_unique']))
    g2=dict(status='PASS' if energy_ok and first_ok and sensor_ok else 'FAIL',
        event_rows=N,first_encounter_rows=checks['rows'],sipm_event_rows=len(sensor_keys),
        energy_fields_valid=energy_ok,first_encounter_fields_valid=first_ok,sipm_fields_valid=sensor_ok,
        active_channels=86,channel_grid_complete=bool(channel_shape_ok))
    result.update(G_I_2=g2,V1=v1,V2=v2,V5=v5)
    # A malformed observation cannot be accepted merely because a ratio looks good.
    for name,valid in [('V1',energy_ok),('V2',first_ok),('V5',sensor_ok)]:
        if not valid:
            result[name]['numerical_hypothesis_status']=result[name]['status']
            result[name]['status']='NOT_VALIDATED_OBSERVATION_FAILURE'
    dump(a.out/'results.json',result)
    rows=[('G-I.1',g1['status']),('G-I.2',g2['status']),('G-I.3','PROPOSAL_REQUIRED' if g3['proposal_required'] else 'MEASURED')]
    rows += [(name,result[name]['status']) for name in ('V1','V2','V5')]
    sidecar(a.out,'validation_status',dict(test=[r[0] for r in rows],status=[r[1] for r in rows]),meta)
    print(json.dumps({k:v['status'] for k,v in result.items() if isinstance(v,dict) and 'status' in v}),flush=True)
    return 0 if g2['status']=='PASS' else 3  # hypothesis FAIL is a result, not a software test failure

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cell',type=Path,required=True)
    p.add_argument('--d1',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    return run(p.parse_args())

if __name__=='__main__':
    raise SystemExit(main())
