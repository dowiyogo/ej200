#!/usr/bin/env python3
"""Read-only transport reanalysis. Never launches a transport executable.

Run with the existing EXEC35 Python environment; see README.md.
All numerical acceptance choices are frozen in EXEC38_PREREGISTRATION.md.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.linalg import solve_triangular
from scipy.stats import chi2
import uproot

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE = Path('/home/rrios')
SEED, B = 38091301, 500
VREF = 299.792 / 1.58
BULK = {'EJ-200': 3800., 'EJ-204': 1600., 'EJ-230': 1200.}
YIELD = {'EJ-200': 10000., 'EJ-204': 10400., 'EJ-230': 9700.}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for x in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(x)
    return h.hexdigest()


def clean(x):
    if isinstance(x, dict): return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [clean(v) for v in x]
    if isinstance(x, np.ndarray): return clean(x.tolist())
    if isinstance(x, np.generic): return clean(x.item())
    if isinstance(x, float) and not math.isfinite(x): return None
    if isinstance(x, Path): return str(x)
    return x


def dump(p, x):
    Path(p).write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')


def latest(p, status='COMPLETE'):
    records = {}
    for line in Path(p).read_text().splitlines():
        r = json.loads(line)
        if r['status'] == status: records[r['cell_id']] = r
    return records


def verified(record):
    p = Path(record['output'])
    for name, h in record['sidecar_sha256'].items():
        if sha(p / name) != h: raise ValueError(f'hash mismatch: {p / name}')
    return p, json.loads((p / 'results.meta.json').read_text())


def bootstrap_indices(n):
    # A fresh identical generator deliberately pairs all equal-N populations.
    return np.random.default_rng(SEED).integers(0, n, size=(B, n), dtype=np.int32)


def sidecars(out, stem, rows, metadata, trees=None):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (out / (stem + '.csv')).open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for row in clean(rows):
            w.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})
    with uproot.recreate(out / (stem + '.root')) as f:
        f['rows_json'] = json.dumps(clean(rows), allow_nan=False)
        # Numeric columns are directly usable without parsing the JSON payload.
        numeric = {k: np.array([r.get(k, np.nan) if r.get(k) is not None else np.nan for r in rows], dtype=float)
                   for k in fields if any(isinstance(r.get(k), (int, float, np.number)) for r in rows)
                   and all(r.get(k) is None or isinstance(r.get(k), (int, float, np.number)) for r in rows)}
        if numeric: f['summary'] = numeric
        for key, values in (trees or {}).items(): f[key] = values
    dump(out / (stem + '.meta.json'), dict(metadata, rows=rows,
        files={str(out / (stem + ext)): sha(out / (stem + ext)) for ext in ['.csv', '.root']}))


def below(value, se, bound):
    if not np.isfinite(value + se) or value <= 0: return 'INVALID'
    if value + 3 * se < bound: return 'PASS'
    if value - 3 * se > bound: return 'FAIL'
    return 'INDETERMINATE'


def moments(t1, t2):
    valid = np.isfinite(t1) & np.isfinite(t2)
    a, b = t1[valid], t2[valid]
    if len(a) < 2: return np.full(10, np.nan)
    s1, s2 = np.std(a, ddof=1), np.std(b, ddof=1)
    sd = np.std(a - b, ddof=1)
    rho = np.corrcoef(a, b)[0, 1]
    q, factor = s1 / (sd / np.sqrt(2)), 1 / np.sqrt(1 - rho)
    exact = np.sqrt(2 * s1**2 / (s1**2 + s2**2 - 2 * rho * s1 * s2))
    return np.array([rho, s1, s2, sd, sd / np.sqrt(2), q, factor, q - factor, exact, s1 / s2])


V7_NAMES = ['rho', 'sigma_T1_ps', 'sigma_T2_ps', 'sigma_delta_ps', 'sigma_delta_over_sqrt2_ps',
            'Q_direct_over_normalized', 'F_equal_variance', 'Q_minus_F', 'Q_exact_unequal_variance', 'sigma_T1_over_T2']


def v7_cell(cell, root, meta, sim, n, indices, endonly=False):
    rows, trees = [], {}
    with uproot.open(root) as f:
        for end in ['left', 'right']:
            for variant in (['V2'] if endonly else ['V2', 'V1']):
                key = end + '_events' if endonly else f'{variant}_{end}_r0p5_f5_th4_events'
                a = f[key].arrays(library='np')
                assert np.array_equal(a['event_id'], np.arange(n))
                t1 = a['T1_ns'] * 1000 if endonly else a['T1_ps']
                t2 = a['T2_ns'] * 1000 if endonly else a['T2_ps']
                accepted = np.isfinite(t1) & np.isfinite(t2)
                assert np.array_equal(a['accepted'].astype(bool), accepted)
                old = next(r for r in meta['rows'] if r['end'] == end and r['variant'] == variant)
                assert old['rise_ns'] == .5 and old['fall_ns'] == 5 and old['threshold_PE'] == 4
                value = moments(t1, t2)
                boot = np.array([moments(t1[i], t2[i]) for i in indices])
                se = np.nanstd(boot, axis=0, ddof=1)
                frac = np.isfinite(boot).all(axis=1).mean()
                status = 'PASS' if value[0] - 3 * se[0] > 0 else ('FAIL' if value[0] + 3 * se[0] < 0 else 'INDETERMINATE')
                if frac < .95 or accepted.sum() < 200: status = 'INVALID'
                assert np.isclose(value[5], value[8], rtol=1e-10, atol=0)
                row = dict(cell_id=cell, material=sim['material'], x_mm=sim['x_mm'], end=end,
                    variant=variant, primary=variant == 'V2', group_ids=old['group_ids'],
                    N_generated=n, n_accepted=int(accepted.sum()), efficiency=accepted.mean(),
                    bootstrap_finite_fraction=frac, status=status, bootstrap_seed=SEED,
                    equal_variance_approximation='PASS' if abs(value[7]) <= 3 * se[7] else 'FAIL')
                for k, v, s in zip(V7_NAMES, value, se): row[k], row[k + '_se'] = v, s
                rows.append(row)
                trees[f'{cell}_{variant}_{end}_bootstrap'] = dict(zip(V7_NAMES, boot.T))
                trees[f'{cell}_{variant}_{end}_events'] = dict(event_id=np.arange(n), T1_ps=t1, T2_ps=t2, accepted=accepted)
    return rows, trees


def fit_exp(d, y, covariance, double=False, initial=None):
    chol = np.linalg.cholesky(covariance)
    whiten = lambda v: solve_triangular(chol, v, lower=True, check_finite=False)
    def model(p):
        p = np.exp(p)
        return p[0] * np.exp(-d / p[1]) + (p[2] * np.exp(-d / p[3]) if double else 0)
    bounds = np.log(([1e-6, 1.] * (2 if double else 1), [1e7, 1e6] * (2 if double else 1)))
    if initial is not None: starts = [initial]
    elif double: starts = [np.log([max(y), s, max(y), l]) for s, l in [(100, 1000), (300, 3000), (500, 10000)]]
    else: starts = [np.log([max(y) * 1.3, 500.])]
    fits = [least_squares(lambda p: whiten(model(p) - y), x, bounds=bounds,
                         max_nfev=2000, ftol=1e-10, xtol=1e-10, gtol=1e-10) for x in starts]
    fit = min(fits, key=lambda r: np.dot(r.fun, r.fun))
    p = np.exp(fit.x)
    if double and p[1] > p[3]: p = p[[2, 3, 0, 1]]
    k = len(p); score = float(np.dot(fit.fun, fit.fun))
    return dict(parameters=p, log_parameters=np.log(p), chi2=score, dof=7-k,
                p_value=chi2.sf(score, 7-k), AICc=score+2*k+2*k*(k+1)/(7-k-1),
                success=fit.success, boundary=bool(np.any(np.abs(fit.active_mask) > 0)))


def fits(profiles, out, provenance):
    v3, v4, consistency, trees = [], [], [], {}
    ratio_boot, ratios = {}, {}
    for material in BULK:
        cells = sorted([p for p in profiles if p['sim']['material'] == material], key=lambda p: p['sim']['x_mm'])
        x = np.array([p['sim']['x_mm'] for p in cells], dtype=float)
        assert len(x) == 7 and np.array_equal(x, [-650, -500, -200, 0, 200, 500, 650])
        for end_index, end in enumerate(['left', 'right']):
            print(f'FITS {material} {end}', flush=True)
            d = 700 + (x if end == 'left' else -x)
            counts = np.array([p['counts'][end_index] for p in cells]).T
            y = counts.mean(axis=0)
            cov = np.cov(counts, rowvar=False, ddof=1) / len(counts)
            yboot = np.array([p['count_boot'][:, end_index] for p in cells]).T
            r = fit_exp(d, y, cov)
            dr = fit_exp(d, y, cov, double=True)
            rb = [fit_exp(d, yy, cov, initial=r['log_parameters']) for yy in yboot]
            db = [fit_exp(d, yy, cov, double=True, initial=dr['log_parameters']) for yy in yboot]
            b1 = np.array([rr['parameters'] if rr['success'] else [np.nan]*2 for rr in rb])
            b2 = np.array([rr['parameters'] if rr['success'] else [np.nan]*4 for rr in db])
            se1, se2 = np.nanstd(b1, axis=0, ddof=1), np.nanstd(b2, axis=0, ddof=1)
            fraction1 = np.isfinite(b1).all(axis=1).mean(); fraction2 = np.isfinite(b2).all(axis=1).mean()
            lam = r['parameters'][1]; bound = below(lam, se1[1], BULK[material])
            if fraction1 < .95 or not r['success']: bound = 'INVALID'
            ratio_boot[material, end] = BULK[material] / b1[:, 1]
            ratios[material, end] = BULK[material] / lam
            pref = ('TWO_COMPONENT' if dr['p_value'] >= .01 and r['AICc']-dr['AICc'] >= 6
                    else 'SINGLE' if r['p_value'] >= .01 else 'NEITHER_ADEQUATE')
            row = dict(material=material, end=end, cell_ids=[p['sim']['cell_id'] for p in cells],
                nominal_attenuation_mm=BULK[material], amplitude=r['parameters'][0], amplitude_se=se1[0],
                lambda_mm=lam, lambda_se_mm=se1[1], nominal_over_effective=BULK[material]/lam,
                nominal_over_effective_se=np.nanstd(ratio_boot[material,end], ddof=1),
                chi2=r['chi2'], dof=r['dof'], p_value=r['p_value'], AICc=r['AICc'],
                bound_status=bound, model_status='PASS' if r['p_value']>=.01 else 'FAIL',
                preferred_model=pref, double_parameters=dr['parameters'], double_parameter_se=se2,
                double_chi2=dr['chi2'], double_dof=dr['dof'], double_p_value=dr['p_value'], double_AICc=dr['AICc'],
                double_at_boundary=dr['boundary'], bootstrap_finite_fraction=fraction1,
                double_bootstrap_finite_fraction=fraction2)
            v3.append(row)
            prefix = material.replace('-', '') + '_' + end
            trees[prefix+'_single_bootstrap'] = dict(amplitude=b1[:,0], lambda_mm=b1[:,1])
            trees[prefix+'_double_bootstrap'] = dict(A1=b2[:,0], lambda1_mm=b2[:,1], A2=b2[:,2], lambda2_mm=b2[:,3])
            trees[prefix+'_count_covariance'] = dict(row=np.repeat(np.arange(7),7), col=np.tile(np.arange(7),7), covariance=cov.ravel())
            mean = np.array([p['mean_time'][end_index] for p in cells])
            tboot = np.array([p['time_boot'][:,end_index] for p in cells]).T
            tcov = np.cov(tboot, rowvar=False, ddof=1)
            design = np.c_[np.ones(7),d]; inv = np.linalg.inv(tcov)
            estimator = np.linalg.solve(design.T @ inv @ design, design.T @ inv)
            par = estimator @ mean; pb = tboot @ estimator.T; vb = 1/pb[:,1]
            se = pb.std(axis=0,ddof=1); v = 1/par[1]; vse = vb.std(ddof=1)
            residual = mean - design @ par; stat = residual @ inv @ residual; pval = chi2.sf(stat,5)
            v4.append(dict(material=material,end=end,cell_ids=row['cell_ids'], intercept_ns=par[0], intercept_se_ns=se[0],
                slope_ns_per_mm=par[1], slope_se_ns_per_mm=se[1], signed_slope_vs_x=par[1]*(1 if end=='left' else -1),
                speed_mm_per_ns=v,speed_se_mm_per_ns=vse,speed_over_c_n=v/VREF,speed_over_c_n_se=vse/VREF,
                chi2=stat,dof=5,p_value=pval, bound_status=below(v,vse,VREF) if par[1]>0 else 'FAIL',
                model_status='PASS' if pval>=.01 else 'FAIL'))
            trees[prefix+'_time_bootstrap'] = dict(intercept_ns=pb[:,0],slope_ns_per_mm=pb[:,1],speed_mm_per_ns=vb)
            trees[prefix+'_time_covariance'] = dict(row=np.repeat(np.arange(7),7), col=np.tile(np.arange(7),7), covariance=tcov.ravel())
    for end in ['left','right']:
        y = np.array([ratios[m,end] for m in BULK]); boot = np.array([ratio_boot[m,end] for m in BULK]).T
        cov = np.cov(boot,rowvar=False); inv = np.linalg.inv(cov); ones = np.ones(3)
        avg = (ones@inv@y)/(ones@inv@ones); stat = (y-avg)@inv@(y-avg); pv=chi2.sf(stat,2)
        consistency.append(dict(end=end,common_ratio=avg,chi2=stat,dof=2,p_value=pv,status='PASS' if pv>=.01 else 'FAIL',
                                covariance=cov,scope='single-exponential descriptive ratios; model adequacy reported separately'))
    sidecars(out,'V3_attenuation',v3,provenance,trees)
    sidecars(out,'V3_material_consistency',consistency,provenance)
    sidecars(out,'V4_velocity',v4,provenance)
    return v3,v4,consistency


def analyze(out):
    out.mkdir(parents=True, exist_ok=False)
    provenance = dict(schema='EXEC38 analysis v1',start_utc=now(),command=shlex.join([sys.executable,*sys.argv]),
        analysis_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        code_sha256=sha(__file__),preregistration=str(HERE/'EXEC38_PREREGISTRATION.md'),
        preregistration_sha256=sha(HERE/'EXEC38_PREREGISTRATION.md'),preregistration_commit='4e08426',
        bootstrap=dict(seed=SEED,replicates=B,unit='generated event, paired within equal-N populations'),
        numpy=np.__version__,scipy=scipy.__version__,uproot=uproot.__version__,python=sys.version,
        statement='Existing data only. No transport executable is invoked; no electronics is injected.')
    simrecs=latest(BASE/'exec34r_20260912/manifest.jsonl','SIMULATION_COMPLETE')
    singles=latest(BASE/'exec36_20260913/single_manifest.jsonl'); groups=latest(BASE/'exec36_20260913/groups_manifest.jsonl')
    assert len(simrecs)==len(singles)==len(groups)==21
    indices=bootstrap_indices(10000); small=bootstrap_indices(2000)
    profiles=[]; v7=[]; v7trees={}; inputs=[]; observations=[]; schemas=[]; configs={}
    for cell in sorted(simrecs):
        print('READ',cell,flush=True)
        sim=simrecs[cell]; root=Path(sim['output'])/'photon_hits_run000.root'
        cp,cm=verified(singles[cell]); gp,gm=verified(groups[cell])
        assert cm['simulation']['root_sha256']==gm['simulation']['root_sha256']==sim['root_sha256']
        assert sim['N_generated']==10000 and sim['workers']==4 and sim['eventModulo']==1
        assert sim['seeds']==[26092601,8349041] and sim['simulation_commit']=='420addf0fd6029d5b2f0e235f472a8ae47f31fac'
        assert sim['N_TOP']==70 and sim['sptr_ns']==0 and sim['diagnostics'] is False
        assert root.is_file() and root.stat().st_size==sim['root_size_bytes']
        assert sha(Path(sim['output'])/'run.mac')==sim['macro_sha256']
        with uproot.open(root) as f:
            schemas.append(dict(cell_id=cell,path=str(root),keys=f.keys(),branches=f['sipm_hits'].keys()))
            assert f['sipm_hits'].num_entries==sim['root_entries']
        # The exact native END stream was already validated in EXEC36; rehash
        # the cache now, retaining that chain rather than rehashing ~80 GB twice.
        with uproot.open(cp/'results.root') as f:
            counts=np.zeros((2,10000)); sums=np.zeros((2,10000))
            for a in f['end_hits'].iterate(['event_id','face_type','global_id','time_ns'],step_size='80 MB',library='np'):
                ev=a['event_id']; face=a['face_type']; gid=a['global_id']; t=a['time_ns']
                assert np.all((ev>=0)&(ev<10000)) and np.isfinite(t).all()
                assert np.all((face==0)|(face==1)) and np.all((gid<8)==(face==0))
                for i in range(2):
                    mask=face==i; counts[i]+=np.bincount(ev[mask],minlength=10000)
                    sums[i]+=np.bincount(ev[mask],weights=t[mask],minlength=10000)
        assert counts[0].sum()==sim['left_total'] and counts[1].sum()==sim['right_total']
        cb=np.array([counts[:,ix].mean(axis=1) for ix in indices])
        tb=np.array([sums[:,ix].sum(axis=1)/counts[:,ix].sum(axis=1) for ix in indices])
        mean=sums.sum(axis=1)/counts.sum(axis=1)
        profiles.append(dict(sim=sim,counts=counts,time_sums=sums,count_boot=cb,time_boot=tb,mean_time=mean))
        config={k:sim[k] for k in ['cell_id','material','x_mm','N_generated','seeds','simulation_commit','workers','eventModulo',
            'configuration','N_TOP','sptr_ns','command','root_sha256','macro_sha256','PDE_sha256','opsc_code','parallelism_provenance']}
        configs[cell]=config
        inputs.append(dict(cell_id=cell,simulation=sim,cache=singles[cell],timestamps=groups[cell],native_recheck='size/schema/entries and macro hash; ROOT full hash inherited from verified EXEC36 cache chain'))
        log=Path(sim['output'])/'stdout.log'; text=log.read_text()
        generated=int(re.findall(r'Scint photons generated:\s*(\d+)',text)[-1])
        entering=int(re.findall(r'Bar -> SiPM \(entering\)\s*:\s*(\d+)',text)[-1])
        observations.append(dict(cell_id=cell,material=sim['material'],x_mm=sim['x_mm'],N_generated=10000,
            scintillation_photons_generated=generated,scintillation_photons_per_event=generated/10000,
            deposited_energy_MeV=None,production_ratio=None,incidents_all_sensors_log=entering,
            detected_left=sim['left_total'],detected_right=sim['right_total'],detected_top=sim['top_total'],
            incidents_left=None,incidents_right=None,log=str(log),log_sha256=sha(log),
            missing='No deposited-energy ledger, photon identities/first encounter, or independently identified per-end incident stream'))
        rows,trees=v7_cell(cell,gp/'results.root',gm,sim,10000,indices)
        v7.extend(rows);v7trees.update(trees)
    for cell,record in latest(BASE/'exec37_20260913/manifest.jsonl').items():
        print('READ',cell,flush=True)
        p,m=verified(record);src=m['source'];sim=src['simulation'];n=sim['N']
        assert n==2000 and sim['N_TOP']==0 and sim['seeds']==[26092601,8349041]
        root=Path(src['raw_ROOT']); assert sha(root)==src['root_sha256']
        with uproot.open(root) as f:
            schemas.append(dict(cell_id=cell,path=str(root),keys=f.keys(),branches=f['sipm_hits'].keys()))
        configs[cell]=dict(cell_id=cell,material=sim['material'],x_mm=sim['x_mm'],N_generated=n,seeds=sim['seeds'],
            simulation_commit=sim['source_commit'],workers=sim['threads'],eventModulo=sim.get('eventModulo'),
            configuration=sim['readout'],N_TOP=0,sptr_ns=0,command=sim['command'],root_sha256=src['root_sha256'],
            reflector=src['reflector'])
        inputs.append(dict(cell_id=cell,source=src,timestamps=record))
        rows,trees=v7_cell(cell,p/'results.root',m,sim,n,small,endonly=True)
        v7.extend(rows);v7trees.update(trees)
        log=Path(sim['cwd'])/'run.log';txt=log.read_text()
        generated=int(re.findall(r'Scint photons generated:\s*(\d+)',txt)[-1])
        observations.append(dict(cell_id=cell,material='EJ-204',x_mm=0,N_generated=n,
            scintillation_photons_generated=generated,scintillation_photons_per_event=generated/n,
            deposited_energy_MeV=None,production_ratio=None,npe_end_mean=m['npe_end_mean'],npe_end_sem=m['npe_end_sem'],
            incidents_left=None,incidents_right=None,log=str(log),log_sha256=sha(log),
            missing='No deposited-energy or first-encounter stream. Per-end archived boundary-state tallies lack photon/event identity and incident spectrum.'))
    dump(out/'inputs.json',inputs);dump(out/'configurations.json',configs);dump(out/'schema_inventory.json',schemas)
    provenance['inputs_sha256']=sha(out/'inputs.json');provenance['configurations']=configs
    sidecars(out,'V1_V2_V5_inventory',observations,provenance)
    sidecars(out,'V7_correlations',v7,provenance,v7trees)
    points=[];pointtrees={}
    for p in profiles:
        sim=p['sim'];cell=sim['cell_id']
        for i,end in enumerate(['left','right']):
            points.append(dict(cell_id=cell,material=sim['material'],x_mm=sim['x_mm'],end=end,
                distance_mm=700+(sim['x_mm'] if i==0 else -sim['x_mm']),N_generated=10000,
                npe=p['counts'][i].mean(),npe_sem=p['counts'][i].std(ddof=1)/100,
                mean_arrival_ns=p['mean_time'][i],mean_arrival_se_ns=p['time_boot'][:,i].std(ddof=1)))
        pointtrees[cell+'_events']=dict(event_id=np.arange(10000),npe_left=p['counts'][0],npe_right=p['counts'][1],
            time_sum_left_ns=p['time_sums'][0],time_sum_right_ns=p['time_sums'][1])
    sidecars(out,'profile_points',points,provenance,pointtrees)
    v3,v4,consistency=fits(profiles,out,provenance)
    # PDE prediction only, using configured tables, independent of hit data.
    pdepath=REPO/'data/sipm/AFBR-S4N66P024M_pde.txt';pde=np.loadtxt(pdepath);pderows=[]
    for material,code in [('EJ-200','100'),('EJ-204','101'),('EJ-230','106')]:
        sim=next(r for r in simrecs.values() if r['material']==material)
        sp=Path(sim['output'])/'sslg4/data/oscnt'/('opsc-'+code)/'scntComp1.txt'
        tab=np.loadtxt(sp);energy=1239.84193/tab[:,0];order=np.argsort(energy)
        energy=energy[order];density=tab[order,1];grid=np.linspace(energy.min(),energy.max(),100001)
        s=np.interp(grid,energy,density);q=np.interp(1239.84193/grid,pde[:,0],pde[:,1])
        avg=np.trapz(s*q,grid)/np.trapz(s,grid)
        wl=np.linspace(tab[:,0].min(),tab[:,0].max(),100001);s_w=np.interp(wl,tab[:,0],tab[:,1])
        wavelength_avg=np.trapz(s_w*np.interp(wl,pde[:,0],pde[:,1]),wl)/np.trapz(s_w,wl)
        pderows.append(dict(material=material,PDE_emission_energy_density=avg,PDE_emission_wavelength_density=wavelength_avg,
            prediction_uncertainty=None,measured_per_end_ratio=None,status='NOT_EVALUABLE',
            spectrum_path=str(sp),spectrum_sha256=sha(sp),PDE_path=str(pdepath),PDE_sha256=sha(pdepath),
            interpretation='Energy-density convention matches property-vector sampling; wavelength-density value is a separate convention, not an uncertainty. Incident spectrum absent.'))
    sidecars(out,'V5_PDE_prediction',pderows,provenance)
    provenance['end_utc']=now();dump(out/'analysis.meta.json',provenance)
    dump(out/'results.json',dict(V3=v3,V3_consistency=consistency,V4=v4,V7=v7,V5_prediction=pderows,inventory=observations))
    print('COMPLETE',out,flush=True)


def self_test():
    # Exact variance decomposition, including unequal variances and missing events.
    a=np.array([1.,2.,4.,8.,np.nan]);b=np.array([3.,1.,7.,5.,2.])
    m=moments(a,b);assert np.isclose(m[5],m[8],rtol=1e-12)
    assert not np.isclose(m[5],m[6]);assert np.isclose(m[1],np.std(a[:4],ddof=1))
    assert np.array_equal(bootstrap_indices(17),bootstrap_indices(17))
    assert below(5,.1,6)=='PASS' and below(7,.1,6)=='FAIL' and below(6,.1,6)=='INDETERMINATE'
    d=np.arange(7)*100.;y=200*np.exp(-d/450)
    fit=fit_exp(d,y,np.eye(7));assert fit['success'] and np.allclose(fit['parameters'],[200,450],rtol=1e-6)
    assert fit['chi2']<1e-12
    print('PASS: exact unequal-variance identity, paired population, bounds, analytical exponential recovery')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path);ap.add_argument('--self-test',action='store_true')
    args=ap.parse_args()
    if args.self_test:self_test()
    elif args.out:analyze(args.out)
    else:ap.error('--out or --self-test required')
