#!/usr/bin/env python
"""Re-analyze the existing EXEC34R ROOTs with preregistered EXEC35 widths."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import datetime
from functools import lru_cache
import hashlib
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
from robust_widths import measure, fixed_timestamps, NAMES, CORE_SOURCE, IQR_SOURCE
from validity_exec35 import verdict, RULE, RULE_PATH

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BASE = Path('/home/rrios/exec35_20260912')
SIM_MANIFEST = Path('/home/rrios/exec34r_20260912/manifest.jsonl')
OLD_MANIFEST = Path('/home/rrios/exec34c_20260912/analysis_manifest.jsonl')


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def clean(x):
    if isinstance(x, dict): return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)): return [clean(v) for v in x]
    if isinstance(x, np.generic): return clean(x.item())
    if isinstance(x, float) and not math.isfinite(x): return None
    return x


def write_json(path, x):
    Path(path).write_text(json.dumps(clean(x), indent=2, allow_nan=False)+'\n')


def latest(path, status):
    result = {}
    for line in path.read_text().splitlines():
        row = json.loads(line)
        if row['status'] == status:
            if row['cell_id'] not in result or row['utc'] > result[row['cell_id']]['utc']:
                result[row['cell_id']] = row
    return result


def csv_rows(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w') as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for row in clean(rows):
            writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k,v in row.items()})


@lru_cache(maxsize=1)
def end_callback():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    header = HERE/'../upstream/related/420addf/analysis/tb_mirror_sigma_vs_x.C'
    assert ROOT.gInterpreter.Declare('#include '+json.dumps(str(header.resolve())))
    # FitResult lives in an anonymous C++ namespace; expose only numeric values.
    assert ROOT.gInterpreter.Declare('''
    namespace exec35 {
    std::vector<double> fit(const std::vector<double>& values, int replica) {
        const auto core = FitCore(values, "exec35_core_"+std::to_string(replica));
        const auto peak = tbmirror::FitPeakSeeded(values, "exec35_peak_"+std::to_string(replica));
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return {core.usedFit ? core.sigmaPs : nan, peak.sigmaPs > 0 ? peak.sigmaPs : nan};
    }
    }
    ''')
    def callback(values_ps, replica):
        values = ROOT.std.vector('double')((values_ps/1000).tolist())
        return list(ROOT.exec35.fit(values,replica))
    return callback, ROOT.gROOT.GetVersion()


def annotate(row, old, arm, procedure, sample, n, gaussian_bootstrap=None):
    row.update(arm=arm, procedure=procedure, sample=sample, N=n, provenance_valid=True,
               sigma_gauss=old['sigma_ps'], sigma_gauss_old_uncertainty=old['uncertainty_ps'],
               sigma_gauss_se=old.get('bootstrap_error_ps'), sigma_gauss_bootstrap_fraction=None,
               gauss_fit_status=old.get('fit_status'), gauss_chi2_ndf=old['chi2_ndf'],
               old_gaussian_record=old)
    if gaussian_bootstrap is not None:
        finite = gaussian_bootstrap[np.isfinite(gaussian_bootstrap)]
        row.update(sigma_gauss_se=float(np.std(finite, ddof=1)),
                   sigma_gauss_bootstrap_fraction=len(finite)/len(gaussian_bootstrap),
                   gauss_fit_status=0 if old.get('fit_used_end', old['sigma_ps']>0) else 1)
    row['core_over_gauss'] = row['sigma_core']/row['sigma_gauss'] if row['sigma_gauss'] else np.nan
    row['core_over_gauss_conditional_se'] = row['sigma_core_se']/abs(row['sigma_gauss']) if row['sigma_gauss'] else np.nan
    row['primary_validity'] = verdict(row)
    row['gaussian_validity'] = verdict(row, gaussian=True)
    if arm == 'END':
        row['gaussian_compatible_diagnostic'] = bool(row['gaussian_validity']['status']=='VALID'
            and abs(row['asymmetry']) <= 3*row['asymmetry_se']
            and abs(row['sigma_core']-row['sigma_gauss']) <= 3*math.hypot(row['sigma_core_se'], row['sigma_gauss_se']))
    return row


def cell(job):
    cell_id, old_record, simulation, provenance = job
    started = now(); clock = time.monotonic()
    out = BASE/'cells'/cell_id
    out.mkdir(parents=True, exist_ok=False)
    old_dir = Path(old_record['output'])
    old_meta_path = old_dir/'analysis.meta.json'
    old = json.loads(old_meta_path.read_text())
    sim = old['simulation']
    assert sim['cell_id'] == cell_id
    assert sim['root_sha256'] == simulation['root_sha256']
    assert sim['root_path'] == str(Path(simulation['output'])/'photon_hits_run000.root')
    assert sim['N_generated'] == 10000 and sim['seeds'] == [26092601,8349041]
    assert sim['workers']==4 and sim['eventModulo']==1 and sim['diagnostics'] is False and sim['sptr_ns']==0
    assert sha(sim['root_path']) == sim['root_sha256']
    assert sha(sim['PDE_path']) == sim['PDE_sha256']
    old_hashes = {str(p): sha(p) for p in [old_meta_path, old_dir/'analysis.root', old_dir/'analysis.csv']}
    for name, record in old['files'].items():
        assert sha(record['path']) == record['sha256']
    # Match frozen model digest to the previous analysis's serialized model.
    canonical = old['top']['0']
    model = canonical['model']; biased = old['biased_full_sample']['model']
    assert hashlib.sha256(json.dumps(clean(model),sort_keys=True).encode()).hexdigest()==canonical['frozen_state_sha256']
    assert model['train_parity']==0 and biased['train_parity'] is None
    channels = sorted(set(model['channels']+biased['channels']))
    pieces = {key: [] for key in ['event_id','global_id','time_ns']}
    with uproot.open(sim['root_path']) as f:
        for arr in f['sipm_hits'].iterate([*pieces,'face_type'], step_size='64 MB', library='np'):
            mask = (arr['face_type']==2) & np.isin(arr['global_id'], channels)
            for key in pieces: pieces[key].append(arr[key][mask])
    arr = {key: np.concatenate(parts) for key,parts in pieces.items()}
    assert np.all((arr['event_id']>=0)&(arr['event_id']<10000)) and np.all(np.isfinite(arr['time_ns']))
    canonical_times = fixed_timestamps(arr['event_id'], arr['global_id'], arr['time_ns'], model['channels'])*1000
    biased_times = fixed_timestamps(arr['event_id'], arr['global_id'], arr['time_ns'], biased['channels'])*1000
    del pieces, arr
    rows = []
    with uproot.recreate(out/'widths.root') as output:
        for name, matrix in [('canonical_times', canonical_times), ('biased_times', biased_times)]:
            output[name] = dict(event_id=np.arange(10000,dtype=np.int32), **{'N'+str(n):matrix[:,n-1] for n in range(1,21)})
        for sample, values, curve, proc in [
            ('TRAIN',canonical_times[::2],model['train_curve'],'a'),
            ('EVAL',canonical_times[1::2],canonical['evaluation']['curve'],'a'),
            ('ALL',biased_times,biased['train_curve'],'b_scan')]:
            for n in range(1,21):
                point = next(r for r in curve if r['N']==n)
                r, bootstrap, edges, _ = measure(values[:,n-1])
                assert r['n_eff']==point['n_eff'] and r['N_generated_partition']==point['N_generated_partition']
                if 'fit_window' in point:
                    v = values[:,n-1]/1000
                    assert int(((v>=point['fit_window'][0])&(v<=point['fit_window'][1])).sum())==point['n_in_fit_window']
                row = annotate(r, point, 'TOP', proc, sample, n)
                row.update(cell_id=cell_id, distribution_id=f'TOP_{sample}_N{n}')
                rows.append(row)
                key = row['distribution_id']
                output[key+'_bootstrap'] = {name:bootstrap[:,j] for j,name in enumerate(NAMES)}
                finite = values[:,n-1][np.isfinite(values[:,n-1])]
                output[key+'_hist'] = (np.histogram(finite,edges)[0].astype(float),edges)
            print(cell_id, sample, '20 fixed indices complete', flush=True)
        for proc, sample, winner in [('b','ALL',biased['winner_N']),('c','EVAL',model['winner_N'])]:
            r = dict(next(r for r in rows if r['sample']==sample and r['N']==winner))
            r['procedure']=proc; rows.append(r)
        with uproot.open(old_dir/'analysis.root') as f:
            end = f['end_events'].arrays(library='np')
        assert np.array_equal(end['event_id'],np.arange(10000))
        delta = end['delta_ns']*1000
        assert np.array_equal(np.isfinite(delta),end['accepted'])
        output['end_events'] = dict(event_id=end['event_id'], delta_ps=delta)
        callback, root_version = end_callback()
        r, bootstrap, edges, extras = measure(delta, callback=callback)
        extras = np.asarray(extras)
        output['END_bootstrap'] = dict(**{name:bootstrap[:,j] for j,name in enumerate(NAMES)},
                                       gauss_FitCore=extras[:,0],gauss_FitPeakSeeded=extras[:,1])
        output['END_hist'] = (np.histogram(delta[np.isfinite(delta)],edges)[0].astype(float),edges)
        for j,name in enumerate(['FitCore','FitPeakSeeded']):
            assert r['n_eff']==old['end'][name]['n_eff']
            row = annotate(dict(r), old['end'][name], 'END', name, 'ALL', None, extras[:,j])
            row.update(cell_id=cell_id, distribution_id='END_DeltaT_LR')
            rows.append(row)
        output['rows_json'] = json.dumps(clean(rows),allow_nan=False)
    csv_rows(out/'widths.csv', rows)
    metadata = dict(schema='EXEC35 widths v1', cell_id=cell_id, simulation=sim,
        inherited_EXEC34C=old, previous_sidecar_sha256=old_hashes,
        preregistration=RULE, preregistration_sha256=sha(RULE_PATH),
        canonical_width='sigma_core', units='ps', provenance=provenance,
        start_utc=started,end_utc=now(),wall_s=time.monotonic()-clock, root_version=root_version,
        numpy_version=np.__version__,uproot_version=uproot.__version__,python_version=sys.version,
        conditional_uncertainty='Frozen existing channels and Gaussian-selected indices. No re-selection in bootstrap.',
        ratio_error='core bootstrap SE / fixed previous Gaussian point; denominator uncertainty not propagated',
        END='Widths of DeltaT_LR, not divided by sqrt(2); electronics excluded',
        paired_design='All cells share simulation seeds; no independent-cell or sigma(x) inference',
        rows=rows,files={name:sha(out/name) for name in ['widths.root','widths.csv']})
    write_json(out/'widths.meta.json', metadata)
    # Recheck the small inherited sidecars; raw ROOT was hashed before use and never opened writable.
    assert all(sha(path)==digest for path,digest in old_hashes.items())
    return dict(utc=now(),cell_id=cell_id,status='COMPLETE',output=str(out),
                sidecar_sha256={name:sha(out/name) for name in ['widths.root','widths.csv','widths.meta.json']})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs',type=int,default=4)
    args = parser.parse_args()
    simulations = latest(SIM_MANIFEST,'SIMULATION_COMPLETE')
    analyses = latest(OLD_MANIFEST,'COMPLETE')
    assert len(simulations)==len(analyses)==21 and simulations.keys()==analyses.keys()
    sources = [Path(__file__),HERE/'robust_widths.py',HERE/'validity_exec35.py',RULE_PATH,CORE_SOURCE,IQR_SOURCE,
        HERE/'../upstream/related/420addf/analysis/congruent_sum4_timing.C',
        HERE/'../upstream/related/420addf/analysis/tb_mirror_sigma_vs_x.C']
    provenance = dict(command=shlex.join([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]]),
        cwd=str(Path.cwd()),pipeline_commit=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip(),
        script_sha256={str(p.resolve()):sha(p) for p in sources},
        manifests={str(p):sha(p) for p in [SIM_MANIFEST,OLD_MANIFEST]},
        preregistration_copy_sha256=sha(BASE/'preregistration.md'),
        environment_repair='venv --system-site-packages; pip --ignore-installed --no-deps numpy==1.23.5: complete same-version wheel, global packages unchanged')
    journal = BASE/'analysis_manifest.jsonl'
    completed = latest(journal,'COMPLETE') if journal.exists() else {}
    jobs = []
    for cell_id in sorted(analyses):
        if cell_id in completed:
            previous = completed[cell_id]
            assert all(sha(Path(previous['output'])/name)==digest for name,digest in previous['sidecar_sha256'].items())
            assert json.loads((Path(previous['output'])/'widths.meta.json').read_text())['provenance']['script_sha256']==provenance['script_sha256']
        else:
            jobs.append((cell_id,analyses[cell_id],simulations[cell_id],provenance))
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(cell,job):job[0] for job in jobs}
        for future in as_completed(futures):
            record = future.result()
            with journal.open('a') as f:
                f.write(json.dumps(record)+'\n');f.flush();os.fsync(f.fileno())
            print(record['cell_id'],'COMPLETE',flush=True)


if __name__=='__main__':
    main()
