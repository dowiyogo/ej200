#!/usr/bin/env python
"""Analyze a native simulation ROOT using the closed EXEC34 TOP/END decisions."""
import argparse
import csv
import datetime
import json
import math
import resource
from pathlib import Path
import shlex
import subprocess
import sys
import time

import numpy as np
import uproot
import ROOT
import top_split as top
from run_simulation import now, sha

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
LABEL = 'INTRINSIC timing resolution — electronics not included'


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value): return None
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False)+'\n')


def git(*args):
    return subprocess.check_output(['git','-C',str(REPO),*args], text=True).strip()


def difference(first, second):
    if first is None or second is None: return dict(difference_ps=None, uncertainty_ps=None, significance=None)
    diff = first['sigma_ps']-second['sigma_ps']
    err = math.hypot(first['uncertainty_ps'],second['uncertainty_ps'])
    return dict(difference_ps=diff, uncertainty_ps=err,
                significance=diff/err if err > 0 else None,
                uncertainty_assumption='quadrature of conditional errors; cross-fit/selection covariance unavailable; diagnostic significance only')


def curvature(curve, index):
    if index is None: return dict(points=[], curve_flat=False, reason='No valid winner')
    points = [r for r in curve if abs(r['N']-index) <= 1]
    winner = next(r for r in points if r['N'] == index)
    comparisons = [dict(N=r['N'], **difference(r,winner)) for r in points if r['N'] != index]
    # Preregistered diagnostic only: adjacent difference <= one quadrature error.
    flat = bool(comparisons and all(abs(r['difference_ps']) <= r['uncertainty_ps'] for r in comparisons))
    return dict(points=[{k:r[k] for k in ('N','sigma_ps','uncertainty_ps','fit_status','efficiency')} for r in points],
                comparisons=comparisons, curve_flat=flat,
                finite_difference_curvature_ps=(points[0]['sigma_ps']-2*points[1]['sigma_ps']+points[2]['sigma_ps']) if len(points)==3 else None,
                boundary_minimum=len(points)<3,
                rule='Adjacent sigma differences <= one quadrature uncertainty: noise-dominated index; no out-of-range index invented',
                interpretation='Winning N is not interpretable as a physical optimum' if flat else 'This diagnostic does not establish a physical optimum')


def analyze(args):
    start = now(); wall = time.monotonic()
    simdir = Path(args.simulation).resolve()
    sim = json.loads((simdir/'simulation.meta.json').read_text())
    raw = Path(sim['root_path'])
    assert raw.name == 'photon_hits_run000.root' and not raw.is_symlink()
    assert sim['exit_code'] == 0 and sim['N_generated'] == 10000
    assert sha(raw) == sim['root_sha256']
    assert sha(sim['PDE_path']) == sim['PDE_sha256']
    assert not git('diff','--','analysis/sigma_t/upstream')
    n_generated = sim['N_generated']
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    rootpath = out/'analysis.root'
    root = ROOT.TFile(str(rootpath), 'RECREATE')
    stages = [dict(stage='simulation', command=sim['command'], cwd=sim['cwd'],
                   start_utc=sim['start_utc'], end_utc=sim['end_utc'], exit_code=sim['exit_code'])]
    read_start = now()
    pieces = {name:[] for name in ('event_id','face_type','global_id','time_ns')}
    counts = np.zeros((3,n_generated), dtype=np.int64)
    with uproot.open(raw) as f:
        for arr in f['sipm_hits'].iterate(list(pieces), library='np', step_size='100 MB'):
            assert np.all((arr['event_id'] >= 0) & (arr['event_id'] < n_generated))
            assert np.all(np.isfinite(arr['time_ns']))
            assert np.all(np.isin(arr['face_type'],[0,1,2]))
            for face in range(3):
                mask = arr['face_type'] == face
                counts[face] += np.bincount(arr['event_id'][mask], minlength=n_generated)
            mask = arr['face_type'] == 2
            assert np.all((arr['global_id'][mask]>=16)&(arr['global_id'][mask]<=85))
            for name in pieces: pieces[name].append(arr[name][mask])
    arr = {name:np.concatenate(value) for name,value in pieces.items()}
    del pieces
    assert int(counts[0].sum()) == sim['left_total'] and int(counts[1].sum()) == sim['right_total']
    npe = (counts[0]+counts[1])/2
    stages.append(dict(stage='adapter_read_native_ROOT', command='analyze.py: uproot sipm_hits iterator; explicit generated IDs 0..9999', start_utc=read_start, end_utc=now(), exit_code=0))
    results, rows = {}, []
    for parity in (0,1):
        phase_start = now()
        print(f'TOP TRAIN parity={parity}', flush=True)
        model = top.learn(arr,n_generated,parity,root_file=root)
        fingerprint = sha_text(model)
        evaluation = top.evaluate(arr,n_generated,model,root_file=root)
        assert fingerprint == sha_text(model), 'EVAL mutated frozen TRAIN model'
        train_winner = next((r for r in model['train_curve'] if r['N']==model['winner_N']),None)
        results[str(parity)] = dict(model=model, evaluation=evaluation,
            train_at_winner=train_winner,
            internal_selection_bias=difference(train_winner,evaluation['primary']),
            curvature_check=curvature(model['train_curve'],model['winner_N']),
            frozen_state_sha256=fingerprint)
        for row in model['train_curve']+evaluation['curve']:
            rows.append(dict(arm='TOP', procedure='a', train_parity=parity, **row))
        if evaluation['primary']:
            rows.append(dict(arm='TOP',procedure='c',train_parity=parity,**evaluation['primary']))
        stages.append(dict(stage=f'TOP_c_train_{parity}',command='top_split.learn -> freeze -> top_split.evaluate; no EVAL optimization',start_utc=phase_start,end_utc=now(),exit_code=0))
    print('TOP biased full-sample comparator', flush=True)
    biased_model = top.learn(arr,n_generated,None,root_file=root)
    biased = top.winner(biased_model['train_curve'])
    for row in biased_model['train_curve']:
        rows.append(dict(arm='TOP',procedure='b_scan',selection_label='biased by selection',**row))
    if biased:
        rows.append(dict(arm='TOP',procedure='b',selection_label='biased by selection',**biased))
    canonical = results['0']
    primary = canonical['evaluation']['primary']
    bmc = difference(biased,primary)
    internal = canonical['internal_selection_bias']
    b_diff, i_diff = bmc['difference_ps'], internal['difference_ps']
    same_sign = (b_diff*i_diff > 0) if b_diff is not None and i_diff is not None else False
    ratio = abs(b_diff/i_diff) if i_diff and b_diff is not None else None
    diagnostic = dict(same_sign=same_sign, magnitude_ratio=ratio,
        same_order_of_magnitude=bool(ratio is not None and .1 <= ratio <= 10),
        convention='TRAIN-selected TRAIN sigma minus EVAL sigma; b is all-generated-events same-sample argmin',
        correlation_note='Full-sample b and split fits overlap; these are not statistically independent bias estimates',
        interpretation='consistent by sign and factor-ten diagnostic' if same_sign and ratio is not None and .1<=ratio<=10 else 'EXEC34 requested diagnostic mismatch: possible implementation defect; not a physical result; not an additional G-P condition')
    symmetry = difference(results['0']['evaluation']['primary'],results['1']['evaluation']['primary'])
    symmetry['convention'] = 'sigma(c) TRAIN even - sigma(c) TRAIN odd'
    symmetry['compatible_at_2sigma'] = bool(symmetry['significance'] is not None and abs(symmetry['significance']) <= 2)
    symmetry['note'] = 'Opposite EVAL samples are disjoint, but trained selections cross-use the halves. No averaging. Two-sigma diagnostic threshold is not a G-P condition.'
    del arr
    print('END preserved FitCore and FitPeakSeeded', flush=True)
    phase_start = now()
    assert ROOT.gInterpreter.Declare('#include '+json.dumps(str(HERE/'end_bridge.h')))
    values = list(ROOT.exec34.analyze_end(str(raw),n_generated,root))
    core = dict(sigma_ps=values[0],uncertainty_ps=values[1],chi2_ndf=values[2],
                n_eff=int(values[3]),fit_used_end=bool(values[4]),mean_ns=values[5],rms_ps=values[6])
    peak = dict(sigma_ps=values[7],uncertainty_ps=values[8],chi2_ndf=values[9],
                n_eff=int(values[10]),mean_ns=values[11],fwhm_ps=values[12],
                fit_status_note='Imported MirrorFit does not expose status; positive sigma indicates its status==0 branch')
    for name, result in [('FitCore',core),('FitPeakSeeded',peak)]:
        result.update(efficiency=result['n_eff']/n_generated,
            discarded_nonfinite_end=n_generated-result['n_eff'],
            nonfinite_left=int(values[13]),nonfinite_right=int(values[14]),
            sigma_delta_over_sqrt2_ps=result['sigma_ps']/math.sqrt(2),
            sigma_delta_over_sqrt2_error_ps=result['uncertainty_ps']/math.sqrt(2),
            normalization='Primary sigma(DeltaT_LR); derived sigma(DeltaT_LR)/sqrt(2)',
            derived_assumption='equal and statistically independent timing contributions from the two ends; NOT verified here',
            measurement_label=LABEL)
        rows.append(dict(arm='END',procedure=name,**result))
    stages.append(dict(stage='END',command='exec34::analyze_end -> LeadingEdgeTime/Earliest -> FitCore and tbmirror::FitPeakSeeded',start_utc=phase_start,end_utc=now(),exit_code=0))
    end_difference = difference(peak,core)
    end_difference['convention'] = 'FitPeakSeeded - FitCore, both sigma(DeltaT_LR)'
    end_difference['uncertainty_assumption'] = 'Same-event fits correlated; quadrature is diagnostic only, covariance not estimated'
    ROOT.TObjString(json.dumps(clean(dict(rows=rows)))).Write('result_rows_json')
    root.Close()
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (out/'analysis.csv').open('w') as f:
        writer = csv.DictWriter(f,fieldnames=fields); writer.writeheader()
        for row in clean(rows):
            writer.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    with uproot.update(rootpath) as f:
        f['result_summary'] = {'row_id':np.arange(len(rows),dtype=np.int32),
            'sigma_ps':np.array([r.get('sigma_ps',float('nan')) for r in rows]),
            'uncertainty_ps':np.array([r.get('uncertainty_ps',float('nan')) for r in rows])}
    command = shlex.join([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]])
    stages.append(dict(stage='analysis',command=command,cwd=str(Path.cwd()),start_utc=start,end_utc=now(),exit_code=0))
    files = {name:dict(path=str(out/name),sha256=sha(out/name)) for name in ('analysis.root','analysis.csv')}
    data = dict(schema='EXEC34 analysis v1',label=LABEL, simulation=sim,
        pipeline_commit=git('rev-parse','HEAD'),
        orchestration_commit=git('log','-1','--format=%H','--','analysis/sigma_t/orchestration/top_split.py'),
        primitive_commits=['30eba3b','fde610f'], root_version=ROOT.gROOT.GetVersion(),
        python_version=sys.version, numpy_version=np.__version__, uproot_version=uproot.__version__,
        config=top.CONFIG, canonical_train_parity='even',canonical_eval_parity='odd',
        split_rule='TRAIN event_id % 2 == 0; EVAL event_id % 2 == 1; reversed orientation as diagnostic; denominator 5000 each including zero hits',
        top=results, biased_full_sample=dict(model=biased_model,result=biased,label='biased by selection'),
        b_minus_c=bmc,partition_symmetry=symmetry,bias_diagnostic=diagnostic,
        end=dict(FitCore=core,FitPeakSeeded=peak,systematic_difference=end_difference,
            map=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]],
            reduction='first finite SUM4 cluster crossing at each end; require both finite',
            rise_ns=.5,fall_ns=5,threshold_amplitude_PE=4,sptr=False,walk=False,ToT_cut=False,
            FitCore_identity='upstream/related/420addf/analysis/congruent_sum4_timing.C:110; four iterative +/-2 sigma fits, n>=20',
            FitPeakSeeded_identity='upstream/related/420addf/analysis/tb_mirror_sigma_vs_x.C:105; peak +/-2 ns, n>=20'),
        electronics_excluded_reason='At least four live electronics variants; unresolved FWHM/sigma ambiguity; SPTR_PROVENANCE.md does not establish sqrt(kN) propagation for order statistics',
        npe_end_mean=float(npe.mean()),npe_end_sem=float(npe.std(ddof=1)/math.sqrt(n_generated)),
        N_generated=n_generated,zero_hit_events=int(np.count_nonzero(counts.sum(axis=0)==0)),
        stages=stages,start_utc=start,end_utc=now(),analysis_wall_s=time.monotonic()-wall,exit_code=0,
        analysis_max_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        analysis_rss_method='getrusage(RUSAGE_SELF).ru_maxrss, Linux KiB',
        files=files,script_sha256={str(p):sha(p) for p in [HERE/'analyze.py',HERE/'top_split.py',HERE/'end_bridge.h']},
        rng=dict(simulation_seeds=sim['seeds'],bootstrap_generator='numpy.default_rng PCG64',
                 bootstrap_seed=top.CONFIG['RANDOM_SEED'],replicates=top.CONFIG['N_BOOTSTRAP'],
                 seed_selection='fixed configuration, not data-derived; reinitialized for each fixed-estimator bootstrap',
                 root_fit='deterministic default ROOT fitter; no injected jitter'),
        interpretation='TOP uncertainty is conditional on the TRAIN-learned estimator; not unconditional training-selection uncertainty')
    write_json(out/'analysis.meta.json',data)
    print(json.dumps(clean(dict(primary=primary,biased=biased,b_minus_c=bmc,end=data['end'],symmetry=symmetry,bias_diagnostic=diagnostic)),indent=2),flush=True)
    return 0


def sha_text(value):
    import hashlib
    return hashlib.sha256(json.dumps(clean(value),sort_keys=True).encode()).hexdigest()


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulation',required=True)
    p.add_argument('--output',required=True)
    raise SystemExit(analyze(p.parse_args()))
