#!/usr/bin/env python
"""Evaluate precisely G-P.1..5; diagnostic biases do not add a sixth gate."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys
import numpy as np
import uproot
from run_simulation import now, sha

CRITERION = '''El piloto PASA si y solo si se cumplen las cinco condiciones:
1. sigma_END finito, con `fit_used_end` verdadero, chi2/ndf registrado
   y n_eff por encima del mínimo del ajuste.
2. sigma_TOP finito con `fit_status` válido en al menos 10 de los 20
   índices N, en AMBAS mitades de la partición.
3. Eficiencia por encima del `EFFICIENCY_FLOOR` configurado para el N
   ganador de (c).
4. Los tres sidecars completos y el comando de extremo a extremo
   reproducible desde el README.
5. La diferencia (b) − (c) es finita y calculable.

NO hay criterio sobre el VALOR de sigma_t: no existe referencia válida
contra la cual juzgarlo. Un valor inesperado no es un fallo.'''


def finite(value):
    return value is not None and math.isfinite(value)


def evaluate(directory):
    out = Path(directory).resolve()
    data = json.loads((out/'analysis.meta.json').read_text())
    end = data['end']['FitCore']; split = data['top']['0']
    train, evals = split['model']['train_curve'],split['evaluation']['curve']
    conditions = []
    def add(index, valid, evidence):
        conditions.append(dict(condition=f'G-P.{index}',status='PASS' if valid else 'FAIL',evidence=evidence))
    add(1, finite(end['sigma_ps']) and end['fit_used_end'] and finite(end['chi2_ndf']) and end['n_eff']>20,
        dict(sigma_END_ps=end['sigma_ps'],fit_used_end=end['fit_used_end'],chi2_ndf=end['chi2_ndf'],n_eff=end['n_eff'],minimum=20))
    valid_counts = {name:sum(finite(r['sigma_ps']) and r['fit_status']==0 for r in curve)
                    for name,curve in [('TRAIN_even',train),('EVAL_odd',evals)]}
    add(2, all(n>=10 for n in valid_counts.values()),valid_counts)
    primary = split['evaluation']['primary']
    efficiency = primary['efficiency'] if primary else None
    floor = data['config']['EFFICIENCY_FLOOR']
    add(3, efficiency is not None and efficiency>floor,dict(efficiency=efficiency,efficiency_floor=floor,winning_N=split['model']['winner_N']))
    checks = {}
    try:
        for filename in ('analysis.root','analysis.csv'):
            checks[filename+'_hash'] = sha(out/filename)==data['files'][filename]['sha256']
        rows = list(csv.DictReader((out/'analysis.csv').open()))
        with uproot.open(out/'analysis.root') as f:
            checks['root_csv_row_count'] = f['result_summary'].num_entries==len(rows)
            saved = f['result_summary']['sigma_ps'].array(library='np')
            expected = np.array([float(r['sigma_ps']) if r['sigma_ps'] else float('nan') for r in rows])
            checks['root_csv_values'] = bool(np.array_equal(saved,expected,equal_nan=True))
            checks['all_generated_end_events'] = f['end_events'].num_entries==data['N_generated']
        checks['metadata_complete'] = all(key in data for key in ('simulation','pipeline_commit','orchestration_commit','split_rule','top','end','stages','rng','script_sha256'))
        checks['simulation_provenance'] = all(key in data['simulation'] for key in ('simulation_commit','geant4_version','seeds','N_generated','material','opsc_code','configuration','N_TOP','x_mm','workers','eventModulo','PDE_path','PDE_sha256'))
        checks['stage_exit_codes'] = all(r['exit_code']==0 and r['command'] and r['start_utc'] and r['end_utc'] for r in data['stages'])
        checks['analysis_scripts_match'] = all(sha(path)==digest for path,digest in data['script_sha256'].items())
        package = Path(__file__).resolve().parents[1]
        readme = (package/'README.md').read_text()
        checks['readme_end_to_end_command'] = 'orchestration/run_cell.py' in readme and '--binary' in readme and '--output' in readme
        checks['native_input_name'] = Path(data['simulation']['root_path']).name=='photon_hits_run000.root'
        checks['separate_arm_results'] = {'TOP','END'}.issubset({r['arm'] for r in rows})
        checks['finite_uncertainties'] = primary is not None and finite(primary['uncertainty_ps']) and finite(end['uncertainty_ps'])
    except Exception as error:
        checks['exception'] = str(error)
    add(4, bool(checks) and all(value is True for value in checks.values()),checks)
    add(5,finite(data['b_minus_c']['difference_ps']),dict(b_minus_c_ps=data['b_minus_c']['difference_ps']))
    result = dict(gate_status='PASS' if all(c['status']=='PASS' for c in conditions) else 'FAIL',
                  evaluated_utc=now(),criterion_verbatim=CRITERION,conditions=conditions,
                  sidecar_sha256={str(out/name):sha(out/name) for name in ('analysis.root','analysis.csv','analysis.meta.json')},
                  diagnostics_not_additional_gates=dict(partition_symmetry=data['partition_symmetry'],bias_diagnostic=data['bias_diagnostic']),
                  grid_started=False)
    (out/'gate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--analysis',required=True)
    raise SystemExit(0 if evaluate(p.parse_args().analysis)['gate_status']=='PASS' else 34)
