#!/usr/bin/env python3
"""Contrastes adicionales de D1–D4 a partir de los estadísticos ya guardados."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammaincc
from diagnose_decay import OUTPUT, TIME_EDGES, model_counts, csv_write, json_write

SCINT_COMPONENT = 1
TAIL_REFERENCE_MULTIPLIER = 5
FIXED_FRACTION_FLOOR = 1.e-8

def main():
    data = np.load(OUTPUT/'decay_sufficient_statistics.npz')
    summary = json.loads((OUTPUT/'d1_d4_summary.json').read_text())
    config = json.loads((OUTPUT/'effective_mpt.json').read_text())
    tau = config['constants_internal_units']['SCINTILLATIONTIMECONSTANT1']
    # Diferencias pareadas de los estimadores de exceso medio.
    n, s = data['tail_n'], data['tail_sum']
    m = s.sum(axis=0)/n.sum(axis=0)
    influence = (s-n*m)/n.sum(axis=0)
    gaps = []
    for j in range(1, len(m)):
        se = np.sqrt(len(n)/(len(n)-1)*np.sum((influence[:,j]-influence[:,j-1])**2))
        gaps.append(dict(lower_cut_ns=summary['cut_scan'][j-1]['cut_ns'],
                         upper_cut_ns=summary['cut_scan'][j]['cut_ns'],
                         tau_increase_ns=float(m[j]-m[j-1]),se_paired_event_ns=float(se)))
    fixed=[]
    for fit in summary['fits']:
        lo, hi = fit['window_low_ns'], fit['window_high_ns']
        sel = (TIME_EDGES[:-1] >= lo-1.e-8) & (TIME_EDGES[1:] <= hi+1.e-8)
        low, high = TIME_EDGES[:-1][sel], TIME_EDGES[1:][sel]
        y=data['event_hist'][:,sel].sum(axis=0).astype(float)
        def mu(p): return np.maximum(model_counts([p[0],tau,p[1],p[2]],low,high,lo),1.e-200)
        def residual(p):
            expected=mu(p); d=expected-y; nz=y>0
            d[nz]+=y[nz]*np.log(y[nz]/expected[nz])
            return np.sign(expected-y)*np.sqrt(2*np.maximum(d,0))
        result=least_squares(residual,[fit['amplitude_extrapolated_above_cut'],fit['tau2_ns']-tau,fit['f1']],
            bounds=([0,.05*tau,FIXED_FRACTION_FLOOR],[np.inf,500*tau,1-FIXED_FRACTION_FLOOR]),
            x_scale='jac',max_nfev=2000)
        if not result.success: raise RuntimeError(result.message)
        fixed.append(dict(window_low_ns=lo,window_high_ns=hi,tau1_fixed_ns=tau,
                          tau2_ns=float(tau+result.x[1]),f1=float(result.x[2]),
                          chi2=float(np.sum((y-mu(result.x))**2/mu(result.x))),ndf=len(y)-3,
                          delta_poisson_deviance=float(np.sum(result.fun**2)-fit['poisson_deviance'])))
        fit['chi2_asymptotic_p']=float(gammaincc(fit['ndf']/2,fit['chi2']/2))
        fit['bins_expected_below5']=int(np.sum(model_counts(fit['parameters'],low,high,lo)<5))
        cut=summary['cut_scan'][-1]['cut_ns']
        fit['extrapolated_count_above42ns']=float(fit['amplitude_extrapolated_above_cut']*(
            fit['f1']*np.exp(-(cut-lo)/fit['tau1_ns'])+fit['f2']*np.exp(-(cut-lo)/fit['tau2_ns'])))
    # El número efectivo por momento conserva el origen global t=0; no es una cola.
    summary.update(paired_cut_increases=gaps,fixed_fast_component_checks=fixed,
        tau_pool_moment_ns=summary['full_scint_pool_mean_creation_ns'],
        tau_pool_moment_se_event_ns=summary['full_scint_pool_se_event_ns'],
        tau_pool_moment_definition='Mean global creation time of all detected source_type=1 photons; t_ref=0 ns',
        tau_pool_moment_warning='Moment-matched exponential surrogate only; not a decay constant or a validated first-order-statistic scale')
    for row in summary['local_cut_scan']:
        if row['count']<2 or row['events_with_tail']<2:
            row['se_event_ns']=None
            row['status']='INSUFFICIENT_STATISTICS'
        else: row['status']='MEASURED'
    csv_write(OUTPUT/'d3_local_tail_stability.csv',summary['local_cut_scan'])
    csv_write(OUTPUT/'d2_paired_cut_increases.csv',gaps)
    csv_write(OUTPUT/'d4_fixed_fast_component.csv',fixed)
    json_write(OUTPUT/'d1_d4_summary.json',summary)
    print(json.dumps({'cut_increases':gaps,'fixed_fast_component':fixed},indent=2))

if __name__=='__main__': main()
