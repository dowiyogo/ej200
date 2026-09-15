#!/usr/bin/env python3
"""D2–D4: diagnóstico exclusivo de EJ200_xp0, sin modificar la compuerta."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shlex
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.optimize import least_squares
import uproot

from exec46_schema import CAMPAIGN_DIR, TREE_NAME

CELL_ID = 'EJ200_xp0'
OUTPUT = Path('/home/rrios/exec46_20260915/decay_diagnostics')
CHUNK = '128 MB'
CUT_MULTIPLIERS = np.array([3, 4, 5, 8, 12, 20])
TIME_BIN_NS = 0.1
TIME_MAX_NS = 200.0
TIME_EDGES = np.arange(0., TIME_MAX_NS + TIME_BIN_NS / 2, TIME_BIN_NS)
X_EDGES_MM = np.array([0., .1, .5, 1., 2., 5., 10., 20., 50., 100., 200., 400., 700., 1600.])
MAP_X_EDGES = np.r_[0., np.geomspace(1.e-6, 1600., 241)]
MAP_T_EDGES = np.r_[0., np.geomspace(1.e-3, 1.e6, 361)]
FIT_WINDOWS_TAU = [(3, 20), (3, 12), (5, 20)]
FIT_STARTS = [(1., 2., .8), (1., 5., .95), (.8, 10., .9)]
FIT_LOW_BOUND_TAU = .05
FIT_HIGH_BOUND_TAU = 500.
FRACTION_FLOOR = 1.e-8
UNRELIABLE_CHI2_NDF = 5.
FINITE_DIFFERENCE_STEP = 1.e-5
SCINTILLATION = 1
LOCALITY_CUT_MM = 0.1
DISPLACED_CUT_MM = 5.0


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def csv_write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def collect(root_path, events, tau):
    cuts = CUT_MULTIPLIERS * tau
    tail_n = np.zeros((events, len(cuts)))
    tail_sum = np.zeros_like(tail_n)
    tail_sq = np.zeros_like(tail_n)
    local_tail_n = np.zeros((2, events, len(cuts)))
    local_tail_sum = np.zeros_like(local_tail_n)
    event_hist = np.zeros((events, len(TIME_EDGES) - 1), dtype=np.int32)
    profile_n = np.zeros((events, len(X_EDGES_MM) - 1))
    profile_sum = np.zeros_like(profile_n)
    profile_sq = np.zeros_like(profile_n)
    profile_late = np.zeros_like(profile_n)
    x_time = np.zeros((len(MAP_X_EDGES) - 1, len(MAP_T_EDGES) - 1))
    profile_hist = np.zeros((len(X_EDGES_MM) - 1, len(TIME_EDGES) - 1))
    moments = np.zeros((4, 6))  # n, sum(x), sum(t), sum(x²), sum(t²), sum(xt)
    pool_moments = np.zeros((events, 3))  # n, sum(t), sum(t²), centelleo detectado
    extremes = [np.inf, -np.inf, 0.]
    overflow = 0
    rows = 0
    with uproot.open(root_path) as f:
        for a in f[TREE_NAME].iterate(['event_id', 'source_type', 't_creation_ns', 'x_creation_mm'],
                                     library='np', step_size=CHUNK):
            rows += len(a['event_id'])
            for source in range(4):
                sel = a['source_type'] == source
                t = a['t_creation_ns'][sel]
                x = np.abs(a['x_creation_mm'][sel])
                moments[source] += [len(t), x.sum(), t.sum(), (x*x).sum(), (t*t).sum(), (x*t).sum()]
            sel = a['source_type'] == SCINTILLATION
            t = a['t_creation_ns'][sel]
            x = np.abs(a['x_creation_mm'][sel])
            e = a['event_id'][sel]
            if not np.all(np.isfinite(t) & np.isfinite(x)):
                raise RuntimeError('Non-finite creation coordinate/time')
            extremes = [min(extremes[0], t.min()), max(extremes[1], t.max()), max(extremes[2], x.max())]
            for j, cut in enumerate(cuts):
                late = t >= cut
                excess = t[late] - cut
                tail_n[:, j] += np.bincount(e[late], minlength=events)
                tail_sum[:, j] += np.bincount(e[late], weights=excess, minlength=events)
                tail_sq[:, j] += np.bincount(e[late], weights=excess**2, minlength=events)
                for scope, spatial in enumerate([x < LOCALITY_CUT_MM, x >= DISPLACED_CUT_MM]):
                    chosen = late & spatial
                    local_tail_n[scope, :, j] += np.bincount(e[chosen], minlength=events)
                    local_tail_sum[scope, :, j] += np.bincount(e[chosen], weights=t[chosen]-cut, minlength=events)
            for j, values in enumerate([np.ones_like(t), t, t*t]):
                pool_moments[:, j] += np.bincount(e, weights=values, minlength=events)
            b = np.searchsorted(TIME_EDGES, t, side='right') - 1
            ok = (b >= 0) & (b < event_hist.shape[1])
            flat = e[ok].astype(np.int64) * event_hist.shape[1] + b[ok]
            # Acumulación exacta por evento, para errores que conservan sus correlaciones.
            np.add.at(event_hist.ravel(), flat, 1)
            overflow += int((~ok).sum())
            xb = np.searchsorted(X_EDGES_MM, x, side='right') - 1
            if np.any((xb < 0) | (xb >= profile_n.shape[1])):
                raise RuntimeError('Creation position outside registered profile range')
            pflat = e.astype(np.int64) * profile_n.shape[1] + xb
            for target, weights in [(profile_n, None), (profile_sum, t), (profile_sq, t*t),
                                    (profile_late, (t >= cuts[2]).astype(float))]:
                target += np.bincount(pflat, weights=weights, minlength=target.size).reshape(target.shape)
            profile_hist += np.histogram2d(x, t, bins=(X_EDGES_MM, TIME_EDGES))[0]
            x_time += np.histogram2d(x, t, bins=(MAP_X_EDGES, MAP_T_EDGES))[0]
        assert rows == f[TREE_NAME].num_entries
    return dict(tail_n=tail_n, tail_sum=tail_sum, tail_sq=tail_sq, local_tail_n=local_tail_n,
                local_tail_sum=local_tail_sum, event_hist=event_hist,
                profile_n=profile_n, profile_sum=profile_sum, profile_sq=profile_sq,
                profile_late=profile_late, profile_hist=profile_hist, x_time=x_time,
                moments=moments, pool_moments=pool_moments, extremes=np.array(extremes),
                overflow=np.array(overflow), rows=np.array(rows))


def ratio_cluster_se(sums, counts):
    mean = sums.sum() / counts.sum()
    influence = sums - mean * counts
    se = np.sqrt(len(counts)/(len(counts)-1) * (influence**2).sum()) / counts.sum()
    return float(mean), float(se)


def model_counts(p, low, high, origin):
    amplitude, tau1, gap, f1 = p
    tau2 = tau1 + gap
    return amplitude * (f1 * (np.exp(-(low-origin)/tau1) - np.exp(-(high-origin)/tau1))
                        + (1-f1) * (np.exp(-(low-origin)/tau2) - np.exp(-(high-origin)/tau2)))


def fit_tail(hist_by_event, tau, window):
    origin, upper = np.array(window) * tau
    selected = (TIME_EDGES[:-1] >= origin - 1.e-8) & (TIME_EDGES[1:] <= upper + 1.e-8)
    low, high = TIME_EDGES[:-1][selected], TIME_EDGES[1:][selected]
    origin, upper = float(low[0]), float(high[-1])
    y_event = hist_by_event[:, selected].astype(float)
    y = y_event.sum(axis=0)
    def residual(p):
        mu = np.maximum(model_counts(p, low, high, origin), 1.e-200)
        dev = mu-y
        nonzero = y > 0
        dev[nonzero] += y[nonzero] * np.log(y[nonzero]/mu[nonzero])
        return np.sign(mu-y) * np.sqrt(2*np.maximum(dev, 0.))
    trials = []
    for fast, slow, fraction in FIT_STARTS:
        start = [y.sum(), fast*tau, (slow-fast)*tau, fraction]
        trials.append(least_squares(residual, start,
            bounds=([0, FIT_LOW_BOUND_TAU*tau, FIT_LOW_BOUND_TAU*tau, FRACTION_FLOOR],
                    [np.inf, FIT_HIGH_BOUND_TAU*tau, FIT_HIGH_BOUND_TAU*tau, 1-FRACTION_FLOOR]),
            x_scale='jac', max_nfev=2000, ftol=1.e-11, xtol=1.e-11, gtol=1.e-9))
    fit = min(trials, key=lambda f: np.sum(f.fun**2))
    if not fit.success: raise RuntimeError(fit.message)
    p = fit.x
    mu = model_counts(p, low, high, origin)
    jac = np.empty((len(mu), len(p)))
    for j in range(len(p)):
        step = FINITE_DIFFERENCE_STEP * max(abs(p[j]), 1.)
        plus, minus = p.copy(), p.copy()
        plus[j] += step; minus[j] -= step
        jac[:,j] = (model_counts(plus,low,high,origin)-model_counts(minus,low,high,origin))/(2*step)
    bread = np.linalg.inv(jac.T @ (jac/mu[:,None]))
    scores = (y_event-y_event.mean(axis=0)) @ (jac/mu[:,None])
    cluster_cov = bread @ (scores.T@scores * len(scores)/(len(scores)-1)) @ bread
    transform = np.array([[0,1,0,0],[0,0,0,1],[0,1,1,0],[0,0,0,-1]])
    errors = np.sqrt(np.maximum(np.diag(transform@cluster_cov@transform.T), 0))
    ndf = len(y)-len(p)
    chi2 = float(np.sum((y-mu)**2/mu))
    result = dict(window_low_ns=origin, window_high_ns=upper, bin_width_ns=TIME_BIN_NS,
                  tau1_ns=float(p[1]), f1=float(p[3]), tau2_ns=float(p[1]+p[2]), f2=float(1-p[3]),
                  tau1_se_event_ns=float(errors[0]), f1_se_event=float(errors[1]),
                  tau2_se_event_ns=float(errors[2]), f2_se_event=float(errors[3]),
                  amplitude_extrapolated_above_cut=float(p[0]), chi2=chi2, ndf=ndf,
                  chi2_ndf=chi2/ndf, poisson_deviance=float(np.sum(fit.fun**2)),
                  unreliable=bool(chi2/ndf > UNRELIABLE_CHI2_NDF),
                  fraction_definition='Integrated mixture weights above window_low, extrapolated to infinity',
                  event_cluster_covariance=cluster_cov.tolist(), parameters=p.tolist())
    return result, low, high, y, mu


def figure_sidecars(stem, rows, root_items, meta):
    csv_write(OUTPUT/f'{stem}.csv', rows)
    with uproot.recreate(OUTPUT/f'{stem}.root') as f:
        for key, value in root_items.items(): f[key] = value
    json_write(OUTPUT/f'{stem}.meta.json', meta)


def main():
    start=time.monotonic()
    OUTPUT.mkdir(exist_ok=True)
    done=json.loads((CAMPAIGN_DIR/'cells'/CELL_ID/'.DONE').read_text())
    root_path=Path(done['root_path'])
    mpt=json.loads((OUTPUT/'effective_mpt.json').read_text())
    tau=mpt['constants_internal_units']['SCINTILLATIONTIMECONSTANT1']
    events=done['N_generated']
    digest=sha256(root_path)
    if digest != done['root_sha256']: raise RuntimeError('Source ROOT SHA mismatch')
    cache=OUTPUT/'decay_sufficient_statistics.npz'
    # Caché ligado por hash al ROOT, configuración y versión del recolector.
    cache_meta=OUTPUT/'decay_cache.meta.json'
    signature=dict(root_sha256=digest, collector_sha256=sha256(Path(__file__)), tau=tau,
                   time_edges=TIME_EDGES.tolist(), x_edges=X_EDGES_MM.tolist())
    if cache.exists() and cache_meta.exists() and json.loads(cache_meta.read_text())==signature:
        data=dict(np.load(cache))
    else:
        data=collect(root_path, events, tau)
        np.savez_compressed(cache, **data)
        json_write(cache_meta, signature)
    meta=dict(utc=datetime.now(timezone.utc).isoformat(), cell_id=CELL_ID, material=done['material'],
              root_source=str(root_path), root_sha256=digest, events=events, seeds=done['seeds'],
              simulation_command=done['command'], command=shlex.join([sys.executable,*sys.argv]),
              optical_model=dict(simulation_commit=done['simulation_commit'], effective_mpt=str(OUTPUT/'effective_mpt.json'),
                                 effective_mpt_sha256=sha256(OUTPUT/'effective_mpt.json')),
              source_selection='source_type == 1', other_cells_scanned=0)
    cuts=[]
    for j,k in enumerate(CUT_MULTIPLIERS):
        n=data['tail_n'][:,j]; s=data['tail_sum'][:,j]; ss=data['tail_sq'][:,j]
        mean,se=ratio_cluster_se(s,n)
        naive=np.sqrt(max((ss.sum()-n.sum()*mean**2)/(n.sum()-1),0)/n.sum())
        cuts.append(dict(cut_multiplier=int(k),cut_ns=float(k*tau),n_tail=int(n.sum()),
                         events_with_tail=int((n>0).sum()),tau_fit_ns=mean,se_event_ns=se,
                         se_iid_photon_ns=float(naive)))
    fig,ax=plt.subplots(figsize=(7,4))
    ax.errorbar([r['cut_ns'] for r in cuts],[r['tau_fit_ns'] for r in cuts],
                yerr=[r['se_event_ns'] for r in cuts],fmt='o-',label='Mean excess; event-cluster SE')
    ax.axhline(tau,ls='--',c='black',label='Effective MPT decay constant')
    ax.set(xlabel='Creation-time cut [ns]',ylabel='Mean excess above cut [ns]',yscale='log',title='EJ-200, x=0: decay-tail stability')
    ax.legend();fig.tight_layout();fig.savefig(OUTPUT/'d2_cut_stability.pdf');plt.close(fig)
    figure_sidecars('d2_cut_stability',cuts,{'cut_scan':{k:np.array([r[k] for r in cuts]) for k in cuts[0]}},
                    dict(meta,binning='unbinned, no upper time truncation',scale='linear x; log y',error='event-cluster ratio influence'))
    profile=[]
    for j,(lo,hi) in enumerate(zip(X_EDGES_MM[:-1],X_EDGES_MM[1:])):
        n=data['profile_n'][:,j];s=data['profile_sum'][:,j]
        if n.sum()==0: continue
        mean,se=ratio_cluster_se(s,n)
        profile.append(dict(abs_x_low_mm=float(lo),abs_x_high_mm=float(hi),count=int(n.sum()),
                            mean_t_creation_ns=mean,se_event_ns=se,late_cut_ns=float(5*tau),
                            late_count=int(data['profile_late'][:,j].sum()),
                            late_fraction=float(data['profile_late'][:,j].sum()/n.sum())))
    local_cuts=[]
    for scope, label in enumerate(['abs_x_lt_0p1mm', 'abs_x_ge_5mm']):
        for j, k in enumerate(CUT_MULTIPLIERS):
            n=data['local_tail_n'][scope,:,j]; sums=data['local_tail_sum'][scope,:,j]
            if n.sum() == 0: continue
            mean, se = ratio_cluster_se(sums,n)
            local_cuts.append(dict(scope=label,cut_ns=float(k*tau),count=int(n.sum()),
                                   events_with_tail=int((n>0).sum()),tau_fit_ns=mean,se_event_ns=se))
    csv_write(OUTPUT/'d3_local_tail_stability.csv',local_cuts)
    correlations=[]
    for scope,m in [('all_hits',data['moments'].sum(axis=0)),('scintillation',data['moments'][1]),('Cherenkov',data['moments'][2])]:
        n,sx,st,sxx,stt,sxt=m
        rho=(sxt-sx*st/n)/np.sqrt((sxx-sx*sx/n)*(stt-st*st/n))
        correlations.append(dict(scope=scope,count=int(n),pearson_absx_t=float(rho),mean_absx_mm=float(sx/n),mean_creation_ns=float(st/n)))
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    mesh=axes[0].pcolormesh(MAP_X_EDGES,MAP_T_EDGES,np.ma.masked_equal(data['x_time'].T,0),norm=LogNorm(),shading='flat')
    axes[0].set_xscale('symlog',linthresh=LOCALITY_CUT_MM)
    axes[0].set(yscale='log',xlim=(0,800),ylim=(.1,max(100,float(data['extremes'][1])*1.1)),xlabel='|Creation x| [mm]',ylabel='Creation time [ns]')
    fig.colorbar(mesh,ax=axes[0],label='Detected scintillation photons / bin')
    centers=np.array([(p['abs_x_low_mm']+p['abs_x_high_mm'])/2 for p in profile])
    axes[1].errorbar(centers,[p['mean_t_creation_ns'] for p in profile],yerr=[p['se_event_ns'] for p in profile],fmt='o-')
    axes[1].set(xscale='log',xlabel='|Creation x| bin center [mm]',ylabel='Mean creation time [ns]')
    fig.suptitle('EJ-200, x=0: source_type=1, exact untruncated profile means');fig.tight_layout();fig.savefig(OUTPUT/'d3_creation_locality.pdf');plt.close(fig)
    map_rows=[dict(abs_x_low_mm=float(MAP_X_EDGES[i]),abs_x_high_mm=float(MAP_X_EDGES[i+1]),
                   time_low_ns=float(MAP_T_EDGES[j]),time_high_ns=float(MAP_T_EDGES[j+1]),count=int(data['x_time'][i,j]))
              for i,j in zip(*np.nonzero(data['x_time']))]
    figure_sidecars('d3_creation_locality',map_rows,{'creation_x_time':(data['x_time'],MAP_X_EDGES,MAP_T_EDGES),
        'profile':{k:np.array([r[k] for r in profile]) for k in profile[0]},
        'time_by_x':(data['profile_hist'],X_EDGES_MM,TIME_EDGES)},dict(meta,binning=dict(x=MAP_X_EDGES.tolist(),t=MAP_T_EDGES.tolist(),profile_x=X_EDGES_MM.tolist()),scale='symlog x (0.1 mm linear core); log t/count color; profile log x',profile_csv='d3_creation_profile.csv'))
    csv_write(OUTPUT/'d3_creation_profile.csv',profile)
    csv_write(OUTPUT/'d3_correlations.csv',correlations)
    fits=[]; fit_rows=[]
    fig,axes=plt.subplots(2,1,figsize=(8,6),sharex=True,gridspec_kw={'height_ratios':[3,1]})
    root_items={'creation_time_scint':(data['event_hist'].sum(axis=0).astype(float),TIME_EDGES)}
    for index,window in enumerate(FIT_WINDOWS_TAU):
        result,low,high,y,mu=fit_tail(data['event_hist'],tau,window)
        fits.append(result)
        for l,h,n,m in zip(low,high,y,mu):
            fit_rows.append(dict(window_low_ns=result['window_low_ns'],window_high_ns=result['window_high_ns'],
                                 bin_low_ns=float(l),bin_high_ns=float(h),count=int(n),model=float(m)))
        root_items[f'fit_{index}']={'low_ns':low,'high_ns':high,'observed':y,'expected':mu}
        root_items[f'parameters_{index}']={k:np.array([result[k]]) for k in ['tau1_ns','tau2_ns','f1','f2','chi2','ndf']}
        if index==0:
            mid=(low+high)/2
            axes[0].errorbar(mid,y,yerr=np.sqrt(y),fmt='.',ms=2,label='Data')
            axes[0].plot(mid,mu,label='Two-exponential mixture')
            axes[0].set(yscale='log',ylabel='Photons / 0.1 ns',title=f"EJ-200 x=0: double exponential; chi2/ndf={result['chi2_ndf']:.2f}"+(' *' if result['unreliable'] else ''))
            axes[0].legend();axes[1].plot(mid,(y-mu)/np.sqrt(mu),'.',ms=2)
            axes[1].axhline(0,c='black',lw=.5);axes[1].set(xlabel='Creation time [ns]',ylabel='Pearson residual')
    fig.tight_layout();fig.savefig(OUTPUT/'d4_double_exponential.pdf');plt.close(fig)
    figure_sidecars('d4_double_exponential',fit_rows,root_items,dict(meta,binning=TIME_EDGES.tolist(),scale='log counts; linear residual',fit_windows_tau=FIT_WINDOWS_TAU,model='A*[f1/tau1*exp(-(t-cut)/tau1)+(1-f1)/tau2*exp(-(t-cut)/tau2)] integrated per bin',fraction_definition=fits[0]['fraction_definition'],fits=fits))
    pool_mean,pool_se=ratio_cluster_se(data['pool_moments'][:,1],data['pool_moments'][:,0])
    summary=dict(meta,cut_scan=cuts,local_cut_scan=local_cuts,correlations=correlations,profile=profile,fits=fits,
                 full_scint_pool_mean_creation_ns=pool_mean,full_scint_pool_se_event_ns=pool_se,
                 scint_min_time_ns=float(data['extremes'][0]),scint_max_time_ns=float(data['extremes'][1]),
                 scint_max_abs_x_mm=float(data['extremes'][2]),time_hist_overflow=int(data['overflow']),
                 source_counts=data['moments'][:,0].astype(int).tolist(),wall_s=time.monotonic()-start,
                 gate_changed=False,grid_resumed=False)
    json_write(OUTPUT/'d1_d4_summary.json',summary)
    print(json.dumps({k:summary[k] for k in ['cut_scan','correlations','full_scint_pool_mean_creation_ns','scint_max_time_ns','time_hist_overflow','wall_s']},indent=2))
    print(json.dumps([{k:v for k,v in f.items() if k not in ['event_cluster_covariance','parameters']} for f in fits],indent=2))

if __name__=='__main__': main()
