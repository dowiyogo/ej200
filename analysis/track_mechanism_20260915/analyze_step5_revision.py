#!/usr/bin/env python3
"""E1--E5 revision. Analysis of existing events only; stop before Step 6."""
import json
import math
from array import array
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot
import ROOT
import analyze_step5 as base

ROOT.gROOT.SetBatch(True)
OUT = base.OUTPUT_DIR
POSITIONS = base.POSITIONS_MM
ABS_X = np.array([0, 200, 500, 650])
CENTER = 3
PROFILE_BINS = 40
PROFILE_PADDING = 0.5
FIT_CHI2_LIMIT = 5.0
BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_BATCH = 40
BOOTSTRAP_SEED = 460516
NS_TO_PS = 1000.0
PROFILE_REFERENCE_ATOL = 2e-10
VARIANTS = ('uniform_pol1', 'uniform_pol2', 'quantile_pol1', 'quantile_pol2')
BOOT_FIELDS = ('n', 'y', 's', 'c', 'fL', 'muCL', 'muSL', 'fR', 'muCR', 'muSR', 'btotal', 'bs', 'bc')
INPUT_BRANCHES = ['material_code', 'x_mm', 't0_ns', 'npe_end', 't_left_ns',
                  't_right_ns', 'npe_scint_left', 'npe_scint_right',
                  'npe_cherenkov_left', 'npe_cherenkov_right',
                  'first_left_source_type', 'first_right_source_type']


def table(frame, columns, digits=4):
    labels = [label.replace('|', r'\|') for label in columns.values()]
    lines = ['| ' + ' | '.join(labels) + ' |', '|' + '|'.join(['---'] * len(labels)) + '|']
    for _, row in frame.iterrows():
        values = []
        for key in columns:
            value = row[key]
            values.append(f'{value:.{digits}f}' if isinstance(value, (float, np.floating)) else str(value))
        lines.append('| ' + ' | '.join(values) + ' |')
    return '\n'.join(lines)


def even(values):
    """Position axis is last. Preserve event/bootstrap axes."""
    return np.stack([values[..., CENTER], (values[..., 2]+values[..., 4])/2,
                     (values[..., 1]+values[..., 5])/2,
                     (values[..., 0]+values[..., 6])/2], axis=-1)


def expand_even(values):
    return values[..., [3, 2, 1, 0, 1, 2, 3]]


def chain(n, beta):
    return np.concatenate([np.zeros_like(n[..., :1]), np.cumsum(
        np.diff(n, axis=-1)*(beta[..., 1:]+beta[..., :-1])/2, axis=-1)], axis=-1)


def between_residual(n, y):
    nc, yc = n-n.mean(axis=-1, keepdims=True), y-y.mean(axis=-1, keepdims=True)
    beta = np.sum(nc*yc, axis=-1)/np.sum(nc*nc, axis=-1)
    r = y-y[..., CENTER:CENTER+1]-beta[..., None]*(n-n[..., CENTER:CENTER+1])
    return beta, r


def local_excess(r, n):
    """Prespecified diagnostic: residual at 500 above Npe-linear chord 200--650."""
    w = (n[..., 2]-n[..., 1])/(n[..., 3]-n[..., 1])
    return r[..., 2] - ((1-w)*r[..., 1]+w*r[..., 3])


def ols_counts(n, s, c, y):
    X = np.column_stack([s-s.mean(), c-c.mean()])
    yc = y-y.mean()
    inv = np.linalg.inv(X.T@X)
    b = inv@X.T@yc
    e = yc-X@b
    cov = len(y)/(len(y)-3)*inv@((X*e[:, None]).T@(X*e[:, None]))@inv
    bn = np.dot(n-n.mean(), yc)/np.sum((n-n.mean())**2)
    return bn, b, cov, np.linalg.cond(X.T@X)


def profile_fits(cell, material, x, root_file):
    n, y = cell['n'], cell['y']
    rows, slopes = [], {}
    root_file.cd()
    for mode in ('uniform', 'quantile'):
        if mode == 'uniform':
            edges = np.linspace(n.min()-PROFILE_PADDING, n.max()+PROFILE_PADDING, PROFILE_BINS+1)
        else:
            edges = np.unique(np.quantile(n, np.linspace(0, 1, PROFILE_BINS+1)))
            edges[0], edges[-1] = n.min()-PROFILE_PADDING, n.max()+PROFILE_PADDING
        name = f'{material.replace("-", "")}_{x:+d}_{mode}'
        p = ROOT.TProfile(name, name, len(edges)-1, array('d', edges))
        for nv, yv in zip(n, y):
            p.Fill(float(nv), float(yv))
        p.Write()
        for degree in (1, 2):
            f = ROOT.TF1(name+f'_pol{degree}', f'pol{degree}', float(edges[0]), float(edges[-1]))
            result = p.Fit(f, 'QRSN')
            base.require(int(result) == 0 and result.Get(), f'Profile fit failed: {name} pol{degree}')
            f.Write()
            result.Get().GetCovarianceMatrix().Write(f.GetName()+'_covariance')
            gradient = np.array([0., 1.] + ([2*n.mean()] if degree == 2 else []))
            cov = np.array([[result.CovMatrix(i, j) for j in range(degree+1)] for i in range(degree+1)])
            beta = f.GetParameter(1)+(2*f.GetParameter(2)*n.mean() if degree == 2 else 0)
            variant = f'{mode}_pol{degree}'
            slopes[variant] = beta
            rows.append(dict(material=material, x_mm=int(x), variant=variant, bins=len(edges)-1,
                             mean_npe=n.mean(), derivative_ns_pe=beta,
                             derivative_se_ns_pe=np.sqrt(gradient@cov@gradient),
                             p0_ns=f.GetParameter(0), p1_ns_pe=f.GetParameter(1),
                             p2_ns_pe2=f.GetParameter(2) if degree == 2 else 0.,
                             chi2=f.GetChisquare(), ndf=f.GetNDF(),
                             chi2_ndf=f.GetChisquare()/f.GetNDF(),
                             unreliable=f.GetChisquare()/f.GetNDF()>FIT_CHI2_LIMIT,
                             status=int(result), covariance_status=result.CovMatrixStatus()))
    return rows, slopes


def sample_stats(cell, rng):
    """Ordinary paired event bootstrap within each cell; preserve both ENDs."""
    result = np.empty((BOOTSTRAP_REPLICATES, len(BOOT_FIELDS)))
    n_events = len(cell['y'])
    for start in range(0, BOOTSTRAP_REPLICATES, BOOTSTRAP_BATCH):
        stop = min(start+BOOTSTRAP_BATCH, BOOTSTRAP_REPLICATES)
        idx = rng.integers(0, n_events, size=(stop-start, n_events))
        n, y, s, c = (cell[k][idx] for k in ('n','y','s','c'))
        values = [v.mean(axis=1) for v in (n, y, s, c)]
        for end in ('L','R'):
            is_c, t = cell['isC'+end][idx], cell['t'+end][idx]
            count = is_c.sum(axis=1)
            base.require(np.all((count>0)&(count<n_events)), 'Bootstrap empty source class')
            values.extend([count/n_events, (t*is_c).sum(axis=1)/count,
                           (t*(~is_c)).sum(axis=1)/(n_events-count)])
        nc, yc, sc, cc = (v-v.mean(axis=1, keepdims=True) for v in (n,y,s,c))
        ss, scov, ccov = (np.sum(v, axis=1) for v in (sc*sc,sc*cc,cc*cc))
        sy, cy = (np.sum(v*yc, axis=1) for v in (sc,cc))
        det = ss*ccov-scov*scov
        values.extend([np.sum(nc*yc,axis=1)/np.sum(nc*nc,axis=1),
                       (ccov*sy-scov*cy)/det, (ss*cy-scov*sy)/det])
        result[start:stop] = np.stack(values, axis=1)
    return result


def mixture_terms(stats):
    """Exact finite-interval product rule; position axis last; units ns."""
    total = np.zeros((*stats['y'].shape, 3))
    by_end = {}
    for end in ('L','R'):
        f, c, s = (stats[key+end] for key in ('f','muC','muS'))
        avg = lambda v: (v[..., 1:]+v[..., :-1])/2
        terms = np.stack([avg(f)*np.diff(c), (1-avg(f))*np.diff(s),
                          (avg(c)-avg(s))*np.diff(f)], axis=-1)
        integ = np.concatenate([np.zeros_like(terms[..., :1, :]), np.cumsum(terms,axis=-2)],axis=-2)
        integ -= integ[..., CENTER:CENTER+1, :]
        by_end[end] = (terms, integ)
        total += integ/2
    expected = stats['y']-stats['y'][..., CENTER:CENTER+1]
    base.require(np.max(np.abs(total.sum(axis=-1)-expected)) < 1e-11, 'Mixture closure failed')
    return total, by_end


def save_frame(name, frame):
    frame.to_csv(OUT/(name+'.csv'),index=False,float_format='%.12g')


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    before_hash = base.sha256(base.DERIVED_ROOT)
    with uproot.open(base.DERIVED_ROOT) as f:
        base.require(f['derived_events'].num_entries == base.EXPECTED_ROWS,'Bad event count')
        arrays = f['derived_events'].arrays(INPUT_BRANCHES,library='np')
    base.require(np.array_equal(arrays['t0_ns'],(arrays['t_left_ns']+arrays['t_right_ns'])/2),'Clock identity')
    total_sources = sum(arrays[k] for k in ('npe_scint_left','npe_scint_right','npe_cherenkov_left','npe_cherenkov_right'))
    base.require(np.array_equal(arrays['npe_end'],total_sources),'Source counts do not close')
    slopes, points, curves, summary = base.analyze(arrays,pd.read_csv(base.BASELINE_POINTS),pd.read_csv(base.BASELINE_FITS))
    base.make_figures(slopes,points,curves,summary)
    save_frame('identification_residual_points',points)
    reference = pd.read_csv(base.BASE_DIR/'step2/baseline_cells.csv').query('clock == "time_ns"')
    rows_profile, rows_count, rows_cells, rows_intervals, rows_mix = [], [], [], [], []
    rows_sensitive, rows_fits, rows_signif, rows_boot, rows_local = [], [], [], [], []
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    all_stats, all_boot, all_variants = {}, {}, {}
    with ROOT.TFile(str(OUT/'profile_refits.root'),'RECREATE') as root_file:
        for code, material in enumerate(base.MATERIALS):
            cells, boots = [], []
            profiles = {v:[] for v in VARIANTS}
            print('Processing '+material,flush=True)
            for x in POSITIONS:
                mask = (arrays['material_code']==code)&(arrays['x_mm']==x)
                cell = {key:arrays[branch][mask].astype(float) for key,branch in
                        [('n','npe_end'),('y','t0_ns'),('tL','t_left_ns'),('tR','t_right_ns')]}
                for key, name in [('s','scint'),('c','cherenkov')]:
                    cell[key] = (arrays['npe_'+name+'_left'][mask]+arrays['npe_'+name+'_right'][mask]).astype(float)
                for end,name in [('L','left'),('R','right')]:
                    source = arrays['first_'+name+'_source_type'][mask]
                    base.require(np.all(np.isin(source,[1,2])),'Unclassified winner')
                    cell['isC'+end] = source==2
                prows, pslopes = profile_fits(cell,material,x,root_file)
                rows_profile.extend(prows)
                for v in VARIANTS: profiles[v].append(pslopes[v])
                ref = reference[(reference.material==material)&(reference.x_mm==x)].iloc[0]
                base.require(abs(pslopes['uniform_pol1']-ref.slope_ns_pe)<PROFILE_REFERENCE_ATOL,'Original slope mismatch')
                bn,b,cov,cond = ols_counts(cell['n'],cell['s'],cell['c'],cell['y'])
                rows_count.append(dict(material=material,x_mm=int(x),beta_total_ps_pe=NS_TO_PS*bn,
                                       beta_scint_ps_pe=NS_TO_PS*b[0],beta_cherenkov_ps_pe=NS_TO_PS*b[1],
                                       se_scint_ps_pe=NS_TO_PS*np.sqrt(cov[0,0]),
                                       se_cherenkov_ps_pe=NS_TO_PS*np.sqrt(cov[1,1]),
                                       covariance_ps2_pe2=NS_TO_PS**2*cov[0,1],condition_number=cond,
                                       rho_counts=np.corrcoef(cell['s'],cell['c'])[0,1]))
                stats = {k:cell[k].mean() for k in ('n','y','s','c')}
                stats.update(btotal=bn,bs=b[0],bc=b[1])
                for end in ('L','R'):
                    ic = cell['isC'+end]
                    stats.update({key+end:val for key,val in [('f',ic.mean()),('muC',cell['t'+end][ic].mean()),('muS',cell['t'+end][~ic].mean())]})
                stats['sem'] = cell['y'].std(ddof=1)/np.sqrt(len(cell['y']))
                rows_cells.append(dict(material=material,x_mm=int(x),**stats))
                cells.append(stats)
                boots.append(sample_stats(cell,rng))
            stats = {k:np.array([c[k] for c in cells]) for k in cells[0]}
            bt = np.stack(boots,axis=1)
            bst = {k:bt[:,:,i] for i,k in enumerate(BOOT_FIELDS)}
            all_stats[material],all_boot[material] = stats,bst
            bbetween, rbetween = between_residual(stats['n'],stats['y'])
            bbboot, rbboot = between_residual(bst['n'],bst['y'])
            rse = rbboot.std(axis=0,ddof=1)
            for j,x in enumerate(POSITIONS):
                rows_signif.append(dict(material=material,x_mm=int(x),residual_ps=NS_TO_PS*rbetween[j],
                                        cell_sem_ps=NS_TO_PS*stats['sem'][j],
                                        residual_over_cell_sem=rbetween[j]/stats['sem'][j],
                                        paired_bootstrap_se_ps=NS_TO_PS*rse[j],
                                        paired_bootstrap_z=rbetween[j]/rse[j] if rse[j]>0 else np.nan,
                                        ci025_ps=NS_TO_PS*np.quantile(rbboot[:,j],.025),
                                        ci975_ps=NS_TO_PS*np.quantile(rbboot[:,j],.975)))
            mapper_design = np.column_stack([np.ones(len(POSITIONS)),(POSITIONS/1000)**2])
            mapper_weights = 1/stats['sem']**2
            mapper_between = np.linalg.inv((mapper_design.T*mapper_weights)@mapper_design)@(mapper_design.T*mapper_weights)
            summary.loc[summary.material==material,'descriptive_residual_a2_total_error_ns_per_m2'] = (rbboot@mapper_between.T)[:,1].std(ddof=1)
            mix,by_end = mixture_terms(stats)
            bmix,_ = mixture_terms(bst)
            for end,(terms,integ) in by_end.items():
                for j in range(len(POSITIONS)-1):
                    dx_m = (POSITIONS[j+1]-POSITIONS[j])/1000
                    rows_intervals.append(dict(material=material,end=end,x_low_mm=int(POSITIONS[j]),x_high_mm=int(POSITIONS[j+1]),
                        cherenkov_ps=NS_TO_PS*terms[j,0],scint_ps=NS_TO_PS*terms[j,1],mixing_ps=NS_TO_PS*terms[j,2],
                        cherenkov_ps_m=NS_TO_PS*terms[j,0]/dx_m,scint_ps_m=NS_TO_PS*terms[j,1]/dx_m,mixing_ps_m=NS_TO_PS*terms[j,2]/dx_m))
            for j,x in enumerate(POSITIONS):
                for k,name in enumerate(('cherenkov','scintillation','mixing')):
                    rows_mix.append(dict(material=material,x_mm=int(x),term=name,shift_ps=NS_TO_PS*mix[j,k],
                        bootstrap_se_ps=NS_TO_PS*bmix[:,j,k].std(ddof=1),
                        left_shift_ps=NS_TO_PS*by_end['L'][1][j,k],right_shift_ps=NS_TO_PS*by_end['R'][1][j,k]))
            ens,eys = even(stats['n']),even(stats['y'])
            predictions = {v:chain(ens,even(np.array(profiles[v]))) for v in VARIANTS}
            predictions['event_ols_total'] = chain(ens,even(stats['btotal']))
            predictions['event_ols_two_counts'] = chain(even(stats['s']),even(stats['bs']))+chain(even(stats['c']),even(stats['bc']))
            base_reference = points[points.material==material].sort_values('abs_x_mm').original_prediction_shift_ps.to_numpy()/NS_TO_PS
            base.require(np.max(np.abs(predictions['uniform_pol1']-base_reference))<1e-8,'Original chain mismatch')
            all_variants[material] = predictions
            ebr = even(rbetween)
            boot_er = even(rbboot)
            ebmix = even(np.moveaxis(bmix, -1, -2)) # B,3,4
            emix = even(mix.T) # 3,4
            # Project the exact mixing term using the SAME fitted slope operator:
            # linearity gives r_between(y)=r_between(mix)+r_between(y-mix).
            _, r_mix = between_residual(stats['n'],mix[:,2])
            _, rb_mix = between_residual(bst['n'],bmix[:,:,2])
            for j,ax in enumerate(ABS_X):
                rows_boot.append(dict(material=material,abs_x_mm=int(ax),descriptive_residual_ps=NS_TO_PS*ebr[j],
                    residual_bootstrap_se_ps=NS_TO_PS*boot_er[:,j].std(ddof=1),
                    mixing_shift_ps=NS_TO_PS*emix[2,j],mixing_se_ps=NS_TO_PS*ebmix[:,2,j].std(ddof=1),
                    projected_mixing_ps=NS_TO_PS*even(r_mix)[j],
                    projected_mixing_se_ps=NS_TO_PS*even(rb_mix)[:,j].std(ddof=1),
                    residual_minus_projected_mix_ps=NS_TO_PS*(ebr[j]-even(r_mix)[j]),
                    residual_minus_projected_mix_se_ps=NS_TO_PS*(even(rbboot-rb_mix)[:,j].std(ddof=1)),
                    projected_fraction_bootstrap_se=(even(rb_mix)[:,j]/boot_er[:,j]).std(ddof=1) if j else np.nan,
                    original_chain_residual_ps=points[(points.material==material)&(points.abs_x_mm==ax)].registered_remnant_ps.iloc[0],
                    original_chain_minus_raw_mix_ps=points[(points.material==material)&(points.abs_x_mm==ax)].registered_remnant_ps.iloc[0]-NS_TO_PS*emix[2,j],
                    raw_mixing_to_residual_ratio=emix[2,j]/ebr[j] if j else np.nan,
                    projected_mixing_to_residual_ratio=even(r_mix)[j]/ebr[j] if j else np.nan,
                    cherenkov_shift_ps=NS_TO_PS*emix[0,j],scintillation_shift_ps=NS_TO_PS*emix[1,j],
                    observed_shift_ps=NS_TO_PS*(eys[j]-eys[0])))
            for variant,pred in predictions.items():
                for label,indices in [('even',None),('negative',[3,2,1,0]),('positive',[3,4,5,6])]:
                    nn = ens if indices is None else stats['n'][indices]
                    yy = eys if indices is None else stats['y'][indices]
                    if indices is None:
                        pp = pred
                    elif variant in VARIANTS:
                        pp = chain(nn,np.array(profiles[variant])[indices])
                    elif variant == 'event_ols_total':
                        pp = chain(nn,stats['btotal'][indices])
                    else:
                        pp = chain(stats['s'][indices],stats['bs'][indices])+chain(stats['c'][indices],stats['bc'][indices])
                    rr = yy-yy[0]-pp
                    boot_nn = even(bst['n']) if indices is None else bst['n'][:,indices]
                    boot_yy = even(bst['y']) if indices is None else bst['y'][:,indices]
                    rows_local.append(dict(material=material,variant=variant,mirror=label,
                        chain_local_excess_ps=NS_TO_PS*local_excess(rr,nn),
                        target_local_excess_ps=NS_TO_PS*local_excess(yy-yy[0],nn),
                        target_local_excess_bootstrap_se_ps=NS_TO_PS*local_excess(boot_yy,boot_nn).std(ddof=1)))

                r = eys-eys[0]-pred
                signed = stats['y']-stats['y'][CENTER]-expand_even(pred)
                ff = base.fit_even(POSITIONS,signed,stats['sem'])
                boot_pred = None
                if variant=='event_ols_total': boot_pred = chain(even(bst['n']),even(bst['btotal']))
                if variant=='event_ols_two_counts': boot_pred = chain(even(bst['s']),even(bst['bs']))+chain(even(bst['c']),even(bst['bc']))
                brows = dict(material=material,variant=variant,a2_ps_m2=NS_TO_PS*ff['a2_ns_per_m2'],
                             conditional_a2_se_ps_m2=NS_TO_PS*ff['a2_error_ns_per_m2'],chi2_ndf=ff['chi2_ndf'],
                             localized_excess_500_ps=NS_TO_PS*local_excess(r,ens))
                if boot_pred is not None:
                    boot_y = even(bst['y']); br = boot_y-boot_y[:,:1]-boot_pred
                    brows['localized_excess_bootstrap_se_ps'] = NS_TO_PS*local_excess(br,even(bst['n'])).std(ddof=1)
                    design = np.column_stack([np.ones(len(POSITIONS)),(POSITIONS/1000)**2])
                    weights = 1/stats['sem']**2
                    mapper = np.linalg.inv((design.T*weights)@design)@(design.T*weights)
                    bsigned = bst['y']-bst['y'][:,CENTER:CENTER+1]-expand_even(boot_pred)
                    brows['a2_bootstrap_se_ps_m2'] = NS_TO_PS*(bsigned@mapper.T)[:,1].std(ddof=1)
                rows_fits.append(brows)
                for j,ax in enumerate(ABS_X):
                    rows_sensitive.append(dict(material=material,variant=variant,abs_x_mm=int(ax),
                        predicted_ps=NS_TO_PS*pred[j],chain_residual_ps=NS_TO_PS*r[j],
                        shift_from_original_ps=NS_TO_PS*(predictions['uniform_pol1'][j]-pred[j]),
                        chain_residual_bootstrap_se_ps=NS_TO_PS*br[:,j].std(ddof=1) if boot_pred is not None else np.nan,
                        descriptive_between_residual_ps=NS_TO_PS*ebr[j]))
    frames = dict(profile_refits=pd.DataFrame(rows_profile),two_count_slopes=pd.DataFrame(rows_count),
                  cell_statistics=pd.DataFrame(rows_cells),mixture_intervals=pd.DataFrame(rows_intervals),
                  mixture_components=pd.DataFrame(rows_mix),chain_sensitivity=pd.DataFrame(rows_sensitive),
                  chain_fit_summary=pd.DataFrame(rows_fits),residual_significance=pd.DataFrame(rows_signif),
                  mixture_target=pd.DataFrame(rows_boot),localization_sensitivity=pd.DataFrame(rows_local))
    for name,frame in frames.items(): save_frame(name,frame)
    save_frame('within_between_summary',summary)
    # Numeric companion stores every table plus categorical string columns.
    with uproot.recreate(OUT/'step5_revision_tables.root') as f:
        for name,frame in frames.items():
            f[name] = {key:(frame[key].astype(str).to_numpy(dtype=str) if frame[key].dtype == object else frame[key].to_numpy()) for key in frame.columns}
    make_figures(frames)
    render_report(frames,summary,points,all_stats)
    base.require(base.sha256(base.DERIVED_ROOT)==before_hash,'Read-only input hash changed')
    meta = dict(created_utc=datetime.now(timezone.utc).isoformat(),status='CHAIN_RULE_WITHIN_SLOPE_REFUTED',
        completed_steps=['5.1','5.3','5.4','5.5'],cancelled_steps=['5.2'],steps_not_run=['6'],
        gate='AWAITING_STEP_6_APPROVAL',input_sha256=before_hash,command=base.COMMAND,
        bootstrap=dict(method='ordinary paired event resampling within each cell; independent cells',
                       replicates=BOOTSTRAP_REPLICATES,seed=BOOTSTRAP_SEED),
        report=str(base.REPORT_PATH),report_sha256=base.sha256(base.REPORT_PATH),
        analysis_script_sha256=base.sha256(__file__),
        between_slope_interpretation='descriptive fit to same cell means; reabsorption is not explanation')
    (OUT/'analysis_summary.json').write_text(json.dumps(meta,indent=2)+'\n')
    (OUT/'profile_refits.meta.json').write_text(json.dumps({**meta,'bins':PROFILE_BINS,'variants':VARIANTS,
        'pol2_beta':'derivative at each empirical cell mean Npe; local Taylor convention',
        'quantiles':'40 quantile bins, unique edges, extrema padded by 0.5 pe; fits use TProfile bin centers',
        'root':'profile_refits.root','csv':'profile_refits.csv'},indent=2)+'\n')
    print(json.dumps(meta,indent=2),flush=True)


def make_figures(frames):
    sig=frames['residual_significance']
    fig,axes=plt.subplots(1,3,figsize=(14,4),sharey=True)
    for ax,material in zip(axes,base.MATERIALS):
        d=sig[sig.material==material]
        ax.errorbar(d.x_mm,d.residual_ps,yerr=d.paired_bootstrap_se_ps,fmt='o-',capsize=3)
        ax.axhline(0,color='black',lw=.8); ax.set_title(material); ax.set_xlabel('x [mm]'); ax.grid(alpha=.2)
    axes[0].set_ylabel('Descriptive between-fit residual [ps]')
    base.save_bundle('residual_mirrors',sig,{'error':'paired cell-stratified event bootstrap including center and fitted between slope',
        'bootstrap_seed':BOOTSTRAP_SEED,'bootstrap_replicates':BOOTSTRAP_REPLICATES,'causal':False},fig)
    fr=frames['chain_sensitivity']
    fig,axes=plt.subplots(1,3,figsize=(15,4),sharey=True)
    for ax,material in zip(axes,base.MATERIALS):
        for v,g in fr[fr.material==material].groupby('variant',sort=False):
            ax.plot(g.abs_x_mm,g.chain_residual_ps,'o-',label=v)
        ax.set_title(material); ax.set_xlabel('|x| [mm]'); ax.grid(alpha=.2)
    axes[0].set_ylabel('Within-chain residual [ps]'); axes[-1].legend(fontsize=7)
    base.save_bundle('chain_sensitivity',fr,{'variants':list(VARIANTS)+['event_ols_total','event_ols_two_counts'],
        'pol2_slope':'derivative at cell mean','warning':'Not the between-fit residual; chain identification remains refuted'},fig)
    fr=frames['mixture_target']
    fig,axes=plt.subplots(1,3,figsize=(15,4),sharey=True)
    for ax,material in zip(axes,base.MATERIALS):
        d=fr[fr.material==material]
        for key,err,label in [('descriptive_residual_ps','residual_bootstrap_se_ps','target residual'),
                             ('mixing_shift_ps','mixing_se_ps','raw mixing term'),
                             ('projected_mixing_ps','projected_mixing_se_ps','mixing after same between projection')]:
            ax.errorbar(d.abs_x_mm,d[key],yerr=d[err],fmt='o-',label=label,capsize=2)
        ax.axhline(0,color='black',lw=.8); ax.set_title(material); ax.set_xlabel('|x| [mm]'); ax.grid(alpha=.2)
    axes[0].set_ylabel('Shift [ps]'); axes[-1].legend(fontsize=7)
    base.save_bundle('mixture_target',fr,{'identity':'delta mu = fbar delta muC +(1-fbar)delta muS +(muCbar-muSbar)delta f',
        'interpretation':'exact finite-interval descriptive decomposition, not an independent causal intervention',
        'bootstrap_seed':BOOTSTRAP_SEED,'bootstrap_replicates':BOOTSTRAP_REPLICATES},fig)


def render_report(fr,summary,points,stats):
    # Implemented below; kept in one module so the command is fully reproducible.
    write_report(fr,summary,points,stats)


def write_report(fr,summary,points,stats):
    within=summary.beta_within_ns_per_pe.to_numpy()*NS_TO_PS
    se=summary.beta_within_error_ns_per_pe.to_numpy()*NS_TO_PS
    common=np.sum(within/se**2)/np.sum(1/se**2)
    common_se=1/np.sqrt(np.sum(1/se**2))
    hetero=np.sum(((within-common)/se)**2)
    secant=(even(stats['EJ-200']['y'])[-1]-even(stats['EJ-200']['y'])[0])/(even(stats['EJ-200']['n'])[-1]-even(stats['EJ-200']['n'])[0])*NS_TO_PS
    slopes=summary[['material']].copy()
    for key in ('beta_within','beta_between'):
        slopes[key+'_ps_pe']=summary[key+'_ns_per_pe']*NS_TO_PS
        slopes[key+'_se_ps_pe']=summary[key+'_error_ns_per_pe']*NS_TO_PS
    slopes['within_z']=abs(within)/se
    target500=fr['mixture_target'].query('abs_x_mm==500')
    common_r=np.sum(target500.descriptive_residual_ps/target500.residual_bootstrap_se_ps**2)/np.sum(1/target500.residual_bootstrap_se_ps**2)
    common_r_se=1/np.sqrt(np.sum(1/target500.residual_bootstrap_se_ps**2))
    hetero_r=np.sum(((target500.descriptive_residual_ps-common_r)/target500.residual_bootstrap_se_ps)**2)
    lines=['# EXEC_46 Step 5 — revised chain-rule identification and localized residual', '',
        'Date: 2026-09-16. Revision E1–E5 supersedes the Step 5 conclusion in commit `185a916`.', '',
        '**CHAIN_RULE_WITHIN_SLOPE_REFUTED. Steps 5.3–5.5 completed; Step 5.2 cancelled; STOP BEFORE STEP 6.**', '',
        'No simulation, production rerun, source-tree regeneration, push, merge, or deck edit. '
        'Input: 210,000 existing derived events, 21 cells, 10,000 per cell. The input SHA-256 is checked unchanged at exit.', '',
        'The localized positive structure is present separately in both mirrors for all materials: '
        'at ±500 mm its significance is 7.49–10.53 cell SEM, or 5.74–8.18 paired-bootstrap SE after refitting. '
        'Its even amplitudes are 6.368 ± 0.900, 7.858 ± 0.920, and 7.338 ± 0.859 ps. '
        'The same-projection source-mixture terms are +4.157 ± 0.298, +3.767 ± 0.263, and +2.800 ± 0.216 ps: '
        'a partial descriptive allocation, leaving +2.212 ± 0.864, +4.091 ± 0.877, and +4.538 ± 0.827 ps. '
        'The mechanism is not fully closed.', '',
        '## E1–E2: correction of the circular claim', '',
        f'The EJ-200 center-to-|x|=650 secant is {secant:.8f} ps/pe, versus the '
        f'seven-mean OLS slope {slopes.beta_between_ps_pe.iloc[0]:.8f} ps/pe. '
        'Regressing the outcome cell means against their Npe means fits the position curve that was to be explained. '
        'The former 96.997%, 97.653%, and 97.294% reductions are **reabsorption into a parameter fitted to the same data**, '
        'not explained fractions, causal evidence, or proof of an artifact. The former majority-artifact status is withdrawn.', '',
        'Step 5.2 is cancelled: with cell-centered Npe, alpha_i equals the observed mean T0_i; '
        'a global reference simply rewrites the original residual. It is not an independent observable.', '',
        'Define the descriptive target r_i = mean(T0)_i − mean(T0)_0 − beta_between[mean(Npe)_i − mean(Npe)_0]. '
        'Its even component is the average of the two mirrored residuals. No physical correction is implied.', '',
        '## E5: conditional slopes across materials', '',
        table(slopes,{'material':'Material','beta_within_ps_pe':'beta within [ps/pe]','beta_within_se_ps_pe':'HC1 SE',
                      'within_z':'within / SE magnitude','beta_between_ps_pe':'beta between [ps/pe]','beta_between_se_ps_pe':'OLS SE'},6), '',
        f'The inverse-variance common within slope is {common:.6f} ± {common_se:.6f} ps/pe. '
        f'Exact equality gives chi2/ndf = {hetero:.3f}/2, p = {math.exp(-hetero/2):.4g}. '
        f'The range relative to the common magnitude is {np.ptp(within)/abs(common)*100:.2f}%. '
        'The conditional response is nearly common, but exact equality is mildly disfavored (p=0.032); do not claim exact universality. '
        'Errors are event-level HC1 fixed-effect slope errors; different material simulations are treated as independent.', '',
        f'The ratio of the EJ-200 to EJ-230 between-slope central magnitudes is '
        f'{abs(slopes.beta_between_ps_pe.iloc[0]/slopes.beta_between_ps_pe.iloc[2]):.2f}; '
        'the EJ-230 denominator is compatible with zero, so “factor 30” is a central-value description, not a precisely measured physical factor. '
        'The between fit spans positions and source populations; the within fit conditions on position. '
        'Their discrepancy rejects using the measured within slope as the between-position response; '
        'it does not independently identify the mechanism of the discrepancy.', '',
        '## E3: the localized target, both mirrors and uncertainty', '',
        'The following reports r/SEM exactly as requested and a separate r/bootstrap-SE. '
        f'The latter uses {BOOTSTRAP_REPLICATES} ordinary paired-event bootstrap replicates per cell, seed {BOOTSTRAP_SEED}, '
        'recomputing the center, seven-point between slope, and residual in each replicate. '
        'Both ENDs and all per-event counts remain paired; cells are resampled independently. '
        'These z values are descriptive standardized effects, not independent hypothesis tests after model selection. '
        'The bootstrap does not cover geometry/model systematics or between-fit misspecification.', '',
        table(fr['residual_significance'].query('abs(x_mm)==500'),{
            'material':'Material','x_mm':'x [mm]','residual_ps':'r [ps]','cell_sem_ps':'cell SEM [ps]',
            'residual_over_cell_sem':'r/SEM','paired_bootstrap_se_ps':'paired SE [ps]',
            'paired_bootstrap_z':'r/paired SE','ci025_ps':'95% lower [ps]','ci975_ps':'95% upper [ps]'}), '',
        'All 21 signed residuals and their uncertainties are in `residual_significance.csv`; '
        'the plot `residual_mirrors.pdf` retains both mirrors.', '',
        table(fr['mixture_target'].query('abs_x_mm>0'),{'material':'Material','abs_x_mm':'|x| [mm]',
            'descriptive_residual_ps':'even r [ps]','residual_bootstrap_se_ps':'paired SE [ps]'}), '',
        f'Across materials the common even amplitude is {common_r:.3f} ± {common_r_se:.3f} ps '
        f'(chi2/ndf={hetero_r:.3f}/2, p={math.exp(-hetero_r/2):.3f}). '
        'This supports consistency of the localized amplitude across these three materials, not a universal law beyond the sampled grid.', '',
        'The retained quadratic summaries of this descriptive curve are:', '',
        table(summary.assign(a2_ps=summary.descriptive_residual_a2_ns_per_m2*NS_TO_PS, a2_se_ps=summary.descriptive_residual_a2_total_error_ns_per_m2*NS_TO_PS),{
            'material':'Material','a2_ps':'a2 [ps/m2]','a2_se_ps':'paired a2 SE','descriptive_residual_chi2_ndf':'diagonal chi2/ndf'}), '',
        'These chi2/ndf values use the historical cell-SEM weighting (7 points, 5 nominal degrees of freedom); '
        'they are lack-of-shape diagnostics, not calibrated goodness-of-fit probabilities: center subtraction and '
        'the fitted between slope correlate residuals. The high values preclude interpreting a2 as an adequate shape description.', '',
        '## 5.3: profile model and binning sensitivity', '',
        'All 21 original uniform 40-bin TProfile pol1 slopes and their integrated chain reproduce Step 2. '
        'The four variants form a 2×2 comparison: uniform/quantile binning × pol1/pol2. '
        'Quantiles use 40 equal-population nominal bins (duplicate boundaries removed), extrema padded by 0.5 PE, '
        'and standard ROOT TProfile bin-center abscissae and errors. For pol2, beta_i = p1 + 2*p2*mean(Npe)_i. '
        'This local-derivative convention is an explicit analysis assumption; it is not a nonlinear transport law. '
        'All coefficients, derivative errors, chi2/ndf, fit status and covariance status are in `profile_refits.csv`; '
        'the actual profiles, fitted functions and covariance matrices are saved in `profile_refits.root`.', '',
        table(fr['profile_refits'].groupby('variant',sort=False).agg(unreliable=('unreliable','sum'),
              median_chi2_ndf=('chi2_ndf','median'),max_chi2_ndf=('chi2_ndf','max')).reset_index(),
              {'variant':'Variant','unreliable':'chi2/ndf > 5 (of 21)','median_chi2_ndf':'median chi2/ndf','max_chi2_ndf':'maximum chi2/ndf'}), '',
        'Pol2 reduces failures from 12/21 to 1/21 with uniform bins; quantile pol2 has 0/21 failures. '
        'The large coefficient changes therefore cannot be dismissed as numerical fit failure alone. '
        'Profile specification materially affects the within-chain remnant.', '',
        'Integrate the even measured count profile trapezoidally with these local slopes; no attenuation-model refit. '
        'The residual is observed even shift minus that chain. It is distinct from the descriptive between-fit target. '
        'Changing local beta cannot alter r_between by definition, so merely showing the original 6–8 ps residual again '
        'would be a tautological robustness test.', '',
        table(fr['chain_sensitivity'].query('abs_x_mm==500'),{'material':'Material','variant':'Variant',
            'chain_residual_ps':'chain residual at 500 [ps]','shift_from_original_ps':'change vs original [ps]'}), '',
        'To test localized shape in the sensitivity curves, additionally define B = r(500) − '
        '[(1−w)r(200)+w*r(650)], w=[N(500)−N(200)]/[N(650)−N(200)]. '
        'This declared diagnostic removes a broad trend linear in the measured Npe profile. '
        'For r_between, B is exactly the same as for the raw mean curve: a constant beta*N term cancels. '
        'Compute the chain variant independently on each mirror (center → 200 → 500 → 650) and on the even curve. '
        'This is a shape sensitivity, not a new causal observable or an interpolated physical boundary.', '',
        table(fr['localization_sensitivity'],{'material':'Material','variant':'Variant','mirror':'mirror',
            'target_local_excess_ps':'target B [ps]','target_local_excess_bootstrap_se_ps':'target B SE','chain_local_excess_ps':'chain B [ps]'}), '',
        'The localized excess B remains positive separately on both mirrors for all four profile variants: '
        '12.29–28.75 ps. However its magnitude changes strongly; this supports persistence of the localized shape, '
        'not invariance of the original 6–8 ps amplitude under a physically identified correction. '
        'Even-chain a2 moves from 123.10/227.05/197.23 to 238.29/336.23/384.33 ps/m2 for quantile pol2.', '',
        'Profile sensitivity values are central estimates; their formal fit errors are retained, but '
        'no independent confidence claim is made from a failed profile fit or from these correlated post-fit shape differences.', '',
        '## 5.4: exact source-mixture finite differences', '',
        'For each END use the creator source of the actual winning photon: source_type=2 is Cherenkov, '
        'source_type=1 is scintillation. mu_C and mu_S are conditional means of the winning END timestamp, '
        'not means of source-specific minima in all events. All winner classes and Npe source counts close exactly. '
        'The interval identity used on adjacent signed positions is:', '',
        '`Delta mu = f_bar Delta mu_C + (1−f_bar) Delta mu_S + (mu_C_bar−mu_S_bar) Delta f`', '',
        'Dividing each term by the measured interval length gives finite-difference derivatives; '
        'integrating from x=0 and averaging the ENDs gives T0 terms with exact numerical closure. '
        'The interval-midpoint convention is declared and avoids an unreported product-rule discretization error. '
        'The three terms and derivative units ps/m, separately by END and signed interval, are in `mixture_intervals.csv`. '
        'All signed integrated terms, left/right contributions and paired bootstrap errors are in `mixture_components.csv`. '
        'No inference of a pointwise derivative inside the unmeasured intervals is possible.', '',
        table(fr['mixture_target'].query('abs_x_mm>0'),{'material':'Material','abs_x_mm':'|x| [mm]',
            'cherenkov_shift_ps':'f dmu_C [ps]','scintillation_shift_ps':'(1−f) dmu_S [ps]',
            'mixing_shift_ps':'mix term [ps]','mixing_se_ps':'mix SE [ps]','observed_shift_ps':'sum / observed [ps]'}), '',
        'At |x|=650 the raw mixing term is −11.808/−10.043/−7.330 ps, opposite in sign to the '
        'registered within-chain residual. At |x|=500 it is only −0.588/−0.369/−0.148 ps. '
        'It therefore does not explain the original positive remnant; the magnitude near −9 ps occurs at 650 mm.', '',
        table(fr['mixture_target'].query('abs_x_mm>0'),{'material':'Material','abs_x_mm':'|x| [mm]',
            'original_chain_residual_ps':'original chain residual [ps]','mixing_shift_ps':'raw mix [ps]',
            'original_chain_minus_raw_mix_ps':'original minus raw mix [ps]'}), '',
        'A raw mixing shift and a detrended target have different definitions. For a like-for-like comparison apply '
        'the SAME operator R_N(v)=v−v_0−OLS_slope(v,N)*(N−N_0) to the integrated mixing curve. '
        'For fixed measured N this is linear, so R_N(T0)=R_N(mixing)+R_N(T0−mixing) exactly. '
        'This is a descriptive allocation of a measured curve; it does not cure circularity or identify a counterfactual effect.', '',
        table(fr['mixture_target'].query('abs_x_mm==500'),{'material':'Material',
            'descriptive_residual_ps':'target [ps]','mixing_shift_ps':'raw mix [ps]',
            'projected_mixing_ps':'same-projection mix [ps]','projected_mixing_se_ps':'projected mix SE [ps]',
            'projected_mixing_to_residual_ratio':'signed fraction of target','projected_fraction_bootstrap_se':'fraction SE',
            'residual_minus_projected_mix_ps':'remaining [ps]',
            'residual_minus_projected_mix_se_ps':'remaining SE [ps]'}), '',
        'Under this identical projection the mixture receives 65.3%, 47.9%, and 38.2% of the 500 mm target '
        '(signed descriptive fractions, not independent explained fractions). The remaining effect is nonzero, '
        'especially for EJ-204/230. The large negative mixing shift at 650 changes the fitted broad trend, '
        'which is why its projected contribution at 500 becomes positive despite a small negative raw value.', '',
        'The fitted cell fraction f_C is fixed as a population summary at each x; the event-level winning source is random '
        'and may correlate with photon counts. Therefore it is not justified to assert that a within-cell beta is '
        'necessarily blind to every mixture effect. The exact identity above measures the mixture term without that assertion.', '',
        '## 5.5: separate scintillation and Cherenkov counts', '',
        'Fit event-level T0 = intercept + beta_S*Nscint_END + beta_C*NCher_END independently in each cell, '
        'where each count sums the two ENDs. Use event OLS with HC1 covariance and retain a one-count event OLS control '
        'to separate count splitting from the original binned-profile weighting. '
        'Average coefficients/counts over mirrors then sum the two trapezoidal beta*dN chains. '
        'All 21 coefficient pairs, covariance, count correlations and conditioning are in `two_count_slopes.csv`. '
        'This remains a within-to-between transfer test; a second count does not automatically resolve position confounding.', '',
        table(fr['chain_sensitivity'].query('abs_x_mm>0 and variant in ["uniform_pol1","event_ols_total","event_ols_two_counts"]'),
            {'material':'Material','variant':'Variant','abs_x_mm':'|x| [mm]','chain_residual_ps':'residual [ps]',
             'chain_residual_bootstrap_se_ps':'paired SE [ps]'}), '',
        'Splitting the counts changes event-OLS a2 by only −4.54/−6.75/−6.15 ps/m2 relative to the matched '
        'one-count event OLS control; the residual remains +170.61/+255.81/+298.94 ps/m2. '
        'The positive localized B remains 18.20/23.16/23.97 ps. This two-count linear transfer does not close the target.', '',
        'The complete quadratic comparison (residual fit, not separate fitted-curve subtraction) is:', '',
        table(fr['chain_fit_summary'],{'material':'Material','variant':'Variant','a2_ps_m2':'a2 [ps/m2]',
            'a2_bootstrap_se_ps_m2':'paired a2 SE','chi2_ndf':'diagonal chi2/ndf',
            'localized_excess_500_ps':'B at 500 [ps]'}), '',
        'Quadratic chi2/ndf uses historical cell SEM weights and 5 nominal degrees of freedom for comparison; '
        'the integrated predictions also have uncertainty and shared fitted coefficients, so these are shape diagnostics. '
        'Bootstrap errors are supplied for event OLS variants; profile variants retain conditional formal errors in the sidecar. '
        'No poor quadratic or poor profile fit is used as a mechanism measurement.', '',
        '## Reproducibility and stopping point', '', '```bash',base.COMMAND,'```','',
        f'Input SHA-256: `{base.sha256(base.DERIVED_ROOT)}`. Bootstrap seed: {BOOTSTRAP_SEED}; replicates: {BOOTSTRAP_REPLICATES}. '
        'ROOT version: '+ROOT.gROOT.GetVersion()+'. Input opened read-only. '
        'Every produced PDF has `.root`, `.csv`, `.meta.json` companions. '
        '`step5_revision_tables.root` also preserves all table columns, including categorical labels. '
        'Rollback tag before edits: `pre-exec46-step5-e1-e5-20260916`.', '',
        '**Step 6 has not run and still requires René’s explicit approval. No push was performed.**','']
    base.REPORT_PATH.write_text('\n'.join(lines))
