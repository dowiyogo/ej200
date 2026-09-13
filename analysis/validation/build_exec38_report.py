#!/usr/bin/env python3
"""Publish EXEC38 tables/figures/contract from saved analysis; never transports photons."""
import csv
import json
from pathlib import Path
import math
import collections
import uproot
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_exec38 import HERE, REPO, BASE, BULK, YIELD, VREF, sha, now, dump, sidecars

OUT=BASE/'exec38_20260913'
REPORT=BASE/'REPORT_EXEC38_20260913.md'
GOLDEN=HERE/'GOLDEN_REFERENCE_20260913.json'
ELJEN='https://eljentechnology.com/products/plastic-scintillators/ej-200-ej-204-ej-208-ej-212'
EJ230='https://eljentechnology.com/images/products/data_sheets/EJ-228_EJ-230.pdf'


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+
                     ['| '+' | '.join(str(v) for v in row)+' |' for row in rows])+'\n'


def pm(r,k,s=None,digits=3):
    return f"{r[k]:.{digits}f} ± {r[s or k+'_se']:.{digits}f}"


def historical(configs):
    refs=[BASE/'REPORT_EXEC26_PHASE1B_20260909.md',BASE/'REPORT_EXEC29_20260910.md',
          BASE/'REPORT_EXEC30_20260910.md',BASE/'REPORT_EXEC33_20260911.md',
          BASE/'REPORT_EXEC34C_20260912.md',BASE/'REPORT_EXEC36_20260913.md']
    b30=BASE/'ej200_exec26_20260909/build_exec30_20260910'
    m=json.loads((b30/'manifest.meta.json').read_text())['cells']['V1']
    cfg=dict(cell_id='EXEC30_V1',material=m['material'],x_mm=m['x_mm'],N_generated=m['N'],
             seeds=m['seeds'],simulation_commit=m['source_commit'],workers=m['threads'],eventModulo=m['eventModulo'],
             configuration=m['readout'],N_TOP=m['N_TOP'],command=m['command'])
    configs['EXEC30_V1']=cfg
    p26=BASE/'ej200_exec26_20260909/build_exec26_phase1b_20260909/run500_ready/boundary_census_run0.meta.json'
    m26=json.loads(p26.read_text())
    configs['EXEC26_R1']=dict(cell_id='EXEC26_R1',material=m26['material'],x_mm=0,N_generated=m26['N'],
        seeds=m26['seeds'],simulation_commit=m26['source_commit'],workers=m26['threads'],eventModulo=None,
        configuration=m26['readout'],N_TOP=m26['N_TOP'],command=m26['execution_command'])
    for cell in ['S1','S2','S3','S4']:
        p=BASE/'exec33_20260911'/cell/'invocation.meta.json';m=json.loads(p.read_text());refs.append(p)
        configs['EXEC33_'+cell]=dict(cell_id='EXEC33_'+cell,material=m['material'],x_mm=m['x_mm'],N_generated=m['N'],
            seeds=m['seeds'],simulation_commit=m['source_commit'],workers=m['workers'],eventModulo=m['eventModulo'],
            configuration=m['mode'],N_TOP=m['N_TOP'],command=m['command'])
    reflection_path=b30/'cells/V1/reflection_panels.csv'
    reflection=next(r for r in csv.DictReader(reflection_path.open()) if r['panel']=='ALL')
    mirror_path=BASE/'exec36_20260913/symmetry.csv';mirrors=list(csv.DictReader(mirror_path.open()))
    reproducibility=BASE/'exec33_20260911/event_reproducibility.json'
    refs += [p26,b30/'manifest.meta.json',reflection_path,mirror_path,reproducibility]
    rows=[
        dict(id='V6_balance',prediction='N_source equals the sum of disjoint terminal physical destinations',
             measured=0.0,unit='percent residual',uncertainty=0.0,status='ARCHIVED_PASS',config_ids=['D0','EXEC30_V1'],
             provenance='EXEC29 report:71 and EXEC30 report:54; D0 scintillation 40219965 produced and terminal, optical 42399341 produced and terminal.',
             claimed_source='EXEC34C',audit='NOT substantiated in EXEC34C: its diagnostics-OFF hit-only files reconcile detections, not a complete terminal photon budget.'),
        dict(id='V6_reflectivity',prediction='Per-encounter optical return probability 0.98 at the physical ESR reflector',
             measured=float(reflection['reflection_pct']),unit='percent',uncertainty=None,status='ARCHIVED_CONSISTENT',
             config_ids=['D0','EXEC30_V1'],provenance='EXEC30 report:219–229, copied ALL row; same value documented in EXEC29:67.',
             archived_row=reflection,audit='98.00010377747189%; not the D3 original dielectric result. No uncertainty was published; no new fit or recount.'),
        dict(id='V6_no_TIR_up_index',prediction='Zero critical-angle total internal reflection when n2>n1',
             measured=0,unit='encounters',uncertainty=0,status='ARCHIVED_PASS',config_ids=['EXEC26_R1'],
             provenance='EXEC26 Phase1b report:47–59, zero over 3623938 eligible forward air-to-wrap encounters.',
             audit='Historical acquisition used an engine-specific classifier; future contract requires physical incidence/refraction geometry, not that classifier.'),
        dict(id='V6_symmetry',prediction='Equal mirrored physical distributions for mirrored configurations',
             measured=2.6638244072,unit='maximum absolute relative core-width difference percent',uncertainty=None,
             status='DOCUMENTED_NO_ARCHIVED_GATE',config_ids=[k for k in configs if k.startswith('EJ')],
             provenance='EXEC36 report:328–375 and symmetry.csv; value transcribed from published maximum, no width recalculation.',
             audit='0.5%/0.3% are not universal archived results. EJ204 x=-650 V2 versus right x=+650 is -0.290%; EJ204 x=0 V2 is +0.556% (specific comparisons). Archived report explicitly introduced no equality threshold.'),
        dict(id='V6_workers',prediction='Identical event-wise detected destination counts for the specified common-RNG CPU worker test',
             measured=0,unit='different events and delta Npe/end',uncertainty=0,status='ARCHIVED_PASS',
             config_ids=['EXEC33_'+c for c in ['S1','S2','S3','S4']],
             provenance='EXEC33 report:21–43; event_reproducibility.json copied unchanged.',
             audit='N=2000, EJ204 x=0, EndTop N_TOP=70; tested 1,4,12,24 workers. Not a 21-cell or cross-RNG GPU identity claim.')]
    sidecars(OUT,'V6_archived_evidence',rows,dict(source_hashes={str(p):sha(p) for p in refs},
             configurations=configs,method='Published results transcribed; no old physical statistic recalculated.'))
    # Copy the source rows losslessly into a dedicated provenance sidecar.
    dump(OUT/'V6_source_rows.json',dict(reflection=reflection,mirror_rows=mirrors,
          worker_rows=json.loads(reproducibility.read_text()),source_hashes={str(p):sha(p) for p in refs}))
    return rows


def legacy_incidence_audit():
    """Explain why legacy classifier totals are not an independent V5 stream.

    These rows are provenance diagnostics only, never acceptance observables.
    In particular, do not silently equate the archived classifier's Detection
    tally with the independently exported detected-hit count.
    """
    rows=[]
    for cell,task in [('D0','29'),('D3','30')]:
        base=BASE/'ej200_exec26_20260909'/('build_exec'+task+'_20260910')/'cells'/cell
        census=base/'boundary_census_run0.csv';raw=list(csv.DictReader(census.open()))
        timing=BASE/'exec37_20260913/cells'/cell/'results.root'
        with uproot.open(timing) as f: counts=f['event_yields'].arrays(library='np')
        for end in ['left','right']:
            tally=collections.Counter()
            for p in raw:
                if p['pre']=='BarPV' and p['post']=='EndSiPM'+end.capitalize()+'_PV':
                    tally[p['status_name']]+=int(p['count'])
            rows.append(dict(cell_id=cell,end=end,legacy_absorption_tally=tally['Absorption'],
                legacy_detection_tally=tally['Detection'],legacy_sum=tally['Absorption']+tally['Detection'],
                exported_detected_hits=int(counts['npe_'+end].sum()),census=str(census),census_sha256=sha(census),
                scope='Schema audit only; internal classifier counts are not a V5 acceptance criterion or an independently identified incident stream.'))
    sidecars(OUT,'V5_legacy_schema_audit',rows,dict(method='Read-only provenance diagnostic; no closure PASS/FAIL inferred from internal status labels.'))
    return rows


def figures(r,points,meta):
    for test in ['V3','V4']:
        fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
        for ax,mat in zip(axes,BULK):
            for end,color in [('left','tab:blue'),('right','tab:orange')]:
                a=sorted([p for p in points if p['material']==mat and p['end']==end],key=lambda p:float(p['distance_mm']))
                d=np.array([float(p['distance_mm']) for p in a]);fit=next(p for p in r[test] if p['material']==mat and p['end']==end)
                ykey,ekey=('npe','npe_sem') if test=='V3' else ('mean_arrival_ns','mean_arrival_se_ns')
                ax.errorbar(d,[float(p[ykey]) for p in a],yerr=[float(p[ekey]) for p in a],fmt='o',ms=4,color=color,label=end)
                grid=np.linspace(50,1350,200)
                if test=='V3':
                    ax.plot(grid,fit['amplitude']*np.exp(-grid/fit['lambda_mm']),color=color,ls='--',alpha=.75)
                    p=fit['double_parameters'];ax.plot(grid,p[0]*np.exp(-grid/p[1])+p[2]*np.exp(-grid/p[3]),color=color,alpha=.65)
                    ax.set_yscale('log');ax.set_ylabel('Detected photons / generated event')
                else:
                    ax.plot(grid,fit['intercept_ns']+fit['slope_ns_per_mm']*grid,color=color,ls='--');ax.set_ylabel('Mean photon arrival [ns]')
                ax.set_title(mat);ax.set_xlabel('Distance to observed end [mm]');ax.grid(alpha=.2)
            if test=='V3':
                vals=[float(p['npe']) for p in points if p['material']==mat]
                ax.set_ylim(min(vals)*.7,max(vals)*1.4)
            ax.legend()
        fig.suptitle('EXEC38: '+('single (dashed) and double (solid) exponentials; both rejected' if test=='V3' else 'global linear fits (dashed) rejected; speed-bound checks pass'))
        stem=test+'_profiles';fig.savefig(OUT/(stem+'.png'),dpi=180);fig.savefig(OUT/(stem+'.pdf'));plt.close(fig)
        sidecars(OUT,stem,points,dict(meta,figure_files={ext:sha(OUT/(stem+ext)) for ext in ['.png','.pdf']},
                                      fit_results=r[test],statistical_only=True))
    prim=[p for p in r['V7'] if p['primary']]
    fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for ax,mat in zip(axes,BULK):
        for end,color in [('left','tab:blue'),('right','tab:orange')]:
            a=sorted([p for p in prim if p['material']==mat and p['end']==end],key=lambda p:p['x_mm'])
            ax.errorbar([p['x_mm'] for p in a],[p['rho'] for p in a],yerr=[p['rho_se'] for p in a],fmt='o-',color=color,label='EndTop '+end)
            if mat=='EJ-204':
                for cell,offset in [('D0',-25),('D3',25)]:
                    p=next(p for p in prim if p['cell_id']==cell and p['end']==end)
                    ax.errorbar([offset],[p['rho']],yerr=[p['rho_se']],fmt='s' if cell=='D0' else '^',color=color,label=cell+' '+end)
        ax.axhline(0,color='black',lw=1);ax.set_xlabel('Gun x [mm]');ax.set_ylabel('Same-end Pearson rho');ax.set_title(mat);ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('EXEC38 V7: fixed V2 groups; errors are event-bootstrap SE; D0/D3 markers offset for visibility')
    stem='V7_rho';fig.savefig(OUT/(stem+'.png'),dpi=180);fig.savefig(OUT/(stem+'.pdf'));plt.close(fig)
    sidecars(OUT,stem,prim,dict(meta,figure_files={ext:sha(OUT/(stem+ext)) for ext in ['.png','.pdf']},
        plotting_offset_only='D0 at -25 and D3 at +25 mm for legibility; both physical x=0, preserved in sidecars'))


def main():
    r=json.loads((OUT/'results.json').read_text());meta=json.loads((OUT/'analysis.meta.json').read_text())
    configs=json.loads((OUT/'configurations.json').read_text());v6=historical(configs)
    legacy=legacy_incidence_audit()
    points=list(csv.DictReader((OUT/'profile_points.csv').open()))
    for row in points:
        for key in ['x_mm','distance_mm','npe','npe_sem','mean_arrival_ns','mean_arrival_se_ns']:
            row[key]=float(row[key])
        row['N_generated']=int(row['N_generated'])
    figures(r,points,meta)
    uncertainty='Per-row *_se fields, 500 generated-event bootstrap replicates, seed 38091301; statistical only.'
    grid=[k for k in configs if k.startswith('EJ')];allcells=grid+['D0','D3']
    tests=[
        dict(id='V1',observable='Scintillation source photons and deposited energy per generated event',
            prediction=dict(expression='sum(N_scint)/(yield*sum(dE_MeV)) = 1',yield_photons_per_MeV=YIELD,source=[ELJEN,EJ230]),
            measured=[dict(material=m,ratio=None,uncertainty=None,status='NOT_EVALUABLE') for m in BULK],uncertainty=None,
            tolerance=dict(rule='abs(ratio-1)<=3*event_paired_SE',justification='Tests specified yield without assuming unavailable datasheet systematic precision'),
            status='NOT_EVALUABLE',config_ids=allcells,required_observables=['event_id','deposited_energy_MeV','N_source_scintillation'],
            limitation='No deposited-energy ledger in existing ROOT/logs. Recorded total production cannot infer its own expected denominator.'),
        dict(id='V2',observable='First bar-air encounter outward escape fraction among scintillation photons; one encounter per photon',
            prediction=dict(expression='sin(theta)*cos(theta) flux measure; P(cone)=1/n^2; Fresnel escape band',n=1.58,
                cone_fraction=1/1.58**2,escape_interval=[.36,.40],source=['Snell law + flux angular integral',ELJEN]),
            measured=None,uncertainty=None,tolerance=dict(rule='0.36-3SE <= escape_fraction <= 0.40+3SE',
                justification='User analytical band plus sampling uncertainty; no extra percentage padding'),status='NOT_EVALUABLE',config_ids=allcells,
            required_observables=['event_id','photon_id','source_kind','surface_id','surface_encounter_ordinal','incident_direction',
                'surface_normal','incident_wavelength_nm','outgoing_destination','incident_n','transmitted_n'],
            applicability='Verify flux-weighted first-hit angular population and bar-air interface; localized isotropic emission alone does not establish it.'),
        dict(id='V3',observable='Mean detected photon counts per end per generated event versus distance to that end',
            prediction=dict(expression='lambda_eff < nominal attenuation; common nominal/lambda ratio under common trapped path factor',
                nominal_mm=BULK,source=[ELJEN,EJ230,'exp(-path_length/attenuation_length)']),measured=r['V3'],
            material_consistency=r['V3_consistency'],uncertainty=uncertainty,
            tolerance=dict(rule='lambda+3SE<nominal for PASS; lambda-3SE>nominal for FAIL; otherwise INDETERMINATE; model GOF p>=0.01; common-ratio GOF p>=0.01',
                model_selection='Two exponentials preferred only if p>=0.01 and AICc improvement>=6',
                justification='3SE separates sampling overlap; 1% descriptive GOF tests stated profile/path-factor approximations'),
            status='FAIL',subtests=dict(bound='PASS',single_model='FAIL',double_model='FAIL',common_material_factor='FAIL'),
            config_ids=grid,required_observables=['event_id','N_generated','gun_position_mm','N_detected_left','N_detected_right'],
            limitation='Descriptive lengths from rejected models are not calibrated physical attenuation measurements. Failure concerns those hypotheses, not proof of incorrect photon transport.'),
        dict(id='V4',observable='Photon-weighted mean detected arrival time relative to source time versus distance to each end',
            prediction=dict(expression='v_eff=1/(dt_mean/dd) < c/n under stable selected path population',c_mm_ns=299.792,n=1.58,
                c_over_n_mm_ns=VREF,source=['t=path_length/(c/n) for constant n',ELJEN]),measured=r['V4'],uncertainty=uncertainty,
            tolerance=dict(rule='Positive slope; v+3SE<c/n PASS; v-3SE>c/n FAIL; overlap INDETERMINATE; linear GOF p>=0.01',
                justification='3SE bound check plus separately reported whole-profile model validity'),
            status='FAIL',subtests=dict(speed_bound='PASS',global_linearity='FAIL'),config_ids=grid,
            required_observables=['event_id','gun_position_mm','source_time_ns','detected_destination','arrival_time_ns'],
            limitation='Rejected linearity prevents interpreting a single material speed; selected mean slope is not an individual causal speed.'),
        dict(id='V5',observable='Independently counted incident and detected photons at one specified end, with spectral accounting',
            prediction=dict(expression='integral S(E)*PDE(hc/E)dE / integral S(E)dE',values=r['V5_prediction'],
                source=[str(REPO/'data/sipm/AFBR-S4N66P024M_pde.txt'),'Configured SSLG4 emission tables, hashes in prediction rows','Bernoulli detection expectation']),
            measured=None,uncertainty=None,tolerance=dict(rule='abs(N_detected/N_incident-PDE_average)<=3*combined_statistical_SE',
                justification='Sampling uncertainty, no invented sensor systematic error'),status='NOT_EVALUABLE',config_ids=allcells,
            required_observables=['event_id','photon_id','detector_end','incident_wavelength_nm','detected_boolean','source_kind'],
            limitation='Grid logs only all-sensor incidence aggregate; D0/D3 have legacy per-end status tallies but no independent photon identities/incident spectra. Pure scintillation emission is not necessarily the incident optical mixture.'),
        dict(id='V6',observable='Physical conservation, reflector return probability, index-order optics, mirror symmetry and worker reproducibility',
            prediction=dict(expression='Exact conservation; ESR return probability .98; zero TIR for increasing index; reciprocal mirrored distributions; CPU common-RNG invariance',
                source=['Conservation of photons','Bernoulli reflection with specified R=.98','Snell law','Mirror symmetry','Deterministic common-RNG event mapping']),
            measured=v6,uncertainty='Copied published uncertainty or null when not published; no historical measurement recalculated.',
            tolerance=dict(conservation='0 missing or duplicate source photons',reflectivity='3*sqrt(.98*.02/N_incident), replace with event-cluster SE when available',
                TIR_increasing_index='0',mirror='abs(paired_difference)<=3*paired_event_bootstrap_SE',CPU_workers='0 event destination-count differences in specified common-RNG fixture',
                justification='Integer identities where exact; 3SE for stochastic surface returns and reciprocal distributions; historical gates are not rewritten'),
            status='PARTIALLY_VERIFIED',config_ids=allcells+['EXEC30_V1','EXEC26_R1']+['EXEC33_'+c for c in ['S1','S2','S3','S4']],
            limitation='EXEC34C full photon-balance attribution unverified. EXEC36 did not assert a universal .5%/.3% or preregister a symmetry gate. GPU need not duplicate CPU random numbers.'),
        dict(id='V7',observable='Event-paired same-end leading-edge group timestamps relative to a fixed source time; standard deviations ddof=1',
            prediction=dict(expression='rho>0; Var(T1-T2)=s1^2+s2^2-2*rho*s1*s2; Q=1/sqrt(1-rho) only if s1=s2',
                source=['Covariance identity and shared event fluctuations, fixed same-end grouping']),measured=r['V7'],uncertainty=uncertainty,
            tolerance=dict(rule='rho-3SE>0 PASS; rho+3SE<0 FAIL; otherwise INDETERMINATE; >=200 accepted events and >=95% finite replicas',
                equal_variance_approximation='abs(Q-F)<=3*SE(Q-F)',exact_variance_identity_relative=1e-10,
                justification='Strict positivity distinguished from compatibility with zero; identity tolerance is floating-point validation only'),
            status='PASS',config_ids=allcells,required_observables=['event_id','N_generated','source_time_ns','detected_sensor_global_id','arrival_time_ns'],
            timestamp_postprocessor=dict(rise_ns=.5,fall_ns=5.,threshold_PE=4.,primary_variant='V2',
                V2_groups=dict(left=[[0,1,2,3],[4,5,6,7]],right=[[8,9,10,11],[12,13,14,15]]),
                V1_secondary_groups=dict(left=[[0,2],[4,6]],right=[[8,10],[12,14]]),
                implementation='Preserved EXEC36/37 LeadingEdgeTime; no added noise, same deterministic pulse and crossing algorithm'),
            limitation='30/88 equal-variance approximations fail their separate comparison; all 88 rho positivity tests pass. Direct sigmas require jitter-free source time; quantile core widths are different observables.')]
    contract=dict(schema='ej200.physical_transport_acceptance.v1',created_utc=now(),
        title='EXEC38 golden reference: measured baseline and preregistered physical hypotheses',
        ready_for_acceptance=False,overall_status='INCOMPLETE_WITH_FAILED_PROFILE_HYPOTHESES',
        fail_closed_rule='Every required test must be evaluable and pass its applicable physical hypothesis; null or unresolved applicability cannot count as PASS.',
        interpretation='This freezes the attempted contract, not a certification of Geant4 or an assertion that these hypotheses universally test every geometry. Do not retune tolerances using this baseline.',
        preregistration=dict(path=str(HERE/'EXEC38_PREREGISTRATION.md'),sha256=sha(HERE/'EXEC38_PREREGISTRATION.md'),commit='4e08426'),
        analysis=dict(path=str(OUT),command=meta['command'],source_commit=meta['analysis_commit'],code_sha256=meta['code_sha256'],
                      bootstrap=meta['bootstrap']),configurations=configs,tests=tests,
        artifact_hashes={str(p):sha(p) for p in OUT.glob('*') if p.is_file()
            and p.name not in ['delivery.meta.json','verification.meta.json']
            and p.suffix in ['.json','.csv','.root']})
    dump(GOLDEN,contract)
    lines=[table(['Test','Prediction','Prediction source','Measured','PASS/FAIL'],[
        ['V1','Nscint/(yield × deposited energy)=1','[Eljen]('+ELJEN+'), [EJ230]('+EJ230+')','Production totals available; deposited energy absent','NOT EVALUABLE'],
        ['V2','First bar-air escape fraction 0.36–0.40','Snell law; specified flux angular distribution','No first-encounter photon identities or angles','NOT EVALUABLE'],
        ['V3 bound','Effective length below 3800/1600/1200 mm','Eljen; exponential path survival','Descriptive single lengths 422.25–601.67 mm','PASS bound only'],
        ['V3 models / material factor','Adequate exponential profile; common nominal/effective ratio','Stated exponential and common-path-factor hypotheses','Single χ²=24936.7–28249.3 (5 dof); double χ²=620.2–740.4 (3 dof); ratios ≈6.3,3.4,2.8','FAIL'],
        ['V4 bound','Effective slope speed below 189.742 mm/ns','c/n with n=1.58','146.655–154.916 mm/ns','PASS bound only'],
        ['V4 global linearity','One arrival-mean slope describes seven positions','Stated global linear model','χ²=47615.4–69135.8 (5 dof)','FAIL'],
        ['V5','Per-end detection/incident ratio equals spectral PDE','Fixed PDE table; Bernoulli expectation','Predicted 0.606655 / 0.619397 / 0.610613; no independent per-end incident stream','NOT EVALUABLE'],
        ['V6','Conservation, optical return, reciprocity and worker invariance','Physical identities and archived measurements','Several archived checks confirmed; balance attribution and universal symmetry percentages unsupported','PARTIALLY VERIFIED'],
        ['V7','Same-end Pearson rho strictly positive','Shared fluctuations; covariance identity','88/88 pass, including D0/D3; 30/88 equal-variance approximations fail separately','PASS rho']]),
        '# EXEC38 — optical transport validation and golden acceptance reference\n',
        'Date: 2026-09-13. **The available evidence does not certify the complete transport contract.** The attenuation and speed bounds pass as descriptive checks, but their whole-profile models fail decisively. V1, V2 and V5 remain unevaluable from the required physical streams. These outcomes are retained; a rejected profile approximation alone is not proof that individual photon transport is incorrect. V7 establishes positive same-end correlation and a bias in inferring direct group standard deviations from a normalized difference.\n',
        f'The machine-readable [golden reference]({GOLDEN}) is explicitly `ready_for_acceptance=false`. Missing inputs and failed hypotheses cannot be silently accepted by a future GPU engine. No new transport simulation, electronics injection, BLUE, push, merge or deck edit was performed. No numerical comparison against experimental data is made in this report.\n',
        '## Preregistration, inputs and uncertainty\n',
        f'[Preregistration]({HERE}/EXEC38_PREREGISTRATION.md), commit `4e08426`, preceded all EXEC38 calculations. Rollback tag `pre-exec38-20260913` points to `35d6bcc`; analysis code at `3197a55`. The preregistration openly records that prior published values were already known, so V6 is retrospective. Tolerances are 3 bootstrap standard errors for one-sided physical-bound/positive-correlation decisions, p≥0.01 for profile and common-ratio goodness of fit, and exact zero for applicable counting identities. No tolerances were broadened after results.\n',
        'The grid has 21 successful EndTop cells: EJ200/OPSC100, EJ204/OPSC101, EJ230/OPSC106; x=0, ±200, ±500, ±650 mm; N=10000; seeds 26092601 and 8349041; 4 workers, eventModulo=1, N_TOP=70, source `420addf0fd6029d5b2f0e235f472a8ae47f31fac`. D0/D3 are EJ204, END-only, x=0, N=2000 with the same seeds, sources `2ddf20592406058c0ad741bd7cdd0f2ce7780acf` and `391fa4e867717088296f25c24efe15f216637817`; 1 and 12 workers respectively. D0 eventModulo was not recorded; D3 is 1.\n',
        f'[Input chain]({OUT}/inputs.json) preserves each simulation command, ROOT hash and archived timestamp/cache hashes. Every consumed EXEC36/37 sidecar was rehashed; native grid files were reopened for schema/size/entry validation and their macro hashes checked. Full grid ROOT hashes are inherited through the previously verified exact END-cache chain, not claimed as newly rehashed here. D0/D3 raw ROOTs were rehashed now. All 42 grid end totals match their native recorded totals exactly. One additional EJ204 x=0 native-stream check compares event counts and time sums against the cache.\n',
        'All new uncertainties use 500 generated-event bootstrap resamples, NumPy PCG64 seed 38091301. The same event index resample is used across equal-N cells/ends/materials, preserving their observed paired covariance. The N=2000 and N=10000 ensembles are separate; no unmeasured cross-population covariance is claimed. Fits use full covariance matrices. All quoted errors are statistical, conditional on geometry, PDE, source and estimator. They exclude model inadequacy and datasheet systematic uncertainty.\n',
        '## V1 — production: required denominator absent\n',
        f'Eljen specifies yields 10000, 10400 and 9700 photons/MeV for EJ200, EJ204 and EJ230. [EJ200/EJ204]({ELJEN}); [EJ230 datasheet, August 2023]({EJ230}). Expected production is yield times the **recorded deposited energy**, not incoming gun energy or an assumed MIP energy loss.\n',
        table(['Material','Production ratio measured/expected','Reason'],[[m,'NOT EVALUABLE','No deposited-energy ledger'] for m in BULK]),
        'The native files contain only a detected-hit TTree. Master logs retain total generated scintillation photons but no deposited energy or per-event production ledger. Source/output inspection found no exported deposited-energy score. Recovering dE from the photon count would make V1 circular. All 23 recorded production totals/per-event means are preserved below and in the inventory sidecars; they do not pass the yield test.\n',
        table(['Cell','N events','Scintillation photons produced','Produced/event'],[[p['cell_id'],p['N_generated'],p['scintillation_photons_generated'],f"{p['scintillation_photons_per_event']:.4f}"] for p in r['inventory']]),
        '## V2 — first encounter and applicability\n',
        f'The stipulated flux integral is ∫₀^θc sinθ cosθ dθ / ∫₀^(π/2) sinθ cosθ dθ = sin²θc = 1/1.58² = {1/1.58**2:.8f}; θc={math.degrees(math.asin(1/1.58)):.4f} degrees. The prescribed Fresnel escape acceptance band is 0.36–0.40, expanded only by 3 event-bootstrap SE. This is restricted to scintillation photons at the **first bar-air surface encounter**. Coupled sensor surfaces have different index ordering.\n',
        'The existing detected-only schema has no photon identity, creator classification, surface-encounter ordinal, surface normal, or incidence/outgoing direction. D0/D3 terminal fates and cumulative boundary records cannot identify the first encounter. No aggregate escape percentage is relabelled as V2. An isotropic localized source does not itself guarantee the flux-weighted incidence population at its first hit on a finite bar; the stated angular assumption must be verified before treating the band as universal.\n',
        'For an engine-neutral future recording, export event/photon identity and source class, the first surface ID, surface normal, incoming wavelength/direction, the two refractive indices, and outgoing physical destination/direction. Classify outward refraction from those physical fields; no stepping callback name or status enum belongs to the acceptance API. This is a proposed output contract only; no instrumentation or simulation was added.\n',
        '## V3 — attenuation and profile adequacy\n',
        f'The nominal light attenuation lengths are 3800/1600/1200 mm. EJ230’s 120 cm is explicitly listed in the [manufacturer datasheet]({EJ230}); its prose separately describes a typical 100 cm emission mean free path. The table value requested here is 120 cm. Eljen calls the quantity **light attenuation length**; equating it with a microscopic bulk absorption length is an additional model assumption used by the simulation.\n',
        'For a specified end the independent variable is distance to that end: dL=700+x and dR=700−x mm. Folding to |x| would merge distinct near/far observations and destroy the per-end attenuation profile. All seven positions remain in each fit. Npe uses the full generated-event denominator, including zeros. Single fits are A exp(−d/λ); double fits have two positive amplitudes and two lengths. The registered multistarts are used for the central fit and its solution initializes each bootstrap refit.\n',
        table(['Material','End','λ single ± SE [mm]','Nominal/λ ± SE','χ² / dof','Bound','Model'],[[p['material'],p['end'],pm(p,'lambda_mm','lambda_se_mm'),pm(p,'nominal_over_effective'),f"{p['chi2']:.2f} / {p['dof']}",p['bound_status'],p['model_status']] for p in r['V3']]),
        table(['Material','End','λ short ± SE [mm]','λ long ± SE [mm]','Double χ² / dof','p','ΔAICc single−double'],[[p['material'],p['end'],f"{p['double_parameters'][1]:.3f} ± {p['double_parameter_se'][1]:.3f}",f"{p['double_parameters'][3]:.3f} ± {p['double_parameter_se'][3]:.3f}",f"{p['double_chi2']:.2f} / 3",f"{p['double_p_value']:.3g}",f"{p['AICc']-p['double_AICc']:.2f}"] for p in r['V3']]),
        'Two components improve AICc substantially but **still fail goodness of fit in every profile**. It would be incorrect to conclude that a second exponential is sufficient. All length bootstrap fits have finite successful replicas and no central double-fit bound is active; the failure is not explained away as a failed optimizer. The small errors are conditional descriptive-fit errors, not credible uncertainty on an adequate physical attenuation model. Single-fit p-values underflow numerically to zero; read them as extremely small, not an exact mathematical probability of zero.\n',
        table(['End','Common-ratio χ² / dof','Consistency'],[[p['end'],f"{p['chi2']:.2f} / 2",p['status']] for p in r['V3_consistency']]),
        'The common zigzag-factor hypothesis is not supported by these descriptive ratios. Different spectra, distance-dependent selection, surface escape and finite sensor acceptance can alter the observed profile, so a universal nominal/axial factor does not follow from equal phase index. A useful next diagnostic is to separate optical path length and destination acceptance by source position; it requires future physical path/source records, not a retuned attenuation fit. No alternate range or extra component was selected to obtain a pass.\n',
        f'![Attenuation profiles]({OUT}/V3_profiles.png)\n',
        '## V4 — effective slope versus microscopic velocity\n',
        f'The stated nondispersive reference is c/n={VREF:.6f} mm/ns for n=1.58. For each end/position, the observable is the sum of all detected photon arrival times divided by its detected photon count, relative to gun time zero. It is not the first photon or the leading-edge timestamp.\n',
        table(['Material','End','dt/dd [ns/mm] ± SE','v ± SE [mm/ns]','v/(c/n) ± SE','χ² / 5','Bound / linearity'],[[p['material'],p['end'],f"{p['slope_ns_per_mm']:.8f} ± {p['slope_se_ns_per_mm']:.8f}",pm(p,'speed_mm_per_ns','speed_se_mm_per_ns'),pm(p,'speed_over_c_n',digits=5),f"{p['chi2']:.2f}",p['bound_status']+' / '+p['model_status']] for p in r['V4']]),
        'The slope versus x is positive on the left and negative on the right; magnitudes are listed versus d. Every descriptive speed passes the registered upper-bound check, but every global straight-line fit fails. Thus these numbers summarize a profile and do not establish a unique bulk group velocity. The observed mean can change with distance through selection of surviving paths; its derivative is not an individual photon’s causal propagation speed. For a dispersive material, the actual group index also involves the derivative of phase index.\n',
        'Historical provenance: commit `e4110ce5f83bbf8bbb7e3e1731805e3631054c09` added a material GROUPVEL diagnostic at 3.17 eV. Commit `55766876ffa66e44e8a42461029d6276ff14d164` introduced an air GROUPVEL override to c/1.58 to address the reported post-reflection velocity-aliasing issue. The current source still contains that override at `src/DetectorConstruction.cc:260–267`. Its original claim of transparency relied on immediately killing escaped photons, a premise later changed by the guard correction. This is a source-history caveat, not a new experimental diagnosis or a reason to edit physics here. Future path/time records can test time increments by physical medium without depending on any engine callback.\n',
        f'![Arrival profiles]({OUT}/V4_profiles.png)\n',
        '## V5 — spectral prediction; independent per-end closure unavailable\n',
        table(['Material','PDE, configured energy-density convention','PDE, wavelength-density convention','Measured per-end ratio'],[[p['material'],f"{p['PDE_emission_energy_density']:.8f}",f"{p['PDE_emission_wavelength_density']:.8f}",'NOT EVALUABLE'] for p in r['V5_prediction']]),
        f'The calculation integrates the frozen [Broadcom PDE table]({REPO}/data/sipm/AFBR-S4N66P024M_pde.txt) against each actual run-local SSLG4 `scntComp1.txt`, with file hashes in [V5 metadata]({OUT}/V5_PDE_prediction.meta.json). Energy-domain amplitudes are linearly interpolated versus E=hc/λ; the tabulated PDE is linearly interpolated versus wavelength. The wavelength-density alternative is explicitly a different spectral convention, not a fitted uncertainty. No sensor curve or normalization changes.\n',
        'Detected END-left/right counts are available and tabulated with profile means; master logs also give an all-sensor entry aggregate. The latter includes TOP and cannot supply a per-end denominator for the grid. D0/D3 have per-end legacy boundary-state tallies, but these lack independent photon/event identity and incident spectra; they are not promoted into an engine-neutral validated incidence sample. The detected-hit `pde` branch samples already detected photons and cannot measure missed incidence. Therefore no measured detector-closure ratio is invented.\n',
        table(['Cell','End','Legacy absorption tally','Legacy detection tally','Legacy sum','Exported detected hits'],[[p['cell_id'],p['end'],p['legacy_absorption_tally'],p['legacy_detection_tally'],p['legacy_sum'],p['exported_detected_hits']] for p in legacy]),
        'This additional schema audit demonstrates a concrete obstacle: the legacy detection tally and the exported detected-hit count are different. Mixing the exported numerator with a classifier-derived denominator would not establish a closed, independently identified photon population. The table is an archival provenance diagnostic, not a new acceptance test based on internal instrumentation. The discrepancy warrants an event/photon-identity reconciliation if those streams are recorded in future; its cause is not inferred here.\n',
        'The pure-emission-weighted PDE equality also requires the incident spectrum to match that emission measure. Wavelength-dependent transport and any other optical source can change it; the hit schema has no source-class branch. A complete independent check records all incident photons at each end, their wavelengths and detection outcomes, then compares the detected count with the sum of incident-photon PDE probabilities. This distinguishes detector response from transport-induced spectral selection.\n',
        '## V6 — archival consolidation, preserving attribution and scope\n',
        table(['Check','Archived measurement','Source and limitation'],[[p['id'],str(p['measured'])+' '+p['unit'],p['provenance']+' '+p['audit']] for p in v6]),
        f'[V6 source rows and hashes]({OUT}/V6_source_rows.json) preserve the archived reflection table, all 42 symmetry rows and the four worker-reproducibility rows. Published physics measurements were not recalculated. In particular, EXEC34C verified detected-hit counts against logs; that is not conservation over all generated photons and their terminal destinations. EXEC29 D0 / EXEC30 V1 do have the complete published zero-residual ledger. The quoted reflectivity is confirmed in EXEC30 **V1**, not D3.\n',
        'The physical wrap is **Vikuiti 3M ESR, R=0.98**. “Mylar” remains a legacy factory/material name documented by commit `522a1d0`; it is not a hardware discrepancy. No reflectivity was changed. The preserved sensor/PDE, timing-reference and quantization provenance issues do not become numerical experimental comparisons in this task.\n',
        '## V7 — same-end correlation and direct source-relative spread\n',
        'Primary V2 uses left groups {0,1,2,3}/{4,5,6,7} and right groups {8,9,10,11}/{12,13,14,15}. These are fixed global IDs for both N_TOP=0 and N_TOP=70. The unchanged peak-normalized difference-of-exponentials response has rise 0.5 ns, fall 5 ns and threshold 4 PE. Archived T1/T2 were reused exactly: no new pulse fitting or electronics smearing. Secondary EndTop V1 retains {0,2}/{4,6}, plus 8 on the right; it is not substituted for V2.\n',
        'All 88 distributions (42 EndTop V2, 42 EndTop V1, four END-only V2) pass rho−3SE>0. Paired crossing efficiencies span '+f"{min(p['efficiency'] for p in r['V7']):.6f}–{max(p['efficiency'] for p in r['V7']):.6f}"+'; missing crossings remain in the generated-event bootstrap population and sigmas condition on paired finite crossings. All bootstrap replicas are finite. The primary rho range is '+f"{min(p['rho'] for p in r['V7'] if p['primary']):.5f}–{max(p['rho'] for p in r['V7'] if p['primary']):.5f}"+'. The V7 positivity decision is separate from the approximation that both group variances are equal.\n',
        table(['Cell','End','ρ ± SE','σ(T1) ± SE [ps]','σ(T2) ± SE [ps]','σ(Δ)/√2 ± SE [ps]','Q ± SE','F ± SE'],[[p['cell_id'],p['end'],pm(p,'rho'),pm(p,'sigma_T1_ps'),pm(p,'sigma_T2_ps'),pm(p,'sigma_delta_over_sqrt2_ps'),pm(p,'Q_direct_over_normalized'),pm(p,'F_equal_variance')] for p in r['V7'] if p['primary'] and p['cell_id'] in ['D0','D3','EJ204_xp0']]),
        'Here σ denotes ordinary sample standard deviation, ddof=1, on paired accepted events. It is **not** a Gaussian core fit or the canonical EXEC35 quantile half-width. The earlier robust results remain unchanged. Pearson covariance has an exact second-moment identity; replacing standard deviations by quantile widths would invalidate that identity.\n',
        'With s1=σ(T1), s2=σ(T2), Var(T1−T2)=s1²+s2²−2ρs1s2. The measured normalization bias Q=√2 s1/σ(T1−T2) equals √[2s1²/(s1²+s2²−2ρs1s2)] without assuming equal variances. It reduces to F=1/√(1−ρ) only for s1=s2. The exact unequal-variance expression agrees to relative 1e−10 in every case; this arithmetic identity is not independent physics evidence. The approximation Q=F fails its |Q−F|≤3SE check in 30/88 cases (14/46 primary V2 cases). Those discrepancies are retained, not absorbed into the correlation error.\n',
        'The direct σ(T1) and σ(T2) relative to a jitter-free gun time are unavailable experimentally without such a zero-time reference. Simulation supplies this quantity and thereby quantifies how a same-end difference can suppress a common timing contribution. This is an interpretive contribution, not a measurement of electronics or an experimental performance comparison. These timestamps still depend on the fixed deterministic pulse/threshold postprocessor; they are not raw single-photon transit widths.\n',
        f'![Same-end correlation]({OUT}/V7_rho.png)\n',
        '### Complete V7 results (ps for all standard deviations)\n',
        table(['Cell','Variant','End','ρ ± SE','s1 ± SE','s2 ± SE','sΔ ± SE','sΔ/√2 ± SE','Q ± SE','F ± SE','Q−F ± SE','ρ / equal-var'],[[p['cell_id'],p['variant'],p['end'],pm(p,'rho'),pm(p,'sigma_T1_ps'),pm(p,'sigma_T2_ps'),pm(p,'sigma_delta_ps'),pm(p,'sigma_delta_over_sqrt2_ps'),pm(p,'Q_direct_over_normalized'),pm(p,'F_equal_variance'),pm(p,'Q_minus_F'),p['status']+' / '+p['equal_variance_approximation']] for p in r['V7']]),
        '## Reproduction and artifacts\n',
        'Analysis command actually executed (does not launch a simulation):\n\n```bash\ncd /home/rrios/ej200_exec33_20260911\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/rrios/exec35_20260912/venv/bin/python -u analysis/validation/run_exec38.py --out /home/rrios/exec38_20260913\n```\n',
        'The output directory must be new; the runner refuses to overwrite an existing analysis. To reproduce, select another output path; the report publisher uses the archived EXEC38 path. Analysis checks: analytical exponential parameter recovery; exact unequal-variance covariance identity including a missing event; generated-event bootstrap reproducibility; and PASS/FAIL/overlap boundaries. No Geant4 build or simulation tests were run.\n',
        f'All figures have `.root`, `.csv` and `.meta.json` sidecars, plus PNG/PDF. [Profile points]({OUT}/profile_points.meta.json), [V3 fits and bootstrap/covariance ROOT]({OUT}/V3_attenuation.meta.json), [V4]({OUT}/V4_velocity.meta.json), [V7 full events and bootstrap]({OUT}/V7_correlations.meta.json), [production/availability inventory]({OUT}/V1_V2_V5_inventory.meta.json), [analysis provenance]({OUT}/analysis.meta.json). The golden JSON contains definitions, sources, measured/error fields, tolerances and exact per-cell configurations. Missing observations are null with explicit requirements.\n',
        '### Original simulation invocations (provenance only; not executed in EXEC38)\n',
        table(['Cell','Material; x; N','Seeds; workers; eventModulo','Simulation commit'],[[k,f"{v['material']}; {v['x_mm']} mm; {v['N_generated']}",f"{v['seeds']}; {v['workers']}; {v['eventModulo']}",v['simulation_commit']] for k,v in configs.items()]),
        '\n'.join(f"- `{k}`: `{v['command']}`" for k,v in configs.items()),
        '\n\nD0/D3 supply only centre-bar information; neither attenuation nor velocity fits are possible for END-only from these two single-position files. No END-only position scan is implied. Missing source/deposition, first-encounter and incident-spectrum observables remain explicit prerequisites for a complete engine-independent acceptance exercise.\n']
    REPORT.write_text('\n'.join(lines))
    dump(OUT/'delivery.meta.json',dict(created_utc=now(),report=str(REPORT),report_sha256=sha(REPORT),
         golden=str(GOLDEN),golden_sha256=sha(GOLDEN),publisher_code_sha256=sha(__file__),
         analysis_code_sha256=sha(HERE/'run_exec38.py'),no_simulations=True))
    print(REPORT);print(GOLDEN)


if __name__=='__main__':main()
