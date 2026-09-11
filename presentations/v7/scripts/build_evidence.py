from pathlib import Path
import json,csv,hashlib,shutil,sys,collections,subprocess,math
import numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'/home/rrios/ej200_exec26_20260909/build_nightly_20260910/audit')
from campaign import table
V=Path(__file__).resolve().parents[1];W=Path('/home/rrios/ej200_exec26_20260909');SRC=V/'sources';FIG=V/'figs'
def rows(p):return list(csv.DictReader(Path(p).open()))
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def copytrio(name,base):
 for ext in ['csv','root','meta.json']:shutil.copy2(str(base)+'.'+ext,SRC/(name+'.'+ext))
 return {'original_base':str(base),'copied_base':'sources/'+name,'sha256':{e:digest(SRC/(name+'.'+e)) for e in ['csv','root','meta.json']}}
runs={}
m=json.loads((W/'build_exec27_20260910/run500/tables/metrics_comparison.meta.json').read_text())
for name in ['R1','R2']:
 runs[name]=m['runs'][name];runs[name]['x_mm']=0
 runs[name]['command']=str(W/('build_exec26_phase1b_20260909/ej200_bar_sim' if name=='R1' else 'build_exec27_20260910/ej200_bar_sim'))+' -m '+str(W/('build_exec26_phase1b_20260909/run500_ready/run500.mac' if name=='R1' else 'build_exec27_20260910/run500/run500.mac'))
 runs[name]['cwd']=str(Path(runs[name]['command'].split(' -m ')[1]).parent)
for name,build in [('A0','build_nightly_20260910'),('B0','build_nightly_20260910'),('D0','build_exec29_20260910'),('D1','build_exec29_20260910'),('D2','build_exec29_20260910'),('V1','build_exec30_20260910'),('D3','build_exec30_20260910')]:
 j=json.loads((W/build/'cells'/name/'cell_summary.json').read_text());runs[name]={**j['metadata'],'x_mm':j['x_mm']};runs[name]['summary']=str(W/build/'cells'/name/'cell_summary.json')
for p in sorted((W/'build_nightly_20260910/cells').glob('*/cell_summary.json')):
 j=json.loads(p.read_text())
 if j['N']==1000:runs[j['cell_id']]={**j['metadata'],'x_mm':j['x_mm'],'summary':str(p)}
refs={}
for n,p in [('guard',W/'build_exec27_20260910/run500/tables/metrics_comparison'),('bar_world',W/'build_exec27_20260910/run500/tables/bar_world_states'),('directional',W/'build_exec26_phase1b_20260909/run500_ready/tables/directional_states'),('metal_reflection',W/'build_exec30_20260910/cells/V1/reflection_panels'),('factorial',W/'build_exec30_20260910/factorial_gains'),('positions',W/'build_nightly_20260910/ab_comparison'),('design',W/'build_exec29_20260910/design_gain'),('multiplicity_dielectric',W/'build_exec30_20260910/cells/D3/encounter_multiplicity'),('multiplicity_metal',W/'build_exec30_20260910/cells/V1/encounter_multiplicity'),('terminal_scint',W/'build_exec29_20260910/cells/D1/terminal_fates_run0_scintillation'),('terminal_audit',W/'build_exec29_20260910/cells/D1/terminal_fates_run0_audit'),('claims_audit',Path('/home/rrios/ej200_deck_20260910/build_exec29_docs/claims_audit'))]:refs[n]=copytrio(n,p)
# Verified source quotations are audit evidence, never replacement physics data.
base='3dc2c37ef5eb363d3dbea8079a6dbec921d9a53d';source_files=['presentations/napkin_first_principles/napkin.py','presentations/napkin_first_principles/values/napkin_values.csv','presentations/v6/results_macros.tex','analysis/optim/phase_ab_optimal.csv','docs/DEFECTS_AND_CORRECTIONS.md']
source_refs={}
for f in source_files:
 result=subprocess.run(['git','-C',str(V.parents[1]),'show',base+':'+f],capture_output=True)
 data=result.stdout if result.returncode==0 else (Path('/home/rrios/ej200')/f).read_bytes()
 target=SRC/('historical_'+Path(f).name);target.write_bytes(data);source_refs[f]={'commit':base if result.returncode==0 else None,'origin':'verified git blob' if result.returncode==0 else 'pre-existing untracked historical artifact, hash recorded; no commit attribution','snapshot':str(target.relative_to(V)),'sha256':digest(target)}
values=[];macros={}
def add(name,value,display,source,run_ids=(),interpretation='measurement',derivation=''):
 macros[name]=display;values.append({'macro':name,'value':value,'display':display,'source':source,'run_ids':list(run_ids),'interpretation':interpretation,'derivation':derivation})
def num(name,value,source,run_ids=(),digits=None,interpretation='measurement',derivation=''):
 display=f'{value:.{digits}f}' if digits is not None else str(value)
 add(name,value,display,source,run_ids,interpretation,derivation)
g={ (r['run'],r['metric']):float(r['value']) for r in rows(SRC/'guard.csv')}
for label,n in [('Before','R1'),('After','R2')]:
 num('Guard'+label,g[n,'Npe_per_end'],'guard',[n],3);num('Guard'+label+'Sem',g[n,'event_SEM'],'guard',[n],3)
bw=[r for r in rows(SRC/'bar_world.csv') if r['run']=='R1'];total=sum(int(r['count']) for r in bw);tir=next(int(r['count']) for r in bw if r['status_name']=='TotalInternalReflection')
num('BarWorldCount',total,'bar_world',['R1']);num('KilledTirCount',tir,'bar_world',['R1']);num('TirFraction',100*tir/total,'bar_world',['R1'],1,derivation='100*TIR / all BarPV-to-WorldPV encounters')
dr=[r for r in rows(SRC/'directional.csv') if r['direction']=='air_to_wrap'];reflected=next(r for r in dr if r['status_name']=='FresnelReflection');metal=next(r for r in rows(SRC/'metal_reflection.csv') if r['panel']=='ALL')
num('OldReflection',float(reflected['percent']),'directional',['R1'],6);num('MetalReflection',float(metal['reflection_pct']),'metal_reflection',['V1'],6)
for name,key in [('OldReflectCount','count'),('OldEligible','eligible_denominator')]:num(name,int(reflected[key]),'directional',['R1'])
for name,key in [('MetalReflectCount','reflections'),('MetalEligible','eligible')]:num(name,int(metal[key]),'metal_reflection',['V1'])
factor=rows(SRC/'factorial.csv')
for i,name,ids in [(0,'EndGain',['D0','D3']),(1,'EndTopGain',['A0','B0']),(2,'GainRatio',['D0','D3','A0','B0'])]:
 num(name,float(factor[i]['value']),'factorial',ids,5);num(name+'Sem',float(factor[i]['SEM_approx']),'factorial',ids,5)
md=next(r for r in rows(SRC/'multiplicity_dielectric.csv') if r['panel']=='ALL');mm=next(r for r in rows(SRC/'multiplicity_metal.csv') if r['panel']=='ALL')
for label,row,rid,src in [('Dielectric',md,'D3','multiplicity_dielectric'),('Metal',mm,'V1','multiplicity_metal')]:
 num('Multiplicity'+label,float(row['encounters_per_started_optical']),src,[rid],8)
num('EncounterRatio',float(mm['encounters_per_started_optical'])/float(md['encounters_per_started_optical']),'multiplicity_metal + multiplicity_dielectric',['D3','V1'],1,derivation='ratio of consistent all-optical means')
design=rows(SRC/'design.csv')[0];num('DesignGain',float(design['value']),'design',['D1','D2'],5);num('DesignGainSem',float(design['SEM_independent_arms']),'design',['D1','D2'],5)
done=json.loads((W/'build_exec29_20260910/cells/D1/cell_summary.json').read_text());num('TopHits',done['metrics']['top'],'D1 summary',['D1']);num('GeneratedScint',done['metrics']['scint_generated'],'D1 summary',['D1']);num('TopHitRatio',100*done['metrics']['top']/done['metrics']['scint_generated'],'D1 summary',['D1'],2,interpretation='mixed-population ratio, not unique interception')
for name,key in [('DOneYield','npe_end_mean'),('DOneSem','npe_end_sem')]:num(name,done[key],'D1 summary',['D1'],4)
ledger=collections.Counter()
for row in rows(SRC/'terminal_scint.csv'):
 if row['kill_reason']!='none':category='Explicit world guard'
 elif row['process']=='OpAbsorption':category='Bulk absorption in bar' if row['volume']=='BarLV' else 'Bulk absorption in wrap'
 else:category='Other terminal transport/boundary'
 ledger[category]+=int(row['count'])
assert sum(ledger.values())==done['metrics']['started_scintillation']==done['metrics']['terminal_scintillation']
for macro,key in [('LedgerBar','Bulk absorption in bar'),('LedgerWrap','Bulk absorption in wrap'),('LedgerGuard','Explicit world guard'),('LedgerOther','Other terminal transport/boundary')]:num(macro,ledger[key],'terminal_scint',['D1'])
num('LedgerStart',done['metrics']['started_scintillation'],'terminal_audit',['D1']);num('LedgerEnd',done['metrics']['terminal_scintillation'],'terminal_audit',['D1']);num('LedgerResidual',0,'terminal_audit',['D1'],0,derivation='100*(starts - terminals)/starts; duplicates and unknowns are zero')
audit=rows(SRC/'claims_audit.csv');counts=collections.Counter(r['status'] for r in audit);num('AuditTotal',len(audit),'claims_audit',interpretation='historical source-audit count')
for status in ['valid','unverified','superseded','refuted']:num('Audit'+status.title(),counts[status],'claims_audit',interpretation='historical source-audit count, not corrected-physics validation')
npkin=(SRC/'historical_napkin.py').read_text();assert '"EJ-200": 1310.8, "EJ-204": 941.0, "EJ-230": 701.3' in npkin
for name,val in [('HistNpeA',1310.8),('HistNpeB',941.0),('HistNpeC',701.3),('HistSigmaLabel',52.1),('HistSigmaCenter',53.36),('HistSigmaMinimum',52.07)]:num(name,val,'historical source snapshots @ '+base,digits=2 if name in ['HistSigmaCenter','HistSigmaMinimum'] else 1,interpretation='AUDIT ONLY: quoted historical literal; not a new or accepted physics result')
# Configuration/code metadata are sourced, not event measurements.
for name,val,display,src in [('Zero',0,'0','run macros'),('Half',0.5,r'\tfrac{1}{2}','END pe/end definition'),('One',1,'1','run macros / ratio unity'),('TopCount',70,'70','D1 metadata'),('NoTopCount',0,'0','D0 metadata'),('GuardN',500,'500','R1/R2 metadata'),('MainN',2000,'2000','D0/D1/D2/D3/V1/A0/B0 metadata'),('ScanN',1000,'1000','position metadata'),('Energy',1,'1','run macros, GeV'),('SeedA',26092601,'26092601','all measured runs'),('SeedB',8349041,'8349041','all measured runs'),('MetalR',0.98,'0.98','f4d90a6 DetectorConstruction.cc:313'),('NominalPercent',98,'98','f4d90a6 R input times 100'),('MetalWorkers',12,'12','V1/D3 metadata'),('XNear',200,'200','position metadata'),('XMiddle',500,'500','position metadata'),('XFar',650,'650','position metadata')]:add(name,val,display,src,interpretation='configuration or definition')
macros.update({'MaterialName':'EJ-204/OPSC-101','GeantVersion':'11.4.0','Seeds':r'\SeedA\ \SeedB','DeckDate':'September 2026','DeckOld':'v6','DeckNew':'v7','MainSha':'3dc2c37','GuardSha':'218241a','ReflectorSha':'f4d90a6','AuditBaseSha':'8349041','AuditFileName':r'CLAIMS\_AUDIT\_20260910.md','ExecCurrent':r'EXEC\_32','ExecMtValidation':r'EXEC\_30 G0','BoundaryFile':r'\texttt{G4OpBoundaryProcess.cc}','NumericFile':r'\texttt{results\_macros\_v7.tex}','HistoricalLine':'71','HistoricalCsvLine':'36','OneWorker':'one worker','ShortAuditSource':r'\texttt{sources/claims\_audit}','HistDictName':r'\texttt{sim\_npe}','HistSigmaName':r'\texttt{sigma\_sim\_ej204}'})
for key,m in runs.items():
 if key in ['R1','R2','A0','B0','D0','D1','D2','D3','V1']:
  name={'R1':'ROne','R2':'RTwo','A0':'AZero','B0':'BZero','D0':'DZero','D1':'DOne','D2':'DTwo','D3':'DThree','V1':'VOne'}[key];macros['Sha'+name]=m['source_commit'][:7];macros['Id'+name]=key
for name,sha in [('GeneAl','51ca4ff'),('GeneTir','af5ddb7'),('GeneRefl','9b8361f'),('GenePort','fc5d8b0'),('GenePortOther','76b582c'),('GeneGap','bd78211'),('GeneGapPort','4cb2c23'),('GeneRegression','f3a8062'),('GeneMerge','87795b4'),('GeneParameter','c7acb7a')]:macros[name]=sha
# Numeric configuration metadata not counted as scientific estimates is still centralized.
macros.update({'HalfWidth':r'0.48\linewidth','WideWidth':r'0.60\linewidth','NarrowWidth':r'0.36\linewidth','FullFigureWidth':r'0.90\linewidth','FigureHeight':r'0.60\textheight','Gap':r'4pt','SmallGap':r'2pt'})
meta={'runs':runs,'source_triplets':refs,'source_snapshots':source_refs,'values':values,'method':'Historical measurements only. All new tables/figures derive from archived sidecars. No timing-resolution simulation or fit. Ratio SEMs retain the historical zero-covariance approximation.'}
(SRC/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
with (SRC/'value_provenance.csv').open('w') as f:
 w=csv.DictWriter(f,['macro','value','display','source','run_ids','interpretation','derivation']);w.writeheader();w.writerows([{**x,'run_ids':';'.join(x['run_ids'])} for x in values])
(V/'results_macros_v7.tex').write_text('% Generated from sources/provenance.json by scripts/build_evidence.py\n'+'\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros.items())+'\n')
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160,'pdf.fonttype':42})
blue='#0066b2';orange='#cc6600';green='#00994c'
def savefig(name,plotrows,description,run_ids):
 p=FIG/name;table(p,plotrows,{'runs':{k:runs[k] for k in run_ids},'sources':refs,'description':description},description)
 plt.tight_layout();plt.savefig(str(p)+'.pdf',bbox_inches='tight');plt.close()
 mp=p.with_suffix('.meta.json');j=json.loads(mp.read_text());j['figure']={'path':str(p)+'.pdf','sha256':digest(str(p)+'.pdf')};mp.write_text(json.dumps(j,indent=2)+'\n')
fig,ax=plt.subplots(figsize=(6.7,3.7));xx=['Original guard','Corrected guard'];yy=[g['R1','Npe_per_end'],g['R2','Npe_per_end']];ee=[g['R1','event_SEM'],g['R2','event_SEM']];ax.bar(xx,yy,yerr=ee,color=[orange,blue],capsize=5,width=.55);ax.set_ylabel('Photoelectrons per end');ax.set_ylim(0,470)
for x,y,e in zip(xx,yy,ee):ax.text(x,y+20,f'{y:.3f} ± {e:.3f}',ha='center')
savefig('guard_yields',[{'run':rid,'N':500,'npe_end':y,'SEM':e} for rid,y,e in zip(['R1','R2'],yy,ee)],'Guard-only comparison, EndTop70, N=500 per arm, same seed pair.', ['R1','R2'])
fig,ax=plt.subplots(figsize=(6.7,3.7));refvals=[float(reflected['percent']),float(metal['reflection_pct'])];ax.bar(['Dielectric (R1, EndTop)','Metal (V1, END-only)'],refvals,color=[orange,blue],width=.55);ax.axhline(98,color='gray',ls='--',lw=1);ax.set_ylim(0,115);ax.set_ylabel('Reflected / eligible encounters (%)')
for i,x in enumerate(refvals):ax.text(i,x+3,f'{x:.6f}%',ha='center')
savefig('boundary_reflectance',[{'run':'R1','reflections':int(reflected['count']),'eligible':int(reflected['eligible_denominator']),'percent':refvals[0]},{'run':'V1','reflections':int(metal['reflections']),'eligible':int(metal['eligible']),'percent':refvals[1]}],'Directional encounter fractions; different geometries explicitly labeled, not a matched yield comparison. Nominal R=0.98 is an input.', ['R1','V1'])
fig,ax=plt.subplots(figsize=(7.2,3.7));keys=list(ledger);vals=[ledger[k] for k in keys];ax.barh(keys,vals,color=[blue,orange,green,'#777777']);ax.set_xlabel('Terminal scintillation photons');ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
savefig('terminal_ledger',[{'category':k,'count':ledger[k],'fraction':ledger[k]/sum(vals)} for k in keys],'Disjoint terminal process/volume/explicit-kill grouping. Other transport/boundary is not equated to detected hits.', ['D1'])
pos=[r for r in rows(SRC/'positions.csv') if int(r['N'])==1000];pos.sort(key=lambda r:int(r['x_mm']));fig,ax=plt.subplots(figsize=(7.4,3.8));ax.errorbar([int(r['x_mm']) for r in pos],[float(r['gain']) for r in pos],yerr=[float(r['gain_sem']) for r in pos],fmt='o',color=blue,capsize=4);ax.axhline(1,color='gray',ls='--');ax.set_xlabel('Muon position x (mm)');ax.set_ylabel('Metal / dielectric END yield');ax.set_ylim(.97,1.34);ax.set_xticks([int(r['x_mm']) for r in pos]);ids=[k for k,v in runs.items() if v['N']==1000]
savefig('position_gain',[{'x_mm':int(r['x_mm']),'N_per_arm':1000,'gain':float(r['gain']),'SEM':float(r['gain_sem'])} for r in pos],'Measured EndTop70 position dependence. Error bars use the original zero-covariance ratio approximation.',ids)
print('Generated',len(values),'value records,',len(runs),'run records and four figure triplets.')
