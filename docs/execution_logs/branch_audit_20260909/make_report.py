import json,re,shlex,posixpath,hashlib
from pathlib import Path
D=Path('/home/rrios/branch_audit_20260909'); R=Path('/home/rrios/REPORT_branch_audit_20260909.md')
rs=[json.loads(l) for l in (D/'evidence.jsonl').read_text().splitlines()]
for i,r in enumerate(rs,1):r['id']=f'E{i:03}'
def get(*args,p='/home/rrios/ej200'):
 c=shlex.join(['git','-C',p,*args]);return next(r for r in rs if r['cmd']==c)
def find(s):return next(r for r in rs if s in r['cmd'])
def cite(r):return '['+r['id']+'](#'+r['id'].lower()+')'
def txt(*args,**kw):return get(*args,**kw)['stdout']
def esc(s):return str(s).replace('|','&#124;').replace('\n','<br>')
def block(s):return '\n```text\n'+s.rstrip()+'\n```\n'
def filetext(b,f):return get('show',b+':'+f)
def tree(b):return {l.split('\t',1)[1]:l.split('\t',1)[0] for l in txt('ls-tree','-r',b).splitlines()}
def excerpt(r,pattern,limit=12,context=0):
 ls=r['stdout'].splitlines();inds=set()
 for i,l in enumerate(ls):
  if re.search(pattern,l):inds.update(range(max(0,i-context),min(len(ls),i+context+1)))
 selected=sorted(inds)
 return '\n'.join(f'{i+1}: {ls[i]}' for i in selected[:limit])
local=txt('for-each-ref','--format=%(refname:short)','refs/heads').splitlines()
origin={l.split('refs/heads/')[1]:l.split()[0] for l in txt('ls-remote','--heads','origin').splitlines()}
localmeta={l.split('|')[0]:l.split('|') for l in txt('for-each-ref','refs/heads','--format=%(refname:short)|%(objectname)|%(upstream:short)|%(upstream:track)').splitlines()}
refmeta={l.split('|')[0]:l.split('|') for l in find('for-each-ref --sort=-committerdate')['stdout'].splitlines()}
models={
'diag/phase7-delta-2026-08-31':('Air-gap + polished dielectric_dielectric reflector border','0.95 active','Historical Phase-7 geometry; pre-EXEC_25 deck'),
'docs/branch-diagnosis-2026-08-31':('Air-gap + polished dielectric_dielectric reflector border','0.95 active','Historical geometry; diagnosis work'),
'exp/pair-scan-2026-06-11':('dielectric_metal polished bar skin','0.98 active','Pre-Phase-7 target optics'),
'feat/bar-end-vikuiti':('Air-gap + polished dielectric_dielectric reflector border','0.95 active','Merged ancestor; pre-EXEC_25 deck'),
'feat/bar-vikuiti':('dielectric_metal polished bar skin','0.98 active','Pre-Phase-7 target optics'),
'feat/ej204-bar-tir-only':('polished dielectric_dielectric TIR-only bar skin','0.98 factory retained, not used by bar skin','Alternative TIR configuration; scientific suitability UNVERIFIED'),
'feat/ej204-event-display-tracks':('dielectric_metal polished bar skin','0.98 active','Pre-Phase-7 target optics'),
'feat/ej228-cylinder':('Cylinder air-gap + dielectric_metal polished outer border','0.98 active; 0.95 unused bar factory retained','Alternative cylindrical geometry; suitability UNVERIFIED'),
'feat/ej228-tir-only':('Cylinder polished dielectric_dielectric TIR boundary','0.98/0.95 factories retained, not used by cylinder reflector','Alternative TIR cylinder; suitability UNVERIFIED'),
'feat/ej230-bar-tir-only':('polished dielectric_dielectric TIR-only bar skin','0.98 factory retained, not used by bar skin','Alternative TIR configuration; suitability UNVERIFIED'),
'feat/ej230-endonly-mylar':('dielectric_metal ground Mylar surface','0.90 active; 0.98 unused bar factory retained','Pre-Phase-7 target optics'),
'feat/ej230-sslg4':('dielectric_metal polished bar skin','0.98 active','Pre-Phase-7 target optics'),
'feat/endonly-mylar':('dielectric_metal ground/polished configurable bar skin','0.90 Mylar default; 0.98 fallback; runtime overrides UNVERIFIED','Pre-Phase-7 target optics'),
'feat/endtop-sslg4':('Air-gap + dielectric_metal polished outer reflector border','0.95 active default','Different pre-Phase-7 reflector implementation'),
'feature/sipm-electronics-response':('Passive Mylar volume; polished dielectric boundaries','No explicit 0.95/0.98 REFLECTIVITY assignment','Older alternative electronics/geometry; suitability UNVERIFIED'),
'main':('Air-gap + polished dielectric_dielectric reflector border','0.98 active; historical deck data 0.95','Current specified artifact, with audit gaps'),
'wip/host-stash-endtop-junio':('dielectric_metal polished explicit sibling-panel border','0.98 active','Pre-Phase-7 target optics'),
'wip/host-uncommitted-2026-08-31':('Air-gap + polished dielectric_dielectric reflector border','0.95 active','Historical geometry; unmerged work')}
branches=[]
for n in origin:
 b=n if n in local else 'origin/'+n; t=tree(b); decks=[f for f in t if f.endswith('talk_v6.tex')]
 behind,ahead=map(int,txt('rev-list','--left-right','--count','main...'+b).split())
 cm=get('log','--oneline','main..'+b); df=get('diff','--stat','main...'+b)
 figs=[]
 for deck in decks:
  dr=filetext(b,deck); source=dr['stdout']; stems=set(re.findall(r'\\anafig\{([^{}]+)\}',source))
  stems={s for s in stems if '#' not in s}
  stems.update(x.removeprefix('figs/') for x in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{(figs/[^{}]+)\}',source) if '#' not in x)
  for stem in sorted(stems):
   base=posixpath.join(posixpath.dirname(deck),'figs',stem)
   figs.append({'deck':deck,'stem':stem,'base':base,'pdf':base+'.pdf' in t,'root':base+'.root' in t,'csv':base+'.csv' in t,'meta':base+'.meta.json' in t,'altmeta':base+'_meta.json' in t,
   'elsewhere':[f for f in t if posixpath.basename(f) in [stem+'.root',stem+'.csv',stem+'.meta.json',stem+'_meta.json'] and not f.startswith(posixpath.dirname(deck)+'/figs/')]})
 branches.append(dict(name=n,ref=b,local=n in local,sha=refmeta[b][2],behind=behind,ahead=ahead,commits=cm,diff=df,tree=t,decks=decks,figs=figs))
(D/'derived_summary.json').write_text(json.dumps([{k:v for k,v in b.items() if k not in ['tree','commits','diff']} for b in branches],indent=2))
o=[]
def w(s=''):o.append(s+'\n')
w('# ej200 forensic branch audit — 2026-09-09')
w('Scope: HOST `t0minidaq`, the discovered additional clone, MSI through the existing reverse tunnel, and live GitHub origin. Audit time: '+find('date --iso')['stdout'].strip()+'. Report and evidence are outside both repositories.')
w('## 1. Verdict')
w('**`main` at `8349041140958226a0ac1cb3bb3e30aff2303435` is the branch containing the specified current deliverable, and matches live `origin/main`. No inspected branch satisfies all requested scientific checks.** The exact PDF passes criterion E, but main fails the strict sidecar completeness check and the requested `69.2 / TOP_SUM4_N1` marker check; design documentation still contains a suspect `701.3` value. A `polishedbackpainted` active reflector was not found. These are observed gaps, not grounds to silently substitute an older branch.')
w('Numbered evidence against criteria A–E:')
w('1. **A — available clones:** both local working trees are clean, HOST is on main, and MSI authentication failed. W1/W2 were not identified in the bounded filesystem search. '+cite(get('status','--porcelain=v2','-b'))+' '+cite(find('ssh -p 9022'))+' '+cite(find('find /home /mnt')))
w('2. **B — topology:** HOST has **11 local branches**, while live origin has **18 branch names**. Only local main and origin/main contain `8349041`; origin/HEAD is a symbolic alias. Only `feat/bar-end-vikuiti` is another local branch fully contained in main. '+cite(get('for-each-ref','--format=%(refname:short)','refs/heads'))+' '+cite(get('ls-remote','--heads','origin'))+' '+cite(get('branch','-a','--contains','8349041'))+' '+cite(get('branch','--merged','main')))
w('3. **C — optical configuration:** main uses `CreateBarSkinReflector()` on the air–wrap border (DetectorConstruction.cc:313,391–396), with `dielectric_dielectric`, `polished`, and `REFLECTIVITY={0.98,0.98}` (Materials.cc:347–357). This is **not** `polishedbackpainted`. `dielectric_metal` also appears in SiPM and another reflector factory; a global text hit alone cannot identify the active reflector. '+cite(filetext('main','src/Materials.cc'))+' '+cite(filetext('main','src/DetectorConstruction.cc')))
w('4. **C — EXEC_25:** main deck contains hybrid **514.9** at line 787 and **~1300** pooled photons at line 533. Its only literal `0.95` is the `\\Rvikuiti` definition at line 33. No `TOP_SUM4_N1` or contextually relevant `69.2` was found in the inspected tracked trees. However, `presentations/v6/FINAL_NUMBERS.md` declares `N_TOP=20 unless noted` at line 6 and still gives **701.3** at lines 24 and 38 without an END-only exception on those rows. **SUSPICIOUS** configuration inconsistency. Explicit END-only 701.3 uses in the deck are retained as historical values, not automatically treated as errors. '+cite(filetext('main','presentations/v6/talk_v6.tex'))+' '+cite(filetext('main','presentations/v6/FINAL_NUMBERS.md'))+' '+cite(get('grep','-n','-I','-E',r'(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted','main','--',':(exclude)src/external/**')))
mb=next(b for b in branches if b['name']=='main')
w(f'5. **C — sidecars:** the main TeX references {len(mb["figs"])} unique figure stems. **0/{len(mb["figs"])}** have all three adjacent `.root`, `.csv`, `.meta.json` files. Two have CSV and `_meta.json` (a different suffix), but no adjacent ROOT. A separate `blue_wscan_x0.root` exists under analysis; it does not complete the figure-specific trios. '+cite(get('ls-tree','-r','main'))+' '+cite(filetext('main','presentations/v6/talk_v6.tex')))
w('6. **D — rollback:** `pre-exec23-260902`, `pre-exec24-260902`, and `pre-exec25-260902` are contained in main. The first two are absent from the live origin tag listing. June `v-exec2*` tags are separate campaigns and many tag tips have no containing branch; keep them. '+cite(get('tag','--list','--sort=-creatordate','--format=%(creatordate:short) %(refname:short) %(objectname:short) %(subject)'))+' '+cite(get('ls-remote','--tags','origin')))
w('7. **E — exact artifact:** `presentations/v6/talk_v6.pdf` at main, origin/main, and 8349041 is **472931 bytes**, SHA-256 **`45784ae6a46ab2306004ade963b9e897611a21bd6f51877975b4347dc11a55de`**. Exact size and supplied hash-prefix match. No full expected SHA-256 was supplied, so only the supplied prefix can be compared with an expectation. '+cite(find('cat-file -p main:presentations/v6/talk_v6.pdf | sha256sum'))+' '+cite(find('cat-file -p main:presentations/v6/talk_v6.pdf | wc -c')))
w('## 2. Clone inventory and divergence')
w('| Clone | Path / endpoint | Branch and full HEAD | HEAD commit date | Working tree / worktrees / stash | Submodules / LFS |\n|---|---|---|---|---|---|')
for p,label in [('/home/rrios/ej200','HOST'),('/home/rrios/ej200_end','Additional local clone (not identified as W1/W2)')]:
 head=get('show','-s','--format=%H|%h|%cI|%aI|%s','HEAD',p=p); h=head['stdout'].strip().split('|');st=get('stash','list',p=p);wt=get('worktree','list','--porcelain',p=p);sub=get('submodule','status','--recursive',p=p)
 w('| '+label+' | `'+p+'` | `'+txt('symbolic-ref','--short','HEAD',p=p).strip()+'`<br>`'+h[0]+'` (`'+h[1]+'`) '+cite(head)+' | '+h[2]+' | Clean; one worktree at this path; '+str(len(st['stdout'].splitlines()))+' stash(es). '+cite(get('status','--porcelain=v2','-b',p=p))+' '+cite(wt)+' '+cite(st)+' | No submodule status entries / no index gitlinks; no tracked HEAD LFS pointers or LFS attributes found. LFS executable unavailable; LFS object-store state UNVERIFIED. '+cite(sub)+' '+cite(get('ls-files','--stage',p=p))+' '+cite(get('grep','-n','-I','-E',r'filter=lfs|version https://git-lfs.github.com/spec/v1','HEAD','--',p=p))+' |')
w('| MSI | Requested `/mnt/d/SHiP/ej200`, `rrios@localhost:9022` | UNVERIFIED | UNVERIFIED | UNVERIFIED — SSH exit 255, permission denied | UNVERIFIED |')
w('| W1 / W2 | Not found in searched roots | UNVERIFIED | UNVERIFIED | Not reachable / identity not established | UNVERIFIED |')
w('| origin | `git@github.com:dowiyogo/ej200.git` | default main; `8349041140958226a0ac1cb3bb3e30aff2303435` | 2026-09-03T00:37:35+02:00 (same commit as HOST) | Working tree, worktrees, stashes: N/A to Git remote | Server LFS storage UNVERIFIED |')
w('MSI was tried using the observed local username `rrios`; whether MSI expects a different user is UNVERIFIED. No new tunnel or alternate user login was attempted. Discovery searched `/home` to depth 4 and `/home /mnt /media /opt /srv` to depth 5. The broader find exited 1 with stderr suppressed; it does not prove absence beyond accessible searched paths. The original requested third physical clone is not substituted with an invented W1/W2 identity. '+cite(find('id -un'))+' '+cite(find('find /home /mnt'))+' '+cite(find('ssh -p 9022'))+' '+cite(get('ls-remote','--symref','origin','HEAD')))
w('| Left tip | Right tip | Left-only / right-only commits | Evidence |\n|---|---|---:|---|')
for left,right,args in [('HOST HEAD','origin/main',('rev-list','--left-right','--count','HEAD...origin/main')),('HOST HEAD','additional clone HEAD',('rev-list','--left-right','--count','HEAD...fb3749def29716dc84a33fcad53a21086bc96822')),('additional clone HEAD','origin/main',('rev-list','--left-right','--count','fb3749def29716dc84a33fcad53a21086bc96822...origin/main'))]:
 rr=get(*args);w(f'| {left} | {right} | '+rr['stdout'].strip().replace('\t',' / ')+' | '+cite(rr)+' |')
w('| MSI HEAD | HOST / origin/main / additional clone | UNVERIFIED / UNVERIFIED | Authentication failure |')
w('**Origin reachability:** all inspected local branch commits on HOST are reachable from live origin branch tips; the additional clone’s HEAD is exactly origin/feat/endonly-mylar. Its second local branch, main at `84e902c522d8af54ace11e5bb38202bf70e8c94d`, also has no commits absent from origin remote-tracking branches after HOST fetch. Being unmerged into main does not mean a commit is absent from origin. '+cite(get('for-each-ref','refs/heads','--format=%(refname:short)|%(objectname)|%(upstream:short)|%(upstream:track)',p='/home/rrios/ej200_end'))+' '+cite(get('log','--oneline','84e902c522d8af54ace11e5bb38202bf70e8c94d','--not','--remotes=origin')))
stash_live=next(r for r in reversed(rs) if r['cmd'].startswith("git -C /home/rrios/ej200 log --oneline 'stash@{0}' --not 000"))
w('The stash reachability test was also repeated against **every live origin-advertised ref**, including tags: the same two commits remain absent. '+cite(stash_live))
w('**Preserve the HOST stash.** The stash and index-parent commits `927c006` and `74d0cd1` are not reachable from origin remote-tracking branches. The stash patch adds 57 lines to `docs/branch_diagnosis/DATA_AUDIT_2026-08-31.md`. No claim is made about private/unadvertised GitHub refs or inaccessible clones. '+cite(get('log','--oneline','stash@{0}','--not','--remotes=origin'))+' '+cite(get('stash','show','-p','stash@{0}')))
w('`feat/bar-end-vikuiti` is 13 commits ahead of its same-named origin branch, but all 13 are already reachable from origin/main; they are not origin-missing commits. '+cite(get('rev-list','--left-right','--count','origin/feat/bar-end-vikuiti...feat/bar-end-vikuiti'))+' '+cite(get('branch','-a','--contains','50cf02e34cd7dddfb847739f6249afd1a600478f')))
w('## 3. Branch table — 18 observed branch names')
w('There are 11 HOST-local heads and seven names represented only by origin refs. Each row uses the local tip when present; otherwise the verified origin tip. Ahead / behind is **branch-only / main-only** (the reverse order of the raw `main...branch` command). “Exclusive” means reachable from branch but not main; it does not establish patch uniqueness or uniqueness against every other feature branch. File stats are the requested merge-base-to-branch diff. Full commits, file lists, markers, and evidence are in section 7.')
w('| Branch / scope / tip | Ahead / behind main | Merged / contains 8349041 | Origin exists / same tip | Exclusive commits / files | Optical surface / reflectivity | EXEC_25: 514.9 / 1300 / 69.2 | Deck sidecars |\n|---|---:|---|---|---|---|---|---|')
for ix,b in enumerate(branches,1):
 n=b['name'];model,rval,cl=models[n];stats=b['diff']['stdout'].strip().splitlines();stat=stats[-1].strip() if stats else '0 files';h='HIT / HIT / MISS' if n=='main' else 'MISS / MISS / MISS';side=(f'0/{len(b["figs"])} complete trios' if b['decks'] else 'MISS: no talk_v6 source')
 w('| `'+n+'`<br>'+('HOST local' if b['local'] else 'origin only')+'; `'+b['sha'][:7]+'` | '+str(b['ahead'])+' / '+str(b['behind'])+' | '+('Yes' if b['ahead']==0 else 'No')+' / '+('Yes' if n=='main' else 'No')+' | Yes / '+('Yes' if b['sha']==origin[n] else 'No: '+origin[n][:7])+' | '+str(b['ahead'])+' commits; '+esc(stat)+'; [details](#b'+str(ix)+') | '+esc(model+'; '+rval)+' | '+h+' | '+side+' |')
w('All rows have **MISS for active `polishedbackpainted`**. All except the passive-wrap electronics branch have `dielectric_metal` factory/source hits, but their active reflector differs as described. TIR-only and cylinder branches are not automatically called invalid because another unused factory is metallic. For branches without a talk_v6 tree, EXEC_25 marker MISS means the requested deliverable context is absent, not that a coincidental raw-data number was impossible.')
w('The differing origin/feat/bar-end-vikuiti tip is `1bfb82743270dcc3ea9b17be03a6af5704cc6850` and retains `analysis/presentation_v6/talk_v6.tex`; its local tip has the relocated `presentations/v6/talk_v6.tex`. Both retain R=0.95 code and lack the final EXEC_25 corrections. '+cite(get('ls-tree','-r','origin/feat/bar-end-vikuiti'))+' '+cite(filetext('origin/feat/bar-end-vikuiti','src/Materials.cc'))+' '+cite(filetext('origin/feat/bar-end-vikuiti','analysis/presentation_v6/talk_v6.tex')))
w('## 4. Taxonomy')
w('These groups overlap: scientific obsolescence is not permission to discard unmerged work.')
w('**(a) Merged, eligible for local branch-label cleanup:** `feat/bar-end-vikuiti` only, at `50cf02e...`. It has 0 commits exclusive to main, and HOST worktree inventory shows main checked out. Main itself is retained. Existing same-named origin is older, so create the proposed exact-tip backup tag before any manual deletion. Remote deletion is deferred while MSI/W1/W2 remain unverified.')
w('**(b) Exclusive, unmerged work — retain:** '+', '.join('`'+b['name']+'` ('+str(b['ahead'])+' commits)' for b in branches if b['ahead'])+'. The complete commit lists are printed for each branch in section 7; none is shortened to the first 20. Shared commits may appear in multiple lists. No patch-equivalence test was used to justify deletion.')
w('**(c) Obsolete for the target Phase-7 bar-reflector baseline:** '+', '.join('`'+n+'`' for n,(_,_,cl) in models.items() if 'Pre-Phase-7' in cl or 'pre-Phase-7' in cl)+'. Their inspected geometry uses a metallic bar skin, older sibling-panel boundary, or the different metallic air-gap reflector implementation. See source excerpts and call sites in section 7. “Obsolete” here means different from the target main configuration; this audit did not rerun physics or invalidate every historical study.')
w('**(d) Indeterminate or incomplete scientifically:** TIR-only EJ-204/EJ-230, both EJ-228 cylinder studies, and the electronics-response branch have different purposes; equivalence to the requested scientific deliverable is UNVERIFIED. Diagnosis/WIP branches with historical dielectric air-gap code lack the final deliverable. Main remains the current specified artifact but is incomplete against the full audit checklist. MSI/W1/W2 are UNVERIFIED.')
w('## 5. Tags and rollback mapping')
w('Dates below are tag creator dates (commit dates for lightweight tags), not dates inferred from tag-name strings. Campaign labels are explicitly taken from names and stored tag/commit messages. June optical EXEC_25 and September presentation EXEC_25 are different contexts. A containing branch is a reachability result, not evidence of the branch on which the tag was originally created.')
taglist=get('tag','--list','--sort=-creatordate','--format=%(creatordate:short) %(refname:short) %(objectname:short) %(subject)')
remote_tags=get('ls-remote','--tags','origin')['stdout']
groups={};tagrows=[]
for l in taglist['stdout'].splitlines():
 date,tag,obj,subject=l.split(' ',3);rr=get('branch','-a','--contains',tag)
 containing=tuple(x.strip().lstrip('* ').removeprefix('remotes/') for x in rr['stdout'].splitlines() if ' -> ' not in x)
 if containing not in groups:groups[containing]='G'+str(len(groups)+1)
 commit=txt('rev-parse',tag+'^{commit}').strip()
 campaign={'checkpoint/pre-endtop-sslg4-2026-06-10':'Before EndTop fork; target commit says CODEX_EXEC_09','checkpoint/pre-pairscan-2026-06-11':'Before pair scan; preserves EXEC_12b','checkpoint/pre-physics-baseline-2026-05-08':'Before physics baseline; EXEC_N UNVERIFIED','diag-photon-budget-v1':'Photon-budget diagnostic; EXEC_N UNVERIFIED','physics-baseline-v1':'Physics baseline; EXEC_N UNVERIFIED','exec07-09-analysis':'EXEC_07–09 (tag name)','v-exec21-optfix':'EXEC_21 label; stored commit message says exec22 follow-up'}.get(tag)
 if campaign is None:
  m=re.search(r'exec(\d+[a-z]?)',tag);campaign=('Before EXEC_'+m.group(1) if tag.startswith(('pre-','checkpoint/pre-')) else 'EXEC_'+m.group(1)) if m else 'UNVERIFIED'
 tagrows.append((date,tag,obj,commit,campaign,groups[containing],rr,subject))
w('| Date | Tag | Tag object / peeled commit | Campaign mapping | Containing branches | On origin |\n|---|---|---|---|---|---|')
for date,tag,obj,commit,campaign,g,rr,subject in tagrows:w('| '+date+' | `'+tag+'` | `'+obj+'` / `'+commit+'` | '+campaign+' | '+g+' '+cite(rr)+' | '+('Yes' if '\trefs/tags/'+tag+'\n' in remote_tags else 'No')+' |')
for containing,g in groups.items():w('- **'+g+'**: '+(', '.join('`'+x+'`' for x in containing) if containing else '**No local or origin-tracking branch contains this tag tip. Preserve the tag itself.**'))
w('Tag dates and subjects: '+cite(taglist)+'. Tag messages: '+cite(get('for-each-ref','refs/tags','--format=%(refname:short)|%(objectname)|%(*objectname)|%(contents)'))+'. Live remote existence: '+cite(get('ls-remote','--tags','origin'))+'. Do not treat `pre-exec25-260902` as a backup of the final 8349041 deliverable: it precedes those corrections.')
w('## 6. Proposed cleanup — NOT EXECUTED')
w('Only one local branch-label deletion is proposed. No source edits, remote deletions, tag deletions, stash drops, or clone removals are proposed. All commands in the following block are for René to review and execute manually. Recheck live main, origin/main, worktrees and the exact branch tip before use; this is a dated snapshot.')
w('```bash\n# Preserve the final deliverable; proposed new tag, NOT created by this audit.\ngit -C /home/rrios/ej200 tag audit/20260909/current-main 8349041140958226a0ac1cb3bb3e30aff2303435\n\n# Preserve the exact merged local branch tip; origin of the same name is 13 commits older.\ngit -C /home/rrios/ej200 tag audit/20260909/feat-bar-end-vikuiti 50cf02e34cd7dddfb847739f6249afd1a600478f\n\n# Preserve the stash commit and its parents; keep the stash itself.\ngit -C /home/rrios/ej200 tag audit/20260909/host-data-audit-stash '+txt('show','-s','--format=%H|%P|%cI|%s','stash@{0}').split('|')[0]+'\n\n# DESTRUCTIVE, MANUAL ONLY: remove this local branch label after verifying the backup tag.\n# Justification: exact tip is contained in main and origin/main; zero main-exclusive commits.\n# Backup: audit/20260909/feat-bar-end-vikuiti. No worktree here has this branch checked out.\ngit -C /home/rrios/ej200 branch -D feat/bar-end-vikuiti\n```')
w('The force spelling is explicit because the branch tracks an older origin counterpart even though its tip is contained in main. This is not permission to use the same command for any unmerged branch. Backup tags above are local proposals, not published remote backups. Tag creation and branch deletion were **not executed**.')
w('## 7. Per-branch evidence: science, complete exclusive commits, and file stats')
for ix,b in enumerate(branches,1):
 n=b['name'];ref=b['ref'];w(f'<a id="b{ix}"></a>\n\n### {ix}. {n}')
 w('Inspected ref `'+ref+'`, full SHA `'+b['sha']+'`. Classification: **'+models[n][2]+'**. '+cite(find('for-each-ref --sort=-committerdate')))
 w('**B:** ahead '+str(b['ahead'])+', behind '+str(b['behind'])+'; merged '+('YES' if b['ahead']==0 else 'NO')+'; contains 8349041 '+('YES' if n=='main' else 'NO')+'. '+cite(get('rev-list','--left-right','--count','main...'+ref))+' '+cite(get('branch','-a','--merged','main'))+' '+cite(get('branch','-a','--contains','8349041')))
 if b['local']:
  meta=localmeta[n];w('Upstream: '+('`'+meta[2]+'` '+(meta[3] or '(no ahead/behind annotation)') if meta[2] else 'none configured')+'. Origin counterpart: `'+origin[n]+'`. '+cite(get('ls-remote','--heads','origin')))
 w('Tip date / author / subject: '+esc(' / '.join([refmeta[ref][1],refmeta[ref][5],refmeta[ref][6]]))+'. '+cite(find('for-each-ref --sort=-committerdate')))
 mats=filetext(ref,'src/Materials.cc');geom=filetext(ref,'src/DetectorConstruction.cc');hdr=filetext(ref,'include/Materials.hh')
 w('**C1/C2:** active configuration: '+models[n][0]+'. Reflectivity: '+models[n][1]+'. Active `polishedbackpainted`: **MISS**. `dielectric_metal` in surface source: **'+('HIT' if 'dielectric_metal' in mats['stdout'] else 'MISS')+'**. Source excerpts below identify factory and binding; comments and unrelated emission-spectrum values are not reflectivity assignments. '+cite(mats)+' '+cite(geom)+' '+cite(hdr))
 w(block(excerpt(mats,r'^G4OpticalSurface\*|SetType\(|SetFinish\(|(?:refl|Reflectivity).*\{.*(?:0\.9|reflectivity)|AddProperty\("REFLECTIVITY"',45)))
 r98={'main','exp/pair-scan-2026-06-11','feat/bar-vikuiti','feat/ej204-event-display-tracks','feat/ej228-cylinder','feat/ej230-sslg4','wip/host-stash-endtop-junio'}
 r95={'diag/phase7-delta-2026-08-31','docs/branch-diagnosis-2026-08-31','feat/bar-end-vikuiti','feat/endtop-sslg4','wip/host-uncommitted-2026-08-31'}
 w('Active/default reflector assignment markers: **R=0.98 '+('HIT' if n in r98 else 'MISS')+'; R=0.95 '+('HIT' if n in r95 else 'MISS')+'**. MISS includes no active reflector assignment for TIR-only/passive-wrap variants; runtime overrides are not certified. Factory assignments and binding lines above/below provide the evidence.')
 if n=='feat/endonly-mylar':
  defaults=filetext(ref,'include/DetectorConstruction.hh');w('Default Mylar mode and R=0.90, with configurable parameters: '+cite(defaults)+block(excerpt(defaults,r'fTopSurface =|fMylarReflectivity =|fMylarSpecularLobe =|fMylarSigmaAlpha =',8)))
 broad=get('grep','-n','-I','-E',r'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1',ref,'--',':(exclude)src/external/**')
 loose=[x for x in broad['stdout'].splitlines() if re.search(r'(?<!\d)0\.95(?!\d)',x)]
 w('All loose `0.95` project-tree hit lines (including documents and data, excluding bundled external libraries): **'+str(len(loose))+'**. Full paths, line numbers and contents are preserved in '+cite(broad)+' and the complete outputs appendix. No fixed line-number assumption was used.')
 w('Geometry call sites: '+cite(geom)+block(excerpt(geom,r'CreateBarSkinReflector\(|CreateMylarReflector\(|CreateVikuitiSurface\(|CreateBarSurface\(|new G4LogicalSkinSurface|TIR.only|passive reflector',35,1)))
 w('Loose `0.95` in project surface/header sources (including historical comments; not automatically active): '+cite(mats)+' '+cite(hdr)+block(('src/Materials.cc\n'+excerpt(mats,r'(?<!\d)0\.95(?!\d)',14)+'\ninclude/Materials.hh\n'+excerpt(hdr,r'(?<!\d)0\.95(?!\d)',8)).strip()))
 nr=get('grep','-n','-I','-E',r'(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted',ref,'--',':(exclude)src/external/**')
 w('**C3:** `69.2` associated with `TOP_SUM4_N1`: **MISS**. Exact contextual search '+cite(nr)+'. '+('Raw matches below are historical finish-test documentation or percentages, not the requested configuration/result marker.'+block(nr['stdout']) if nr['stdout'] else 'No matches outside the bundled external libraries.'))
 if b['decks']:
  for deck in b['decks']:
   dr=filetext(ref,deck);ls=dr['stdout'].splitlines();l95=[f'{i+1}: {x}' for i,x in enumerate(ls) if re.search(r'(?<!\d)0\.95(?!\d)',x)];hits=[('514.9',any('514.9' in x for x in ls)),('1300',any('1300' in x for x in ls))]
   w('Deck `'+deck+'`: '+', '.join('**'+v+' '+('HIT' if h else 'MISS')+'**' for v,h in hits)+'. '+cite(dr))
   w(block(excerpt(dr,r'514\.9|1300|701\.3|Rvikuiti|NpeG',45)))
   w('Exact loose `0.95` occurrences in the deck: **'+str(len(l95))+'**. '+('Only the macro definition; no unapplied literal remains in this deck.' if n=='main' else '**OPEN-07-style macro gap:** historical literals remain; the requested final macro cleanup is absent. This label describes the requested audit check, not a verified issue identifier.')+block('\n'.join(l95)))
   w('701.3 contextual caution: '+('**SUSPICIOUS** in companion FINAL_NUMBERS.md under its N_TOP=20 default; deck S20 now correctly uses 514.9. END-only-labeled 701.3 occurrences alone are not errors.' if n=='main' else '**SUSPICIOUS**: the old deck retains 701.3 before the hybrid correction; inspect the design frame below.'))
   w(block(excerpt(dr,r'701\.3|NpeG|Design Decision',55,3)))
 else:w('**514.9 hybrid / 1300 S13: MISS / MISS.** No tracked talk_v6 source exists in this ref; the whole-tree marker search is retained in evidence. '+cite(get('ls-tree','-r',ref)))
 w('**C4:** '+('No talk_v6 source; deck figure provenance **MISS**, not a vacuous pass.' if not b['figs'] else f'**0/{len(b["figs"])} complete exact-suffix figure trios.**')+' '+cite(get('ls-tree','-r',ref)))
 if b['figs']:
  w('| Figure stem (relative to deck figs/) | Figure PDF | .root | .csv | .meta.json | _meta.json alternative | Matching sidecars elsewhere in tree |\n|---|---|---|---|---|---|---|')
  for f in b['figs']:w('| `'+f['stem']+'` | '+' | '.join('HIT' if f[k] else 'MISS' for k in ['pdf','root','csv','meta','altmeta'])+' | '+(', '.join('`'+v+'`' for v in f['elsewhere']) or 'none')+' |')
  w('Figure list is derived from literal `\\anafig{...}` and `\\includegraphics{figs/...}` calls in the tracked TeX. No dynamic figure paths were assumed; TeX build execution and contents inside figure PDFs were not tested.')
 w('**Complete exclusive commits** (`main..ref`; the first 20 entries also fulfill the requested preview): '+cite(b['commits'])+block(b['commits']['stdout'] or '(none)'))
 w('**Files touched on the branch since its merge base with main:** '+cite(b['diff'])+block(b['diff']['stdout'] or '(none)'))
 live=next(r for r in rs if r['cmd'].startswith(shlex.join(['git','-C','/home/rrios/ej200','log','--oneline',ref,'--not'])+' ') and '--remotes' not in r['cmd'])
 w('Commits absent from all 18 advertised origin branch histories: **'+str(len(live['stdout'].splitlines()))+'**. '+cite(live))
w('## 8. Additional main evidence and uncertainties')
w('The tracked `analysis/top_npe_diag/top_npe_diag.csv` and its metadata provide the observable data behind the hybrid correction. This audit read their bytes; it did not rerun the raw ROOT event analysis. '+cite(filetext('main','analysis/top_npe_diag/top_npe_diag.csv'))+' '+cite(filetext('main','analysis/top_npe_diag/top_npe_diag_meta.json')))
w('Main still has historical R=0.95 prose in Materials.cc comments, DetectorConstruction.cc comments, CONFIGURATION_AUDIT.md, FINAL_NUMBERS.md, and README/REVISION_NOTES. The single remaining literal in the current deck is the deliberately defined macro. The report does not equate every historical 0.95 statement to an active 0.95 assignment. All marker hits, including loose literals outside the deck, are preserved in the linked evidence outputs.')
w('`OPEN-07` in the current EXEC_25 report actually labels a superseded PowerPoint artifact and is marked resolved there; it is not the same issue description as the user’s OPEN-07 shorthand. '+cite(get('grep','-n','-I','-E',r'TOP_SUM4_N1|69[.,]2|69\.1|sigmaTOP|TOP.*N1|OPEN-07','main','--',':(exclude)src/external/**')))
w('- MSI HEAD, branches, worktrees, dirty files, stash, submodules, LFS and divergence: **UNVERIFIED** because SSH authentication failed. W1/W2 identity and reachability: **UNVERIFIED** within the bounded search. Historical reports describing MSI are not substituted for a live inspection.\n- Full repository history outside refs, dangling objects, ignored/generated worktree data and unadvertised origin refs were not comprehensively audited. “Clean” means the requested porcelain status has no tracked/untracked changes; it does not mean there are no ignored datasets.\n- LFS executable is unavailable. No HEAD pointer or filter=lfs marker was found in either local clone; actual LFS server/object-store content remains **UNVERIFIED**. No submodule status entries or index gitlinks were found.\n- No simulation, ROOT regeneration, TeX build, timing recalculation, or PDF rendering was run. Scientific correctness beyond the requested source/marker/artifact checks is **UNVERIFIED**. The old R=0.95 dataset cannot be certified as an R=0.98 result from code text.\n- No valid `69.2 / TOP_SUM4_N1` pairing was found. Incidental values in percentages, CSV numbers, or external material spectra were excluded as scientific evidence.\n- Sidecar absence is verified in the committed trees. Existence of ignored copies, MSI-only files, renamed provenance mappings not explicitly encoded in these file names, or reproducible regeneration is **UNVERIFIED**. `_meta.json` is reported separately from the specified `.meta.json`.\n- No branch count was forced to 18 local heads. The 18-row table is the observed branch-name union; the seven remote-only rows are explicit.\n- Patch equivalence between unmerged histories was not tested. All main-exclusive commits remain preservation candidates.\n- Original branch of tag creation and EXEC_N mappings not stated in tag names/messages are **UNVERIFIED**.\n- Proposed backup tags do not exist as a result of this audit; no cleanup command was executed.')
w('## Appendix A. Exact command ledger and complete outputs')
w('Every reference E### below identifies the exact recorded command, UTC timestamp and exit status. Complete untruncated stdout/stderr are available in [evidence_outputs.md](branch_audit_20260909/evidence_outputs.md), keyed by the same identifiers, and machine-readable [evidence.jsonl](branch_audit_20260909/evidence.jsonl). `git grep` exit 1 means no matches; missing `git lfs` also exited 1 and is interpreted from stderr. The deliberately attempted `main:src/Materials.hh` path failed (128); the real header was subsequently located and read as `include/Materials.hh`. No missing-path result is treated as source evidence.')
w('Git read commands used environment `GIT_OPTIONAL_LOCKS=0`, `GIT_TERMINAL_PROMPT=0`, `GIT_SSH_COMMAND="ssh -o BatchMode=yes -o ConnectTimeout=10"`. The only authorized repository update was HOST fetch, with `fetch.prune=false`, `fetch.pruneTags=false`, `gc.auto=0`, `maintenance.auto=false`, and `--no-prune`. No checkout, merge, rebase, push, destructive reset, cleanup, branch/tag deletion, or stash mutation was executed. Initial and final porcelain outputs for both local clones are identical. Fetch success is verified by its exit code, not inferred from stale tracking refs.')
for r in rs:
 w('<a id="'+r['id'].lower()+'"></a>\n\n**'+r['id']+'** — '+r['time']+' — exit `'+str(r['exit'])+'`\n\n```bash\n'+r['cmd']+'\n```')
w('## Appendix B. Method and preservation')
w('All Git commands specified `git -C <clone-path>`. Branch content was read with git show, git grep, git cat-file and git ls-tree; no branch checkout was used. The report, collection scripts, derived table and evidence were written exclusively under `/home/rrios` outside ej200 and ej200_end. Discovery and prerequisite instruction reads were filesystem reads. The first discovery/probe was repeated in the timestamped ledger for reproducibility.\n\nFigure counts, source line numbers, branch counts and group membership are deterministic aggregations of the recorded command outputs. Line numbers in source excerpts are enumerated from the exact git-show output; discovery used patterns rather than assumed line positions. Figures use tree membership checks; no renderer or uncommitted worktree state was used to certify a sidecar. The table-generation implementation is available at [make_report.py](branch_audit_20260909/make_report.py).')
R.write_text(''.join(o))
with (D/'evidence_outputs.md').open('w') as f:
 f.write('# Complete command evidence — ej200 audit 2026-09-09\n\n')
 for r in rs:
  fence='`' * max(3, 1+max([len(v) for v in re.findall(r'`+',r['stdout']+r['stderr'])] or [0]))
  f.write('## '+r['id']+'\n\n'+r['time']+'; exit '+str(r['exit'])+'\n\n```bash\n'+r['cmd']+'\n```\n\nstdout:\n'+fence+'text\n'+r['stdout']+'\n'+fence+'\n\nstderr:\n'+fence+'text\n'+r['stderr']+'\n'+fence+'\n\n')
print('Report',R,'bytes',R.stat().st_size,'records',len(rs),'branches',len(branches),'main figures',len(mb['figs']))
