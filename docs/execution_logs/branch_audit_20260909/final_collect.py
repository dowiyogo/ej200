from audit import *
import re,shlex
p='/home/rrios/ej200'
rs=[json.loads(l) for l in DB.read_text().splitlines()]
trees={r['cmd'].split('ls-tree -r ',1)[1]:r['stdout'] for r in rs if 'ls-tree -r ' in r['cmd'] and len(shlex.split(r['cmd']))==6}
# Only the 18 audited tips, plus differing remote counterpart.
for b,t in trees.items():
 files=[l.split('\t',1)[1] for l in t.splitlines()]
 for f in files:
  if f.endswith('talk_v6.tex') or f in ['include/Materials.hh','src/Materials.hh']:
   git(p,'show',b+':'+f)
 git(p,'grep','-n','-I','-E',r'(^|[^0-9])69\.2([^0-9]|$)|TOP_SUM4_N1|polishedbackpainted',b,'--',':(exclude)src/external/**')
for f in ['analysis/top_npe_diag/top_npe_diag.csv','analysis/top_npe_diag/top_npe_diag_meta.json','presentations/v6/figs/v5_veff_fit_meta.json','presentations/v6/figs/v5_top_position_loo_meta.json','docs/branch_diagnosis/REFLECTIVITY_CHANGE.md']:
 git(p,'show','main:'+f)
git(p,'ls-tree','-r','origin/feat/bar-end-vikuiti');git(p,'show','origin/feat/bar-end-vikuiti:src/Materials.cc')
git(p,'rev-list','--left-right','--count','origin/feat/bar-end-vikuiti...feat/bar-end-vikuiti')
git(p,'log','--oneline','origin/feat/bar-end-vikuiti..feat/bar-end-vikuiti')
git(p,'branch','-a','--contains','50cf02e34cd7dddfb847739f6249afd1a600478f')
git(p,'log','--oneline','84e902c522d8af54ace11e5bb38202bf70e8c94d','--not','--remotes=origin')
git(p,'stash','show','-p','stash@{0}')
git(p,'ls-remote','origin')
for clone in [p,'/home/rrios/ej200_end']:
 git(clone,'status','--porcelain=v2','-b');git(clone,'for-each-ref','refs/heads','refs/tags','refs/stash','--format=%(refname)|%(objectname)')
# Resolve tags to commit IDs, rather than confusing annotated tag objects with commits.
for tag in git(p,'tag','--list').splitlines():git(p,'rev-parse',tag+'^{commit}')
print('Final collection done')
