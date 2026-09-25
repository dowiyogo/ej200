from audit import *
import re,shlex
path='/home/rrios/ej200'
rs=[json.loads(l) for l in DB.read_text().splitlines()]
local=next(r['stdout'].splitlines() for r in rs if r['cmd'].endswith("'--format=%(refname:short)' refs/heads"))
remote=next(r['stdout'] for r in rs if r['cmd'].endswith('ls-remote --heads origin'))
refs=[b if b in local else 'origin/'+b for b in [l.split('refs/heads/')[1] for l in remote.splitlines()]]
for b in refs:
 if b not in local:
  for args in [('rev-list','--left-right','--count',f'main...{b}'),('log','--oneline',f'main..{b}'),('diff','--stat',f'main...{b}'),('ls-tree','-r',b),('grep','-n','-I','-E',r'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti',b,'--')]:git(path,*args)
 for f in ['src/Materials.cc','src/DetectorConstruction.cc']:
  git(path,'show',f'{b}:{f}')
 git(path,'grep','-n','-I','-E',r'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti|TOP_SUM4_N1',b,'--',':(exclude)src/external/**')
 git(path,'log','--oneline',b,'--not',*[l.split()[0] for l in remote.splitlines()])
for args in [('branch','-a','--merged','main'),('for-each-ref','--format=%(refname)|%(objectname)','refs/heads','refs/remotes'),('stash','show','--stat','stash@{0}'),('show','-s','--format=%H|%P|%cI|%s','stash@{0}'),('log','--oneline','stash@{0}','--not','--remotes=origin'),('rev-list','--left-right','--count','HEAD...origin/main'),('rev-list','--left-right','--count','HEAD...fb3749def29716dc84a33fcad53a21086bc96822'),('rev-list','--left-right','--count','fb3749def29716dc84a33fcad53a21086bc96822...origin/main'),('ls-remote','--symref','origin','HEAD')]:git(path,*args)
for p in [path,'/home/rrios/ej200_end']:
 git(p,'for-each-ref','refs/heads','--format=%(refname:short)|%(objectname)|%(upstream:short)|%(upstream:track)')
 git(p,'grep','-n','-I','-E',r'filter=lfs|version https://git-lfs.github.com/spec/v1','HEAD','--')
 git(p,'ls-tree','-r','HEAD','.gitmodules','.gitattributes')
for b in ['main','origin/main','8349041']:
 for op in ['sha256sum','wc -c']:
  run(['bash','-o','pipefail','-c',f'git -C {path} cat-file -p {shlex.quote(b+":presentations/v6/talk_v6.pdf")} | {op}'])
for f in ['presentations/v6/talk_v6.tex','presentations/v6/results_macros.tex','presentations/v6/FINAL_NUMBERS.md','presentations/v6/README.md','presentations/v6/EXEC_25_REPORT.md','src/Materials.hh']:
 git(path,'show','main:'+f)
git(path,'grep','-n','-I','-E',r'TOP_SUM4_N1|69[.,]2|69\.1|sigmaTOP|TOP.*N1|OPEN-07','main','--',':(exclude)src/external/**')
print('Extended',len(refs),'refs')
