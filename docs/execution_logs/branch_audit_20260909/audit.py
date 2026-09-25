import subprocess,json,os,datetime
from pathlib import Path
OUT=Path('/home/rrios/branch_audit_20260909')
DB=OUT/'evidence.jsonl'
def run(args,timeout=90):
 env=os.environ.copy();env['GIT_OPTIONAL_LOCKS']='0';env['GIT_TERMINAL_PROMPT']='0';env['GIT_SSH_COMMAND']='ssh -o BatchMode=yes -o ConnectTimeout=10'
 try:
  p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env,timeout=timeout)
  r={'cmd':__import__('shlex').join(args),'exit':p.returncode,'stdout':p.stdout.decode(errors='replace'),'stderr':p.stderr.decode(errors='replace')}
 except subprocess.TimeoutExpired as e:r={'cmd':__import__('shlex').join(args),'exit':'TIMEOUT','stdout':(e.stdout or b'').decode(errors='replace'),'stderr':(e.stderr or b'').decode(errors='replace')}
 r['time']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 with DB.open('a') as f:f.write(json.dumps(r)+'\n')
 return r['stdout']
def git(path,*args):return run(['git','-C',path,*args])
if __name__=='__main__':
 run(['date','--iso-8601=seconds']);run(['hostname']);run(['id','-un'])
 run(['bash','-c','find /home -maxdepth 4 -type d -name .git 2>/dev/null | grep -i ej200'])
 run(['bash','-c','find /home /mnt /media /opt /srv -maxdepth 5 \\( -type d -name .git -o -type f -name .git \\) -print 2>/dev/null'],timeout=60)
 run(['ssh','-p','9022','-o','BatchMode=yes','-o','ConnectTimeout=5','rrios@localhost','echo OK'])
 for path in ['/home/rrios/ej200','/home/rrios/ej200_end']:
  for args in [('remote','-v'),('status','--porcelain=v2','-b'),('symbolic-ref','--short','HEAD'),('show','-s','--format=%H|%h|%cI|%aI|%s','HEAD'),('worktree','list','--porcelain'),('stash','list'),('submodule','status','--recursive'),('ls-files','--stage'),('lfs','version'),('lfs','ls-files'),('config','--get-regexp','^(fetch\.|remote\..*\.(prune|pruneTags)|maintenance\.|gc\.)')]:git(path,*args)
 path='/home/rrios/ej200'
 git(path,'-c','fetch.prune=false','-c','fetch.pruneTags=false','-c','gc.auto=0','-c','maintenance.auto=false','fetch','--all','--no-prune')
 git(path,'ls-remote','--heads','origin');git(path,'ls-remote','--tags','origin')
 git(path,'for-each-ref','--sort=-committerdate','refs/heads','refs/remotes','--format=%(refname:short)|%(committerdate:iso8601)|%(objectname)|%(upstream:short)|%(upstream:track)|%(authorname)|%(contents:subject)')
 git(path,'branch','--merged','main');git(path,'branch','--no-merged','main');git(path,'branch','-a','--contains','8349041')
 branches=git(path,'for-each-ref','--format=%(refname:short)','refs/heads').splitlines()
 for b in branches:
  for args in [('rev-list','--left-right','--count',f'main...{b}'),('log','--oneline',f'main..{b}'),('diff','--stat',f'main...{b}'),('ls-tree','-r',b),('grep','-n','-I','-E',r'polishedbackpainted|dielectric_metal|REFLECTIVITY|0\.98|0\.95|514\.9|1300|69\.2|701\.3|Rvikuiti',b,'--')]:git(path,*args)
 git(path,'tag','--list','--sort=-creatordate','--format=%(creatordate:short) %(refname:short) %(objectname:short) %(subject)')
 git(path,'for-each-ref','refs/tags','--format=%(refname:short)|%(objectname)|%(*objectname)|%(contents)')
 for tag in git(path,'tag','--list').splitlines():
  git(path,'branch','-a','--contains',tag)
 git(path,'show','--stat','8349041')
 print('Collected',len(branches),'branches. Evidence:',DB)
