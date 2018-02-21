import sys
sys.path.append("../../")

from JobEnsemble import JobEnsemble
import Recipes as job
import pickle as pkl
from copy import deepcopy
import os
from QWalkRunner import QWalkRunnerPBS
from PySCF import dm_set_spins

xyz="N 0. 0. 0.; N 0. 0. 2.5"

cas_opts={
    'xyz':"N 0. 0. 0.; N 0. 0. 2.5",
    'method':'ROHF',
    'postHF':True,
    'cas': {
        'ncore':2, # s bonding and antibonding.
        'ncas':6,
        'nelec':(3,3), 
        'tol': 0.02,
        'method': 'CASSCF'
      }
  }

variance_opts={
    'iterations':15
  }

energy_opts={
    'total_nstep':8192
  }

dmc_opts={
    'timesteps':[.03],
     'extra_observables':[ {'name':'region_fluctuation','maxn':6}] 
  }


methods={'UKS':['pbe,pbe','b3lyp','pbe0','lda,vwn'],
         'RKS':['pbe,pbe','b3lyp','pbe0','lda,vwn'],
         'ROHF':[''],
         'UHF':['']
        }

joblist=[]
for m,func in methods.items():
  for f in func:
    pyscf_opts={'xyz':xyz,'method':m,'dft':f}
    if m[0]=='U':
      pyscf_opts['dm_generator']=dm_set_spins(
          atomspins=[-1,1],
          double_occ={'N':[0]}
        )

    joblist.append(job.PySCFQWalk(
                                  'n2'+m+f.replace(',',''),
                                  pyscf_opts=pyscf_opts,
                                  variance_opts=variance_opts,
                                  energy_opts=energy_opts,
                                  dmc_opts=dmc_opts,
                                  qwalkrunner=QWalkRunnerPBS(np=4) 
                                 )
                   )

joblist.append(
    job.PySCFQWalk('n2_cas',
       pyscf_opts=cas_opts,
       variance_opts=variance_opts,
       energy_opts=energy_opts,
       dmc_opts=dmc_opts,
       qwalkrunner=QWalkRunnerPBS(np=4) )
  )

test=JobEnsemble(joblist)
with open('plan.pickle','wb') as outf:
  pkl.dump(test,outf)
