import sys
sys.path.append("../../")

from JobEnsemble import JobEnsemble
import Recipes as job
import pickle as pkl
from copy import deepcopy
import os
from PySCFRunner import PySCFRunnerPBS
from QWalkRunner import QWalkRunnerPBS
from PySCF import dm_set_spins

sicif=open('si.cif','r').read()

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


methods={'RKS':['pbe,pbe','b3lyp','pbe0','lda,vwn']}

joblist=[]
for m,func in methods.items():
  for f in func:
    pyscf_opts={
        'cif':sicif,
        'method':m,
        'dft':f,
        'kpts':[2,2,2],
        'gs':[4,4,4],
        'bfd_library':'../../BFD_Library.xml'}

    joblist.append(job.PySCFQWalk(
                                  'si'+m+f.replace(',',''),
                                  pyscf_opts=pyscf_opts,
                                  variance_opts=variance_opts,
                                  energy_opts=energy_opts,
                                  dmc_opts=dmc_opts,
                                  pyscfrunner=PySCFRunnerPBS(np=8),
                                  qwalkrunner=QWalkRunnerPBS(np=8) 
                                 )
                   )

test=JobEnsemble(joblist)
with open('plan.pickle','wb') as outf:
  pkl.dump(test,outf)
