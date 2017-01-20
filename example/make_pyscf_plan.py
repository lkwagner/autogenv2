import sys
sys.path.append("../")

from JobEnsemble import JobEnsemble
import Recipes as job
import pickle as pkl
from copy import deepcopy
import os
from QWalkRunner import QWalkRunnerPBS

# TODO
#dft_opts={
#    "xyz":"H 0. 0. 0.; H 0. 0. 0.7"
#  }

hf_opts={
    'xyz':"N 0. 0. 0.; N 0. 0. 2.5",
    'method':'ROHF',
  }

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
    'total_nstep':512
  }

dmc_opts={
    'timesteps':[.03],
    'extra_observables':[{
      'name':'average_derivative_dm',
      'nmo':8,
      'orbfile':'qw.orb',
      'basis':'qw.basis',
      'states':[3,4,5,6,7,8]
    }]
  }

test = JobEnsemble([
    job.PySCFQWalk('n2_hf',
       pyscf_opts=hf_opts,
       variance_opts=variance_opts,
       energy_opts=energy_opts,
       dmc_opts=dmc_opts,
       qwalkrunner=QWalkRunnerPBS(np=4) ),
    job.PySCFQWalk('n2_cas',
       pyscf_opts=cas_opts,
       variance_opts=variance_opts,
       energy_opts=energy_opts,
       dmc_opts=dmc_opts,
       qwalkrunner=QWalkRunnerPBS(np=4) )
    ]
  )

with open('plan.pickle','wb') as outf:
  pkl.dump(test,outf)
