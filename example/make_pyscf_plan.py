import sys
sys.path.append("../")

from JobEnsemble import JobEnsemble
import Recipes as job
import pickle as pkl
from copy import deepcopy
import os
from QWalkRunner import LocalQWalkRunner

dft_opts={
    "xyz":"H 0. 0. 0.; H 0. 0. 0.7"
  }

dft_opts_n2={'xyz':"N 0. 0. 0.; N 0. 0. 2.5"
  }
    

variance_opts={
    'iterations':15
    }

energy_opts={
    'total_nstep':512
    }

dmc_opts={
    'timesteps':[.03]
    }

test = JobEnsemble([
#    job.PySCFQWalk('h2',"H 0. 0. 0.; H 0. 0. 0.7",
#      pyscf_opts=dft_opts,
#      variance_opts=variance_opts,
#      energy_opts=energy_opts,
#      dmc_opts=dmc_opts,
#      qwalkrunner=LocalQWalkRunner() ),
  job.PySCFQWalk('n2',"H 0. 0. 0.; H 0. 0. 0.7",
     pyscf_opts=dft_opts_n2,
     variance_opts=variance_opts,
     energy_opts=energy_opts,
     dmc_opts=dmc_opts,
     qwalkrunner=LocalQWalkRunner() )
  ]
  )

with open('plan.pickle','wb') as outf:
  pkl.dump(test,outf)
