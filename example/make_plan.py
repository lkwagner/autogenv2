import sys
sys.path.append("../")

from JobEnsemble import JobEnsemble
import Job as job
import pickle as pkl
from copy import deepcopy
import os

dft_opts={
    'xml_name':os.getcwd()+'/../BFD_Library.xml',
    'basis_params':[0.2,0,3],
    'cutoff':0.0,
    'dftgrid':'LGRID',
    'spin_polarized':False
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

pbc_opts={
    'xml_name':os.getcwd()+'/../BFD_Library.xml',
    'basis_params':[0.4,1,3],
    'cutoff':0.2,
    'dftgrid':'LGRID',
    'spin_polarized':False,
    'kmesh':[2,2,2] 
  }

test = JobEnsemble([

    job.LocalCrystalQWalk('h2',open('h2.xyz','r').read(),
      crystal_opts=dft_opts,
      variance_opts=variance_opts,
      energy_opts=energy_opts,
      dmc_opts=dmc_opts,
      structtype='xyz'),

    job.LocalCrystalQWalk('n2',open('n2.xyz','r').read(),
      crystal_opts=dft_opts,
      variance_opts=variance_opts,  
      energy_opts=energy_opts,
      dmc_opts=dmc_opts,      
      structtype='xyz'),

    job.LocalCrystalQWalk('si',open('si.cif','r').read(),
      crystal_opts=pbc_opts,
      variance_opts=variance_opts,
      energy_opts=energy_opts,
      dmc_opts=dmc_opts,            
      structtype='cif')
    
  ])

with open('plan.pickle','wb') as outf:
  pkl.dump(test,outf)
