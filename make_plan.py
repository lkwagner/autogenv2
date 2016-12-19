from JobEnsemble import JobEnsemble
import Job as job
import pickle as pkl
from copy import deepcopy

dft_opts={
    'xml_name':'/home/brian/programs/autogenv2/BFD_Library.xml',
    'basis_params':[0.2,0,3],
    'cutoff':0.0,
    'dftgrid':'LGRID',
    'spin_polarized':False
  }

test = JobEnsemble([
    job.LocalCrystalDFT('h2',open('h2.xyz','r').read(),dft_opts,'xyz'),
    job.LocalCrystalDFT('n2',open('n2.xyz','r').read(),dft_opts,'xyz')
  ])

with open('plan.pickle','wb') as outf:
  pkl.dump(test,outf)
