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

pbc_opts={
    'xml_name':os.getcwd()+'/../BFD_Library.xml',
    'basis_params':[0.4,1,3],
    'cutoff':0.2,
    'dftgrid':'LGRID',
    'spin_polarized':False
  }

test = JobEnsemble([
    job.LocalCrystalDFT('h2',open('h2.xyz','r').read(),dft_opts,'xyz'),
    job.LocalCrystalDFT('n2',open('n2.xyz','r').read(),dft_opts,'xyz'),
    job.LocalCrystalDFT('si',open('si.cif','r').read(),pbc_opts,'cif')
    
  ])

with open('plan.pickle','wb') as outf:
  pkl.dump(test,outf)
