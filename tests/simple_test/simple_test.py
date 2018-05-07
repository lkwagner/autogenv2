'''
Some simple tests to check that autogen runs without erroring out.
'''

from autogenv2 import manager,autopyscf,autorunner
from autopyscf import PySCFWriter
from manager import PySCFManager,CrystalManager
from autorunner import PySCFRunnerLocal,PySCFRunnerPBS
import sys

h2='\n'.join([
    'H 0.0 0.0 0.0',
    'H 0.0 0.0 0.74'
  ])
h2stretch='\n'.join([
    'H 0.0 0.0 0.0',
    'H 0.0 0.0 1.4'
  ])

def h2_tests():
  ''' Simple tests that check PySCF and queue interaction.'''

  # Most basic possible job.
  eqwriter=PySCFWriter({'xyz':h2})
  eqman=PySCFManager(
      name='scf',
      path='h2equil',
      writer=eqwriter,
      runner=PySCFRunnerLocal()
    )

  # Change some options and run with PBS.
  stwriter=PySCFWriter({
      'xyz':h2stretch,
      'method':'UKS',
      'dm_generator':autopyscf.dm_set_spins([1,-1],[]),
    })
  stman=PySCFManager(
      name='scf',
      path='h2stretch',
      writer=stwriter,
      runner=PySCFRunnerPBS(nn=1,np=2,ppath=sys.path),
    )

  return [eqman,stman]

def run_tests():
  ''' Choose which tests to run and execute `nextstep()`.'''
  jobs=[]
  jobs+=h2_tests()

  for job in jobs:
    job.nextstep()

if __name__=='__main__':
  run_tests()
