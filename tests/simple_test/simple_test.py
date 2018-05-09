'''
Some simple tests to check that autogen runs without erroring out.

Run repeatedly to check if runner is also checking queued items correctly.
'''

from autogenv2 import manager,autopyscf,autorunner,crystal
from crystal import CrystalWriter
from autopyscf import PySCFWriter,PySCFPBCWriter
from manager import PySCFManager,CrystalManager,QWalkManager
from autorunner import PySCFRunnerLocal,PySCFRunnerPBS,RunnerPBS
from variance import VarianceWriter,VarianceReader
from linear import LinearWriter,LinearReader
from dmc import DMCWriter,DMCReader
from trialfunc import SlaterJastrow
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

def si_crystal_test():
  ''' Simple tests that check PBC is working Crystal, and that QMC can be performed on the result.'''
  jobs=[]

  cwriter=CrystalWriter({
      'xml_name':'../BFD_Library.xml',
      'kmesh':(4,4,4),
    })
  cwriter.set_struct_fromcif(open('si.cif','r').read(),primitive=True)
  cwriter.set_options({'symmetry':True})

  cman=CrystalManager(
      name='crys',
      path='sicrys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,
          walltime='0:10:00'
        )
    )
  jobs.append(cman)

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=(0,0,0))
    )

  jobs.append(var)

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,queue='secondary',walltime='1:00:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=(0,0,0))
    )

  jobs.append(lin)

  for kpt in cman.qwfiles['kpts']:
    dmc=QWalkManager(
        name='dmc_%d%d%d'%kpt,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,queue='secondary',walltime='1:00:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kpt)
      )
    jobs.append(dmc)

  return jobs

def si_pyscf_test():
  ''' Simple tests that check PBC is working Crystal, and that QMC can be performed on the result.'''
  jobs=[]

  pwriter=PySCFPBCWriter({
      'cif':open('si.cif','r').read()
    })
  pman=PySCFManager(
      name='scf',
      path='sipyscf',
      writer=pwriter,
      runner=PySCFRunnerPBS(
          queue='secondary',
          nn=1,
          walltime='0:20:00',
          ppath=sys.path
        )
    )
  jobs.append(pman)

  var=QWalkManager(
      name='var',
      path=pman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(pman,kpoint=0)
    )
  jobs.append(var)

  lin=QWalkManager(
      name='lin',
      path=pman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,queue='secondary',walltime='0:30:00'
        ),
      trialfunc=SlaterJastrow(slatman=pman,jastman=var,kpoint=0)
    )
  jobs.append(lin)

  return jobs

def run_tests():
  ''' Choose which tests to run and execute `nextstep()`.'''
  jobs=[]
  #jobs+=h2_tests()
  jobs+=si_crystal_test()
  jobs+=si_pyscf_test()

  for job in jobs:
    job.nextstep()

if __name__=='__main__':
  run_tests()
