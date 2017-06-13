import sys
sys.path.append("../")

from JobEnsemble import JobEnsemble
import Recipes as job
import pickle as pkl
from copy import deepcopy
import os
from PySCFRunner import PySCFRunnerPBS
from QWalkRunner import QWalkRunnerPBS
from PySCF import dm_set_spins,dm_from_chkfile
from time import sleep

def make_quick_test():
  ''' Quick n2 tests can check basic logic for runs. '''
  cas_opts={
      'xyz':"N 0. 0. 0.; N 0. 0. 2.5",
      'method':'RHF',
      'postHF':True,
      'cas': {
          'ncore':2, # s bonding and antibonding.
          'ncas':6,
          'nelec':(3,3), 
          'tol': 0.02,
          'method': 'CASSCF'
        }
    }

  spin_opts={
      'xyz':"N 0. 0. 0.; N 0. 0. 2.5",
      'method':'UHF'
    }
  spin_opts['dm_generator']=dm_set_spins(
      atomspins=[-1,1],
      double_occ={'N':[0]}
    )

  restart_opts={
      'xyz':"N 0. 0. 0.; N 0. 0. 2.5",
      'method':'RHF',
      'max_cycle':5
    }

  variance_opts={
      'iterations':5
    }

  energy_opts={
      'total_nstep':100
    }

  dmc_opts={
      'timesteps':[.05],
      'nblock':5,
      'savetrace':True
    }
  post_opts={'extra_observables':
      [
        {'name':'region_fluctuation','maxn':6}
      ]
    }
  jobs=[
      job.PySCFQWalk('quick_cas_test',
        pyscf_opts=cas_opts,
        variance_opts=variance_opts,
        energy_opts=energy_opts,
        dmc_opts=dmc_opts,
        post_opts=post_opts,
        pyscfrunner=PySCFRunnerPBS(np=2),
        qwalkrunner=QWalkRunnerPBS(np=4) ),
      job.PySCFQWalk('quick_spin_test',
        pyscf_opts=spin_opts,
        variance_opts=variance_opts,
        energy_opts=energy_opts,
        dmc_opts=dmc_opts,
        post_opts=post_opts,
        pyscfrunner=PySCFRunnerPBS(np=2),
        qwalkrunner=QWalkRunnerPBS(np=4) ),
      job.PySCFQWalk('quick_restart_test',
        pyscf_opts=restart_opts,
        variance_opts=variance_opts,
        energy_opts=energy_opts,
        dmc_opts=dmc_opts,
        post_opts=post_opts,
        pyscfrunner=PySCFRunnerPBS(np=2),
        qwalkrunner=QWalkRunnerPBS(np=4) )
    ]
  return jobs

def make_pbc_test():
  ''' Test pbc functionality that quick test misses. '''
  sicif=open('si.cif','r').read()

  pyscf_opts={
      'cif':sicif,
      #'max_cycle':3, # tests restart.
      'method':'RKS',
      'dft':'pbe,pbe',
      'kpts':[2,2,2],
      'gs':[4,4,4],
      'bfd_library':'../BFD_Library.xml'
    }

  variance_opts={
      'iterations':5
    }

  energy_opts={
      'total_nstep':100
    }

  dmc_opts={
      'timesteps':[.03],
      'extra_observables':[ {'name':'region_fluctuation','maxn':6}],
      'savetrace':True
    }
  post_opts={'extra_observables':
      [
        {'name':'region_fluctuation','maxn':6}
      ]
    }

  jobs=[
      job.PySCFQWalk('pbc_test',
        pyscf_opts=pyscf_opts,
        variance_opts=variance_opts,
        energy_opts=energy_opts,
        dmc_opts=dmc_opts,
        post_opts=post_opts,
        qwalkrunner=QWalkRunnerPBS(np=4) ),
    ]
  return jobs

if __name__=='__main__':

  jobs=make_quick_test()+make_pbc_test()
  # For faster turn-around:
  #jobs=make_quick_test()
  ensemble=JobEnsemble(jobs)

  # This continually resubmits to check submission logic.
  for resub in range(1000):
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Calling nextstep (ctrl-c to exit)")
    ensemble.nextstep()
    print()
    sleep(2)

