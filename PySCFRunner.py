from __future__ import print_function
import os
import sys
import numpy as np
import subprocess as sub
import shutil
import submitter
from submitter import LocalSubmitter

class LocalPySCFRunner(LocalSubmitter):
  def __init__(self, BIN=''):
    self.BIN=BIN
    self.np=1
    self.nn=1
    self.jobname='ag_pyscf'
    self._queueid=None
  #-------------------------------------------------      
  def check_status(self):
    # This is only for a queue, so just never return 'running'.
    return 'ok'

  #-------------------------------------------------      
  def run(self,inpfn,outfn):
    # TODO This code organization doesn't match CRYSTAL.
    """ Submits executibles using _qsub. """
    
    exe = self.BIN+"python3  %s"%inpfn[0]
    prep_commands = []
    final_commands = []

    loc = os.getcwd()

    qids=self._qsub(exe,prep_commands,final_commands,self.jobname,outfn[0],loc)
    self._queueid=qids

####################################################

class PySCFRunnerPBS:
  _name_='PySCFRunnerPBS'
  def __init__(self,queue='batch',
                    walltime='12:00:00',
                    np=1,
                    nn=1,
                    jobname=os.getcwd().split('/')[-1]+'_pyscf',
                    prefix="",#for example, load modules
                    postfix=""#for example, remove tmp files.
                    ):
    self.np=np
    if nn!=1: raise NotImplementedError
    self.nn=nn
    self.jobname=jobname
    self.queue=queue
    self.walltime=walltime
    self.prefix=prefix
    self.postfix=postfix
    self.queueid=None

  #-------------------------------------
  def check_status(self):
    return submitter.check_PBS_status(self.queueid)
  #-------------------------------------

  def run(self,inpfns,outfns):
    #just running in serial for now
    exe="python3"
    np_tot=self.np*self.nn
    #"mpirun -np %i %s > %s\n"%(np_tot,exe,crysoutfn) +
    
    for pyscf_driver,pyscf_out in zip(inpfns,outfns):
      jobout=pyscf_driver+".jobout"
      qsublines=[
           "#PBS -q %s"%self.queue,
           "#PBS -l nodes=%i:ppn=%i"%(self.nn,self.np),
           "#PBS -l walltime=%s"%self.walltime,
           "#PBS -j oe",
           "#PBS -N %s"%self.jobname,
           "#PBS -o %s"%jobout,
           self.prefix,
           "export OMP_NUM_THREADS=%i"%(self.np),
           "cd ${PBS_O_WORKDIR}",
           "export PYTHONPATH=%s"%(':'.join(sys.path)),
           "%s %s > %s \n"%(exe,pyscf_driver,pyscf_out),
           self.postfix
         ]
      qsubfile=pyscf_driver+".qsub"
      with open(qsubfile,'w') as f:
        f.write('\n'.join(qsublines))
      result = sub.check_output("qsub %s"%(qsubfile),shell=True)
      self.queueid = result.decode().split()[0]
      print("Submitted as %s"%self.queueid)
