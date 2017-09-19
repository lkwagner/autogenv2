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
                    np='allprocs',
                    ppath=None,
                    nn=1,
                    jobname=os.getcwd().split('/')[-1]+'_pyscf',
                    prefix=None,
                    postfix=None
                    ):
    self.np=np
    if nn!=1: raise NotImplementedError
    self.nn=nn
    self.jobname=jobname
    self.queue=queue
    self.walltime=walltime
    self.prefix=prefix
    self.postfix=postfix
    self.exelines=[]
    if ppath is None: self.ppath=[]
    else:             self.ppath=ppath
    if prefix is None: self.prefix=[]
    else:              self.prefix=prefix
    if postfix is None: self.postfix=[]
    else:               self.postfix=postfix
    self.queueid=None

  #-------------------------------------
  def check_status(self):
    return submitter.check_PBS_status(self.queueid)

  #-------------------------------------
  def run(self,exestr):
    ''' Accumulate executable commands.'''
    self.exelines.append(exestr)

  #-------------------------------------
  def script(self,scriptfile):
    ''' Dump accumulated commands into a script for another job to run.
    Returns true if the runner had lines to actually execute.'''

    if len(self.exelines)==0:
      return False

    # Prepend mpi specs.
    exelines=[]
    for line in self.exelines:
      exelines.append('mpirun -n {} {}'.format(
        self.nn*self.np,line))

    with open(scriptfile,'w') as outf:
      outf.write('\n'.join(self.prefix + exelines + self.postfix))

    # Remove exelines so the runner is ready for the next go.
    self.exelines=[]
    return True

  #-------------------------------------

  def submit(self,jobname=None):
    if len(self.exelines)==0: 
      print("Nothing to run.")
      return
    if jobname is None: jobname=self.jobname
    
    jobout=jobname+".jobout"
    qsublines=[
         "#PBS -q %s"%self.queue,
       ]
    if self.np == 'allprocs' :
      qsublines+=[
           "#PBS -l nodes=%i,flags=allprocs"%self.nn,
         ]
    else:
      qsublines+=[
           "#PBS -l nodes=%i:ppn=%i"%(self.nn,self.np),
         ]
    qsublines+=[
         "#PBS -l walltime=%s"%self.walltime,
         "#PBS -j oe",
         "#PBS -N %s"%self.jobname,
         "#PBS -o %s"%jobout,
         "cd ${PBS_O_WORKDIR}",
         "export OMP_NUM_THREADS=%d"%(self.nn*self.np),
         "export PYTHONPATH=%s"%(':'.join(self.ppath))
       ] + self.prefix + self.exelines + self.postfix
    qsubfile="qsub.in"
    with open(qsubfile,'w') as f:
      f.write('\n'.join(qsublines))
    result = sub.check_output("qsub %s"%(qsubfile),shell=True)
    self.queueid = result.decode().split()[0]
    print("Submitted as %s"%self.queueid)

    # Clear out the lines to set up for the next job.
    self.exelines=[]
