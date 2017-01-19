from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
import submitter

from submitter import LocalSubmitter

####################################################

class LocalQWalkRunner(LocalSubmitter):
  """ Runs a crystal job defined by CrystalWriter. """
  _name_='LocalCrystalRunner'
  def __init__(self, BIN='~/bin/'):
    self.BIN=BIN
    self.np=1
    self.nn=1
    self.jobname='ag_qwalk'
    self._queueid=None
  #-------------------------------------------------      
  def check_status(self):
    # This is only for a queue, so just never return 'running'.
    return 'ok'

  #-------------------------------------------------      
  def run(self,qwinpfns,qwoutfns):
    """ Submits executibles using _qsub. """
    for qwinpfn,qwoutfn in zip(qwinpfns,qwoutfns):
      exe = self.BIN+"qwalk %s"%qwinpfn
      prep_commands=[]
      final_commands = []
      loc = os.getcwd()
      qids=self._qsub(exe,prep_commands,final_commands,self.jobname,qwoutfn,loc)
      self._queueid=qids


####################################################

class QWalkRunnerPBS:
  _name_='QWalkRunnerPBS'
  def __init__(self,exe='~/bin/qwalk',
                    queue='batch',
                    walltime='12:00:00',
                    prefix="",#for example, load modules
                    postfix="",#for any clean-up.
                    np=1,nn=1
                    ):
    self.exe=exe
    self.np=np
    self.nn=nn
    self.jobname='QWalkRunnerPBS'
    self.queue=queue
    self.walltime=walltime
    self.prefix=prefix
    self.postfix=postfix
    self.queueid=[]

  #-------------------------------------
  def check_status(self):
    return submitter.check_PBS_stati(self.queueid)
  #-------------------------------------

  def run(self,qwinps,qwouts):
    #Just supporting one run for the moment
    qwinp=qwinps[0]
    qwout=qwouts[0]
    for qwinp,qwout in zip(qwinps,qwouts):
      jobout=qwinp+".jobout"
      np_tot=self.np*self.nn
    
      qsub="#PBS -q %s \n"%self.queue +\
         "#PBS -l nodes=%i:ppn=%i\n"%(self.nn,self.np) +\
         "#PBS -l walltime=%s\n"%self.walltime +\
         "#PBS -j oe \n" +\
         "#PBS -N %s \n"%self.jobname +\
         "#PBS -o %s \n"%jobout +\
         self.prefix+"\n" +\
         "cd ${PBS_O_WORKDIR}\n" +\
         "cp %s INPUT\n"%(qwinp) +\
         "mpirun -np %i %s %s \n"%(np_tot,self.exe,qwinp) +\
         self.postfix
      qsubfile=qwinp+".qsub"
      with open(qsubfile,'w') as f:
        f.write(qsub)
      result = sub.check_output("qsub %s"%(qsubfile),shell=True)
      self.queueid.append(result.decode().split()[0].split('.')[0])
      print("Submitted as %s"%self.queueid)
      
