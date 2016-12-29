from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
from submitter import LocalSubmitter

####################################################

class LocalCrystalRunner(LocalSubmitter):
  """ Runs a crystal job defined by CrystalWriter. """
  _name_='LocalCrystalRunner'
  def __init__(self, BIN='~/bin/'):
    self.BIN=BIN
    self.np=1
    self.nn=1
    self.jobname='ag_crystal'
    self._queueid=None
  #-------------------------------------------------      
  def check_status(self):
    # This is only for a queue, so just never return 'running'.
    return 'ok'

  #-------------------------------------------------      
  def run(self,crysinpfn,crysoutfn):
    """ Submits executibles using _qsub. """
    
    exe = self.BIN+"crystal < %s"%crysinpfn

    prep_commands=["cp %s INPUT"%crysinpfn]
    # Not needed for nonparallel.
    #final_commands = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]
    final_commands = []

    loc = os.getcwd()

    qids=self._qsub(exe,prep_commands,final_commands,self.jobname,crysoutfn,loc)
    self._queueid=qids


####################################################

class CrystalRunnerPBS:
  _name_='CrystalRunnerPBS'
  def __init__(self,BIN='~/bin/',
                    queue='batch',
                    walltime='12:00:00',
                    prefix="",#for example, load modules
                    postfix="rm fort.*.pe*"
                    ):
    self.BIN=BIN
    self.np=1
    self.nn=1
    self.jobname='CrystalRunnerPBS'
    self.queue=queue
    self.walltime=walltime
    self.prefix=prefix
    self.postfix=postfix
    self.queueid=None

  #-------------------------------------
  def check_status(self):
    try:
      qstat = sub.check_output(
          "qstat %s"%self.queueid, stderr=sub.STDOUT, shell=True
        ).decode().split('\n')[-2].split()[4]
    except sub.CalledProcessError:
      return "unknown"
    if qstat == "R" or qstat == "Q":
      return "running"
    if qstat == "C" or qstat == "E":
      return "finished"
    return status
  #-------------------------------------

  def run(self,crysinpfn,crysoutfn):
    exe=self.BIN+"Pcrystal"
    jobout=crysinpfn+".jobout"
    np_tot=self.np*self.nn
    qsub="#PBS -q %s \n"%self.queue +\
         "#PBS -l nodes=%i:ppn=%i\n"%(self.nn,self.np) +\
         "#PBS -l walltime=%s\n"%self.walltime +\
         "#PBS -j oe \n" +\
         "#PBS -N %s \n"%self.jobname +\
         "#PBS -o %s \n"%jobout +\
         self.prefix+"\n" +\
         "cd ${PBS_O_WORKDIR}\n" +\
         "cp %s INPUT\n"%(crysinpfn) +\
         "mpirun -np %i %s > %s\n"%(np_tot,exe,crysoutfn) +\
         self.postfix
    qsubfile=crysinpfn+".qsub"
    with open(qsubfile,'w') as f:
      f.write(qsub)
    result = sub.check_output("qsub %s"%(qsubfile),shell=True)
    self.queueid = result.decode().split()[0]
    print("Submitted as %s"%self.queueid)
    



