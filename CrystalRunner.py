from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
from submitter import LocalTorqueSubmitter

####################################################
class LocalCrystalRunner(LocalTorqueSubmitter):
  """ Runs a crystal job defined by CrystalWriter. """
  _name_='LocalCrystalRunner'
  def __init__(self,nn=1,np=1,queue='batch',walltime='4:00:00',BIN='~/bin/'):
    self.BIN=BIN
    self.np=np
    self.nn=nn
    self.queue=queue
    self.walltime=walltime
    self.jobname='ag_crystal'
    self._queueid=None
  #-------------------------------------------------      
  def check_status(self):
    # This is only for a queue, so just never return 'running'.
    return 'ok'

  #-------------------------------------------------      
  def run(self,crysinpfn,crysoutfn):
    """ Submits executibles using _qsub. """
    
    exe = self.BIN+"Pcrystal"

    prep_commands=[]
    prep_commands+=[". ~/bin/add_intel"]
    prep_commands+=[". ~/bin/add_ompi"]
    prep_commands+=["cp %s INPUT"%crysinpfn]
    final_commands = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]

    loc = os.getcwd()

    qids=self._qsub(exe,prep_commands,final_commands,self.jobname,crysoutfn,loc)
    self._queueid=qids
