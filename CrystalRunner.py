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
  def __init__(self, BIN='~/bin/',cryinpfn='autogen.d12'):
    self.cryinpfn = cryinpfn
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
  def run(self):
    """ Submits executibles using _qsub. """
    
    exe = self.BIN+"crystal < %s"%self.cryinpfn

    prep_commands=["cp %s INPUT"%self.cryinpfn]
    # Not needed for nonparallel.
    #final_commands = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]
    final_commands = []

    outfn = self.cryinpfn+".o"
    loc = os.getcwd()

    qids=self._qsub(exe,prep_commands,final_commands,self.jobname,outfn,loc)
    self._queueid=qids
