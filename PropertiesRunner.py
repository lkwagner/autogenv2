from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
from submitter import LocalSubmitter

class LocalPropertiesRunner(LocalSubmitter):
  """ Runs a crystal properties job defined by CrystalWriter. """
  _name_='LocalPropertiesRunner'
  def __init__(self, BIN='~/bin/'):
    self.BIN=BIN
    self.np=1
    self.nn=1
    self.jobname='ag_properties'
    self._queueid=None
  #-------------------------------------------------      
  def check_status(self):
    # This is only for a queue, so just never return 'running'.
    return 'ok'

  #-------------------------------------------------      
  def run(self,propinpfn):
    """ Submits executibles using _qsub. """
    
    exe = self.BIN+"properties < %s"%propinpfn

    # Not needed for nonparallel.
    #final_commands = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]
    prep_commands = []
    final_commands = []

    outfn = propinpfn+".o"
    loc = os.getcwd()

    qids=self._qsub(exe,prep_commands,final_commands,self.jobname,outfn,loc)
    self._queueid=qids
