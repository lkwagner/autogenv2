from __future__ import print_function
import os
import sys
import numpy as np
import subprocess as sub
import shutil
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
