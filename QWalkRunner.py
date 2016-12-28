from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
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
