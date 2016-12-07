from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
from submission_tools import LocalSubmitter

####################################################

class LocalCrystalRunner(LocalSubmitter):
  """ Runs a crystal job defined by CrystalWriter. """
  _name_='LocalCrystalRunner'
  def __init__(self, BIN='~/bin/',cryinpfn='autogen.d12'):
    self.cryinpfn = cryinpfn
    self.BIN=BIN
    self.np=1
    self.nn=1
  #-------------------------------------------------      
  def check_status(self):
    # This is only for a queue, so just never return 'running'.
    return 'ok'

  #-------------------------------------------------      
  def run(self,job_record):
    """ Calls submitter with appropriate filenames. """
    self.execute(
        job_record,
        [self.cryinpfn],
        [self.cryinpfn],
        self.cryinpfn+'.o',
        'crystalrun'
      )

    # Local runner will not proceed until finished.
    return 'ok'

  #-------------------------------------------------      
  def _submit_job(self,inpfns,outfn="stdout",jobname="",loc=""):
    """ Sets up the execution line and calls _qsub. """

    assert len(inpfns)==1, "Bundle not implmented yet."
    inpfn=inpfns[0]
    
    exe = self.BIN+"crystal < %s"%inpfn

    prep_commands=["cp %s INPUT"%inpfn]
    # Not needed for nonparallel.
    #final_commands = ["rm *.pe[0-9]","rm *.pe[0-9][0-9]"]
    final_commands = []

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    qid = self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc)
    return qid
