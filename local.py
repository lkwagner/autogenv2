# Submission scripts for running on your computer.
import subprocess as sub
from submission_tools import LocalSubmitter
import os

##########################################################
class LocalRunner(LocalSubmitter):
  """Abstract submission class. Defines interaction with the queuing system."""
#-------------------------------------------------------  
  def __init__(self,BIN="~/bin/"):
    """ Initialize a submitter object. 

    This part should contain all parts of calculation that user should care
    about, and all parts that are common to calculations of a certain type (all
    crystal SCF calculations, for instance)."""
    self.BIN=BIN
    self.np=1
    self.nn=1

### Below is legacy code for now. ###
###############################################################
class LocalProperties(LocalRunner):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfn,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    exe = self.BIN+"properties < %s"%inpfn
    prep_commands = []
    final_commands = []

    if self.nn != 1 or self.np != 1:
      print("Refusing to run properties in parallel!")
      self.nn = 1
      self.np = 1

    if jobname == "":
      jobname = outfn
    if loc == "":
      loc = os.getcwd()

    qid = self._qsub(exe,prep_commands,final_commands,jobname,outfn,loc)
    return qid

#####################################################################
class LocalQWalk(LocalRunner):
  """Fully defined submission class. Defines interaction with specific
  program to be run."""
  def _submit_job(self,inpfns,outfn="stdout",jobname="",loc=""):
    """ Submit a specific job to the queue. 
    
    Should not interact with user, and should receive only information specific
    to instance of a job run."""
    prep_commands=[]
    final_commands=[]
    qid=[]
    
    for f in inpfns:
      exe = self.BIN+"qwalk %s"%f
      fulloutfn = f+outfn # Ugly but more functional.

      if jobname == "":
        jobname = outfn
      if loc == "":
        loc = os.getcwd()

      qid+=self._qsub(exe,prep_commands,final_commands,jobname,fulloutfn,loc)
    return qid

#####################################################################
