import subprocess as sub
import importlib
import os
import shutil
import sys
import time

#####################################################################################
class LocalSubmitter:
  """Abstract submission class. Defines interaction with queuing system (for
  nonqueuing defines how to execute jobs). Child classes must define:
  __init__: 
     can have any parameters, but should set up any variables like queue time and number of
     processor cores
  _submit_job(self,
              inpfns : list of input filenames (list of strings)
              outfn  : where stdout should be sent (string)
              jobname : job name for the queue (string)
              loc :  directory where job should be run (string)
            )
            returns a list of queue ids (list of strings)
  """
  #-------------------------------------------------------
  def _qsub(self,exe,prep_commands=[],final_commands=[],
      name="",stdout="runstdout",loc=""):
    """ Helper function for executable submitters. 
    Should work in most cases to simplify code."""

    if loc=="": loc=os.getcwd()
    if name=="": name=loc
    header = []
    exeline = exe
    commands = header +  prep_commands + [exeline] + final_commands

    outstr = ""
    for c in commands:
      print(c)
      sub.call(c,stdout=open(stdout,'w'),shell=True)
    return []

#-------------------------------------------------------

def check_PBS_status(queueid):
  """Utility function to determine the status of a PBS job."""
  try:
    qstat = sub.check_output(
        "qstat %s"%queueid, stderr=sub.STDOUT, shell=True
      ).decode().split('\n')[-2].split()[4]
  except sub.CalledProcessError:
    return "unknown"
  if qstat == "R" or qstat == "Q":
    return "running"
  if qstat == "C" or qstat == "E":
    return "finished"
  return 'unknown'

#-------------------------------------------------------
def check_PBS_stati(queueid):
  """Utility function to determine the status of a set PBS job.
  Can be done with one qstat call which can improve speed."""
  try:
    qstat = sub.check_output(
        "qstat ", stderr=sub.STDOUT, shell=True
      ).decode()
  except sub.CalledProcessError:
    return "unknown"
  print(qstat)
  qstat=qstat.split('\n')
  for qid in queueid:
    for line in qstat:
      spl=line.split()
      if qid in line and len(spl) > 4:
        stat=line.split()[4]
        if stat == "R" or stat == "Q":
          return "running"
  return 'unknown'

