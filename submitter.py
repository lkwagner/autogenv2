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

    if stdout=="": stdout="stdout"
    if loc=="": loc=os.getcwd()
    if name=="": name=stdout
    header = []
    exeline = exe
    commands = header +  prep_commands + [exeline] + final_commands

    outstr = ""
    for c in commands:
      print(c)
      # Problem: this method doesn't allow you to watch its progress.
      #outstr+=sub.check_output(c,shell=True).decode()
      sub.call(c,stdout=open(stdout,'w'),shell=True)
   # with open(stdout,'w') as outf:
   #   outf.write(outstr)
    return []
      
