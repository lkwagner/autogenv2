import subprocess as sub
import importlib
import os
import shutil
import sys
import time

#####################################################################################
class LocalTorqueSubmitter:
  """Abstract submission class. Defines interaction with standard Torque
  queuing systems that exist on many Linux computer and Campus Cluster.  """
  #-------------------------------------------------------
  def _qsub(self,exe,prep_commands=[],final_commands=[],
      name="",stdout="",loc=""):
    """ Submit a job through qsub command."""

    if stdout=="": stdout="stdout"
    if loc=="": loc=os.getcwd()
    if name=="": name=stdout
    header = []
    header.append("#!/bin/bash")
    if self.np=="allprocs": 
      header.append("#PBS -l nodes=%d,flags=allprocs"%self.nn)
    else:
      header.append("#PBS -l nodes=%d:ppn=%d"%(self.nn,self.np))
    header += [
      "#PBS -q %s"%self.queue,
      "#PBS -l walltime=%s"%self.walltime,
      "#PBS -j oe",
      "#PBS -m n",
      "#PBS -N %s"%name,
      "#PBS -o {0}".format(loc+"/qsub.out"),
    ]
    if self.np=="allprocs":
      exeline = "mpirun %s &> %s"%(exe, stdout)
    elif self.nn*self.np > 1:
      exeline = "mpirun -n %d %s &> %s"%(self.nn*self.np, exe, stdout)
    else:
      exeline = "%s &> %s"%(exe, stdout)
    commands = header + ["cd %s"%loc] + prep_commands + [exeline] + final_commands
    outstr = '\n'.join(commands)
    with open(loc+"/qsub.in",'w') as qsin:
      qsin.write(outstr)
    result = sub.check_output("qsub %s"%(loc+"/qsub.in"),shell=True)
    qid = result.decode().split()[0]
    print("Submitted as %s"%qid)
    return [qid]

#####################################################################################
class LocalCommandLineSubmitter:
  """Abstract submission class. Submits to command line, and is mostly for
  debugging purposes. """
  #-------------------------------------------------------
  def _qsub(self,exe,prep_commands=[],final_commands=[],
      name="",stdout="",loc=""):
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
      # Problem: this method doesn't allow you to watch it's progress.
      outstr+=sub.check_output(c,shell=True).decode()
    with open(stdout,'w') as outf:
      outf.write(outstr)
    return []
