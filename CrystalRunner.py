from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
import CrystalWriter

####################################################

class CrystalRunner:
  _name_="RunCrystal"
  def __init__(self, submitter,cwriter):
    self._submitter = submitter
    self._cwriter=cwriter
#-------------------------------------------------      
  def run(self,job_record):
    with open("autogen.d12",'w') as f:
      f.write(self._cwriter.crystal_input())
      
    self._submitter.execute(job_record, ['autogen.d12'],
          'autogen.d12', 'autogen.d12.o',self._name_)
    return 'running'
#-------------------------------------------------      
  def output(self):
    """ Collect results from output."""
    out={}
    if os.path.isfile('autogen.d12.o'):
      f = open('autogen.d12.o', 'r')
      lines = f.readlines()
      for li,line in enumerate(lines):
        if 'SCF ENDED' in line:
          out['total_energy']=float(line.split()[8])    
        elif 'TOTAL ATOMIC SPINS' in line:
          moms = []
          shift = 1
          while "TTT" not in lines[li+shift]:
            moms += map(float,lines[li+shift].split())
            shift += 1
          out['mag_moments']=moms
      
    return out
#-------------------------------------------------      
  # This can be made more efficient if it's a problem: searches whole file for
  # each query.
  def check_outputfile(self,outfilename,acceptable_scf=10.0):
    """ Check output file. 

    Current return values:
    no_record, not_started, ok, too_many_cycles, finished (fall-back),
    scf_fail, not_enough_decrease, divergence, not_finished
    """
    if os.path.isfile("autogen.d12.o"):
      outf = open("autogen.d12.o",'r')
    else:
      return "not_started"

    outlines = outf.read().split('\n')
    reslines = [line for line in outlines if "ENDED" in line]

    if len(reslines) > 0:
      if "CONVERGENCE" in reslines[0]:
        return "ok"
      elif "TOO MANY CYCLES" in reslines[0]:
        print("CrystalRunner: Too many cycles.")
        return "too_many_cycles"
      else: # What else can happen?
        print("CrystalRunner: Finished, but unknown state.")
        return "finished"
      
    detots = [float(line.split()[5]) for line in outlines if "DETOT" in line]
    if len(detots) == 0:
      print("CrystalRunner: Last run completed no cycles.")
      return "scf_fail"

    detots_net = sum(detots[1:])
    if detots_net > acceptable_scf:
      print("CrystalRunner: Last run performed poorly.")
      return "not_enough_decrease"

    etots = [float(line.split()[3]) for line in outlines if "DETOT" in line]
    if etots[-1] > 0:
      print("CrystalRunner: Energy divergence.")
      return "divergence"
    
    print("CrystalRunner: Not finished.")
    return "not_finished"
  
#-------------------------------------------------      
  # Diagnose routines basically decide 'not_finished' or 'failed'
  def stubborn_diagnose(self,status):
    if status in ['too_many_cycles','not_finished']:
      return 'not_finished'
    else:
      return 'failed'
#-------------------------------------------------      
  def optimistic_diagnose(self,status):
    if status == 'not_finished':
      return 'not_finished'
    else:
      return 'failed'
#-------------------------------------------------      
  def conservative_diagnose(self,status):
    return 'failed'
#-------------------------------------------------      
  def check_status(self,job_record):
    """ Decide status of job (in queue or otherwise). """
    outfilename="autogen.d12.o"
    diagnose_options = {
        'stubborn':self.stubborn_diagnose,
        'optimistic':self.optimistic_diagnose,
        'conservative':self.conservative_diagnose
      }

    status=self._submitter.status(job_record,self._name_)
    if 'running' in status:
      return 'running'

    self._submitter.transfer_output(job_record, ['autogen.d12.o', 'fort.9'])
    status=self.check_outputfile(outfilename)
    print("status",status)
    return status
    
#-------------------------------------------------      
  def resume(self,job_record,maxresume=5):
    """ Continue a crystal run using GUESSP."""
    trynum = 0
    while os.path.isfile("%d.autogen.d12.o"%trynum):
      trynum += 1
      if trynum > maxresume:
        print("Not resuming because resume limit reached ({}>{}).".format(
          trynum,maxresume))
        return 'failed'
    for filename in ["autogen.d12","autogen.d12.o","fort.79"]:
      shutil.copy(filename,"%d.%s"%(trynum,filename))
    shutil.copy("fort.79","fort.20")
    self._cwriter.set_restart(True)
    return self.run(job_record)
