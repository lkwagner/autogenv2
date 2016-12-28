import os
import pickle as pkl
from copy import deepcopy

class JobEnsemble:
  """ List of Jobs that will be run in different directories. 

  * Manages directories and calls Job methods in appropriate places.
  * Performs consistency check between Job objects recorded in directories and
    plan object that exists outside directories."""
  def __init__(self,joblist):
    # Could have this read in a plan pickle from a filename instead.
    self.plan=joblist
    self.executed_plan=[]

  #----------------------------------------------
  def addjob(self,job):
    self.plan.append(job)
    
  #----------------------------------------------
  def nextstep(self):
    self.executed_plan=[]
    for job in self.plan:
      print(" ### jobid %s"%job.jobid)

      if not os.path.exists(job.jobid):
        os.mkdir(job.jobid)
      os.chdir(job.jobid)

      if os.path.exists(job.picklefn):
        with open(job.picklefn,'rb') as inpf:
          rec=pkl.load(inpf)
        if not job.is_consistent(rec):
          raise NotImplementedError("Job not consistent.")
      else: # Nothing done yet. Record starts as plan.
        rec=deepcopy(job)

      rec.nextstep()
      self.executed_plan.append(rec)
      with open(rec.picklefn,'wb') as outf:
        pkl.dump(rec,outf)
     
      os.chdir('..')

  #-----------------------------------------------
  def write_summary(self):
    
    for job in self.executed_plan:
      job.write_summary()
