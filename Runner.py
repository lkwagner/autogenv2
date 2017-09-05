from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import submitter

####################################################

class RunnerPBS:
  _name_='RunnerPBS'
  def __init__(self,queue='batch',
                    walltime='12:00:00',
                    jobname='AGRunner',
                    prefix=[],
                    postfix=[],
                    np=1,nn=1
                    ):
    ''' Note: exelines are prefixed by appropriate mpirun commands.'''
    self.exelines=[]
    self.np=np
    self.nn=nn
    self.queue=queue
    self.walltime=walltime
    self.prefix=prefix
    self.postfix=postfix
    self.queueid=[]

  #-------------------------------------
  def check_status(self):
    return submitter.check_PBS_stati(self.queueid)

  #-------------------------------------
  def run(self,exestr):
    self.exelines.append(exestr)

  #-------------------------------------
  def submit(self,jobname=None):
    if jobname is None:
      jobname=self.jobname

    if len(self.exelines)==0:
      print("RunnerPBS: warning, no jobs to run.")
      return

    # Prepend mpi specs.
    exelines=[]
    for line in self.exelines:
      exelines.append('mpirun -n {} {}'.format(
        self.nn*self.np,line))

    jobout=jobname+'.qsub.out'
    # Submet all jobs.
    qsub=[
        "#PBS -q %s"%self.queue,
        "#PBS -l nodes=%i:ppn=%i"%(self.nn,self.np),
        "#PBS -l walltime=%s"%self.walltime,
        "#PBS -j oe ",
        "#PBS -N %s "%jobname,
        "#PBS -o %s "%jobout,
        "cd ${PBS_O_WORKDIR}",
      ] + self.prefix + exelines + self.postfix
    qsubfile=jobname+".qsub"
    with open(qsubfile,'w') as f:
      f.write('\n'.join(qsub))
    try:
      result = sub.check_output("qsub %s"%(qsubfile),shell=True)
      self.queueid.append(result.decode().split()[0].split('.')[0])
      print("Submitted as %s"%self.queueid)
    except sub.CalledProcessError:
      print("Error submitting job. Check queue settings.")

####################################################

class RunnerBW:
  _name_='RunnerBW'
  def __init__(self,queue='normal',
                    walltime='48:00:00',
                    jobname='AGRunner',
                    prefix=[],
                    postfix=[],
                    np=32,nn=1
                    ):
    ''' Note: exelines are prefixed by appropriate mpirun commands.'''
    self.exelines=[]
    self.np=np
    self.nn=nn
    self.queue=queue
    self.walltime=walltime
    self.queueid=[]
    self.prefix=prefix
    self.postfix=postfix

  #-------------------------------------
  def check_status(self):
    return submitter.check_PBS_stati(self.queueid)

  #-------------------------------------
  def run(self,exestr):
    ''' Accumulate executable commands.'''
    self.exelines.append(exestr)

  #-------------------------------------
  def submit(self,jobname=None):
    ''' Submit series of commands.'''
    if jobname is None:
      jobname=self.jobname

    if len(self.exelines)==0:
      print("RunnerPBS: warning, no jobs to run.")
      return

    # Prepend mpi specs.
    exelines=[]
    for line in self.exelines:
      exelines.append('aprun -n {} {}'.format(
        self.nn*self.np,line))

    jobout=jobname+'.qsub.out'
    # Submit all jobs.
    qsub=[
        "#PBS -q %s"%self.queue,
        "#PBS -l nodes=%i:ppn=%i:xe"%(self.nn,self.np),
        "#PBS -l walltime=%s"%self.walltime,
        "#PBS -j oe ",
        "#PBS -A bahu",
        "#PBS -N %s "%jobname,
        "#PBS -o %s "%jobout,
        "cd %s"%os.getcwd(),
        "module swap PrgEnv-cray PrgEnv-gnu",
        "module load acml",
        "module load cblas",
      ] + self.prefix + exelines + self.postfix
    qsubfile=jobname+".qsub"
    with open(qsubfile,'w') as f:
      f.write('\n'.join(qsub))
    try:
      result = sub.check_output("qsub %s"%(qsubfile),shell=True)
      self.queueid.append(result.decode().split()[0].split('.')[0])
      print("Submitted as %s"%self.queueid)
    except sub.CalledProcessError:
      print("Error submitting job. Check queue settings.")
