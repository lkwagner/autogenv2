from __future__ import print_function
import os
import sys
import numpy as np
import subprocess as sub
import shutil
import submitter

####################################################
class RunnerPBS:
  ''' Object that can accumulate jobs to run and run them together in one submission. '''
  def __init__(self,queue='batch',
                    walltime='48:00:00',
                    jobname='AGRunner',
                    np='allprocs',nn=1,
                    prefix=None,
                    postfix=None
                    ):
    ''' Note: exelines are prefixed by appropriate mpirun commands.'''

    # Good prefix choices (Blue Waters).
    # These are needed for Crystal runs.
        # module swap PrgEnv-cray PrgEnv-intel
    # These are needed for QWalk runs.
        #"module swap PrgEnv-cray PrgEnv-gnu",
        #"module load acml",
        #"module load cblas",
    self.exelines=[]
    self.np=np
    self.nn=nn
    self.jobname=jobname
    self.queue=queue
    self.walltime=walltime
    if prefix is None: self.prefix=[]
    else:              self.prefix=prefix
    if postfix is None: self.postfix=[]
    else:               self.postfix=postfix
    self.queueid=[]

  #-------------------------------------
  def check_status(self):
    return submitter.check_PBS_stati(self.queueid)

  #-------------------------------------
  def add_task(self,exestr,loc):
    ''' Accumulate executable commands.
    Args: 
      exestr (str): executible statement. Will be prepended with appropriate mpirun. 
      loc (str): absolute path to location to execute.
    '''

    for line in self.exelines:
      if self.np=='allprocs':
        self.exelines.append("cd {loc}; mpirun {exe}; cd $cwd".format(loc=loc,exe=exestr))
      else:
        self.exelines.append("cd {loc}; mpirun -n {tnp} {exe}; cd $cwd".format(loc=loc,tnp=self.nn*self.np,exe=exestr))

  #-------------------------------------
  def script(self,scriptfile):
    ''' Dump accumulated commands into a script for another job to run.
    Returns true if the runner had lines to actually execute.'''

    if len(self.exelines)==0:
      return False

    with open(scriptfile,'w') as outf:
      outf.write('\n'.join(self.prefix + exelines + self.postfix))

    # Remove exelines so the runner is ready for the next go.
    self.exelines=[]
    return True

  #-------------------------------------
  def submit(self,jobname=None):
    ''' Submit series of commands.'''
    if jobname is None:
      jobname=self.jobname

    if len(self.exelines)==0:
      print("RunnerPBS: no jobs to run.")
      return

    jobout=jobname+'.qsub.out'
    # Submit all jobs.
    qsub=[
        "#PBS -q %s"%self.queue,
        "#PBS -l nodes=%i%s"%(self.nn,ppnstr),
        "#PBS -l walltime=%s"%self.walltime,
        "#PBS -j oe ",
        "#PBS -N %s "%jobname,
        "#PBS -o %s "%jobout,
        "cd %s"%os.getcwd(),
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

    # Remove exelines so the runner is ready for the next go.
    self.exelines=[]
    return qsubfile

####################################################
class PySCFRunnerPBS(RunnerPBS):
  ''' Specialized Runner for dealing with OMP python commands.'''
  def __init__(self,queue='batch',
                    walltime='12:00:00',
                    np='allprocs',
                    ppath=None,
                    nn=1,
                    jobname=os.getcwd().split('/')[-1]+'_pyscf',
                    prefix=None,
                    postfix=None
                    ):
    self.np=np
    if nn!=1: raise NotImplementedError
    self.nn=nn
    self.jobname=jobname
    self.queue=queue
    self.walltime=walltime
    self.prefix=prefix
    self.postfix=postfix
    self.exelines=[]
    if ppath is None: self.ppath=[]
    else:             self.ppath=ppath
    if prefix is None: self.prefix=[]
    else:              self.prefix=prefix
    if postfix is None: self.postfix=[]
    else:               self.postfix=postfix
    self.queueid=None

  #-------------------------------------
  def script(self,scriptfile):
    ''' Dump accumulated commands into a script for another job to run.
    Returns true if the runner had lines to actually execute.'''

    if len(self.exelines)==0:
      return False

    # Shuffle executions with directory changes.
    actions=[]
    for loc,exe in zip(self.loclines,self.exelines):
      actions.append(loc)
      actions.append(exe)
      actions.append('cd $cwd')

    # Prepend mp specs.
    actions=["export OMP_NUM_THREADS=%d"%(self.nn*self.np)]+actions

    # Dump script.
    with open(scriptfile,'w') as outf:
      outf.write('\n'.join(self.prefix + actions + self.postfix))

    # Remove exelines so the runner is ready for the next go.
    self.loclines=[]
    self.exelines=[]

    return True

  #-------------------------------------
  def submit(self,jobname=None):
    if len(self.exelines)==0: 
      print("Nothing to run.")
      return
    if jobname is None: jobname=self.jobname

    # Shuffle executions with directory changes.
    actions=[]
    for loc,exe in zip(self.loclines,self.exelines):
      actions.append(loc)
      actions.append(exe)
      actions.append('cd $cwd')
    
    jobout=jobname+".jobout"
    qsublines=[
         "#PBS -q %s"%self.queue,
       ]
    if self.np == 'allprocs' :
      qsublines+=[
           "#PBS -l nodes=%i,flags=allprocs"%self.nn,
         ]
    else:
      qsublines+=[
           "#PBS -l nodes=%i:ppn=%i"%(self.nn,self.np),
         ]
    qsublines+=[
         "#PBS -l walltime=%s"%self.walltime,
         "#PBS -j oe",
         "#PBS -N %s"%self.jobname,
         "#PBS -o %s"%jobout,
         "cd ${PBS_O_WORKDIR}",
         "export OMP_NUM_THREADS=%d"%(self.nn*self.np),
         "export PYTHONPATH=%s"%(':'.join(self.ppath)),
         "cwd=`pwd`"
       ] + self.prefix + actoins + self.postfix
    qsubfile="qsub.in"
    with open(qsubfile,'w') as f:
      f.write('\n'.join(qsublines))
    result = sub.check_output("qsub %s"%(qsubfile),shell=True)
    self.queueid = result.decode().split()[0]
    print("Submitted as %s"%self.queueid)

    # Clear out the lines to set up for the next job.
    self.loclines=[]
    self.exelines=[]

# TODO Specialize a runner for running QWalk jobs in the same directory together. 
# Should just have to specialize the run command.
