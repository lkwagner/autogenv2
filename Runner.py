from __future__ import print_function
import os
import sys
import numpy as np
import subprocess as sub
import shutil
import submitter

# TODO organize with inheritance.

####################################################
class RunnerLocal:
  ''' Object that can accumulate jobs to run and run them together locally.'''
  def __init__(self,np='allprocs',nn=1):
    ''' Note: exelines are prefixed by appropriate mpirun commands.'''

    self.exelines=[]
    self.np=np
    self.nn=nn

  #-------------------------------------
  def check_status(self):
    return 'done'

  #-------------------------------------
  def add_task(self,exestr):
    ''' Accumulate executable commands.
    Args: 
      exestr (str): executible statement. Will be prepended with appropriate mpirun. 
    '''

    if self.np=='allprocs':
      self.exelines.append("mpirun {exe}".format(exe=exestr))
    else:
      self.exelines.append("mpirun -n {tnp} {exe}".format(tnp=self.nn*self.np,exe=exestr))

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
      return ''
    
    try:
      for line in self.exelines:
        result = sub.check_output(line,shell=True)
        print(self.__class__.__name__,": executed %s"%line)
    except sub.CalledProcessError:
      print(self.__class__.__name__,": Error executing line.")

    # Remove exelines so the runner is ready for the next go.
    self.exelines=[]
    return ''

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

  def add_command(self,cmdstr):
    ''' Accumulate commands that don't get an MPI command.
    Args: 
      cmdstr (str): executible statement. Will be prepended with appropriate mpirun. 
    '''
    self.exelines.append(cmdstr)

  #-------------------------------------
  def add_task(self,exestr):
    ''' Accumulate executable commands.
    Args: 
      exestr (str): executible statement. Will be prepended with appropriate mpirun. 
    '''

    if self.np=='allprocs':
      self.exelines.append("mpirun {exe}".format(exe=exestr))
    else:
      self.exelines.append("mpirun -n {tnp} {exe}".format(tnp=self.nn*self.np,exe=exestr))

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
      #print(self.__class__.__name__,": All tasks completed or queued.")
      return
    
    if self.np=='allprocs':
      ppnstr=',flags=allprocs'
    else:
      ppnstr=':ppn=%d'%self.np

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
      ] + self.prefix + self.exelines + self.postfix
    qsubfile=jobname+".qsub"
    with open(qsubfile,'w') as f:
      f.write('\n'.join(qsub))
    try:
      result = sub.check_output("qsub %s"%(qsubfile),shell=True)
      self.queueid.append(result.decode().split()[0].split('.')[0])
      print(self.__class__.__name__,": Submitted as %s"%self.queueid)
    except sub.CalledProcessError:
      print(self.__class__.__name__,": Error submitting job. Check queue settings.")

    # Remove exelines so the runner is ready for the next go.
    self.exelines=[]
    return qsubfile

####################################################
class PySCFRunnerLocal:
  ''' Object that can accumulate jobs to run and run them together locally.'''
  def __init__(self,np='allprocs'):
    ''' Note: exelines are prefixed by appropriate mpirun commands.'''
    self.exelines=[]
    self.queueid=[]

  #-------------------------------------
  def check_status(self):
    return 'done'

  #-------------------------------------
  def add_task(self,exestr):
    ''' Accumulate executable commands.
    Args: 
      exestr (str): executible statement. Will be prepended with appropriate mpirun. 
    '''
    self.exelines.append(exestr)

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
      #print(self.__class__.__name__,": All tasks completed or queued.")
      return ''    

    try:
      for line in self.exelines:
        result = sub.check_output(line,shell=True)
        print(self.__class__.__name__,": executed %s"%line)
    except sub.CalledProcessError:
      print(self.__class__.__name__,": Error executing line.")

    # Remove exelines so the runner is ready for the next go.
    self.exelines=[]
    return ''

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
    self.queueid=[]

  #-------------------------------------
  def add_task(self,exestr):
    ''' Accumulate executable commands.
    Args: 
      exestr (str): executible statement. Will be prepended with appropriate mpirun. 
    '''
    self.exelines.append(exestr)

  #-------------------------------------
  def script(self,scriptfile):
    ''' Dump accumulated commands into a script for another job to run.
    Returns true if the runner had lines to actually execute.'''

    if len(self.exelines)==0:
      return False

    # Prepend mp specs.
    actions=["export OMP_NUM_THREADS=%d"%(self.nn*self.np)]+actions

    # Dump script.
    with open(scriptfile,'w') as outf:
      outf.write('\n'.join(self.prefix + actions + self.postfix))

    # Remove exelines so the runner is ready for the next go.
    self.exelines=[]

    return True

  #-------------------------------------
  def submit(self,jobname=None):
    if len(self.exelines)==0: 
      #print(self.__class__.__name__,": All tasks completed or queued.")
      return
    if jobname is None: jobname=self.jobname

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
       ] + self.prefix + self.exelines + self.postfix
    qsubfile=jobname+".qsub"
    with open(qsubfile,'w') as f:
      f.write('\n'.join(qsublines))
    try: 
      result = sub.check_output("qsub %s"%(qsubfile),shell=True)
      self.queueid.append(result.decode().split()[0])
      print(self.__class__.__name__,": Submitted as %s"%self.queueid)
    except sub.CalledProcessError:
      print(self.__class__.__name__,": Error submitting job. Check queue settings.")

    # Clear out the lines to set up for the next job.
    self.exelines=[]

# TODO Specialize a runner for running QWalk jobs in the same directory together. 
# Should just have to specialize the run command.
