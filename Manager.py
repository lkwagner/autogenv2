import os
import shutil as sh
import numpy as np
from copy import deepcopy

######################################################################
def resolve_status(runner,reader,outfiles, method='not_defined'):
  #Check if the reader is done
  if reader.completed:
    return 'done'

  #Check if the job is in the queue
  # or running. If so, we just return that.
  currstat=runner.check_status()
  print("Current %s status:"%runner.__class__.__name__,currstat)
  if currstat=='running':
    return currstat
  
  #Now we are in a state where either there was an error,
  #the job hasn't been run, or we haven't collected the results
  for outfile in outfiles:
    if not os.path.exists(outfile):
      return 'not_started'

  #We are in an error state or we haven't collected 
  #the results. 
  return "ready_for_analysis"

######################################################################
def diff_keys(old,new,skip_keys=[]):
  ''' Check if two objects have different keys, and what those keys are. '''
  issame=True
  diff={'old':[],'new':[]}

  for newkey in new.__dict__.keys():
    if newkey not in old.__dict__.keys():
      issame=False
      diff['new'].append(newkey)
  for oldkey in old.__dict__.keys():
    if oldkey not in new.__dict__.keys():
      issame=False
      diff['old'].append(oldkey)
  for key in old.__dict__.keys():
    if (key not in diff['new']) and (key not in diff['old']) and \
        (old.__dict__[key]!=new.__dict__[key]) and (key not in skip_keys):
      issame=False
      diff['old'].append(key)
      diff['new'].append(key)
  return issame,diff

######################################################################
def update_attributes(old,new,skip_keys=[],safe_keys=[]):
  ''' Replace attributes that do not affect accuracy. 
  Raise an AssertionError if there's a problematic discrepancy. 
  By default, all keys are not safe, so this mainly checks consistency.
  skip_keys are not checked or replaced.'''

  issame,diff=diff_keys(old,new,skip_keys)
  if not issame:
    print("Key update: {} from one doesn't match {} from new."\
        .format(diff['old'],diff['new']))
    for key in diff['new']:
      if key in safe_keys:
        print("Keeping {} from the latter.".format(diff['new']))
        old.__dict__[key]=new.__dict__[key]
      else:
        raise AssertionError("Unsafe update; new setting affects accuracy.")
  return not issame

#######################################################################
class CrystalManager:
  """ Internal class managing process of running a DFT job though Crystal.
  Has authority over file names associated with this task."""
  def __init__(self,writer,runner,crys_reader,prop_reader):
    ''' convert controls if QWalk input files are produced. '''
    self.writer=writer
    self.creader=crys_reader
    self.preader=prop_reader
    self.runner=runner

    # Internal.
    self.scriptfile=None
    self.crysinpfn='crys.in'
    self.crysoutfn='crys.in.o'
    self.cresinpfn='crys_restart.in'
    self.cresoutfn='crys_restart.in.o'
    self.propinpfn='prop.in'
    self.propoutfn='prop.in.o'
    self.location=os.getcwd()
    self._runready=False
    self.completed=False

  #----------------------------------------
  def nextstep(self):
    """ Check write status, then if it's running. If not running check if
    finished. If not finished, attempt to run. 
    
    Note: this will not submit any jobs to the queue, but updates the runner 
    object with the ability to start the run. Call either submit() or use a bundler to 
    actually start the run.""" 

    # Generate input files.
    if not self.writer.completed:
      if self.writer.guess_fort is not None:
        sh.copy(self.writer.guess_fort,'fort.20')
      with open(self.crysinpfn,'w') as f:
        self.writer.write_crys_input(self.crysinpfn)
      with open(self.propinpfn,'w') as f:
        self.writer.write_prop_input(self.propinpfn)

    #Check on the CRYSTAL run
    # TODO while True is not doing anything anymore.
    while True:
      status=resolve_status(self.runner,self.creader,[self.crysoutfn,self.propoutfn],method='crystal')
      print(status)
    
      print("Crystal status",status)
      if status=="running":
        return
      elif status=="not_started":
        self.runner.prefix.append("cp %s INPUT"%self.crysinpfn)
        self.runner.run("Pcrystal &> %s"%self.crysoutfn)
        self.runner.postfix.append("properties < %s &> %s"%(self.propinpfn,self.propoutfn))
        return
      elif status=="ready_for_analysis":
        #This is where we (eventually) do error correction and resubmits
        status=self.creader.collect(self.crysoutfn)
        self.preader.collect(self.propoutfn)
        if status=='killed':
          self.writer.restart=True
          self.writer.guess_fort='./fort.79'
          self.writer.write_crys_input(self.cresinpfn)
          self.runner.run("Pcrystal &> %s"%self.cresoutfn)
        break
      elif status=='done':
        break
      else:
        return

    self.completed=(self.creader.completed and self.preader.completed)

  #------------------------------------------------
  def script(self,jobname=None):
    ''' Script execution lines for a bundler to pick up and run.'''
    if jobname is None: jobname=self.runner.jobname
    self.scriptfile="%s.run"%jobname
    self._runready=self.runner.script(self.scriptfile)

  #------------------------------------------------
  def submit(self,jobname=None):
    ''' Submit the runner's job to the queue. '''
    qsubfile=self.runner.submit(jobname)
    return qsubfile

  #------------------------------------------------
  def update_queueid(self,qid):
    ''' If a bundler handles the submission, it can update the queue info with this.'''
    self.runner.queueid.append(qid)
    self._runready=False # After running, we won't run again without more analysis.

  #------------------------------------------------
  def update_options(self,other):
    ''' Safe copy options from other to self. '''

    updated=update_attributes(old=self.runner,new=other.runner,
        safe_keys=['queue','walltime','np','nn','jobname','prefix','postfix'],
        skip_keys=['queueid'])

    updated=update_attributes(old=self.writer,new=other.writer,
        safe_keys=['maxcycle','edifftol'],
        skip_keys=['completed'])
    if updated:
      self.writer.completed=False

  #----------------------------------------
  def to_json(self):
    raise NotImplementedError

  #----------------------------------------
  def write_summary(self):
    self.creader.write_summary()
    
  #----------------------------------------
  def status(self):
    if self.completed:
      return 'ok'
    else:
      return 'not_finished'

#######################################################################

class QWalkfromCrystalManager:
  """Set up a QWalk job from Crystal. 
  In this we will Convert from a CRYAPI_OUT-ed properties run. 
  """
  #------------------------------------------------
  def __init__(self,convert_runner,convert_checker):
    self.runner=convert_runner
    self.reader=convert_checker
    
  #------------------------------------------------
  def update_options(self,other):
    ''' Safe copy options from other to self. '''

    updated=update_attributes(old=self.runner,new=other.runner,
        safe_keys=['queue','walltime','np','nn','jobname','prefix','postfix'],
        skip_keys=['queueid'])

    if updated:
      self.writer.completed=False

  #------------------------------------------------
  def is_consistent(self,other):
    return self.runner.is_consistent(other.runner)
  
  #------------------------------------------------
  def nextstep(self):
    if not self.reader.completed:
      self.runner.run()
      self.reader.collect()
  #------------------------------------------------

  def write_summary(self):
    print("K-points",len(self.reader.out['basenames']))

  #----------------------------------------
  def status(self):
    if self.reader.completed:
      return 'ok'
    else:
      return 'not_finished'
    
      
#######################################################################

class QWalkRunManager:
  def __init__(self,writer,runner,reader):
    self.writer=writer
    self.runner=runner
    self.reader=reader
    self.infiles=[]
    self.outfiles=[]

  #------------------------------------------------
  def is_consistent(self,other):
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)

  #------------------------------------------------
  def update_options(self,other):
    ''' Safe copy options from other to self. '''

    update_attributes(old=self.runner,new=other.runner,
        safe_keys=['queue','walltime','np','nn','jobname','prefix','postfix'],
        skip_keys=['queueid'])

    update_attributes(old=self.writer,new=other.writer,
        skip_keys=[ 'completed','sysfiles','slaterfiles','jastfiles',
          'basenames','wffiles','tracefiles'])
    
  #------------------------------------------------
  def nextstep(self):
    if not self.writer.completed:
      self.infiles,self.outfiles=self.writer.qwalk_input()
    
    while True:
      status=resolve_status(self.runner,self.reader,self.outfiles)
      print("%s status: %s"%(self.writer.qmc_type,status))
      if status=="running":
        return
      elif status=="not_started":
        stdoutfiles=[x+".stdout" for x in self.infiles]
        exestr="~/bin/qwalk {}".format(*self.infiles)
        # TODO this outfile handling is shitty.
        exestr+=" &> {}".format(self.infiles[0]+'.out')
        self.runner.run(exestr)
        print("%s status: submitted"%(self.writer.qmc_type))
        return
      elif status=="ready_for_analysis":
        #This is where we (eventually) do error correction and resubmits
        self.reader.collect(self.outfiles)
        break
      elif status=='done':
        break
      else:
        return
      
  #------------------------------------------------
  def write_summary(self):
    self.reader.write_summary()

  #------------------------------------------------
  def generate_report(self):
    return {}

  #----------------------------------------
  def status(self):
    if self.reader.completed:
      return 'ok'
    else:
      return 'not_finished'
    

#######################################################################


class PySCFManager:
  def __init__(self,writer,runner,reader):
    self.writer=writer
    self.runner=runner
    self.reader=reader
    self.driverfn='pyscf_driver.py'
    self.completed=False
    self.infile=''
    self.restart_infile=''
    self.outfile=''
  #------------------------------------------------
  # Obsolete with update_options?
  def is_consistent(self,other):
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)

  #------------------------------------------------
  def update_options(self,other):
    ''' Safe copy options from other to self. '''

    updated=update_attributes(old=self.runner,new=other.runner,
        safe_keys=['queue','walltime','np','nn','jobname','prefix','postfix'],
        skip_keys=['queueid'])

    updated=update_attributes(old=self.writer,new=other.writer,
        safe_keys=['max_cycle'],
        skip_keys=['completed','chkfile','dm_generator'])
    if updated:
      self.writer.completed=False
    
  #------------------------------------------------
  def nextstep(self):
    if not self.writer.completed:
      self.infile,self.restart_infile,self.outfile,self.chkfile=self.writer.pyscf_input(self.driverfn)
    
    status=resolve_status(self.runner,self.reader,self.outfile, 'pyscf')
    print("PySCF status",status)
    if status=="running":
      pass
    elif status=="not_started":
      self.runner.run("python %s &> %s"%(self.infile,self.outfile))
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      self.reader.collect(self.outfile,self.chkfile)
    elif status=='done':
      pass
    #If we need to restart the run
    elif status=='retry':
      self.runner.run(self.restart_infile, self.outfile)

    self.completed=self.reader.completed
      
  #------------------------------------------------

  def write_summary(self):
    self.reader.write_summary()
    

  #----------------------------------------
  def status(self):
    current_status = resolve_status(self.runner,self.reader,self.outfiles, 'pyscf')
    if current_status == 'done':
      return 'ok'
    elif current_status == 'retry':
      return 'retry'
    else:
      return 'not_finished'
