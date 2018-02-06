import Crystal
import PropertiesReader
import os
#import PySCF # TODO uncomment, due to library issue 
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

  #We are in an error state or we haven't collected the results. 
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
        print("Keeping {} from the latter.".format(key))
        old.__dict__[key]=new.__dict__[key]
      else:
        print("Problem with update of {}".format(key))
        raise AssertionError("Unsafe update; new setting affects accuracy.")
  return not issame

#######################################################################
class CrystalManager:
  """ Internal class managing process of running a DFT job though Crystal.
  Has authority over file names associated with this task."""
  def __init__(self,writer,runner,crys_reader=None,prop_reader=None,trylev=False):
    ''' convert controls if QWalk input files are produced. None makes a default instance.'''
    self.writer=writer

    if crys_reader is None:
      self.creader=Crystal.CrystalReader()
    else:
      self.creader=crys_reader
    if prop_reader is None:
      self.preader=PropertiesReader.PropertiesReader()
    else:
      self.preader=prop_reader
    self.runner=runner

    # Internal.
    self.scriptfile=None
    self.crysinpfn='crys.in'
    self.crysoutfn='crys.in.o'
    #self.cresinpfn='crys_restart.in'
    #self.cresoutfn='crys_restart.in.o'
    self.propinpfn='prop.in'
    self.propoutfn='prop.in.o'
    self.location='unset'
    self.restarts=0
    self._runready=False
    self.completed=False

    # Smart error detection.
    self.trylev=trylev
    self.savebroy=[]
    self.lev=False


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
      self.location=os.getcwd()

    #Check on the CRYSTAL run
    # TODO while True is not doing anything anymore.
    status=resolve_status(self.runner,self.creader,[self.crysoutfn],method='crystal')
    print(status)
  
    print("Crystal status",status)
    if status=="not_started":
      sh.copy(self.crysinpfn,'INPUT')
      self.runner.run("Pcrystal &> %s"%self.crysoutfn)
      self.runner.postfix += ["cp %s INPUT"%self.propinpfn,"mpirun Pproperties &> %s.o"%self.propinpfn]
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.creader.collect(self.crysoutfn)
      self.preader.collect(self.propoutfn)
      if status=='killed':
        self.writer.restart=True
        if self.trylev:
          print("Trying LEVSHIFT.")
          self.writer.levshift=[10,1] # No mercy.
          self.savebroy=deepcopy(self.writer.broyden)
          self.writer.broyden=[]
          self.lev=True
        sh.copy(self.crysinpfn,"%d.%s"%(self.restarts,self.crysinpfn))
        sh.copy(self.crysoutfn,"%d.%s"%(self.restarts,self.crysoutfn))
        sh.copy('fort.79',"%d.fort.79"%(self.restarts))
        self.writer.guess_fort='./fort.79'
        sh.copy(self.writer.guess_fort,'fort.20')
        self.writer.write_crys_input(self.crysinpfn)
        sh.copy(self.crysinpfn,'INPUT')
        self.runner.run("Pcrystal &> %s"%self.crysoutfn)
        self.restarts+=1
    elif status=='done' and self.lev:
      # We used levshift to converge. Now let's restart to be sure.
      print("Recovering from LEVSHIFTer.")
      self.writer.restart=True
      self.writer.levshift=[]
      self.creader.completed=False
      self.lev=False
      sh.copy(self.crysinpfn,"%d.%s"%(self.restarts,self.crysinpfn))
      sh.copy(self.crysoutfn,"%d.%s"%(self.restarts,self.crysoutfn))
      sh.copy('fort.79',"%d.fort.79"%(self.restarts))
      self.writer.guess_fort='./fort.79'
      sh.copy(self.writer.guess_fort,'fort.20')
      self.writer.write_crys_input(self.crysinpfn)
      sh.copy(self.crysinpfn,'INPUT')
      self.runner.run("Pcrystal &> %s"%self.crysoutfn)
      self.restarts+=1

    # Sometimes an error puts it in a state where the crystal is done but not properties.
    elif status=='done' and not self.preader.completed:
      print("Crystal is done, but properties not...trying to read properties.")
      self.preader.collect(self.propoutfn)

    self.completed=(self.creader.completed and self.preader.completed)

  #------------------------------------------------
  def script(self,jobname=None):
    ''' Script execution lines for a bundler to pick up and run.'''
    if jobname is None: jobname=self.runner.jobname
    self.scriptfile="%s.run"%jobname
    self._runready=self.runner.script(self.scriptfile)
    self.location=os.getcwd()

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
        safe_keys=['queue','walltime','np','nn','jobname'],
        skip_keys=['queueid','prefix','postfix'])

    updated=update_attributes(old=self.writer,new=other.writer,
        safe_keys=['maxcycle','edifftol'],
        skip_keys=['completed','modisymm','restart','guess_fort'])
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
  ''' Manage files involved in a QWalk run. 

  This object should work with all types of qwalk runs, only the writer and
  runner changes. 

  Attributes:
    writer (Writer object): Object to write input files (see DMCWriter).
    runner (Runner object): Object to handle job submission (see RunnerPBC).
    reader (Reader object): Object to handle reading of QMC output (see DMCReader)
    completed (bool): Whether the Manager has finished all tasks.

  Args: 
    writer (Writer object): Set as writer attribute.
    runner (Runner object): Set as runner attribute. 
    reader (Reader object): Set as reader attribute.
    infiles (list of str): List of input files as defined by conversion.
      (TODO fix dependency)
      for each file.
  '''

  def __init__(self,writer,runner,reader,infiles):
    self.writer=writer
    self.runner=runner
    self.reader=reader
    self.completed=False

    self.scriptfile=None
    self._runready=False
    self.infiles=infiles
    self.outfiles=["%s.o"%f for f in infiles]

  #------------------------------------------------
  def is_consistent(self,other):
    ''' Checks that the input settings that determine accuracy are the same. 
    
    Args: 
      other (QWalkRunManager): Check this instance is the same as other.
    Returns:
      bool: Whether the accuracy of the two should be the same.
    '''
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)

  #------------------------------------------------
  def update_options(self,other):
    ''' Safe copy options from other to self. 
    
    Args: 
      other (QWRunManager): New object to copy attributes from.
    '''

    update_attributes(old=self.runner,new=other.runner,
        safe_keys=['queue','walltime','np','nn','jobname','prefix','postfix'],
        skip_keys=['queueid'])

    update_attributes(old=self.writer,new=other.writer,
        skip_keys=[ 'completed','sysfiles','slaterfiles',
          'jastfiles','wffiles','tracefiles'])
    
  #------------------------------------------------
  def nextstep(self):
    ''' Perform next step in calculation. '''
    if not self.writer.completed:
      self.writer.qwalk_input(self.infiles)
    
    status=resolve_status(self.runner,self.reader,self.outfiles)
    print("%s status: %s"%(self.writer.qmc_type,status))
    if status=="not_started":
      exestr="/home/busemey2/bin/qwalk %s"%' '.join(self.infiles)
      exestr+=" &> %s.out"%self.infiles[-1] # TODO better outfile name.
      self.runner.run(exestr)
      print("%s status: submitted"%(self.writer.qmc_type))
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.reader.collect(self.outfiles)
      if status=='ok':
        self.completed=True
      else:
        print("Status: %s, attempting rerun."%status)
        exestr="/home/busemey2/bin/qwalk %s"%' '.join(self.infiles)
        exestr+=" &> %s.out"%self.infiles[-1]
        self.runner.run(exestr)
        print(self.runner.exelines)
    elif status=='done':
      self.completed=True
    else:
      print("Status: %s"%status)

  #------------------------------------------------
  def script(self,jobname=None):
    ''' Script execution lines for a bundler to pick up and run.
    
    Args:
      jobname (str): name for job in queue.
    '''
    if jobname is None: jobname=self.runner.jobname
    #if jobname is None: jobname="QWManager"
    self.scriptfile="%s.run"%jobname
    self._runready=self.runner.script(self.scriptfile)
    self.location=os.getcwd()
      
  #------------------------------------------------
  def update_queueid(self,qid):
    ''' If a bundler handles the submission, it can update the queue info with this.
    Args:
      qid (str): new queue id from submitting a job. The Manager will check if this is running.
    '''
    self.runner.queueid.append(qid)
    self._runready=False # After running, we won't run again without more analysis.

  #------------------------------------------------
  def write_summary(self):
    ''' Write a summary of the job results.'''
    self.reader.write_summary()

  #------------------------------------------------
  def generate_report(self):
    ''' Generate a useful report of the run inputs and results.'''
    # TODO Not implemented.
    return {}

  #----------------------------------------
  def status(self):
    ''' Check if this Manager has completed all it's tasks.
    Returns:
      str: 'ok' or 'not_finished'.
    '''
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
    self.outfile=self.driverfn+'.o'
    self.chkfile=self.driverfn+'.chkfile'
    self.completed=False
    self.restarts=0
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
      self.writer.pyscf_input(self.driverfn,self.chkfile)
    
    status=resolve_status(self.runner,self.reader,[self.outfile], 'pyscf')
    print("PySCF status",status)
    if status=="running":
      pass
    elif status=="not_started":
      self.runner.run("/usr/bin/python %s &> %s"%(self.driverfn,self.outfile))
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.reader.collect(self.outfile,self.chkfile)
      if status=='killed':
        sh.copy(self.driverfn,"%d.%s"%(self.restarts,self.driverfn))
        sh.copy(self.outfile,"%d.%s"%(self.restarts,self.outfile))
        sh.copy(self.chkfile,"%d.%s"%(self.restarts,self.chkfile))
        self.writer.dm_generator=PySCF.dm_from_chkfile("%d.%s"%(self.restarts,self.chkfile))
        self.writer.pyscf_input(self.driverfn,self.chkfile)
        self.runner.run("/usr/bin/python %s &> %s"%(self.driverfn,self.outfile))
        self.restarts+=1

    self.completed=self.reader.completed

  #------------------------------------------------
  def script(self,jobname=None):
    ''' Script execution lines for a bundler to pick up and run.'''
    if jobname is None: jobname=self.runner.jobname
    self.scriptfile="%s.run"%jobname
    self._runready=self.runner.script(self.scriptfile,self.driverfn)

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

  def write_summary(self):
    self.reader.write_summary()

  #----------------------------------------
  def status(self):
    current_status = resolve_status(self.runner,self.reader,[self.outfile], 'pyscf')
    if current_status == 'done':
      return 'ok'
    elif current_status == 'retry':
      return 'retry'
    else:
      return 'not_finished'
