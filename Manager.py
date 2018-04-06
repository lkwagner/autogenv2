import Crystal
import PropertiesReader
import os
import PySCF 
import shutil as sh
import numpy as np
import pyscf2qwalk
import pickle as pkl
import Runner
from copy import deepcopy

#TODO These should all probably be subclasses, if not just for logical purposes. 

######################################################################
# Misc functions. 
######################################################################
def resolve_status(runner,reader,outfile, method='not_defined'):
  #Check if the reader is done
  if reader.completed:
    return 'done'

  #Check if the job is in the queue or running. If so, we just return that.
  currstat=runner.check_status()
  if currstat=='running':
    return currstat
  
  #Now we are in a state where either there was an error,
  #the job hasn't been run, or we haven't collected the results
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
    print(old.__class__.__name__,": Key update-- {} from one doesn't match {} from new."\
        .format(diff['old'],diff['new']))
    for key in diff['new']:
      if key in safe_keys:
        print("Keeping {} from the latter.".format(key))
        print("old")
        print(old.__dict__[key])
        print("new")
        print(new.__dict__[key])
        old.__dict__[key]=new.__dict__[key]
      else:
        raise AssertionError("Unsafe update; new setting affects accuracy.")
  return not issame

######################################################################
# Manager classes.
#######################################################################
class CrystalManager:
  """ Internal class managing process of running a DFT job though Crystal.
  Has authority over file names associated with this task."""
  def __init__(self,writer,runner,name='crystal_run',path=None,crys_reader=None,prop_reader=None,trylev=False,bundle=False):
    ''' CrystalManager manages the writing of a Crystal input file, it's running, and keeping track of the results.
    Args:
      writer (PySCFWriter): writer for input.
      reader (PySCFReader): to read PySCF output.
      runner (runner object): to run job.
      name (str): identifier for this job. This names the files associated with run.
      path (str): directory where this manager is free to store information.
      bundle (bool): False - submit jobs. True - dump job commands into a script for a bundler to run.
    '''
    # Where to save self.
    self.name=name
    self.pickle="%s.pkl"%self.name

    # Ensure path is set up correctly.
    if path is None:
      path=os.path.getcwd()
    if path[-1]!='/': path+='/'
    self.path=path

    print(self.__class__.__name__,": initializing")
    print(self.__class__.__name__,": name= %s"%self.name)
    print(self.__class__.__name__,": path= %s"%self.path)

    self.writer=writer
    if crys_reader is None: self.creader=Crystal.CrystalReader()
    else: self.creader=crys_reader
    if prop_reader is None: self.preader=PropertiesReader.PropertiesReader()
    else: self.preader=prop_reader
    if runner is None: self.runner=RunnerPBS()
    else: self.runner=runner

    # Internal.
    self.crysinpfn=self.name+'.in'
    self.propinpfn=self.name+'.prop.in'
    self.crysoutfn=self.crysinpfn+'.o'
    self.propoutfn=self.propinpfn+'.o'
    self._runready=False
    self.scriptfile=None
    self.completed=False
    self.bundle=bundle

    # Smart error detection.
    self.trylev=trylev
    self.savebroy=[]
    self.lev=False

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      print(self.__class__.__name__,": rebooting old manager.")
      old=pkl.load(open(self.path+self.pickle,'rb'))
      if False: #not self.is_consistent(old): TODO check consistency.
        raise NotImplementedError("Handling updated input files is not implemented yet.")
      else:
        self.__dict__=old.__dict__

    # Update the file.
    if not os.path.exists(self.path): os.mkdir(self.path)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

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
  def nextstep(self):
    ''' Determine and perform the next step in the calculation.'''

    print(self.__class__.__name__,": next step.")
    cwd=os.getcwd()
    os.chdir(self.path)

    # Generate input files.
    if not self.writer.completed:
      if self.writer.guess_fort is not None:
        sh.copy(self.writer.guess_fort,'fort.20')
      with open(self.crysinpfn,'w') as f:
        self.writer.write_crys_input(self.crysinpfn)
      with open(self.propinpfn,'w') as f:
        self.writer.write_prop_input(self.propinpfn)

    #Check on the CRYSTAL run
    status=resolve_status(self.runner,self.creader,self.crysoutfn,method='crystal')
    print(self.__class__.__name__,": %s status= %s"%(self.name,status))

    if status=="not_started":
      self.runner.add_command("cp %s INPUT"%self.crysinpfn)
      self.runner.add_task("Pcrystal &> %s"%self.crysoutfn)

    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.creader.collect(self.crysoutfn)
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
        self.runner.add_task("Pcrystal &> %s"%self.crysoutfn)
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
      self.runner.add_task("Pcrystal &> %s"%self.crysoutfn)
      self.restarts+=1

    # Ready for bundler or else just submit the jobs as needed.
    if self.bundle:
      self.scriptfile="%s.run"%self.name
      self.bundle_ready=self.runner.script(self.scriptfile,self.driverfn)
    else:
      qsubfile=self.runner.submit(self.name)

    self.completed=self.creader.completed

    # Update the file.
    with open(self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    os.chdir(cwd)

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
class PySCFManager:
  def __init__(self,writer,reader=None,runner=None,name='psycf_run',path=None,bundle=False):
    ''' PySCFManager manages the writing of a PySCF input file, it's running, and keep track of the results.
    Args:
      writer (PySCFWriter): writer for input.
      reader (PySCFReader): to read PySCF output.
      runner (runner object): to run job.
      name (str): identifier for this job. This names the files associated with run.
      path (str): directory where this manager is free to store information.
      bundle (bool): False - submit jobs. True - dump job commands into a script for a bundler to run.
    '''
    # Where to save self.
    self.name=name
    self.pickle="%s.pkl"%self.name

    # Ensure path is set up correctly.
    if path is None:
      path=os.path.getcwd()
    if path[-1]!='/': path+='/'
    self.path=path

    print(self.__class__.__name__,": initializing")
    print(self.__class__.__name__,": name= %s"%self.name)
    print(self.__class__.__name__,": path= %s"%self.path)

    self.writer=writer
    if reader is not None: self.reader=reader
    else: self.reader=PySCF.PySCFReader()
    if runner is not None: self.runner=runner
    else: self.runner=PySCFRunnerPBS()
    self.bundle=bundle

    self.driverfn="%s.py"%name
    self.outfile=self.driverfn+'.o'
    self.chkfile=self.driverfn+'.chkfile'
    self.qwfiles={}
    self.completed=False
    self.bundle_ready=False
    self.restarts=0

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      print(self.__class__.__name__,": rebooting old manager.")
      old=pkl.load(open(self.path+self.pickle,'rb'))
      if False: #not self.is_consistent(old): TODO check consistency.
        raise NotImplementedError("Handling updated input files is not implemented yet.")
      else:
        self.__dict__=old.__dict__

    # Update the file.
    if not os.path.exists(self.path): os.mkdir(self.path)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #------------------------------------------------
  def update_options(self,other):
    ''' Safe copy options from other to self. '''
    updated=update_attributes(old=self.writer,new=other.writer,
        safe_keys=['max_cycle'],
        skip_keys=['completed','chkfile','dm_generator'])
    if updated:
      self.writer.completed=False

    # Update the file.
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    
  #------------------------------------------------
  def nextstep(self):
    ''' Determine and perform the next step in the calculation.'''
    print(self.__class__.__name__,": next step.")
    cwd=os.getcwd()
    os.chdir(self.path)

    if not self.writer.completed:
      self.writer.pyscf_input(self.driverfn,self.chkfile)
    
    status=resolve_status(self.runner,self.reader,self.outfile, 'pyscf')
    print(self.__class__.__name__,": %s status= %s"%(self.name,status))

    if status=="not_started":
      self.runner.add_task("/usr/bin/python3 %s > %s"%(self.driverfn,self.outfile))
    elif status=="ready_for_analysis":
      status=self.reader.collect(self.outfile,self.chkfile)
      if status=='killed':
        sh.copy(self.driverfn,"%d.%s"%(self.restarts,self.driverfn))
        sh.copy(self.outfile,"%d.%s"%(self.restarts,self.outfile))
        sh.copy(self.chkfile,"%d.%s"%(self.restarts,self.chkfile))
        if os.path.exists(self.chkfile):
          self.writer.dm_generator=PySCF.dm_from_chkfile("%d.%s"%(self.restarts,self.chkfile))
        self.writer.pyscf_input(self.driverfn,self.chkfile)
        self.runner.add_task("/usr/bin/python3 %s > %s"%(self.driverfn,self.outfile))
        self.restarts+=1
      elif status=='ok':
        print(self.__class__.__name__,": %s status= %s, task complete."%(self.name,status))

    # Ready for bundler or else just submit the jobs as needed.
    if self.bundle:
      self.scriptfile="%s.run"%self.name
      self.bundle_ready=self.runner.script(self.scriptfile,self.driverfn)
    else:
      qsubfile=self.runner.submit(self.name)

    self.completed=self.reader.completed
    # Update the file.
    with open(self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    os.chdir(cwd)

  #------------------------------------------------
  def update_queueid(self,qid):
    ''' If a bundler handles the submission, it can update the queue info with this.'''
    self.runner.queueid.append(qid)
    # Update the file.
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    self._runready=False # After running, we won't run again without more analysis.
      
  #------------------------------------------------
  def export_qwalk(self):
    ''' Export QWalk input files into current directory.'''
    if len(self.qwfiles)==0:
      self.nextstep()
      if not self.completed:
        return False
      else:
        print(self.__class__.__name__,": %s generating QWalk files."%self.name)
        cwd=os.getcwd()
        os.chdir(self.path)
        self.qwfiles=pyscf2qwalk.print_qwalk_chkfile(self.chkfile)
        os.chdir(cwd)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    return True

  #----------------------------------------
  def status(self):
    ''' Determine the course of action based on info from reader and runner.'''
    current_status = resolve_status(self.runner,self.reader,self.outfile, 'pyscf')
    if current_status == 'done':
      return 'ok'
    elif current_status == 'retry':
      return 'retry'
    else:
      return 'not_finished'

#######################################################################
class QWalkManager:
  def __init__(self,writer,reader,runner=None,trialfunc=None,name='qw_run',path=None,bundle=False,qwalk='~/bin/qwalk'):
    ''' QWalkManager managers the writing of a QWalk input files, it's running, and keeping track of the results.
    Args:
      writer (qwalk writer): writer for input.
      reader (qwalk reader): to read job.
      runner (runner object): to run job.
      trialfunc (TrialFunction): TrialFunction object for generating trail function input. 
        Note: This is only used if write.trailfunc arguement==''. 
      name (str): identifier for this job. This names the files associated with run.
      path (str): directory where this manager is free to store information.
      bundle (bool): False - submit jobs. True - dump job commands into a script for a bundler to run.
      qwalk (str): absolute path to qwalk executible.
    '''
    self.name=name
    self.pickle="%s.pkl"%(self.name)

    # Ensure path is set up correctly.
    if path is None:
      path=os.path.getcwd()
    if path[-1]!='/': path+='/'
    self.path=path

    print(self.__class__.__name__,": initializing")
    print(self.__class__.__name__,": name= %s"%self.name)
    print(self.__class__.__name__,": path= %s"%self.path)

    self.writer=writer
    self.reader=reader
    self.trialfunc=trialfunc
    self.qwalk=qwalk
    if runner is not None: self.runner=runner
    else: self.runner=Runner.RunnerPBS()
    self.bundle=bundle

    self.completed=False
    self.scriptfile=None
    self.bundle_ready=False
    self.infile="%s.%s"%(name,writer.qmc_abr)
    self.outfile="%s.o"%self.infile
    self.stdout="%s.out"%self.infile

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      print(self.__class__.__name__,": Recovering old manager.")
      old=pkl.load(open(self.path+self.pickle,'rb'))
      old.update_options(self)
      self=old

    # Update the file.
    if not os.path.exists(self.path): os.mkdir(self.path)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #------------------------------------------------
  def update_options(self,other):
    ''' Safe copy options from other to self. 
    
    Args: 
      other (QWRunManager): New object to copy attributes from.
    '''

    update_attributes(old=self.runner,new=other.runner,
        safe_keys=['queue','walltime','np','nn','jobname','prefix','postfix'],
        skip_keys=['queueid'])

    # trialfunc gets updated as the manager generating it finishes its calculation.
    update_attributes(old=self.writer,new=other.writer,
        skip_keys=['completed','trialfunc'])

    # Update the file.
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    
  #------------------------------------------------
  def nextstep(self):
    ''' Perform next step in calculation. trialfunc managers are updated if they aren't completed yet.'''

    print(self.__class__.__name__,": next step.")

    # Check dependency is completed first.
    if self.writer.trialfunc=='':
      self.writer.trialfunc=self.trialfunc.export(self.path)

    # Work on this job.
    cwd=os.getcwd()
    os.chdir(self.path)

    dir(self.reader)

    # Write the input file.
    if not self.writer.completed:
      self.writer.qwalk_input(self.infile)
    
    status=resolve_status(self.runner,self.reader,self.outfile)
    print(self.__class__.__name__,": %s status= %s"%(self.name,status))
    if status=="not_started" and self.writer.completed:
      exestr="%s %s &> %s"%(self.qwalk,self.infile,self.stdout)
      self.runner.add_task(exestr)
      print(self.__class__.__name__,": %s status= submitted"%(self.name))
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.reader.collect(self.outfile)
      if status=='ok':
        print(self.__class__.__name__,": %s status= %s, task complete."%(self.name,status))
        self.completed=True
      else:
        print(self.__class__.__name__,": %s status= %s, attempting rerun."%(self.name,status))
        exestr="%s %s"%' '.join(self.qwalk,self.infile)
        exestr+=" &> %s.out"%self.infile[-1]
        self.runner.add_task(exestr)
    elif status=='done':
      self.completed=True
    else:
      print(self.__class__.__name__,": %s status= %s"%(self.name,status))

    # Ready for bundler or else just submit the jobs as needed.
    if self.bundle:
      self.scriptfile="%s.run"%self.name
      self.bundle_ready=self.runner.script(self.scriptfile,self.driverfn)
    else:
      qsubfile=self.runner.submit(self.name)

    # Update the file.
    with open(self.pickle,'wb') as outf:
      pkl.dump(self,outf)

    os.chdir(cwd)

  #------------------------------------------------
  def update_queueid(self,qid):
    ''' If a bundler handles the submission, it can update the queue info with this.
    Args:
      qid (str): new queue id from submitting a job. The Manager will check if this is running.
    '''
    self.runner.queueid.append(qid)
    self._runready=False # After running, we won't run again without more analysis.

    # Update the file.
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

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
