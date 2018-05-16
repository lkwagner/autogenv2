''' Manager classes handle the file maintainance of various runs.

They:
  - Call write, run, and read in the appropriate order.
  - Keep a record of progress on disk.
  - Hold a list of all file paths together in one place.
  - Report progress to the user.

All choices revolving around file names should be handled here.

Tips on avoiding bugs with the write/run/read cycle:
  - Every time data is updated, write the updates to the disk.
  - Every time the data is used, update the Manager from the disk.
'''
import crystal
import propertiesreader
import os
import autopyscf 
import shutil as sh
import numpy as np
import pyscf2qwalk
import crystal2qmc
import pickle as pkl
import autorunner
from autopaths import paths

from copy import deepcopy

#TODO These should all probably be subclasses, if not just for logical purposes. 
#TODO slow queue checking might be made more efficient by saving the qstat results somewhere.

######################################################################
# Misc functions. 
######################################################################
def resolve_status(runner,reader,outfile):
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
def deep_compare(d1,d2):
  '''I have to redo dict comparison because numpy will return a bool array when comparing.'''
  if type(d1)!=type(d2):
    return False
  if type(d1)==dict:
    if d1.keys()!=d2.keys():
      return False
    allsame=True
    for key in d1.keys():
      allsame=allsame and deep_compare(d1[key],d2[key])
    return allsame
  else:
    try:
      return np.array_equal(d1,d2)
    except TypeError:
      return d1==d2

######################################################################
def update_attributes(copyto,copyfrom,skip_keys=[],take_keys=[]):
  ''' Save update of class attributes. If copyfrom has additional attributes, they are ignored.

  Args:
    copyto (obj): class who's attributes are being updated.
    copyfrom (obj): class who's attributes will be copied from.
    skip_keys (list): list of attributes (str) not to update. 
    take_keys (list): list of attributes (str) that are ok to update. Others will raise warning and be skipped.
  Returns:
    bool: Whether any changes were made.
  '''
  updated=False
  for key in copyfrom.__dict__.keys():
    if key in skip_keys: 
      #print("Skipping key (%s)"%key)
      pass
    elif key not in copyto.__dict__.keys():
      print("Warning: Object update. An attribute (%s) was skipped because it doesn't exist in both objects."%key)
    elif not deep_compare(copyto.__dict__[key],copyfrom.__dict__[key]):
      if key not in take_keys:
        print("Warning: update to attribute (%s) cancelled, because it requires job to be rerun."%key)
      else:
        #print("Copy",key)
        copyto.__dict__[key]=copyfrom.__dict__[key]
        updated=True
    else:
      #print("Keys match (%s)"%key)
      pass
  return updated

def seperate_jastrow(wffile,optimizebasis=False):
  ''' Seperate the jastrow section of a QWalk wave function file.'''
  # Copied from utils/seperate_jastrow TODO: no copy, bad
  wff=open(wffile,'r')
  tokens=wff.read().split('\n')
  in_jastrow=False
  nopen=0
  nclose=0
  jastlines=[]
  for line in tokens:
    if 'jastrow2' in line.lower():
      in_jastrow=True
    if in_jastrow:
      if not optimizebasis and 'optimizebasis' in line.lower():
        continue
      nopen+=line.count("{")
      nclose+=line.count("}")
    if in_jastrow and nopen >= nclose:
      jastlines.append(line)
  return '\n'.join(jastlines)

######################################################################
# Manager classes.
#######################################################################
class CrystalManager:
  """ Internal class managing process of running a DFT job though crystal.
  Has authority over file names associated with this task."""
  def __init__(self,writer,runner,creader=None,name='crystal_run',path=None,
      preader=None,prunner=None,
      trylev=False,bundle=False,max_restarts=5):
    ''' CrystalManager manages the writing of a Crystal input file, it's running, and keeping track of the results.
    Args:
      writer (PySCFWriter): writer for input.
      reader (PySCFReader): to read PySCF output.
      runner (runner object): to run job.
      creader (CrystalReader): Reads the crystal results, (None implies use default reader).
      preader (PropertiesReader): Reads properties results, if any (None implies use default reader).
      prunner (runner object): run properties if needed (None implies use same runner as crystal).
      name (str): identifier for this job. This names the files associated with run.
      trylev (bool): When restarting use LEVSHIFT option to encourage convergence, then do a rerun without LEVSHIFT.
      bundle (bool): Whether you'll use a bundling tool to run these jobs.
      max_restarts (int): maximum number of times you'll allow restarting before giving up (and manually intervening).
    '''
    # Where to save self.
    self.name=name
    self.pickle="%s.pkl"%self.name

    # Ensure path is set up correctly.
    if path is None:
      path=os.path.getcwd()
    if path[-1]!='/': path+='/'
    self.path=path

    self.logname="%s@%s"%(self.__class__.__name__,self.path+self.name)

    #print(self.logname,": initializing")

    # Handle reader and runner defaults.
    self.writer=writer
    if creader is None: self.creader=crystal.CrystalReader()
    else: self.creader=creader
    if preader is None: self.preader=propertiesreader.PropertiesReader()
    else: self.preader=preader
    if prunner is None: self.prunner=runner
    else: self.prunner=prunner
    if runner is None: self.runner=RunnerPBS()
    else: self.runner=runner

    # Internal.
    self.crysinpfn=self.name
    self.propinpfn=self.name+'.prop'
    self.crysoutfn=self.crysinpfn+'.o'
    self.propoutfn=self.propinpfn+'.o'
    self.restarts=0
    self._runready=False
    self.scriptfile=None
    self.completed=False
    self.bundle=bundle
    self.qwfiles={ 
        'kpoints':[],
        'basis':'',
        'jastrow2':'',
        'orbplot':{},
        'orb':{},
        'sys':{},
        'slater':{}
      }

    # Smart error detection.
    self.trylev=trylev
    self.max_restarts=max_restarts
    self.savebroy=[]
    self.lev=False

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      #print(self.logname,": rebooting old manager.")
      old=pkl.load(open(self.path+self.pickle,'rb'))
      self.recover(old)

    # Update the file.
    if not os.path.exists(self.path): os.mkdir(self.path)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #------------------------------------------------
  def recover(self,other):
    ''' Recover old class by copying over data. Retain variables from old that may change final answer.'''
    # Practically speaking, the run will preserve old `take_keys` and allow new changes to `skip_keys`.
    # This is because you are taking the attributes from the older instance, and copying into the new instance.

    update_attributes(copyto=self,copyfrom=other,
        skip_keys=['writer','runner','creader','preader','prunner','lev','savebroy',
                   'path','logname','name',
                   'trylev','max_restarts','bundle'],
        take_keys=['restarts','completed','qwfiles'])

    # Update queue settings, but save queue information.
    update_attributes(copyto=self.runner,copyfrom=other.runner,
        skip_keys=['queue','walltime','np','nn','jobname'],
        take_keys=['queueid'])
    update_attributes(copyto=self.prunner,copyfrom=other.prunner,
        skip_keys=['queue','walltime','np','nn','jobname'],
        take_keys=['queueid'])

    update_attributes(copyto=self.creader,copyfrom=other.creader,
        skip_keys=[],
        take_keys=['completed','output'])

    update_attributes(copyto=self.preader,copyfrom=other.preader,
        skip_keys=[],
        take_keys=['completed','output'])

    updated=update_attributes(copyto=self.writer,copyfrom=other.writer,
        skip_keys=['maxcycle','edifftol'],
        take_keys=['completed','modisymm','restart','guess_fort','_elements'])
    if updated:
      self.writer.completed=False

  #----------------------------------------
  def nextstep(self):
    ''' Determine and perform the next step in the calculation.'''
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    print(self.logname,": next step.")
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

    # Check on the CRYSTAL run
    status=resolve_status(self.runner,self.creader,self.crysoutfn)
    print(self.logname,": status= %s"%(status))

    if status=="not_started":
      self.runner.add_command("cp %s INPUT"%self.crysinpfn)
      self.runner.add_task("%s &> %s"%(paths['Pcrystal'],self.crysoutfn))

    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.creader.collect(self.crysoutfn)
      print(self.logname,": status %s"%status)
      if status=='killed':
        if self.restarts >= self.max_restarts:
          print(self.logname,": restarts exhausted (%d previous restarts). Human intervention required."%self.restarts)
        else:
          print(self.logname,": attempting restart (%d previous restarts)."%self.restarts)
          self.writer.restart=True
          if self.trylev:
            print(self.logname,": trying LEVSHIFT.")
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
          self.runner.add_task("%s &> %s"%(paths['Pcrystal'],self.crysoutfn))
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
      self.runner.add_task("%s &> %s"%(paths['Pcrystal'],self.crysoutfn))
      self.restarts+=1

    # Ready for bundler or else just submit the jobs as needed.
    if self.bundle:
      self.scriptfile="%s.run"%self.name
      self.bundle_ready=self.runner.script(self.scriptfile)
    else:
      qsubfile=self.runner.submit(self.name)

    self.completed=self.creader.completed

    # Update the file.
    with open(self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    os.chdir(cwd)

  #----------------------------------------
  def collect(self):
    ''' Call the collect routine for readers.'''
    print(self.logname,": collecting results.")
    self.creader.collect(self.path+self.crysoutfn)

    # Update the file.
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

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
    
  #------------------------------------------------
  def export_qwalk(self):
    ''' Export QWalk input files into current directory.'''
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    ready=False
    if len(self.qwfiles['slater'])==0:
      self.nextstep()

      if not self.completed:
        return False

      cwd=os.getcwd()
      os.chdir(self.path)

      print(self.logname,": %s attempting to generate QWalk files."%self.name)

      # Check on the properties run
      status=resolve_status(self.prunner,self.preader,self.propoutfn)
      print(self.logname,": properties status= %s"%(status))
      if status=='not_started':
        ready=False
        self.prunner.add_command("cp %s INPUT"%self.propinpfn)
        self.prunner.add_task("%s &> %s"%(paths['Pproperties'],self.propoutfn))

        if self.bundle:
          self.scriptfile="%s.run"%self.name
          self.bundle_ready=self.prunner.script(self.scriptfile,self.driverfn)
        else:
          qsubfile=self.runner.submit(self.name)
      elif status=='ready_for_analysis':
        self.preader.collect(self.propoutfn)

      if self.preader.completed:
        ready=True
        print(self.logname,": converting crystal to QWalk input now.")
        self.qwfiles=crystal2qmc.convert_crystal(base=self.name,propoutfn=self.propoutfn)
      else:
        print(self.logname,": conversion failed due to problem with properties run. Check for GRED.DAT and KRED.DAT.")

      os.chdir(cwd)
    else:
      ready=True

    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

    return ready
    
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
      runner (Runner object): to run job.
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

    self.logname="%s@%s"%(self.__class__.__name__,self.path+self.name)

    #print(self.logname,": initializing")

    self.writer=writer
    if reader is not None: self.reader=reader
    else: self.reader=autopyscf.PySCFReader()
    if runner is not None: self.runner=runner
    else: self.runner=autorunner.PySCFRunnerPBS()
    self.bundle=bundle

    self.driverfn="%s.py"%name
    self.outfile=self.driverfn+'.o'
    self.chkfile=self.driverfn+'.chkfile'
    self.qwfiles={ 
        'kpoints':[],
        'basis':'',
        'jastrow2':'',
        'orbplot':{},
        'orb':{},
        'sys':{},
        'slater':{}
      }
    self.completed=False
    self.bundle_ready=False
    self.restarts=0

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      print(self.logname,": rebooting old manager.")
      old=pkl.load(open(self.path+self.pickle,'rb'))
      self.recover(old)

    # Update the file.
    if not os.path.exists(self.path): os.mkdir(self.path)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #------------------------------------------------
  def recover(self,other):
    ''' Safe copy options from other to self. '''
    # Practically speaking, the run will preserve old `take_keys` and allow new changes to `skip_keys`.
    # This is because you are taking the attributes from the older instance, and copying into the new instance.
    updated=update_attributes(copyto=self,copyfrom=other,
        skip_keys=['writer','runner','reader', 'path','logname','name','max_restarts','bundle'],
        take_keys=['restarts','completed','qwfiles'])

    update_attributes(copyto=self.runner,copyfrom=other.runner,
        skip_keys=['queue','walltime','np','nn','jobname'],
        take_keys=['queueid'])

    update_attributes(copyto=self.reader,copyfrom=other.reader,
        skip_keys=[],
        take_keys=['completed','output'])

    updated=update_attributes(copyto=self.writer,copyfrom=other.writer,
        skip_keys=['max_cycle'],
        take_keys=['completed','dm_generator'])

    # Update the file.
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    
  #------------------------------------------------
  def nextstep(self):
    ''' Determine and perform the next step in the calculation.'''
    # Recover old data.
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    print(self.logname,": next step.")
    cwd=os.getcwd()
    os.chdir(self.path)

    if not self.writer.completed:
      self.writer.pyscf_input(self.driverfn,self.chkfile)
    
    status=resolve_status(self.runner,self.reader,self.outfile)
    print(self.logname,": %s status= %s"%(self.name,status))

    if status=="not_started":
      self.runner.add_task("python3 %s > %s"%(self.driverfn,self.outfile))
    elif status=="ready_for_analysis":
      status=self.reader.collect(self.outfile,self.chkfile)
      if status=='killed':
        print(self.logname,": attempting restart (%d previous restarts)."%self.restarts)
        sh.copy(self.driverfn,"%d.%s"%(self.restarts,self.driverfn))
        sh.copy(self.outfile,"%d.%s"%(self.restarts,self.outfile))
        sh.copy(self.chkfile,"%d.%s"%(self.restarts,self.chkfile))
        if os.path.exists(self.chkfile):
          self.writer.dm_generator=PySCF.dm_from_chkfile("%d.%s"%(self.restarts,self.chkfile))
        self.writer.pyscf_input(self.driverfn,self.chkfile)
        self.runner.add_task("/usr/bin/python3 %s > %s"%(self.driverfn,self.outfile))
        self.restarts+=1
      elif status=='done':
        print(self.logname,": %s status= %s, task complete."%(self.name,status))

    # Ready for bundler or else just submit the jobs as needed.
    if self.bundle:
      self.scriptfile="%s.run"%self.name
      self.bundle_ready=self.runner.script(self.scriptfile,self.driverfn)
    else:
      qsubfile=self.runner.submit(jobname=self.name,ppath=[paths['pyscf']])

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
    # Recover old data.
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    if len(self.qwfiles['slater'])==0:
      self.nextstep()
      if not self.completed:
        return False
      print(self.logname,": %s generating QWalk files."%self.name)
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
    current_status = resolve_status(self.runner,self.reader,self.outfile)
    if current_status == 'done':
      return 'ok'
    elif current_status == 'retry':
      return 'retry'
    else:
      return 'not_finished'

#######################################################################
class QWalkManager:
  def __init__(self,writer,reader,runner=None,trialfunc=None,
      name='qw_run',path=None,bundle=False):
    ''' QWalkManager managers the writing of a QWalk input files, it's running, and keeping track of the results.
    Args:
      writer (qwalk writer): writer for input.
      reader (qwalk reader): to read job.
      runner (Runner object): to run job.
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

    self.logname="%s@%s"%(self.__class__.__name__,self.path+self.name)

    #print(self.logname,": initializing")

    self.writer=writer
    self.reader=reader
    self.trialfunc=trialfunc
    if runner is not None: self.runner=runner
    else: self.runner=autorunner.RunnerPBS()
    self.bundle=bundle

    self.completed=False
    self.scriptfile=None
    self.bundle_ready=False
    self.infile=name
    self.outfile="%s.o"%self.infile
    # Note: qwfiles stores file names of results, used for exporting trial wave functions.
    self.qwfiles={ 
        'jastrow2':'',
        'wfout':''
      }
    self.stdout="%s.out"%self.infile

    # Handle old results if present.
    if os.path.exists(self.path+self.pickle):
      print(self.logname,": rebooting old manager.")
      old=pkl.load(open(self.path+self.pickle,'rb'))
      self.recover(old)

    # Update the file.
    if not os.path.exists(self.path): os.mkdir(self.path)
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #------------------------------------------------
  def recover(self,other):
    ''' Recover old class by copying over data. Retain variables from old that may change final answer.'''
    # Practically speaking, the run will preserve old `take_keys` and allow new changes to `skip_keys`.
    # This is because you are taking the attributes from the older instance, and copying into the new instance.

    update_attributes(copyto=self,copyfrom=other,
        skip_keys=['writer','runner','reader','path','logname','name','bundle'],
        take_keys=['restarts','completed','trialfunc','qwfiles'])

    # Update queue settings, but save queue information.
    update_attributes(copyto=self.runner,copyfrom=other.runner,
        skip_keys=['queue','walltime','np','nn','jobname'],
        take_keys=['queueid'])

    update_attributes(copyto=self.reader,copyfrom=other.reader,
        skip_keys=[],
        take_keys=['completed','output'])

    updated=update_attributes(copyto=self.writer,copyfrom=other.writer,
        skip_keys=['maxcycle','errtol','minblocks','nblock','savetrace'],
        take_keys=['completed','tmoves','extra_observables','timestep','trialfunc'])
    if updated:
      self.writer.completed=False

  #------------------------------------------------
  def nextstep(self):
    ''' Perform next step in calculation. trialfunc managers are updated if they aren't completed yet.'''
    # Recover old data.
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    print(self.logname,": next step.")

    # Check dependency is completed first.
    if self.writer.trialfunc=='':
      print(self.logname,": checking trial function.")
      self.writer.trialfunc=self.trialfunc.export(self.path)

    # Work on this job.
    cwd=os.getcwd()
    os.chdir(self.path)

    # Write the input file.
    if not self.writer.completed:
      self.writer.qwalk_input(self.infile)
    
    status=resolve_status(self.runner,self.reader,self.outfile)
    print(self.logname,": %s status= %s"%(self.name,status))
    if status=="not_started" and self.writer.completed:
      exestr="%s %s &> %s"%(paths['qwalk'],self.infile,self.stdout)
      self.runner.add_task(exestr)
      print(self.logname,": %s status= submitted"%(self.name))
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      status=self.reader.collect(self.outfile)
      if status=='ok':
        print(self.logname,": %s status= %s, task complete."%(self.name,status))
        self.completed=True
      else:
        print(self.logname,": %s status= %s, attempting rerun."%(self.name,status))
        exestr="%s %s &> %s"%(paths['qwalk'],self.infile,self.stdout)
        self.runner.add_task(exestr)
    elif status=='done':
      self.completed=True

    # Ready for bundler or else just submit the jobs as needed.
    if self.bundle:
      self.scriptfile="%s.run"%self.name
      self.bundle_ready=self.runner.script(self.scriptfile)
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

  #----------------------------------------
  def collect(self):
    ''' Call the collect routine for readers.'''
    print(self.logname,": collecting results.")
    self.reader.collect(self.path+self.outfile)

    # Update the file.
    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)

  #----------------------------------------
  def export_qwalk(self):
    ''' Extract jastrow from the optimization run and return that file name.'''
    # Theoretically more than just Jastrow can be provided, but practically that's the only type of wavefunction we tend to export.

    # Recover old data.
    self.recover(pkl.load(open(self.path+self.pickle,'rb')))

    assert self.writer.qmc_abr!='dmc',"DMC doesn't provide a wave function."

    if self.qwfiles['wfout']=='':
      self.nextstep()
      if not self.completed:
        return False
      print(self.logname,": %s generating QWalk files."%self.name)
      cwd=os.getcwd()
      os.chdir(self.path)
      self.qwfiles['wfout']="%s.wfout"%self.infile
      newjast=seperate_jastrow(self.qwfiles['wfout'])
      self.qwfiles['jastrow2']="%s.jast"%self.infile
      with open(self.qwfiles['jastrow2'],'w') as outf:
        outf.write(newjast)
      os.chdir(cwd)

    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    return True
