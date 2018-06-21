from manager_tools import resolve_status, update_attributes
from autopyscf import PySCFReader,dm_from_chkfile
from autorunner import PySCFRunnerPBS
import os
import shutil as sh 
import pickle as pkl
import pyscf2qwalk
from autopaths import paths

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
      path=os.getcwd()
    if path[-1]!='/': path+='/'
    self.path=path

    self.logname="%s@%s"%(self.__class__.__name__,self.path+self.name)

    #print(self.logname,": initializing")

    self.writer=writer
    if reader is not None: self.reader=reader
    else: self.reader=PySCFReader()
    if runner is not None: self.runner=runner
    else: self.runner=PySCFRunnerPBS()
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
          self.writer.dm_generator=dm_from_chkfile("%d.%s"%(self.restarts,self.chkfile))
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
      qsubfile=self.runner.submit(jobname=self.path.replace('/','-')+self.name,ppath=[paths['pyscf']])

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
    ''' Export QWalk input files into current directory.
    Returns:
      bool: whether it was successful.'''
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
