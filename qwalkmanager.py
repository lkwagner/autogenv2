from manager_tools import resolve_status, update_attributes, seperate_jastrow
from autorunner import RunnerPBS
import os
import pickle as pkl
from autopaths import paths

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
    else: self.runner=RunnerPBS()
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
      qsubfile=self.runner.submit(self.path.replace('/','-')+self.name)

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
    ''' Store resulting wave function into self.qwfiles['wfout']. Extract Jastrow and store in self.qwfiles['jastrow2']
    Returns:
      bool: Whether it was successful.'''
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
      newjast=separate_jastrow(self.qwfiles['wfout'])
      self.qwfiles['jastrow2']="%s.jast"%self.infile
      with open(self.qwfiles['jastrow2'],'w') as outf:
        outf.write(newjast)
      os.chdir(cwd)

    with open(self.path+self.pickle,'wb') as outf:
      pkl.dump(self,outf)
    return True
