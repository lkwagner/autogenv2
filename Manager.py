# Overall comments:
# - Can we assert that there is only one runner for each manager? That way
#   bundling is handled the same as nonbundling: the details of each behavior is
#   in the runner.

#############################
class CrystalManager:
  """ Internal class managing process of running a DFT job though Crystal.
  Has authority over file names associated with this task."""
  def __init__(self,writer,crys_reader,crys_runner,prop_reader,prop_runner):
    self.writer=writer
    self.creader=crys_reader
    self.crunner=crys_runner
    self.preader=prop_reader
    self.prunner=prop_runner
    self.crysinpfn='crys.in'
    self.crysoutfn='crys.in.o'
    self.propinpfn='prop.in'
    self.propoutfn='prop.in.o'
    self.completed=False

  #----------------------------------------
  def nextstep(self):
    """ Check write status, then if it's running. If not running check if
    finished. If not finished, attempt to run. """ 

    # Generate input files.
    if not self.writer.completed:
      with open(self.crysinpfn,'w') as f:
        self.writer.write_crys_input(self.crysinpfn)
      with open(self.propinpfn,'w') as f:
        self.writer.write_prop_input(self.propinpfn)

    if not self.creader.completed:
      self.crunner.run(self.crysinpfn,self.crysoutfn)
      self.creader.collect(self.crysoutfn)
    print("Crystal done: ",self.creader.completed)

    if self.creader.completed and not self.preader.completed:
      self.prunner.run(self.propinpfn,self.propoutfn)
      self.preader.collect(self.propoutfn)
    print("Crystal properties done: ",self.preader.completed)

    if self.creader.completed and self.preader.completed:
      self.completed=True

  #-----------------------------------
  def is_consistent(self,other):
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)

  #-----------------------------------
  def get_fileinfo(self):
    """ Export information about file name conventions for use in things like
    conversion. """
    return {
        'crysoutfn':self.crysoutfn,
        'propoutfn':self.propoutfn
      }

  #----------------------------------------
  def to_json(self):
    raise NotImplementedError

#############################
class QMCManager:
  """ Internal class managing process of running a set of QMC jobs though QWalk.
  Has authority over file names associated with this task.
  
  * Readers and writers in the list correspond to QMC jobs that can run in
    parallel, which are only grouped because they depend on the same trial
    wave function.

  * Different QWalkManagers can point to the same QWalkConverter, so that it
    doesn't need to be run for each job.
    """
  #-----------------------------------
  def __init__(self,writers,readers,runner,qw_base="qs_{}"):
    """ Needs writers and runners, plus the qw_base string, which should take a
    format argument for runids (labeling variations of QMC parameters). 

    qw_base should be decided by the converter script, which passes this info to
    the Job object, which passes that info here."""
    self.writer=writers
    self.reader=readers
    self.runner=runner
    self.completed=False
    self.current_runs=[]
    self.base=qw_base
    self.runids=list(range(len(writers)))

    assert len(writers)==len(readers),\
        'Need same number of writers and runners.'

  #-----------------------------------
  def nextstep(self):
    """ Check write status, then if it's running. If not running check if
    finished. If not finished, attempt to run. """ 

    self.current_runs=[]
    self.completed=True

    for writer,reader,runid in zip(self.writers,self.runners,self.runids):
      inpfn=self.qw_base.format(runid)
      # Generate input files.
      if not self.writer.completed:
        self.writer.write_input(inpfn)

      if not self.reader.completed:
        self.current_runs.append(inpfn)
        self.completed=False
      else:
        self.reader.collect(self.crysoutfn)

      print("QWalk file {} done: ".format(inpfn),self.reader.completed)
    self.runner.run(self.current_runs,["%s.out"%rfn for rfn in self.current_runs])

  #-----------------------------------
  def is_consistent(self,other):
    consistent=True
    # TODO I'd be nice if the order of the readers and writer didn't matter. Maybe
    # find some canonical sorting routine. 
    for swriter,owriter in zip(self.writers,other.writers):
      consistent = consistent and swriter.is_consistent(owriter)
    return consistent

  #-----------------------------------
  def to_json(self):
    raise NotImplementedError


