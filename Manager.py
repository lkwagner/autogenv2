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

  ####################################
  def nextstep(self):
    """ Check write status, then if it's running. If not running check if
    finished. If not finished, attempt to run. """ 

    # Generate input files.
    if not self.writer.completed:
      with open(self.crysinpfn,'w') as f:
        # Change: Need to have writer do the writing otherwise how does it *know* that
        # it's actually done, if it just returns the string to write?
        self.writer.write_crys_input(self.crysinpfn)
      with open(self.propinpfn,'w') as f:
        self.writer.write_prop_input(self.propinpfn)

    if not self.creader.completed:
      self.crunner.run(self.crysinpfn,self.crysoutfn)
      self.creader.collect(self.crysoutfn)
    print("Crystal done: ",self.creader.completed)

    if not self.preader.completed:
      self.prunner.run(self.propinpfn,self.propoutfn)
      self.preader.collect(self.propoutfn)
    print("Crystal properties done: ",self.preader.completed)

    if self.creader.completed and self.preader.completed:
      self.completed=True

  def check_consistent(self,other):
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)

  ####################################
  def to_json(self):
    raise NotImplementedError
