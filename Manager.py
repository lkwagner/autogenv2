
#######################################################################
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

    if not self.preader.completed:
      self.prunner.run(self.propinpfn,self.propoutfn)
      self.preader.collect(self.propoutfn)
    print("Crystal properties done: ",self.preader.completed)

    if self.creader.completed and self.preader.completed:
      self.completed=True

  def is_consistent(self,other):
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)

  #----------------------------------------
  def to_json(self):
    raise NotImplementedError

#######################################################################

class QWalkfromCrystalManager:
  """Set up and run a QWalk job from Crystal"""
  #------------------------------------------------
  def __init__(self,convert_runner,convert_checker):
    self.convert_runner=convert_runner
    self.convert_checker=convert_checker
  
  #------------------------------------------------
  def nextstep(self):
    if not self.convert_checker.completed:
      self.convert_runner.run()
      self.convert_checker.collect()
