import os
def resolve_status(runner,reader,outfiles):
  #Check if the reader is done
  if reader.completed:
    return 'done'

  #Check if the job is in the queue
  # or running. If so, we just return that.
  currstat=runner.check_status()
  print("Current status:",currstat)
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
  

#######################################################################
class CrystalManager:
  """ Internal class managing process of running a DFT job though Crystal.
  Has authority over file names associated with this task."""
  def __init__(self,writer,crys_runner,crys_reader,prop_runner,prop_reader):
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


    #Check on the CRYSTAL run
    status=resolve_status(self.crunner,self.creader,[self.crysoutfn])
    
    print("status",status)
    if status=="running":
      return
    elif status=="not_started":
      self.crunner.run(self.crysinpfn,self.crysoutfn)
      return
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      self.creader.collect(self.crysoutfn)


    if not self.preader.completed:
      self.prunner.run(self.propinpfn,self.propoutfn)
      self.preader.collect(self.propoutfn)
    print("Crystal properties done: ",self.preader.completed)

    if self.creader.completed and self.preader.completed:
      self.completed=True
  #----------------------------------------
  def is_consistent(self,other):
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)

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
    self.convert_runner=convert_runner
    self.convert_checker=convert_checker
    
  def is_consistent(self,other):
    # This documents what needs to be checked.
    return self.convert_runner.is_consistent(other.convert_runner)
  
  #------------------------------------------------
  def nextstep(self):
    if not self.convert_checker.completed:
      self.convert_runner.run()
      self.convert_checker.collect()
  #------------------------------------------------

  def write_summary(self):
    print("K-points",self.convert_checker.out)

  #----------------------------------------
  def status(self):
    if self.convert_checker.completed:
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
  def nextstep(self):
    if not self.writer.completed:
      self.infiles,self.outfiles=self.writer.qwalk_input()
    
    status=resolve_status(self.runner,self.reader,self.outfiles)
    
    print("status",status)
    if status=="running":
      return
    elif status=="not_started":
      self.runner.run(self.infiles,self.outfiles)
      return
    elif status=="ready_for_analysis":
      #This is where we (eventually) do error correction and resubmits
      self.reader.collect(self.outfiles)
      
  #------------------------------------------------

  def write_summary(self):
    self.reader.write_summary()
    

  #----------------------------------------
  def status(self):
    if self.reader.completed:
      return 'ok'
    else:
      return 'not_finished'
    
    
