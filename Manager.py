import os
import numpy as np
def resolve_status(runner,reader,outfiles):
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

  
  #We are in an error state 
  for outf in outfiles:
    lines = open(outf, 'r').read().split()
    if ('HF_done' in  lines) and ('All_done' not in lines):
      return  "retry"

  #we haven't collected the results. 
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
    while True:
      status=resolve_status(self.crunner,self.creader,[self.crysoutfn])
    
      print("Crystal status",status)
      if status=="running":
        return
      elif status=="not_started":
        self.crunner.run(self.crysinpfn,self.crysoutfn)
      elif status=="ready_for_analysis":
        #This is where we (eventually) do error correction and resubmits
        self.creader.collect(self.crysoutfn)
        break
      elif status=='done':
        break
      else:
        return


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
    self.runner=convert_runner
    self.reader=convert_checker
    
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
  def is_consistent(self,other):
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)
    
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
        self.runner.run(self.infiles,stdoutfiles)
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
    self.infiles=[]
    self.outfiles=[]
  #------------------------------------------------
    
  def is_consistent(self,other):
    # This documents what needs to be checked.
    return self.writer.is_consistent(other.writer)
    
  #------------------------------------------------
  def nextstep(self):
    if not self.writer.completed:
      self.infiles,self.outfiles,self.chkfiles=self.writer.pyscf_input(self.driverfn)
    
    while True:
      status=resolve_status(self.runner,self.reader,self.outfiles)
      print("PySCF status",status)
      if status=="running":
        return
      elif status=="not_started":
        self.runner.run(self.infiles,self.outfiles)
      elif status=="ready_for_analysis":
        #This is where we (eventually) do error correction and resubmits
        self.reader.collect(self.outfiles,self.chkfiles)
        break
      elif status=='done':
        break
      #If we need to restart the run
      elif status=='retry':
        self.writer.update_input(self.driverfn)
        self.runner.run(self.infiles, self.outfiles)
        break
      else:
        return
      
  #------------------------------------------------

  def write_summary(self):
    self.reader.write_summary()
    

  #----------------------------------------
  def status(self):
    current_status = resolve_status(self.runner,self.reader,self.outfiles)
    if current_status == 'done':
      return 'ok'
    elif current_status == 'retry':
      return 'retry'
    else:
      return 'not_finished' 
