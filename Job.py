import Manager as mgmt
from copy import deepcopy
from CrystalWriter import CrystalWriter
from CrystalReader import CrystalReader
from CrystalRunner import LocalCrystalRunner
from PropertiesReader import PropertiesReader
from PropertiesRunner import LocalPropertiesRunner


#########################################################
class Job:
  """ Contains lists DFT and QMC steps in the order they need to be performed.

  * One job per directory. 
  * For most users, a child class of this should be used."""
  
  def __init__(self,jobid,jobplans,managers):
    self.jobid=jobid
    self.managers=managers
    self.picklefn="%s.pickle"%jobid
#---------------------------------------
  def is_consistent(self,other):
    result=True
    if len(other.managers)!=len(self.managers):
      print('You have added or removed tasks for this job.')
      result=False


    for rec_manager,plan_manager in zip(other.managers,self.managers):
      plancheck=plan_manager.is_consistent(rec_manager)
      if plancheck==False:
        print('You have modified a job.')
        result=False
    return result
#---------------------------------------

  def nextstep(self):
    for manager in self.managers:
      manager.nextstep()


##########################################################
class LocalCrystalDFT(Job):
  """ An example of a Job where it simply calculates a crystal DFT job and
  quits. """
  
  #TODO sometimes you don't need properties?
  def __init__(self,jobid,struct,crystal_opts,structtype='cif'):
    # May have it automatically detect file type? Probably wouldn't be too hard.
    inpcopy=deepcopy(crystal_opts)
    self.jobid=jobid

    #TODO primitive option.
    cwriter=CrystalWriter()
    if structtype=='cif':
      cwriter.set_struct_fromcif(struct)
    elif structtype=='xyz':
      cwriter.set_struct_fromxyz(struct)
    else:
      raise ValueError("structtype not recognized.")
    cwriter.set_options(crystal_opts)


    # For this simple case, only one Manager is needed.
    self.managers=[mgmt.CrystalManager(
        cwriter,
        CrystalReader(),
        LocalCrystalRunner(),
        PropertiesReader(),
        LocalPropertiesRunner()
      )]
    self.picklefn="%s.pickle"%jobid

