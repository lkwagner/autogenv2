import os 
from crystal2qmc import read_gred,read_kred,read_outputfile

class PropertiesReader:
  """ Gets results of a properties run. """
  def __init__(self):
    self.completed=False
    self.out={}
#-------------------------------------------------      
  def collect(self,outfilename):
    """ Just check that results are there. The actual data is too large to want to store."""
    if os.path.isfile("GRED.DAT") and os.path.isfile("KRED.DAT"):
      self.completed=True
    else:
      self.completed=False
