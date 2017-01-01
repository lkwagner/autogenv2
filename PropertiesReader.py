import os 

class PropertiesReader:
  """ Gets results of a properties run. """
  def __init__(self):
    self.completed=False
    self.out={}
#-------------------------------------------------      
  def collect(self,outfilename):
    """ Collect results from output."""
    # TODO actually gather results and check if run is successful
    if os.path.isfile(outfilename):
      self.completed=True
    else:
      self.completed=False
