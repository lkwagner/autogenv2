import os 

class PropertiesReader:
  """ Gets results of a properties run. """
  def __init__(self):
    self.completed=False
    self.out={}
#-------------------------------------------------      
  def collect(self):
    """ Collect results from output."""
    if not self.completed and os.path.isfile('autogen.d12.o'):
      f = open('autogen.d12.o', 'r')
      lines = f.readlines()
      for li,line in enumerate(lines):
        if 'SCF ENDED' in line:
          self.out['total_energy']=float(line.split()[8])    
        elif 'TOTAL ATOMIC SPINS' in line:
          moms = []
          shift = 1
          while "TTT" not in lines[li+shift]:
            moms += map(float,lines[li+shift].split())
            shift += 1
          self.out['mag_moments']=moms
    self.completed=True
    return 'ok'
      
#-------------------------------------------------      
  def check_status(self,outfilename):
    """ Decide status of crystal run. """

    if os.path.exists(outfilename):
      return 'ok'
    else:
      return 'not_started'
    
