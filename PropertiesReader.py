import os 
from crystal2qmc import read_gred,read_kred,read_outputfile

class PropertiesReader:
  """ Gets results of a properties run. """
  def __init__(self):
    self.completed=False
    self.out={}
#-------------------------------------------------      
  def collect(self,outfilename):
    """ Collect results from output."""
    if os.path.isfile(outfilename):
      info, lat_parm, ions, basis, pseudo = read_gred()
      eigsys = read_kred(info,basis)

      if eigsys['nspin'] > 1:
        eigsys['totspin'] = read_outputfile(outfilename)
      else:
        eigsys['totspin'] = 0

      # Some of this info is redundant with the input. 
      # But having all here makes conversion simpler.
      # This data can be fed directly into the write_files in crystal2qmc.
      self.out.update({
          'lat_parm':lat_parm,
          'ions':ions,
          'basis':basis,
          'pseudo':pseudo,
          'eigsys':eigsys
        })

      self.completed=True
    else:
      self.completed=False
