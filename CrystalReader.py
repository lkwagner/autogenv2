import os 

class CrystalReader:
  """ Tries to extract properties of crystal run, or else diagnose what's wrong. """
  def __init__(self):
    self.completed=False
    self.out={}
#-------------------------------------------------      
  def collect(self,outfilename):
    """ Collect results from output."""
    if os.path.isfile(outfilename):
      f = open(outfilename, 'r')
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
    else:
      # Just to be sure/clear...
      self.completed=False
      
#-------------------------------------------------      
  # This can be made more efficient if it's a problem: searches whole file for
  # each query.
  def check_outputfile(self,outfilename,acceptable_scf=10.0):
    """ Check output file. 

    Current return values:
    no_record, not_started, ok, too_many_cycles, finished (fall-back),
    scf_fail, not_enough_decrease, divergence, not_finished
    """
    if os.path.isfile(outfilename):
      outf = open(outfilename,'r')
    else:
      return "not_started"

    outlines = outf.read().split('\n')
    reslines = [line for line in outlines if "ENDED" in line]

    if len(reslines) > 0:
      if "CONVERGENCE" in reslines[0]:
        return "ok"
      elif "TOO MANY CYCLES" in reslines[0]:
        print("CrystalRunner: Too many cycles.")
        return "too_many_cycles"
      else: # What else can happen?
        print("CrystalRunner: Finished, but unknown state.")
        return "finished"
      
    detots = [float(line.split()[5]) for line in outlines if "DETOT" in line]
    if len(detots) == 0:
      print("CrystalRunner: Last run completed no cycles.")
      return "scf_fail"

    detots_net = sum(detots[1:])
    if detots_net > acceptable_scf:
      print("CrystalRunner: Last run performed poorly.")
      return "not_enough_decrease"

    etots = [float(line.split()[3]) for line in outlines if "DETOT" in line]
    if etots[-1] > 0:
      print("CrystalRunner: Energy divergence.")
      return "divergence"
    
    print("CrystalRunner: Not finished.")
    return "not_finished"
  
#-------------------------------------------------      
  def check_status(self,outfilename):
    """ Decide status of crystal run. """

    status=self.check_outputfile(outfilename)
    print("status",status)
    return status
    
