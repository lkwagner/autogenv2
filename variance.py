from __future__ import print_function
import os
####################################################
class VarianceWriter:
  def __init__(self,options={}):
    ''' Object for producing input into a variance optimization QWalk run. 
    Args:
      options (dict): editable options are as follows.
        trialfunc (str): system and trial wavefunction section.
        errtol (float): tolerance for the variance. 
        minblocks (int): minimum number of VMC steps to take.
        iterations (int): number of VMC steps to attempt.
        macro_iterations (int): Number of optimize calls to make.
    '''
    self.trialfunc=''
    self.errtol=10
    self.minblocks=0
    self.iterations=10
    self.macro_iterations=3

    self.qmc_type='Variance optimization'
    self.qmc_abr='variance'
    self.completed=False
    self.set_options(options)
    
  #-----------------------------------------------
  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        raise AssertionError("Error:",k,"not a keyword for VarianceWriter.")
      selfdict[k]=d[k]
    
  #-----------------------------------------------
  def qwalk_input(self,infile):
    if self.trialfunc=='':
      print(self.__class__.__name__,": Trial function not ready. Postponing input file generation.")
      self.completed=False
    else:
      with open(infile,'w') as f:
        for j in range(self.macro_iterations):
          f.write("method { optimize iterations %i }\n"%self.iterations)
        f.write(self.trialfunc)
      self.completed=True
     
####################################################
class VarianceReader:
  def __init__(self,vartol=10,vardifftol=0.1,minsteps=2):
    self.output={}
    self.completed=False

    self.vartol=vartol
    self.vardifftol=vardifftol
    self.minsteps=minsteps

  #------------------------------------------------
  def read_outputfile(self,outfile):
    ret={}
    ret['sigma']=[]
    with open(outfile,'r') as f:
      for line in f:
        if 'dispersion' in line:
          ret['sigma'].append(float(line.split()[4]))
    return ret

  #------------------------------------------------
  def check_complete(self):
    ''' Check if a variance optimize run is complete.
    Returns:
      bool: If self.results are within error tolerances.
    '''
    if len(self.output['sigma']) < self.minsteps:
      print(self.__class__.__name__,": Variance optimize incomplete: number of steps (%f) less than minimum (%f)"%\
          (len(self.output['sigma']),self.minsteps))
      return False
    if self.output['sigma'][-1] > self.vartol:
      print(self.__class__.__name__,": Variance optimize incomplete: variance ({}) does not meet tolerance ({})"\
          .format(self.output['sigma'],self.vartol))
      return False
    if (self.output['sigma'][-1]-self.output['sigma'][-2]) > self.vardifftol:
      print(self.__class__.__name__,": Variance optimize incomplete: change in variance (%f) less than tolerance (%f)"%\
          (self.output['sigma'],self.vardifftol))
      return False
    return True
          
  #------------------------------------------------
  def collect(self,outfile,errtol=None,minblocks=None):
    ''' Collect results for output file and resolve if the run needs to be resumed. 

    Args: 
      outfile (str): Output file name to open and read.
    Returns:
      str: status of run = {'ok','restart'}
    '''
    # Gather output from files.
    status='unknown'
    if os.path.exists(outfile):
      self.output=self.read_outputfile(outfile)
      self.output['file']=outfile

    # Check files.
    self.completed=self.check_complete()
    if not self.completed:
      status='restart'
    else:
      status='ok'
    return status

  #------------------------------------------------
  def write_summary(self):
    print("#### Variance optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("dispersion",run['sigma'])
