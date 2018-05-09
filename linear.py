from __future__ import print_function
import os
####################################################
class LinearWriter:
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
    self.total_nstep=2048*4 # 2048 gets stuck pretty often.
    self.total_fit=2048
    self.qmc_abr='energy'

    self.qmc_type='Linear optimization'
    self.qmc_abr='energy'
    self.completed=False
    self.set_options(options)

  #-----------------------------------------------
  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for LinearWriter")
        raise InputError
      selfdict[k]=d[k]

  #-----------------------------------------------
  def is_consistent(self,other):
    #In principle we should check for the files, but 
    #they are often determined *after* the plan has been 
    #written so it's not currently practical to check them.
    skipkeys = ['completed','sysfiles','wffiles','basenames']
    for otherkey in other.__dict__.keys():
      if otherkey not in self.__dict__.keys():
        print('other is missing a key.')
        return False
    for selfkey in self.__dict__.keys():
      if selfkey not in other.__dict__.keys():
        print('self is missing a key.')
        return False
    for key in self.__dict__.keys():
      if self.__dict__[key]!=other.__dict__[key] and key not in skipkeys:
        print("Different keys [{}] = \n{}\n or \n {}"\
            .format(key,self.__dict__[key],other.__dict__[key]))
        return False
    return True
    
  #-----------------------------------------------
  def qwalk_input(self,infile):
    if self.trialfunc=='':
      print(self.__class__.__name__,": Trial function not ready. Postponing input file generation.")
      self.completed=False
    else:
      with open(infile,'w') as f:
        f.write("method { linear \n")
        f.write("total_nstep %i \n"%self.total_nstep)
        f.write("total_fit %i \n"%self.total_fit)
        f.write("}\n")
        f.write(self.trialfunc)

    self.completed=True

     
####################################################
class LinearReader:
  def __init__(self,sigtol=2.0,minsteps=2):
    ''' Object for reading, diagnosing, and storing Linear results.

    The arguements control when the object sends a 'restart' flag.

    Args:
      sigtol (float): How many standard errors away from zero to you consider zero energy change?
      minsteps (int): Minimum number of steps allowed for no restart.
    Attributes:
      output (dict): Results for energy, error, and other information.
      completed (bool): Whether no more runs are needed.
    '''
    self.output={}
    self.completed=False
    self.sigtol=sigtol
    self.minsteps=minsteps

  def read_outputfile(self,outfile):
    ret={}
    ret['energy']=[]
    ret['energy_err']=[]
    with open(outfile) as f:
      for line in f:
        if 'current energy' in line:
          ret['energy'].append(float(line.split()[4]))
          ret['energy_err'].append(float(line.split()[6]))
    return ret

  #------------------------------------------------
  def check_complete(self):
    ''' Check if a variance optimize run is complete.
    Returns:
      bool: If self.results are within error tolerances.
    '''
    if len(self.output['energy']) < self.minsteps:
      print(self.__class__.__name__,"Linear optimize incomplete: number of steps (%f) less than minimum (%f)"%\
          (len(self.output['energy']),self.minsteps))
      return False
    else:
      ediff=self.output['energy'][-1]-self.output['energy'][-2]
      ediff_err=(self.output['energy_err'][-1]**2 + self.output['energy_err'][-2]**2)**0.5
      if ediff > self.sigtol*ediff_err:
        print(self.__class__.__name__,"Linear optimize incomplete: change in energy (%.5f) less than tolerance (%.2f*%.2f=%.5f)"%\
            (ediff,self.sigtol,ediff_err,self.sigtol*ediff_err))
        return False
    return True
          
  #------------------------------------------------
  def collect(self,outfile,errtol=None,minblocks=None):
    ''' Collect results for each output file and resolve if the run needs to be resumed. 

    Args: 
      outfile (str): output file to read.
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
    print("#### Linear optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("energy",run['energy'])
        print("energy_err",run['energy_err'])
      
      
      

