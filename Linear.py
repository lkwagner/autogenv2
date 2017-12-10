from __future__ import print_function
import os
####################################################
class LinearWriter:
  def __init__(self,options={}):
    self.qmc_type='Linear optimization'
    self.sysfiles=['qw_000.sys']
    self.wffiles=[]
    #self.basenames=['qw_000']
    self.completed=False
    self.errtol=10
    self.minblocks=0
    self.total_nstep=2048*4 # 2048 gets stuck pretty often.
    self.total_fit=2048
    self.qmc_abr='energy'
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
  def qwalk_input(self,infiles):
    nfiles=len(infiles)
    assert nfiles==len(self.sysfiles), "Check sysfiles"
    assert nfiles==len(self.wffiles), "Check wffiles"

    for inp,sys,wf in zip(infiles,self.sysfiles,self.wffiles):
      
      with open(inp,'w') as f:
        f.write("method { linear \n")
        f.write("total_nstep %i \n"%self.total_nstep)
        f.write("total_fit %i \n"%self.total_fit)
        f.write("}\n")
        f.write("include "+sys+"\n")
        f.write("trialfunc { include %s\n"%wf)
        f.write("}\n")
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
    complete={}
    for fname,results in self.output.items():
      complete[fname]=True
      print(results['energy'])
      if len(results['energy']) < self.minsteps:
        print("Linear optimize incomplete: number of steps (%f) less than minimum (%f)"%\
            (len(results['energy']),self.minsteps))
        complete[fname]=False
      else:
        ediff=results['energy'][-1]-results['energy'][-2]
        ediff_err=(results['energy_err'][-1]**2 + results['energy_err'][-2]**2)**0.5
        if ediff > self.sigtol*ediff_err:
          print("Linear optimize incomplete: change in energy (%.5f) less than tolerance (%.2f*%.2f=%.5f)"%\
              (ediff,self.sigtol,ediff_err,self.sigtol*ediff_err))
          complete[fname]=False
    return complete
          
  #------------------------------------------------
  def collect(self,outfiles,errtol=None,minblocks=None):
    ''' Collect results for each output file and resolve if the run needs to be resumed. 

    Args: 
      outfiles (list): list of output file names to open and read.
    Returns:
      str: status of run = {'ok','restart'}
    '''
    # Gather output from files.
    self.completed=True
    status='unknown'
    for f in outfiles:
      if os.path.exists(f):
        self.output[f]=self.read_outputfile(f)

    # Check files.
    file_complete=self.check_complete()
    self.completed=all([c for f,c in file_complete.items()])
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
      
      
      

