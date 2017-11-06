from __future__ import print_function
import os
####################################################
class VarianceWriter:
  def __init__(self,options={}):
    self.qmc_type='Variance optimization'
    self.sysfiles=['qw_000.sys']
    self.slaterfiles=['qw_000.slater']
    self.jastfiles=['qw.jast2']
    self.errtol=10
    self.minblocks=0
    #self.basenames=['qw_000']
    self.iterations=10
    self.macro_iterations=3
    self.completed=False
    self.qmc_abr='variance'
    self.set_options(options)
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for VarianceWriter")
        raise AssertionError
      selfdict[k]=d[k]
  #-----------------------------------------------
  def is_consistent(self,other):
    #In principle we should check for the files, but 
    #they are often determined *after* the plan has been 
    #written so it's not currently practical to check them.
    skipkeys = ['completed','sysfiles','slaterfiles','jastfiles']
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
    assert nfiles==len(self.slaterfiles), "Check slaterfiles"
    assert nfiles==len(self.jastfiles), "Check jastfiles"

    for inp,sys,slater,jast in zip(infiles,self.sysfiles,self.slaterfiles,self.jastfiles):
      with open(inp,'w') as f:
        for j in range(self.macro_iterations):
          f.write("method { optimize iterations %i }\n"%self.iterations)
        f.write("include "+sys+"\n")
        f.write("trialfunc { slater-jastrow \n")
        f.write("wf1 { include "+slater + "}\n")
        f.write("wf2 { include " + jast + "}\n")
        f.write("}\n")
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
    complete={}
    for fname,results in self.output.items():
      complete[fname]=True
      if results['sigma'][-1] > self.vartol:
        print("Variance optimize incomplete: variance ({}) does not meet tolerance ({})"\
            .format(results['sigma'],self.vartol))
        complete[fname]=False
      if len(results['sigma']) < self.minsteps:
        print("Variance optimize incomplete: number of steps (%f) less than minimum (%f)"%\
            (len(results['sigma']),self.minsteps))
        complete[fname]=False
      if (results['sigma'][-1]-results['sigma'][-2]) > self.vardifftol:
        print("Variance optimize incomplete: change in variance (%f) less than tolerance (%f)"%\
            (results['sigma'],self.vardifftol))
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
    print("#### Variance optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("dispersion",run['sigma'])
      
      
      

