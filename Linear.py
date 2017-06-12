from __future__ import print_function
####################################################
class LinearWriter:
  def __init__(self,options={}):
    self.qmc_type='Linear optimization'
    self.sysfiles=['qw_000.sys']
    self.wffiles=[]
    self.basenames=['qw_000']
    self.completed=False
    self.total_nstep=2048
    self.total_fit=2048
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
  def qwalk_input(self):
    nfiles=len(self.sysfiles)
    for i in range(nfiles):
      sys=self.sysfiles[i]
      wf=self.wffiles[i]
      base=self.basenames[i]
      
      with open(base+'.energy','w') as f:
        f.write("method { linear \n")
        f.write("total_nstep %i \n"%self.total_nstep)
        f.write("total_fit %i \n"%self.total_fit)
        f.write("}\n")
        f.write("include "+sys+"\n")
        f.write("trialfunc { include %s\n"%wf)
        f.write("}\n")
    infiles=[x+".energy" for x in self.basenames]
    outfiles=[x+".o" for x in infiles]        
    self.completed=True
    return infiles,outfiles

     
####################################################
class LinearReader:
  def __init__(self):
    self.output={}
    self.completed=False

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
  def collect(self,outfiles):
    self.completed=True
    for f in outfiles:
      if f not in self.output.keys():
        self.output[f]=[]
      results=self.read_outputfile(f)
      self.output[f].append(results)
      # minimal error checking.
      self.completed=(self.completed and len(results)>1)
      
  #------------------------------------------------
  def write_summary(self):
    print("#### Linear optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("energy",run['energy'])
        print("energy_err",run['energy_err'])
      
      
      

