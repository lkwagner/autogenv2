from __future__ import print_function
####################################################
class VarianceWriter:
  def __init__(self,options={}):
    self.qmc_type='Variance optimization'
    self.sysfiles=['qw_000.sys']
    self.slaterfiles=['qw_000.slater']
    self.jastfiles=['qw.jast2']
    self.basenames=['qw_000']
    self.iterations=10
    self.macro_iterations=3
    self.completed=False
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
    skipkeys = ['completed','sysfiles','slaterfiles','jastfiles','basenames']
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
    assert nfiles==len(self.slaterfiles)
    assert nfiles==len(self.jastfiles)
    for i in range(nfiles):
      sys=self.sysfiles[i]
      slater=self.slaterfiles[i]
      jast=self.jastfiles[i]
      base=self.basenames[i]
      
      with open(base+'.variance','w') as f:
        for j in range(self.macro_iterations):
          f.write("method { optimize iterations %i }\n"%self.iterations)
        f.write("include "+sys+"\n")
        f.write("trialfunc { slater-jastrow \n")
        f.write("wf1 { include "+slater + "}\n")
        f.write("wf2 { include " + jast + "}\n")
        f.write("}\n")
    infiles=[x+".variance" for x in self.basenames]
    outfiles=[x+".o" for x in infiles]        
    self.completed=True
    return infiles,outfiles

     
####################################################
class VarianceReader:
  def __init__(self):
    self.output={}
    self.completed=False

  def read_outputfile(self,outfile):
    ret={}
    ret['sigma']=[]
    with open(outfile) as f:
      for line in f:
        if 'dispersion' in line:
          ret['sigma'].append(float(line.split()[4]))
    return ret
          
  #------------------------------------------------
  def collect(self,outfiles):
    self.completed=True
    for f in outfiles:
      if f not in self.output.keys():
        self.output[f]=[]
      results=self.read_outputfile(f)
      self.output[f].append(results)
      # Minimal error checking.
      self.completed=(self.completed and len(results)>1)
      
  #------------------------------------------------
  def write_summary(self):
    print("#### Variance optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("dispersion",run['sigma'])
      
      
      

