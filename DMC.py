from __future__ import print_function
####################################################
class DMCWriter:
  def __init__(self):
    self.sysfiles=['qw_000.sys']
    self.wffiles=[]
    self.basenames=['qw_000']
    self.completed=False
    self.timesteps=[0.01]
    self.nblock=20
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for EnergyWriter")
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
    infiles=[]
    for i in range(nfiles):
      sys=self.sysfiles[i]
      wf=self.wffiles[i]
      base=self.basenames[i]
      
      for t in self.timesteps:
        fname=base+'t'+str(t)+".dmc"
        infiles.append(fname)
        with open(fname,'w') as f:
          f.write("method { dmc timestep %g nblock %i }\n"\
              %(t,self.nblock))
          f.write("include "+sys+"\n")
          f.write("trialfunc { include %s\n"%wf)
          f.write("}\n")
    outfiles=[x+".log" for x in infiles]        
    self.completed=True
    return infiles,outfiles

     
####################################################
import subprocess as sub
import json
class DMCReader:
  def __init__(self):
    self.output={}
    self.completed=False
    self.gosling="gosling"

  def read_outputfile(self,outfile):
    return json.loads(sub.check_output([self.gosling,"-json",outfile]).decode())
          
  #------------------------------------------------
  def collect(self,outfiles):
    for f in outfiles:
      self.output[f]=self.read_outputfile(f)
    self.completed=True
      
  #------------------------------------------------
  def write_summary(self):
    print("#### Diffusion Monte Carlo")
    for f,out in self.output.items():
      print(f,out)
      
      
      

