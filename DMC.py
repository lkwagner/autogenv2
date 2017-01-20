from __future__ import print_function
import average_tools as avg
####################################################
class DMCWriter:
  def __init__(self,options={}):
    self.sysfiles=['qw_000.sys']
    self.wffiles=[]
    self.basenames=['qw_000']
    self.completed=False
    self.timesteps=[0.01]
    self.nblock=20
    # For Docs:
    # nmo: (int) number of orbitals needed from orbital file.
    # orbfile: (str) location of QWalk orbital file.
    # basis: (str) basis set the orbitals are expressed in.
    # states: (list) of the orbitals read in, which states are you using?
    self.extra_observables=[]
    self.set_options(options)
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        raise ValueError("Error:",k,"not a keyword for DMCWriter")
      selfdict[k]=d[k]

    # Check completeness of average generator options.
    for avg_generator in self.extra_observables:
      avg.check_opts(avg_generator)

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
      wflines=open(self.wffiles[i],'r').read().split('\n')

      # NOTE: Not needed if default is to optimize determinants?
      ## May need to modify wave function if doing derivatives.
      #if any(['average_derivative_dm'==opts['name'] for opts in self.extra_observables]):
      #  for lidx,line in enumerate(wflines):
      #    if 'slater' in line.lower() and 'jastrow' not in line.lower():
      #      wflines.insert(lidx,'optimize_det')
      #      break

      base=self.basenames[i]
      
      for t in self.timesteps:
        fname=base+'t'+str(t)+".dmc"
        infiles.append(fname)
        outlines=[
            "method { dmc timestep %g nblock %i"%(t,self.nblock)
          ]
        for avg_opts in self.extra_observables:
          outlines+=avg.average_section(avg_opts)
        outlines+=[
            "}",
            "include %s"%sys,
            "trialfunc { ",
          ] + wflines + [
            '}'
          ]

        with open(fname,'w') as f:
          f.write('\n'.join(outlines))
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
      
      
      

