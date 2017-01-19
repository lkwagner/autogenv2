from __future__ import print_function
####################################################
class DMCWriter:
  def __init__(self,options={}):
    self.sysfiles=['qw_000.sys']
    self.wffiles=[]
    self.basenames=['qw_000']
    self.completed=False
    self.timesteps=[0.01]
    self.nblock=20
    self.set_options(options)
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for EnergyWriter")
        raise InputError
      selfdict[k]=d[k]

    # Check completeness of average generator options.
    for avg_generator in self.extra_observables:
      check=['nmo','orbfile','include','orbfile','basis','states']
      for key in check:
        assert key in avg_generator.keys(),"""
          %s missing from 'extra_observables' settings!
          Make sure all of %s are set."""%', '.join(key,check)

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
  def _average_section(self,opts):
    outlines=[]
    if opts['name'].lower()=='average_derivative_dm':
      outlines=[
          "  average { average_derivative_dm",
          "    average { tbdm_basis",
          "      orbitals { ",
          "        cutoff_mo",
          "        magnify 1",
          "        nmo %d"%opt['nmo'],
          "        orbfile %s"%opt['orbfile'],
          "        include %s"%opt['basis'],
          "        centers { useglobal }",
          "      }",
          "      states { %s }"%' '.join(opt['states']),
          "    }",
          "  }"
        ]
    else:
      raise NotImplementedError("""
      '%s' is not implemented in autogen yet: 
      You should implement it, it should be easy!"""%opts['name'])
    return outlines

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
        outlines=[
            "method { dmc timestep %g nblock %i"%(t,self.nblock)
          ]
        for avg_opts in self.extra_observables:
          outlines+=self._average_section(avg_opts)
        outlines+=[
            "}",
            "include %s"%sys,
            "trialfunc { include %s }"%wf
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
      
      
      

