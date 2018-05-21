from __future__ import print_function
from mython import fix_duped_json
import average_tools as avg
####################################################
class PostprocessWriter:
  def __init__(self,options={}):
    self.qmc_type='Postprocess'
    self.sysfiles=['qw_000.sys']
    self.wffiles=[]
    self.nskip=0
    self.tracefiles=[]
    self.qmc_abr='post'
    self.completed=False
    # For Docs:
    # nmo: (int) number of orbitals needed from orbital file.
    # orbfile: (str) location of QWalk orbital file.
    # basis: (str) basis set the orbitals are expressed in.
    # states: (list) of the orbitals read in, which states are you using?
    self.extra_observables=[
        {'name':'region_fluctuation','maxn':20}
      ]
    self.set_options(options)
    
  #-----------------------------------------------

  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        raise ValueError("Error:",k,"not a keyword for PostprocessWriter")
      selfdict[k]=d[k]

    # Check completeness of average generator options.
    for avg_generator in self.extra_observables:
      avg.check_opts(avg_generator)

  #-----------------------------------------------

  def is_consistent(self,other):
    #In principle we should check for the files, but 
    #they are often determined *after* the plan has been 
    #written so it's not currently practical to check them.
    skipkeys = ['completed','sysfiles','wffiles','tracefiles']
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
    assert nfiles==len(self.tracefiles), "Check tracefiles"

    for inp,sys,wf,trace in zip(infiles,self.sysfiles,self.wffiles,self.tracefiles):
      # TODO calculation skip.
      outlines=[
          "method { postprocess",
          "  noenergy",
          "  readconfig %s"%trace,
          "  nskip %s"%self.nskip
        ]
      for avg_opts in self.extra_observables:
        outlines+=avg.average_section(avg_opts)
      outlines+=[
          "}",
          "include %s"%sys,
          "trialfunc { include %s }"%wf
        ]

      with open(inp,'w') as f:
        f.write('\n'.join(outlines))
    self.completed=True
     
####################################################
import subprocess as sub
import json
class PostprocessReader:
  def __init__(self):
    self.output={}
    self.completed=False
    self.gosling="gosling"

  def read_outputfile(self,outfile):
    return json.load(open(outfile,'r'))
          
  #------------------------------------------------
  def collect(self,outfiles):
    for f in outfiles:
      # Hack to fix JSONS when there are multiple tbdm_basis runs.
      fix_duped_json(f.replace('.o','.json'),newfn='fixed.json')
      self.output[f]=self.read_outputfile('fixed.json')
    self.completed=True
      
  #------------------------------------------------
  def write_summary(self):
    print("#### Postprocess")
    for f,out in self.output.items():
      print(f,out)
      
