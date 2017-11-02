from __future__ import print_function
import os
import average_tools as avg
####################################################
class DMCWriter:
  ''' Object to store DMC options and write input files.

  Note: nfiles needs to be consistent across options.

  Attributes:
    qmc_type (str): Name for QMC run.
    sysfiles (list): list of nfiles sys file names for runs.
    wffiles (list): list of nfiles wave function file names for runs. 
    timesteps (list): list of float time steps to evaluate.
    nblock (int): number of blocks to run.
    tmoves (bool): Use T-Moves.
    savetrace (bool): Save trace for later analysis.
    tracefiles (list): list of nfiles trace file names if savetrace is used.
    errtol (float): Tolerance on stochastic error on total energy.
    minblocks (int): Minimum number of blocks needed for good equillibration analysis.
    extra_observables (list of dicts):  TODO Fill in these details.
      nmo: (int) number of orbitals needed from orbital file.
      orbfile: (str) location of QWalk orbital file.
      basis: (str) basis set the orbitals are expressed in.
      states: (list) of the orbitals read in, which states are you using?
  Args:
    options (dict): values to initialize attributes.
  '''
  def __init__(self,options={}):
    self.sysfiles=['qw_000.sys']
    self.wffiles=[]
    #self.basenames=['qw_000']
    self.timesteps=[0.01]
    self.nblock=20
    self.tmoves=True
    self.savetrace=False
    self.tracefiles=[]
    # For Docs:
    # nmo: (int) number of orbitals needed from orbital file.
    # orbfile: (str) location of QWalk orbital file.
    # basis: (str) basis set the orbitals are expressed in.
    # states: (list) of the orbitals read in, which states are you using?
    self.extra_observables=[]

    self.qmc_type='DMC'
    self.qmc_abr='dmc'
    self.completed=False
    self.set_options(options)
    
  #-----------------------------------------------
    
  def set_options(self, d):
    ''' Save setting of options.
    Args: 
      d (dict): attributes to update.
    '''
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
    skipkeys = ['completed','sysfiles','wffiles','basenames','tracefiles']
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
    ''' Produce nfiles input files for each of the wave functions, sys, etc.

    Args:
      infiles (list): List of input file names that this will write to.
    '''
    nfiles=len(infiles)
    assert nfiles==len(self.sysfiles), "Check sysfiles"
    assert nfiles==len(self.wffiles), "Check wffiles"

    for inp,sys,wf in zip(infiles,self.sysfiles,self.wffiles):
      
      for t in self.timesteps:
        outlines=[
            "method { dmc timestep %g nblock %i"%(t,self.nblock)
          ]
        if self.tmoves:
          outlines+=['tmoves']
        if self.savetrace:
          tracename = "%s.trace"%inp
          outlines+=['save_trace %s'%tracename]
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
class DMCReader:
  ''' Reads results from a DMC calculation. 

  Attributes:
    output (dict): results of calculation. 
    completed (bool): whether the run has converged to a final answer.
  '''
  def __init__(self,errtol=0.01,minblocks=15):
    self.output={}
    self.completed=False

    self.errtol=errtol
    self.minblocks=minblocks
    self.gosling="gosling"

  def read_outputfile(self,outfile):
    ''' Read output file results.

    Args:
      outfile (str): output to read.
    '''
    return json.loads(sub.check_output([self.gosling,"-json",outfile.replace('.o','.log')]).decode())

  def check_complete(self):
    ''' Check if a DMC run is complete.
    Returns:
      bool: If self.results are within error tolerances.
    '''
    complete={}
    for fname,results in self.output.items():
      complete[fname]=True
      if results['properties']['total_energy']['error'][0] > self.errtol:
        print("DMC incomplete: (%f) does not meet tolerance (%f)"%\
            (results['properties']['total_energy']['error'][0],self.errtol))
        complete[fname]=False
      if results['total blocks']-results['warmup blocks'] < self.minblocks:
        print("DMC incomplete: Run completed %d blocks, but requires %d."%\
            (results['total blocks']-results['warmup blocks'],self.minblocks))
        complete[fname]=False
    return complete
          
  #------------------------------------------------
  def collect(self,outfiles):
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
    ''' Print out all the items in output. '''
    print("#### Diffusion Monte Carlo")
    for f,out in self.output.items():
      print(f,out)
      
      
      

