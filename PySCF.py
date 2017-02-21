from __future__ import print_function
import os
import shutil as sh

def dm_from_minao():
  return ["init_dm=scf.uhf.init_guess_by_minao(mol)"]

def dm_set_spins(spins,double_occ={}):
  return [
    "atomspins=%r"%spins
    "double_occ=%r"%double_occ
    "init_dm=scf.uhf.init_guess_by_minao(mol)",
    "print(init_dm[0].diagonal())",
    "for atmid, (shl0,shl1,ao0,ao1) in enumerate(mol.offset_nr_by_atom()):",
    "  opp=int((atomspins[atmid]+1)/2)",
    "  s=(opp+1)%2",
    "  sym=mol.atom_pure_symbol(atmid)",
    "  print(sym,atmid,s,opp)",
    "  docc=[]",
    "  if sym in double_occ:",
    "    docc=double_occ[sym]",
    "  for ii,i in enumerate(range(ao0,ao1)):",
    "    if ii not in docc:",
    "      init_dm[opp][i,i]=0.0",
  ]

def dm_from_chkfile(chkfile):
  # It might be nice to catch this error somewhere and have it just continue to
  # the next job. The chkfile may be created by other members of the job
  # ensemble.
  assert os.path.exists(chkfile), """
  chkfile doesn't exist. 

  Supposedly it's

  %s

  I'm currently in 

  %s
  """%(chkfile,os.getcwd())
  return ["init_dm=mol.from_chk(%s)"]

####################################################
class PySCFWriter:
  def __init__(self,options={}):
    self.basis='bfd_vtz'
    self.charge=0
    self.completed=False
    self.dft="" #Any valid input for PySCF. This gets put into the 'xc' variable
    self.diis=True
    self.diis_start_cycle=1
    self.ecp="bfd"
    self.level_shift=0.0
    self.max_cycle=50
    self.method='ROHF' 
    self.postHF=False   
    self.pyscf_path=[]
    self.spin=0
    self.xyz=""
    
    # ncore: %d     -- Number of core states.
    # ncas: %d      -- Number of states in active space. 
    # nelec: (%d,%d)-- Number of up (x) and down (y) electrons in active space.
    # tol: %f       -- tolerance on coefficient for det to be included in QWalk calculations.
    # method: %s    -- CASSCF or CASCI.
    self.cas={}

    self.basename ='qw'

    self.dm_generator=dm_from_minao()
    self.double_occ={}
    self.atomspins=[]

    self.set_options(options)
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__

    # Check important keys are set. 
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for VarianceWriter")
        raise InputError
      selfdict[k]=d[k]

    # If postHF got set, new options are required input.
    if self.postHF==True:
      for key in ['ncore','nelec','ncas','tol','method']:
        assert key in self.cas.keys(),"%s missing from 'cas' settings! "%key+\
            "Make sure all of 'ncore','nelec','ncas','tol','method' are set."
  #-----------------------------------------------
  def is_consistent(self,other):
    skipkeys = ['completed','chkfile']
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
  def pyscf_input(self,fname):
    f=open(fname,'w')
    chkfile=fname+".chkfile"
    add_paths=[]
    for i in self.pyscf_path:
      add_paths.append("sys.path.append('"+i+"')")
    outlines=[
        "import sys",
      ] + add_paths + [
        "from pyscf import gto,scf,mcscf",
        "from pyscf.scf import ROHF, UHF",
        "from pyscf.dft.rks import RKS",
        "from pyscf.dft.uks import UKS",
        "from pyscf2qwalk import print_qwalk",
        "mol=gto.Mole(verbose=4)",
        "mol.build(atom='''"+self.xyz+"''',",
        "basis='%s',"%self.basis,
        "ecp='%s')"%self.ecp,
        "mol.charge=%i"%self.charge,
        "mol.spin=%i"%self.spin,
        "m=%s(mol)"%self.method,
        "m.max_cycle=%d"%self.max_cycle,
        "m.chkfile='%s'"%chkfile,
        "m.diis=%r"%self.diis,
        "m.diis_start_cycle=%d"%self.diis_start_cycle
      ] + dm_generator

    if self.level_shift>0.0:
      outlines+=["m.level_shift=%f"%self.level_shift]
    
    if self.dft!="":
      outlines+=['m.xc="%s"'%self.dft]

    #if self.restart is not None:
    #  if os.path.exists(self.restart):
    #    sh.copy(self.restart,"%s/%s"%(os.getcwd(),chkfile))
    #  else:
    #    print("%s not found."%self.restart)
    #    raise AssertionError("Currently, PySCFWriter can't handle this error")
    #    self.completed=False
    #    return [],[]
    #  outlines+=["m.init_guess='chkfile'"]
    #  outlines+=["print('E(HF) =',m.kernel())"]
    #elif self.special_guess: # Should this be changed to "spins" or something more specific?
    #  outlines+=[self.dm_generator]
    #  outlines+=["print('E(HF) =',m.kernel(init_dm))"]
    #else:
    #  outlines+=["print('E(HF) =',m.kernel())"]

    outlines+=["print('E(HF) =',m.kernel(init_dm))"]
    
    if self.postHF :
      outlines += ["mc=mcscf.%s(m, ncas=%i, nelecas=(%i, %i),ncore= %i)"%( 
                   self.cas['method'], self.cas['ncas'], self.cas['nelec'][0], 
                   self.cas['nelec'][1], self.cas['ncore']), 

                   "mc.kernel()",

                   "print_qwalk(mol, mc, method= 'mcscf', tol = %f , basename = '%s')"%(
                    self.cas['tol'], self.basename)]
    else:
      outlines +=[ "print_qwalk(mol,m)"]
    f.write('\n'.join(outlines))

    self.completed=True
    return [fname],[fname+".o"]

     
####################################################
class PySCFReader:
  def __init__(self):
    self.output={}
    self.completed=False

  def read_outputfile(self,outfile):
    ret={}
    with open(outfile, 'r') as of: 
      lines = of.readlines() 
    for line in lines:
      if 'E(HF)' in line and 'print' not in line:
        ret['HF_Energy'] = float(line.split('=')[1]) 
      if 'CASCI energy' in line and 'print' not in line: 
        ret['CASCI_Energy'] =float(line.split()[3]) 
      if 'CASSCF energy' in line and 'print' not in line:
        ret['CASSCF_Energy'] =float(line.split()[3])
    return ret
          
  #------------------------------------------------
  def collect(self,outfiles):
    problem=False
    for f in outfiles: 
      if f not in self.output.keys():
        self.output[f]={}
      if 'converged' not in open(f,'r').read().split():
        problem=True
   #   self.output[f].append(self.read_outputfile(f))
      self.output[f]['energy'] = self.read_outputfile(f)
    if not problem:
      self.completed=True
    else: 
      print('Problem detected in PySCF run.')
      
  #------------------------------------------------
  def write_summary(self):
    print("#### Variance optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("dispersion",run['sigma'])
      
      
      

