from __future__ import print_function
####################################################
class PySCFWriter:
  def __init__(self,options={}):
    self.basis='bfd_vtz'
    self.xyz=""
    self.ecp="bfd"
    self.spin=0
    self.max_cycle=50
    self.charge=0
    self.pyscf_path=[]
    self.completed=False
    # For docs later:
    # Options include: 'ROHF','UHF'
    # TODO: RKS,UKS (DFT)
    self.method='ROHF' 
    self.postHF=False   
    
    # ncore: %d     -- Number of core states.
    # ncas: %d      -- Number of states in active space. 
    # nelec: (%d,%d)-- Number of up (x) and down (y) electrons in active space.
    # tol: %f       -- tolerance on coefficient for det to be included in QWalk calculations.
    # method: %s    -- CASSCF or CASCI.
    self.cas={}
    self.dft="" #Any valid input for PySCF. This gets put into the 'xc' variable

    self.basename ='qw'


    self.special_guess=False
    self.dm_generator="""
def generate_guess(atomspins,mol,
                   double_occ={} ):
  dm_init_guess=scf.uhf.init_guess_by_minao(mol)
  print(dm_init_guess[0].diagonal())
  for atmid, (shl0,shl1,ao0,ao1) in enumerate(mol.offset_nr_by_atom()):

    opp=int((atomspins[atmid]+1)/2)
    s=(opp+1)%2
    sym=mol.atom_pure_symbol(atmid)
    print(sym,atmid,s,opp)
    docc=[]
    if sym in double_occ:
      docc=double_occ[sym]

    for ii,i in enumerate(range(ao0,ao1)):
      if ii not in docc:
        dm_init_guess[opp][i,i]=0.0
  return dm_init_guess
"""    
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
    skipkeys = ['completed']
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
        "mol=gto.Mole()",
        "mol.build(atom='''"+self.xyz+"''',",
        "basis='%s',"%self.basis,
        "ecp='%s')"%self.ecp,
        "mol.charge=%i"%self.charge,
        "mol.spin=%i"%self.spin,
        "m=%s(mol)"%self.method,
        "m.max_cycle=%d"%self.max_cycle
      ]

    if self.dft!="":
      outlines+=['m.xc="%s"'%self.dft]

    if self.special_guess:
      outlines+=[self.dm_generator]
      outlines+=["init_dm=generate_guess(%s,mol,%s)"%(str(self.atomspins),str(self.double_occ))]
      outlines+=["print('E(HF) =',m.kernel(init_dm))"]
      
    else:
      outlines+=["print('E(HF) =',m.kernel())"]
    
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
      if 'E(HF)' in line:
        ret['HF_Energy'] = float(line.split('=')[1]) 
      if 'CASCI' in line: 
        ret['CASCI_Energy'] =float(line.split()[3]) 
      if 'CASSCF' in line:
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
      
      
      

