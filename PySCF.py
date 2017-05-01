from __future__ import print_function
import os
import shutil as sh
from pymatgen.io.cif import CifParser
from pyscf import lib
from pyscf.scf.uhf import UHF,mulliken_meta


####################################################
class PySCFWriter:
  def __init__(self,options={}):
    self.basis='bfd_vtz'
    self.charge=0
    self.completed=False
    self.dft="" #Any valid input for PySCF. This gets put into the 'xc' variable
    self.diis_start_cycle=1
    self.ecp="bfd"
    self.level_shift=0.0
    self.conv_tol=1e-10
    self.max_cycle=50
    self.method='ROHF' 
    self.postHF=False   
    self.direct_scf_tol=1e-10
    self.pyscf_path=[]
    self.spin=0
    self.xyz=""
    
    # ncore: %d     -- Number of core states.
    # ncas: %d      -- Number of states in active space. 
    # nelec: (%d,%d)-- Number of up (x) and down (y) electrons in active space.
    # tol: %g       -- tolerance on coefficient for det to be included in QWalk calculations.
    # method: %s    -- CASSCF or CASCI.
    self.cas={}

    self.basename ='qw'

    # Default chosen by method at runtime.
    self.dm_generator=None

    self.set_options(options)
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__

    # Check important keys are set. 
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for PySCFWriter")
        raise AssertionError
      selfdict[k]=d[k]

    # If charge and spin should have same parity.
    assert selfdict['charge']%2==selfdict['spin']%2,"""
      Spin and charge should both be even or both be odd.
      Charge=%d, spin=%d."""%(selfdict['charge'],selfdict['spin'])

    # If postHF got set, new options are required input.
    if self.postHF==True:
      for key in ['ncore','nelec','ncas','tol','method']:
        assert key in self.cas.keys(),"%s missing from 'cas' settings! "%key+\
            "Make sure all of 'ncore','nelec','ncas','tol','method' are set."
  #-----------------------------------------------
  def is_consistent(self,other):
    # dm_generator currently gets printed differently because of the dictionaries.
    skipkeys = ['completed','chkfile','dm_generator']
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

    # Figure out correct default initial guess (if not set).
    if self.dm_generator is None:
      if self.method in ['RKS','RHF','ROHF']:
        self.dm_generator=dm_from_rhf_minao()
      elif self.method in ['UKS','UHF']:
        self.dm_generator=dm_from_uhf_minao()
      else:
        print("Warning: default guess not set for method=%s.\n Trying UHF."%self.method)
        self.dm_generator=dm_from_uhf_minao()

    for i in self.pyscf_path:
      add_paths.append("sys.path.append('"+i+"')")
    outlines=[
        "import sys",
      ] + add_paths + [
        "import pyscf",
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
        "m.direct_scf_tol=%g"%self.direct_scf_tol,
        "m.chkfile='%s'"%chkfile,
        "m.conv_tol=%g"%self.conv_tol,
        "m.diis_start_cycle=%d"%self.diis_start_cycle
      ] + self.dm_generator

    if self.level_shift>0.0:
      outlines+=["m.level_shift=%g"%self.level_shift]
    
    if self.dft!="":
      outlines+=['m.xc="%s"'%self.dft]

    outlines+=["print('E(HF) =',m.kernel(init_dm))"]
    
    if self.postHF :
      outlines += ["mc=mcscf.%s(m, ncas=%i, nelecas=(%i, %i),ncore= %i)"%( 
                   self.cas['method'], self.cas['ncas'], self.cas['nelec'][0], 
                   self.cas['nelec'][1], self.cas['ncore']), 
                   "mc.direct_scf_tol=%f"%self.direct_scf_tol,

                   "mc.kernel()",

                   "print_qwalk(mol, mc, method= 'mcscf', tol = %g , basename = '%s')"%(
                    self.cas['tol'], self.basename)]
    else:
      outlines +=[ "print_qwalk(mol,m)"]
    f.write('\n'.join(outlines))

    self.completed=True
    return [fname],[fname+".o"],[chkfile]
     

####################################################


from xml.etree.ElementTree import ElementTree

def generate_pbc_basis(xml_name,symbol,min_exp=0.2,naug=2,alpha=3,
                       cutoff=0.2,basis_name='vtz',
                       nangular={"s":1,"p":1,"d":1,"f":1,"g":0}
                       ):
  transition_metals=["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"]
  if symbol in transition_metals:
    nangular['s']=min(nangular['s'],2)
  tree = ElementTree()
  tree.parse(xml_name)
  element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
  basis_path = './Basis-set[@name="{}"]/Contraction'.format(basis_name)
  allbasis=[]

  # add in the first nangular basis functions.
  found_orbitals = []  
  for contraction in element.findall(basis_path):
    angular = contraction.get('Angular_momentum')
    if found_orbitals.count(angular) >= nangular[angular]:
      continue
    nterms = 0
    basis_sec=symbol+" " + angular +"\n"
    for basis_term in contraction.findall('./Basis-term'):
      exp = basis_term.get('Exp')
      coeff = basis_term.get('Coeff')
      if float(exp) > cutoff:
        basis_sec += '  {} {} \n'.format(exp, coeff)
        nterms+=1
    if nterms > 0:
      allbasis.append(basis_sec)
      found_orbitals.append(angular)

  angular_uncontracted=['s','p']
  if symbol in transition_metals:
    angular_uncontracted.append('d')

  for angular in angular_uncontracted:
    for i in range(0,naug):
      exp=min_exp*alpha**i
      #print(symbol,angular)
      basis_sec=symbol+ " " + angular + "\n"
      basis_sec+='{} {}\n'.format(exp,1.0)
      allbasis.append(basis_sec)
  #print(" ".join(allbasis))
  return " ".join(allbasis)
  

    
           
  
####################################################

class PySCFPBCWriter:
  def __init__(self,options={}):
    self.basis='bfd_vtz'
    self.charge=0
    self.cif=''
    self.completed=False
    self.dft="pbe,pbe" #Any valid input for PySCF. This gets put into the 'xc' variable
    self.diis_start_cycle=1
    self.ecp="bfd"
    self.level_shift=0.0
    self.conv_tol=1e-7
    self.max_cycle=50
    self.method='RKS' 
    self.direct_scf_tol=1e-7
    self.pyscf_path=[]
    self.spin=0
    self.gs=[4,4,4]
    self.xyz=""
    self.latticevec=""
    self.kpts=[2,2,2]
    self.bfd_library="BFD_Library.xml"
    self.basis_parameters={'cutoff':0.2,'basis_name':'vtz',
                          'naug':2,'alpha':3,'min_exp':0.2 } 
    self.special_basis={}
    
    self.basename ='qw'

    # Default chosen by method at runtime.
    self.dm_generator=None

    self.set_options(options)
  #-----------------------------------------------

  def from_cif(self,cifstring,primitive=True,supercell=[[1,0,0],[0,1,0],[0,0,1]]):
    struct=CifParser.from_string(cifstring).get_structures(primitive=primitive)[0]
    struct.make_supercell(supercell)
    struct=struct.as_dict()
    self.latticevec=" "
    for a in struct['lattice']['matrix']:
      for b in a:
        self.latticevec+= str(b)+ " "

    #print(struct['sites'])
    self.xyz=""
    elements=set()
    for s in struct['sites']:
      if len(s['species']) > 1:
        print("More than one species per site.. Taking the first.")
      self.xyz+=s['species'][0]['element']+" " + " ".join(map(str,s['xyz'])) + "\n"
      elements.add(s['species'][0]['element'])

    for e in elements:
      self.special_basis[e]=generate_pbc_basis(self.bfd_library,e,**self.basis_parameters)
    #print(self.xyz)
    
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__

    # Check important keys are set. 
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for PySCFWriter")
        raise AssertionError
      selfdict[k]=d[k]

    # Must be done after bdf_library is set.
    if 'cif' in d.keys():
      self.from_cif(d['cif'])
  #-----------------------------------------------

  def is_consistent(self,other):
    # dm_generator currently gets printed differently because of the dictionaries.
    skipkeys = ['completed','chkfile','dm_generator']
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

    # Figure out correct default initial guess (if not set).
    if self.dm_generator is None:
      if self.method in ['RKS','RHF','ROHF']:
        self.dm_generator=dm_from_rhf_minao()
      elif self.method in ['UKS','UHF']:
        self.dm_generator=dm_from_uhf_minao()
      else:
        print("Warning: default guess not set for method=%s.\n Trying UHF."%self.method)
        self.dm_generator=dm_from_uhf_minao()

    print(self.dm_generator)
    for i in self.pyscf_path:
      add_paths.append("sys.path.append('"+i+"')")
    outlines=[
        "import sys",
      ] + add_paths + [
        "import pyscf",
        "import numpy",
        "from pyscf.pbc import gto,scf",
        "from pyscf.pbc.scf import KRHF as RHF",
        "from pyscf.pbc.scf import KUHF as UHF",
        "from pyscf.pbc.dft import KRKS as RKS",
        "from pyscf.pbc.dft import KUKS as UKS",
        "from pyscf2qwalk import print_qwalk"]

    #The basis
    outlines+=["basis={"]
    for el in self.special_basis.keys():
      outlines+=["'"+el+"':pyscf.gto.basis.parse('''"]
      outlines+=[self.special_basis[el] + "'''),"]
    outlines+=['}']
    # The cell/molecule
    outlines+=[
        "mol=gto.M(verbose=4,",
        "gs="+str(self.gs)+",",
        "atom='''"+self.xyz+"''',",
        "a='''"+str(self.latticevec) +"''',",
        "basis=basis,",
        "spin=%i,"%self.spin,
        "ecp='%s')"%self.ecp,
        "mol.charge=%i"%self.charge
        ]
    #Set up k-points
    outlines+=['kpts=mol.make_kpts('+str(self.kpts) + ')']
    
    #Mean field
    outlines+=[
        "m=%s(mol,kpts)"%self.method,
        "m.max_cycle=%d"%self.max_cycle,
        "m.direct_scf_tol=%g"%self.direct_scf_tol,
        "m.chkfile='%s'"%chkfile,
        "m.conv_tol=%g"%self.conv_tol,
        "m.diis_start_cycle=%d"%self.diis_start_cycle
      ] 
      
    outlines+=self.dm_generator
    if self.method in ['UKS','UHF']:
      outlines+=['dm_kpts= numpy.array([[init_dm[0] for k in range(len(kpts))],' +\
                          '[init_dm[1] for k in range(len(kpts))]])'] 
    else: 
      outlines+=['dm_kpts= [init_dm for k in range(len(kpts))]']

    if self.level_shift>0.0:
      outlines+=["m.level_shift=%g"%self.level_shift]
    
    if self.dft!="":
      outlines+=['m.xc="%s"'%self.dft]

    outlines+=["print('E(HF) =',m.kernel(dm_kpts))"]
    
    outlines +=[ "print_qwalk(mol,m)"]
    f.write('\n'.join(outlines))

    self.completed=True
    return [fname],[fname+".o"],[chkfile]
    
    
####################################################
class PySCFReader:
  def __init__(self):
    self.output={}
    self.completed=False

  def read_outputfile(self,outfile):
    ''' Read energy from outputfile (obsolete).'''
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
  def read_chkfile(self,chkfile):
    ''' Read all data from the chkfile.'''
    ret={}
    mol=lib.chkfile.load_mol(chkfile)

    # TODO density matrix for mcscf parts.
    # I don't think those results are saved in the chkfile, necessarily,
    # unfortunately.
    uhf=UHF(mol)
    dm=uhf.from_chk('pyscf_driver.py.chkfile')
    ret['basis_labels']=mol.spherical_labels()
    ret['density_matrix']=dm

    for key in ('scf','mcscf'):
      ret[key]=lib.chkfile.load(chkfile,key)
    return ret
          
  #------------------------------------------------
  def collect(self,outfiles,chkfiles):
    problem=False
    for outf,chkf in zip(outfiles,chkfiles): 
      if outf not in self.output.keys():
        self.output[outf]={}
      if 'converged' not in open(outf,'r').read().split():
        problem=True
   #   self.output[outf].append(self.read_outputfile(outf))
      self.output[outf] = self.read_chkfile(chkf)
      self.output[outf]['chkfile']=chkf
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

# You should be also able to use the DM read into autogenv2, but I haven't
# thought about how precisely this should work.
def dm_from_rhf_minao():
  return ["init_dm=pyscf.scf.rhf.init_guess_by_minao(mol)"]
      
def dm_from_uhf_minao():
  return ["init_dm=pyscf.scf.uhf.init_guess_by_minao(mol)"]

def dm_set_spins(atomspins,double_occ={},startdm=None):
  ''' startdm should be the location of a chkfile if not None. '''
  if startdm is None:
    dmstarter='pyscf.scf.uhf.init_guess_by_minao(mol)'
  else:
    dmstarter=startdm

  return [
    "atomspins=%r"%atomspins,
    "double_occ=%r"%double_occ,
    "init_dm=%s"%dmstarter,
    "print(init_dm[0].diagonal())",
    "for atmid, (shl0,shl1,ao0,ao1) in enumerate(mol.offset_nr_by_atom()):",
    "  if atmid < len(atomspins) and atomspins[atmid]!=0:",
    "    opp=int((atomspins[atmid]+1)/2)",
    "    s=(opp+1)%2",
    "    sym=mol.atom_pure_symbol(atmid)",
    "    print(sym,atmid,s,opp)",
    "    docc=[]",
    "    if sym in double_occ:",
    "      docc=double_occ[sym]",
    "    for ii,i in enumerate(range(ao0,ao1)):",
    "      if ii not in docc:",
    "        init_dm[opp][i,i]=0.0",
  ]

def dm_from_chkfile(chkfile):
  """ Read a dm from a chkfile produced by a PySCF calculation. 

  It's preferrable to use absolute file paths, because the working directory
  will change to the folder where the driver is executed when this line is read."""
  return ["init_dm=m.from_chk('%s')"%chkfile]
