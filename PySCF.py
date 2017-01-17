from __future__ import print_function
####################################################
class PySCFWriter:
  def __init__(self,options={}):
    self.basis='bfd_vtz'
    self.xyz=""
    self.ecp="bfd"
    self.spin=0
    self.charge=0
    self.pyscf_path=[]
    self.completed=False
    
    
    self.cas={'ncore':0,
              'nelec':(0,0),
              'ncas':4
              }
              

    self.set_options(options)
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for VarianceWriter")
        raise InputError
      selfdict[k]=d[k]
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
        "sys.path.append('../..')" # For pyscf2qwalk.py TODO cleaner?
      ] + add_paths + [
        "from pyscf import gto,scf,mcscf",
        "from pyscf2qwalk import print_qwalk",
        "mol=gto.Mole()",
        "mol.build(atom='''"+self.xyz+"''',",
        "basis='%s',"%self.basis,
        "ecp='%s')\n"%self.ecp,
        "mol.charge=%i\n"%self.charge,
        "mol.spin=%i\n"%self.spin,
        "m=scf.ROHF(mol)",
        "print('E(HF) =',m.kernel())",
        "print_qwalk(mol,m)"
      ]
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
    return ret
          
  #------------------------------------------------
  def collect(self,outfiles):
    for f in outfiles:
      if f not in self.output.keys():
        self.output[f]=[]
      self.output[f].append(self.read_outputfile(f))
    self.completed=True
      
  #------------------------------------------------
  def write_summary(self):
    print("#### Variance optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("dispersion",run['sigma'])
      
      
      

