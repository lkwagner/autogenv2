from crystal2qmc import convert_crystal
import subprocess as sub

class Crystal2QWalkConverter:
  """ Convert a crystal run into the appropriate qwalk input files. Keep track
  of which input files correspond to what QWalk parameters.
  
  Operates at the same level as a Manager, since it manages files within the
  directory. """
  # ------------------
  def __init__(self,crysoutfn,propoutfn=None):
    self.propoutfn=propoutfn
    self.crysoutfn=crysoutfn
    self.kpoints='real' # or 'all'
    self.kpointmap={}
    self.excitation_map={}
    self.completed=False

  # ------------------
  def set_options(self, d):
    selfdict=self.__dict__
    for k in d.keys():
      if not k in selfdict.keys():
        print("Error:",k,"not a keyword for Crystal2QWalkConverter")
        raise InputError
      selfdict[k]=d[k]

  # ------------------
  def nextstep(self,tmpout='patched.o'):
    # TODO This can be replaced with contents of method once new converter is tested.
    self.convert_oldcrystal2qmc(tmpout)
    self.completed=True

  # ------------------
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

  # ------------------
  def convert_oldcrystal2qmc(self,tmpout='patched.o'):
    self.crystal_patch_output(tmpout)
    if self.kpoints=='real':
      stdout=sub.check_output('crystal2qmc -o qw %s'%tmpout,shell=True).decode()
      with open('crystal2qmc.out','w') as outf: outf.write(stdout)
    elif self.kpoints=='all':
      stdout=subprocess.check_output('crystal2qmc -c -o qw %s'%tmpout,shell=True)
      with open('crystal2qmc.out','w') as outf: outf.write(stdout)
    else:
      raise AssertionError('Crystal2QWalkConverter needs kpoints = real or all')
  
  # ------------------
  def crystal_patch_output(self,patchname):
    """ Hack to put together properties and crystal output. 
    Should be replaced by new crystal2qmc.py script soon. """
    prop=open(self.propoutfn,'r')
    shrink=[1,1,1]
    for line in prop:
      if "SHRINK FACTORS(MONK.)" in line:
        spl=line.split()
        shrink[0]=int(spl[2])
        shrink[1]=int(spl[3])
        shrink[2]=int(spl[4])

    patch=open(patchname,'w')

    out=open(self.crysoutfn,'r')
    for line in out:
      if "SHRINK. FACT.(MONKH.)" in line:
        patch.write("SHRINK. FACT.(MONKH.)    %i  %i  %i  NUMBER OF K POINTS IN THE IBZ    XXX\n"%(shrink[0],shrink[1],shrink[2]))
      else:
        patch.write(line)
      
      if "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT EDFT" in line:
        break
    out.close()

    prop=open(self.propoutfn,'r')
    patch.write("NEWK EIGENVECTORS\n \n")
    found_hamil=False
    for line in prop:
      if "HAMILTONIAN EIGENVECTORS" in line:
        found_hamil=True
      if found_hamil:
        patch.write(line)

  # ------------------
  # Let's try to make new converter the new default for autogenv2?
  def convert_crystal2qmc(self):
    raise NotImplementedError
    _,_,_,_,_,eigsys =\
      convert_crystal(
          base="qw",
          kfmt='int',
          kset=job_record['qmc']['kpoints']
        )
    # Geometric k-point weights. 
    self.kpt_weights = eigsys['kpt_weights'].tolist()

  def convert_pyscf2qmc(self):
    raise NotImplementedError('Last update: pyscf2qmc.py still under dev.')
