
import numpy as np
from crystal2qmc import periodic_table,read_gred, read_kred
import pyscf

def crystal2pyscf(propoutfn="prop.in.o"):
  ''' Make a PySCF object with solution from a crystal run.

  Args:
    propoutfn (str): properties or crystal stdout.
  Returns:
    PySCF object with run.
  '''
  #TODO what will this return?

  info, crylat_parm, cryions, crybasis, crypseudo = read_gred()
  cryeigsys = read_kred(info,crybasis)

  vecmat=[np.array(cryeigsys['eigvecs'][(0,0,0)]['real'][s]) for s in [0,1]]
  #print(basis)

  atom=[
      [periodic_table[atnum%200-1],tuple(pos)]\
      for atnum,pos in zip(cryions['atom_nums'],cryions['positions'])
    ]

  basis=make_basis(crybasis,cryions)

  mol=pyscf.gto.Mole(atom=atom,unit='bohr')
  mol.basis=basis
  mol.build()

  return mol

def make_basis(crybasis,cryions):
  ''' Make the basis string for PySCF.

  Args:
    crybasis: basis from read_gred().
    cryions: ions from read_gred().
  Returns:
    dict: basis for pyscf.gto.Mole.basis.
  '''
  # Mostly a copy from crystal2qmc.write_basis
  # PySCF API doesn't seem to have a way to feed the basis in directly as a dict, unfortunately.

  blines=[]

  coefs = crybasis['coef_s'] + crybasis['coef_p'] + crybasis['coef_dfg']

  shell_type = np.tile("Unknown...",crybasis['shell_type'].shape)
  typemap = [0,-1,1,2,3,5,6]
  for i in range(5): shell_type[crybasis['shell_type']==i] = typemap[i]

  cnt = 0
  aidx = 0
  atom_type = periodic_table[cryions['atom_nums'][aidx]%200-1]

  for sidx in range(len(shell_type)):
    new_aidx = crybasis['atom_shell'][sidx]-1

    new_atom_type = periodic_table[cryions['atom_nums'][new_aidx]%200-1]
    if aidx != new_aidx:
      if new_atom_type in basis:
        cnt+=crybasis['prim_shell'][sidx]
        continue
      else:
        atom_type = new_atom_type
        aidx = new_aidx
        basis[periodic_table[atom_type-200-1]]=[]

    nprim = crybasis['prim_shell'][sidx]
    basis[atom_type].append(shell_type[sidx])
    for pidx in range(nprim):
      basis[atom_type].append([crybasis['prim_gaus'][cnt],coefs[cnt]])
      cnt += 1
  return basis

def compute_energy(mol,mf):
  ''' Compute the energy of the solution in mol,mf.

  Used to check consistency between PySCF run and conversion run.

  Args:
    mol (Mole): system.
    mf (SCF object): SCF system with solution inside.
  Returns:
    float: Energy from calcualation
  '''

  dm=mf.make_rdm1()
  mf.energy_tot
  h1e=mf.get_hcore(mol)
  vhf=mf.get_veff(mol,dm)
  return mf.energy_tot(dm,h1e,vhf)

def test_crystal2pyscf(ref_chkfile="pyscf_driver.py.o",propoutfn="prop.in.o"):
  ''' Test converter by reading an identical PySCF run and checking results are consistent.'''
  #TODO Doc

  # Reference for checking.
  check_mol=pyscf.lib.chkfile.load_mol(ref_chkfile)
  check_mf=pyscf.scf.UKS(check_mol)
  check_mf.__dict__.update(pyscf.lib.chkfile.load(ref_chkfile,'scf'))

  check_energy=compute_energy(check_mol,check_mf)
  print(check_energy)

  mol=crystal2pyscf()
  mf=pyscf.scf.UKS(mol)
  mf.__dict__.update(pyscf.lib.chkfile.load(ref_chkfile,'scf'))

  check_energy=compute_energy(mol,mf)




if __name__=='__main__':
  test_crystal2pyscf(ref_chkfile="../h2_pyscf/pyscf_driver.py.chkfile")
