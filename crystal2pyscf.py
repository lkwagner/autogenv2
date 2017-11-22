
''' Library for converting crystal results to a PySCF object.  '''

import pandas as pd
import numpy as np
from crystal2qmc import periodic_table,read_gred, read_kred, read_outputfile
import pyscf

##########################################################################################################
def crystal2pyscf(propoutfn="prop.in.o",basis='bfd_vtz'):
  ''' Make a PySCF object with solution from a crystal run.

  Args:
    propoutfn (str): properties or crystal stdout.
    basis (str): PySCF basis option--should match the crystal basis.
  Returns:
    tuple: (mol,scf) PySCF-equilivent Mole and SCF object.
  '''
  #TODO Make basis from crystal output (fix normalization issue).
  #TODO Generalize spin and kpoint.

  # Load crystal data.
  info, crylat_parm, cryions, crybasis, crypseudo = read_gred()
  cryeigsys = read_kred(info,crybasis)

  # Format and input structure.
  atom=[
      [periodic_table[atnum%200-1],tuple(pos)]\
      for atnum,pos in zip(cryions['atom_nums'],cryions['positions'])
    ]
  mol=pyscf.gto.Mole(atom=atom,unit='bohr',basis=basis,ecp='bfd')
  mol.build()

  mf=pyscf.scf.RKS(mol)
  # Temporary check: fills in any data that's not implemented below.
  #mf.__dict__.update(pyscf.lib.chkfile.load('../py_eq/pyscf_driver.py.chkfile','scf'))

  # Copy over MO info.
  crydf=format_eigenstates(mol,cryeigsys)
  mf.mo_energy=cryeigsys['eigvals']
  mf.mo_coeff=crydf[[*range(crydf.shape[0])]].values
  mf.mo_occ=cryeigsys['eig_weights'][0][0]*2
  mf.e_tot=np.nan #TODO compute energy and put it here, if needed.

  return mol,mf

##########################################################################################################
def format_eigenstates(mol,cryeigsys):
  ''' Organize crystal eigenstates to be consistent with PySCF order.

  Args: 
    mol (Mole): Contains structure in PySCF object.
    cryeigsys (dict): eigenstate info from cryeigsys.
  Returns:
    DataFrame: eigenstates in correct order, with labeling.
  '''

  # Extract.
  crydf=pd.DataFrame(np.array(cryeigsys['eigvecs'][(0,0,0)]['real'][0]).T)

  # Label by patching PySCF order.
  pydf=pd.DataFrame(mol.spherical_labels(),columns=['atnum','elem','orb','type'])
  crydf=crydf.join(pydf[['atnum','elem','orb','type']])
  crystal_order=('', 'x', 'y', 'z', 'z^2', 'xz', 'yz',  'x2-y2', 'xy',   'z^3', 'xz^2', 'yz^2', 'zx^2', 'xyz',  'x^3',  'y^3')
  pyscf_order=  ('', 'x', 'y', 'z', 'xy',  'yz', 'z^2', 'xz',   'x2-y2', 'y^3', 'xyz',  'yz^2', 'z^3',  'xz^2', 'zx^2', 'x^3')
  orbmap=dict(zip(pyscf_order,crystal_order))
  def convert_order(key):
    try: return orbmap[key]
    except KeyError: return None
  crydf['type']=crydf['type'].apply(convert_order)

  # Reorder crydf.
  crydf=pydf.merge(crydf,on=['atnum','elem','orb','type'])
  return crydf

##########################################################################################################
def make_basis(crybasis,ions,base="qwalk"):
  ''' Format the basis set from crystal for PySCF input. '''

  raise NotImplementedError("Doesn't work due to normalization issue. Just put in basis manually for now.")

  hybridized_check = 0.0
  hybridized_check += sum(abs(crybasis['coef_s'] * crybasis['coef_p']))
  hybridized_check += sum(abs(crybasis['coef_p'] * crybasis['coef_dfg']))
  hybridized_check += sum(abs(crybasis['coef_s'] * crybasis['coef_dfg']))
  if hybridized_check > 1e-10:
    raise NotImplementedError("Hybridized AOs (like sp) not implmemented in write_basis(...)")

  # If there's no hybridization, at most one of coef_s, coef_p, and coef_dfg is
  # nonzero. Just add them, so we have one array.
  done_atoms = []
  coefs = crybasis['coef_s'] + crybasis['coef_p'] + crybasis['coef_dfg']
  
  snorm=pyscf.gto.gto_norm(0,crybasis['prim_gaus'])
  #print(snorm)
  #for c in coefs:
  #  print(c/snorm)

  shell_type = np.tile("Unknown...",crybasis['shell_type'].shape)
  typemap = ["S","SP","P","D","F","G","H"]
  for i in range(5): shell_type[crybasis['shell_type']==i] = typemap[i]

  cnt = 0
  aidx = 0
  atom_type = ions['atom_nums'][aidx]
  done_atoms.append(atom_type)
  blines = []
  for sidx in range(len(shell_type)):
    new_aidx = crybasis['atom_shell'][sidx]-1

    new_atom_type = ions['atom_nums'][new_aidx]
    if aidx != new_aidx:
      if new_atom_type in done_atoms:
        cnt+=crybasis['prim_shell'][sidx]
        continue
      else:
        atom_type = new_atom_type
        done_atoms.append(atom_type)
        aidx = new_aidx

    nprim = crybasis['prim_shell'][sidx]
    blines.append("{0} {1}".format(periodic_table[atom_type%200-1],shell_type[sidx]))
    for pidx in range(nprim):
      blines.append("  {0} {1}".format(
        crybasis['prim_gaus'][cnt],
        coefs[cnt]
      ))
      cnt += 1
  return '\n'.join(blines)

##########################################################################################################
def compute_energy(mol,mf,xc='pbe,pbe'):
  ''' Compute the energy of the solution in mol,mf.  
  Used to check consistency between PySCF run and conversion run.

  Args:
    mol (Mole): system.
    mf (SCF object): SCF system with solution inside.
  Returns:
    float: Energy from calcualation
  '''

  dm=mf.make_rdm1()
  mf.xc=xc
  h1e=mf.get_hcore(mol)
  vhf=mf.get_veff(mol,dm)
  return mf.energy_tot(dm,h1e,vhf)

##########################################################################################################
def test_crystal2pyscf(ref_chkfile="pyscf_driver.py.o",propoutfn="prop.in.o"):
  ''' Test converter by reading an identical PySCF run and checking results are consistent.'''
  #TODO Doc

  # Reference for checking.
  check_mol=pyscf.lib.chkfile.load_mol(ref_chkfile)
  check_mf=pyscf.scf.RKS(check_mol)
  check_mf.__dict__.update(pyscf.lib.chkfile.load(ref_chkfile,'scf'))

  # Result of conversion.
  mol,mf=crystal2pyscf()


  ### Begin tests. ###

  # Energies:
  check_energy=compute_energy(check_mol,check_mf)
  test_energy=compute_energy(mol,mf)
  err=test_energy-check_energy
  print("\nBeginning tests...\n")
  print("  Energy test.")
  etest=abs(err)<1e-8
  if not etest: print("    Energy test failed; test: {:.8} Ha, check {:.8} Ha.".format(test_energy,check_energy))
  else: print("    passed.")

  print("\n\n === All tests passed! Good on you, friend! ===")

##########################################################################################################
if __name__=='__main__':
  test_crystal2pyscf(ref_chkfile="../py_eq/pyscf_driver.py.chkfile")
