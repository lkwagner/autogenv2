
''' Library for converting crystal results to a PySCF object.  '''

import pandas as pd
import numpy as np
from functools import reduce
import crystal2qmc
from crystal2qmc import periodic_table,read_gred, read_kred, read_outputfile
import pyscf
import pyscf.lo
import pyscf.pbc
import pyscf.pbc.dft
import downfold_tools as dt
from collections import Counter

##########################################################################################################
def crystal2pyscf_mol(propoutfn="prop.in.o",
    basis='bfd_vtz',
    gs=(8,8,8),
    basis_order=None):
  ''' Make a PySCF object with solution from a crystal run.

  Args:
    propoutfn (str): properties or crystal stdout.
    basis (str): PySCF basis option--should match the crystal basis.
  Returns:
    tuple: (mol,scf) PySCF-equilivent Mole and SCF object.
  '''
  #TODO Make basis from crystal output (fix normalization issue).
  #TODO Automatically detect and use correct spin: currently set for unrestricted calculations.
  nspin=2

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

  mf=pyscf.scf.UKS(mol)

  # Copy over MO info.
  crydfs=format_eigenstates_mol(mol,cryeigsys,basis_order)
  nvals=len(cryeigsys['eigvals'])//nspin
  mf.mo_energy=[cryeigsys['eigvals'][:nvals],cryeigsys['eigvals'][nvals:]]
  mf.mo_coeff=np.array([df[[*range(df.shape[0])]].values for df in crydfs])
  mf.mo_occ=np.array([(cryeigsys['eig_weights'][0][s]>1e-8).astype(float) for s in [0,1]])
  mf.e_tot=np.nan #TODO compute energy and put it here, if needed.

  return mol,mf

##########################################################################################################
def crystal2pyscf_cell(propoutfn="prop.in.o",
    basis='bfd_vtz',
    mesh=(16,16,16),
    basis_order=None):
  ''' Make a PySCF object with solution from a crystal run.

  Args:
    propoutfn (str): properties or crystal stdout.
    basis (str): PySCF basis option--should match the crystal basis.
  Returns:
    tuple: (cell,scf) PySCF-equilivent Mole and SCF object.
  '''
  #TODO Make basis from crystal output (fix normalization issue).
  #TODO Generalize spin and kpoint.

  # Load crystal data.
  info, crylat_parm, cryions, crybasis, crypseudo = read_gred()
  cryeigsys = read_kred(info,crybasis)

  totspin=read_outputfile(propoutfn)
  ntot=int(round(sum(crybasis['charges'])))
  nmo=int(round(sum(crybasis['nao_shell'])))
  nup=int(round(0.5*(ntot + totspin)))
  ndn=int(round(0.5*(ntot - totspin)))

  # Format and input structure.
  atom=[
      [periodic_table[atnum%200-1],tuple(pos)]\
      for atnum,pos in zip(cryions['atom_nums'],cryions['positions'])
    ]

  cell=pyscf.pbc.gto.Cell()
  cell.build(atom=atom,a=crylat_parm['latvecs'],unit='bohr',
      mesh=mesh,basis=basis,ecp='bfd')

  # Get kpoints that PySCF expects.
  # TODO only Gamma for now.
  kpts=cell.make_kpts((1,1,1))
  kpts=cell.get_scaled_kpts(kpts)

  mf=pyscf.pbc.dft.KUKS(cell)

  # Copy over MO info.
  crydfs=format_eigenstates_cell(cell,cryeigsys,basis_order)
  mf.mo_coeff=np.array([[df[[*range(df.shape[0])]].values] for df in crydfs])
  mf.mo_energy=cryeigsys['eigvals'].reshape(cryeigsys['eig_weights'].shape).swapaxes(0,1)
  mf.mo_occ=np.zeros((2,1,nmo))
  mf.mo_occ[0,0,0:nup]+=1
  mf.mo_occ[1,0,0:ndn]+=1
  mf.e_tot=np.nan #TODO compute energy and put it here, if needed.K

  return cell,mf

##########################################################################################################
def fix_basis_order(basis_order):
  ''' Return the order to rearrange the basis to order it by angular momentum (like PySCF).

  Example:
    >>> # for 's','p','d','s','s','p','p','d','d','f'
    >>> fix=_fix_basis_order([0,1,2,0,0,1,1,2,2,3])
    >>> print(fix)

    array([0,9,10,1,2,3,11,12,13,
           14,15,16,4,5,6,7,8,17,
           18,19,20,21,22,23,24,
           25,26,27,28,29,30,31,32,33])
  Args:
    basis_order (list): Current basis order.
  Returns:
    ndarray: index for reordering.
  '''
  # Sort by l
  ams=[]
  for i,l in enumerate(basis_order):
    ams+=[l for i in range(2*l+1)]
  amsorted=np.argsort(ams)
  
  # Sort within l by order of appearance.
  start=0
  counts=dict(Counter(basis_order))
  for l in range(0,max(basis_order)+1):
    counts[l]=(2*l+1)*counts[l]
    end=start+counts[l]
    amsorted[start:end]=sorted(amsorted[start:end])
    start+=counts[l]

  return amsorted

##########################################################################################################
def format_eigenstates_mol(mol,cryeigsys,basis_order=None):
  ''' Organize crystal eigenstates to be consistent with PySCF order.

  Args: 
    mol (Mole): Contains structure in PySCF object.
    cryeigsys (dict): eigenstate info from cryeigsys.
    basis_order (list): order that basis set is entered. 
      None means it's sorted by angular momentum, like PySCF.
      Example: for 's','p','d','s','s','p','p','d','d',
        use [0,1,2,0,0,1,1,2,2,3].
  Returns:
    DataFrame: eigenstates in correct order, with labeling.
  '''

  # Extract.
  crydfs=[pd.DataFrame(np.array(cryeigsys['eigvecs'][(0,0,0)]['real'][s]).T) for s in [0,1]]

  # PySCF basis order (our goal).
  pydf=pd.DataFrame(mol.sph_labels(fmt=False),columns=['atnum','elem','orb','type'])

  # Info about atoms.
  crydfs=[df.join(pydf[['atnum','elem']]) for df in crydfs]

  # Reorder basis.
  def _apply_fix(df):
    elem=df['elem'].values[0]
    fixed_basis_order=fix_basis_order(basis_order[elem])
    return df.reset_index(drop=True).loc[fixed_basis_order]
  if basis_order is not None:
    crydfs=[df.groupby(['atnum','elem']).apply(_apply_fix).reset_index(drop=True) for df in crydfs]

  # Label by patching PySCF order.
  crydfs=[df.join(pydf[['orb','type']]) for df in crydfs]
  crystal_order=('', 'x', 'y', 'z', 'z^2', 'xz', 'yz',  'x2-y2', 'xy',   'z^3', 'xz^2', 'yz^2', 'zx^2', 'xyz',  'x^3',  'y^3')
  pyscf_order=  ('', 'x', 'y', 'z', 'xy',  'yz', 'z^2', 'xz',   'x2-y2', 'y^3', 'xyz',  'yz^2', 'z^3',  'xz^2', 'zx^2', 'x^3')
  orbmap=dict(zip(pyscf_order,crystal_order))
  def convert_order(key):
    try: return orbmap[key]
    except KeyError: return None
  for s in [0,1]:
    crydfs[s]['type']=crydfs[s]['type'].apply(convert_order)

  # Reorder crydf.
  crydfs=[pydf.merge(df,on=['atnum','elem','orb','type']) for df in crydfs]

  ## DEBUG check a vector
  #ref_chkfile="../py_str/pyscf_driver.py.chkfile"
  #check_mol=pyscf.lib.chkfile.load_mol(ref_chkfile)
  #check_mf=pyscf.dft.UKS(check_mol)
  #check_mf.__dict__.update(pyscf.pbc.lib.chkfile.load(ref_chkfile,'scf'))
  #vector=0
  ##print(pd.DataFrame({'energy':check_mf.mo_energy[0]}))
  #crydfs[0]['check']=check_mf.mo_coeff[0][:,vector].real
  #crydfs[0]['diff']=crydfs[0][vector]-crydfs[0]['check']
  #print(crydfs[0][['atnum','elem','orb','type','check',vector,'diff']].round(4))
  return crydfs

##########################################################################################################
def format_eigenstates_cell(cell,cryeigsys,basis_order=None):
  ''' Organize crystal eigenstates to be consistent with PySCF order.

  Args: 
    cell (Cell): Contains structure in PySCF object.
    cryeigsys (dict): eigenstate info from cryeigsys.
    basis_order (list): order that basis set is entered. 
      None means it's sorted by angular momentum, like PySCF.
      Example: for 's','p','d','s','s','p','p','d','d',
        use [0,1,2,0,0,1,1,2,2,3].
  Returns:
    DataFrame: eigenstates in correct order, with labeling.
  '''

  # Extract.
  # TODO non-Gamma points.
  crydfs=[pd.DataFrame(np.array(cryeigsys['eigvecs'][(0,0,0)]['real'][s]).T) for s in [0,1]]

  # PySCF basis order (our goal).
  pydf=pd.DataFrame(cell.sph_labels(fmt=False),columns=['atnum','elem','orb','type'])

  # Info about atoms.
  crydfs=[df.join(pydf[['atnum','elem']]) for df in crydfs]

  # Reorder basis.
  def _apply_fix(df):
    elem=df['elem'].values[0]
    fixed_basis_order=fix_basis_order(basis_order[elem])
    return df.reset_index(drop=True).loc[fixed_basis_order]
  if basis_order is not None:
    crydfs=[df.groupby(['atnum','elem']).apply(_apply_fix).reset_index(drop=True) for df in crydfs]

  # Label by patching PySCF order.
  crydfs=[df.join(pydf[['orb','type']]) for df in crydfs]
  crystal_order=('', 'x', 'y', 'z', 'z^2', 'xz', 'yz',  'x2-y2', 'xy',   'z^3', 'xz^2', 'yz^2', 'zx^2', 'xyz',  'x^3',  'y^3')
  pyscf_order=  ('', 'x', 'y', 'z', 'xy',  'yz', 'z^2', 'xz',   'x2-y2', 'y^3', 'xyz',  'yz^2', 'z^3',  'xz^2', 'zx^2', 'x^3')
  orbmap=dict(zip(pyscf_order,crystal_order))
  def convert_order(key):
    try: return orbmap[key]
    except KeyError: return None
  for s in [0,1]:
    crydfs[s]['type']=crydfs[s]['type'].apply(convert_order)

  # Reorder crydf.
  crydfs=[pydf.merge(df,on=['atnum','elem','orb','type']) for df in crydfs]

  return crydfs

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
if __name__=='__main__':
  # Place calculation-specific info in this script and call test_crystal2pyscf_{mol,cell}.
  import test_convert
  test_convert.run_test()

