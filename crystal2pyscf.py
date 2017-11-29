
''' Library for converting crystal results to a PySCF object.  '''

import pandas as pd
import numpy as np
from crystal2qmc import periodic_table,read_gred, read_kred, read_outputfile
import pyscf
import pyscf.pbc
import pyscf.pbc.dft
from collections import Counter

# PySCF basis hardcoded for now.
basis={
    'Si':pyscf.gto.basis.parse('''
    Si s
      0.282668 0.464290 
      0.614115 -0.002322 
      1.334205 -0.268234 
      2.898645 0.031921 
      6.297493 -0.000106 
      13.681707 -0.000145 
      29.724387 0.000067 
     Si s
      0.6000000000000001 1.0
     Si s
      0.2 1.0
     Si p
      0.330843 0.281676 
      0.689658 0.069876 
      1.437625 -0.056306 
      2.996797 0.000744 
      6.246966 -0.000259 
      13.022097 -0.000022 
     Si p
      0.6000000000000001 1.0
     Si p
      0.2 1.0
     Si d
      0.539756 1.000000 
     Si f
      0.352999 1.000000 
    ''')
  }

##########################################################################################################
def crystal2pyscf_mol(propoutfn="prop.in.o",basis='bfd_vtz'):
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
def crystal2pyscf_cell(propoutfn="prop.in.o",
    basis='bfd_vtz',
    gs=(8,8,8),
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

  ## DEBUG available reference.
  #ref_chkfile="../pyscf/pyscf_driver.py.chkfile"
  #check_cell=pyscf.pbc.lib.chkfile.load_cell(ref_chkfile)
  #check_mf=pyscf.pbc.dft.KRKS(check_cell)
  #check_mf.__dict__.update(pyscf.lib.chkfile.load(ref_chkfile,'scf'))

  # Load crystal data.
  info, crylat_parm, cryions, crybasis, crypseudo = read_gred()
  cryeigsys = read_kred(info,crybasis)

  # Format and input structure.
  atom=[
      [periodic_table[atnum%200-1],tuple(pos)]\
      for atnum,pos in zip(cryions['atom_nums'],cryions['positions'])
    ]

  cell=pyscf.pbc.gto.Cell()
  cell.build(atom=atom,a=crylat_parm['latvecs'],unit='bohr',
      gs=gs,basis=basis,ecp='bfd')

  # Get kpoints that PySCF expects.
  kpts=cell.make_kpts((4,4,4))
  kpts=cell.get_scaled_kpts(kpts)

  #print("CRYSTAL")
  #print(cryeigsys['kpt_coords'])

  # TODO make cell yourself.
  #cell=check_cell

  mf=pyscf.pbc.dft.KRKS(cell)
  # Temporary check: fills in any data that's not implemented below.
  #mf.__dict__.update(pyscf.lib.chkfile.load('../py_eq/pyscf_driver.py.chkfile','scf'))

  # Copy over MO info.
  # TODO only Gamma for now.
  #mf=check_mf
  crydf=format_eigenstates_cell(cell,cryeigsys,basis_order)
  mf.mo_energy=cryeigsys['eigvals']
  mf.mo_coeff=[crydf[[*range(crydf.shape[0])]].values]
  mf.mo_occ=[(cryeigsys['eig_weights'][0][0]>1e-8).astype(float)*2]
  mf.e_tot=np.nan #TODO compute energy and put it here, if needed.

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
def format_eigenstates_mol(mol,cryeigsys):
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
  ## DEBUG available reference.
  #ref_chkfile="../pyscf/pyscf_driver.py.chkfile"
  #check_cell=pyscf.pbc.lib.chkfile.load_cell(ref_chkfile)
  #check_mf=pyscf.pbc.dft.KRKS(check_cell)
  #check_mf.__dict__.update(pyscf.pbc.lib.chkfile.load(ref_chkfile,'scf'))

  # Extract.
  # TODO non-Gamma points.
  crydf=pd.DataFrame(np.array(cryeigsys['eigvecs'][(0,0,0)]['real'][0]).T)

  # PySCF basis order (our goal).
  pydf=pd.DataFrame(cell.spherical_labels(),columns=['atnum','elem','orb','type'])

  # Info about atoms.
  crydf=crydf.join(pydf[['atnum','elem']])

  # Reorder basis.
  def _apply_fix(df):
    elem=df['elem'].values[0]
    fixed_basis_order=fix_basis_order(basis_order[elem])
    return df.reset_index(drop=True).loc[fixed_basis_order]
  if basis_order is not None:
    crydf=crydf.groupby(['atnum','elem']).apply(_apply_fix).reset_index(drop=True)

  # Label by patching PySCF order.
  crydf=crydf.join(pydf[['orb','type']])
  crystal_order=('', 'x', 'y', 'z', 'z^2', 'xz', 'yz',  'x2-y2', 'xy',   'z^3', 'xz^2', 'yz^2', 'zx^2', 'xyz',  'x^3',  'y^3')
  pyscf_order=  ('', 'x', 'y', 'z', 'xy',  'yz', 'z^2', 'xz',   'x2-y2', 'y^3', 'xyz',  'yz^2', 'z^3',  'xz^2', 'zx^2', 'x^3')
  orbmap=dict(zip(pyscf_order,crystal_order))
  def convert_order(key):
    try: return orbmap[key]
    except KeyError: return None
  crydf['type']=crydf['type'].apply(convert_order)

  # Reorder crydf.
  crydf=pydf.merge(crydf,on=['atnum','elem','orb','type'])

  ## DEBUG check a vector
  #vector=0
  ##print(pd.DataFrame({'energy':check_mf.mo_energy[0]}))
  #crydf['check']=check_mf.mo_coeff[0][:,vector].real
  #crydf['diff']=crydf[vector]-crydf['check']
  ##print(crydf[['atnum','elem','orb','type','check',vector,'diff']].round(4))
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
def compute_energy_mol(mol,mf,xc='pbe,pbe'):
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
def compute_energy_cell(cell,mf,xc='pbe,pbe'):
  ''' Compute the energy of the solution in mol,mf.  
  Used to check consistency between PySCF run and conversion run.

  Args:
    cell (Cell): system.
    mf (SCF object): SCF system with solution inside.
  Returns:
    float: Energy from calcualation
  '''

  dm=mf.make_rdm1()
  mf.xc=xc
  mf.direct_scf_tol=1e-5
  h1e=mf.get_hcore(cell)
  vhf=mf.get_veff(cell,dm[:1,:,:])
  return mf.energy_tot(dm[:1,:,:],h1e,vhf)

##########################################################################################################
def test_crystal2pyscf_mol(ref_chkfile="pyscf_driver.py.o",propoutfn="prop.in.o"):
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
  check_energy=compute_energy_mol(check_mol,check_mf)
  test_energy=compute_energy_mol(mol,mf)
  err=test_energy-check_energy
  print("\nBeginning tests...\n")
  print("  Energy test.")
  etest=abs(err)<1e-8
  if not etest: print("    Energy test failed; test: {:.8} Ha, check {:.8} Ha.".format(test_energy,check_energy))
  else: print("    passed.")

  print("\n\n === All tests passed! Good on you, friend! ===")

##########################################################################################################
def test_crystal2pyscf_cell(ref_chkfile="pyscf_driver.py.o",propoutfn="prop.in.o",basis_order=None):
  ''' Test converter by reading an identical PySCF run and checking results are consistent.'''
  #TODO Doc

  # Reference for checking.
  check_cell=pyscf.pbc.lib.chkfile.load_cell(ref_chkfile)
  check_mf=pyscf.pbc.dft.KRKS(check_cell)
  check_mf.__dict__.update(pyscf.lib.chkfile.load(ref_chkfile,'scf'))

  # Result of conversion.
  cell,mf=crystal2pyscf_cell(basis=basis,gs=(8,8,8),basis_order=basis_order)

  #print("DIFF")
  #for key in cell.__dict__.keys():
  #  try:
  #    if cell.__dict__[key]!=check_cell.__dict__[key]:
  #      print("Key:  {}".format(key))
  #      print("  Ref:  {}".format(check_cell.__dict__[key]))
  #      print("  Test: {}".format(cell.__dict__[key]))
  #  except ValueError:
  #    try:
  #      if all(cell.__dict__[key]!=check_cell.__dict__[key]):
  #        print("Key:  {}".format(key))
  #        print("  Ref:  {}".format(check_cell.__dict__[key]))
  #        print("  Test: {}".format(cell.__dict__[key]))
  #    except Exception as e: 
  #      print("Key:  {}".format(key))
  #      print("  Cannot compare.", type(e), e)
  #      print("  Ref:  {}".format(check_cell.__dict__[key]))
  #      print("  Test: {}".format(cell.__dict__[key]))

  ### Begin tests. ###

  #assert 0

  # Energies:
  print("\nBeginning tests...\n")
  print("  Energy test:")
  print("    Test energy...")
  test_energy=compute_energy_cell(cell,mf)
  print("    {} Ha".format(test_energy))
  print("    Reference energy...")
  check_energy=compute_energy_cell(check_cell,check_mf)
  print("    {} Ha".format(check_energy))
  err=test_energy-check_energy
  etest=abs(err)<1e-5
  if not etest: 
    print("    Energy test failed; test: {:.8} Ha, check {:.8} Ha.".format(test_energy,check_energy))
  else: 
    print("    passed.")

  if not all([etest]):
    print("\n\n !!! Some tests failed! Shame on you, friend! !!!")
  else:
    print("\n\n === All tests passed! Good on you, friend! ===")

##########################################################################################################
def test_counter():
  print("These should match:")
  print(np.array([0,9,10,1,2,3,11,12,13,
           14,15,16,4,5,6,7,8,17,
           18,19,20,21,22,23,24,
           25,26,27,28,29,30,31,32]))
  print(fix_basis_order([0,1,2,0,0,1,1,2,2,3]))

##########################################################################################################
if __name__=='__main__':
  pass
  # DEBUG
  #test_crystal2pyscf_cell(ref_chkfile="../pyscf/pyscf_driver.py.chkfile",basis_order={'Si':[0,1,2,3,0,0,1,1]})
  #_test_counter()
