''' 
crystal2qmc

This module can read Crystal output and organize the data into python dictionaries.
It then takes those dictionaries and writes out the appropriate QWalk files.

Example:

  >>> qw_writer=QWalkWriter()
  >>> qw_writer.read_gredkred("GRED.DAT","KRED.DAT")
  >>> qw_writer.write_files()

You can also just use `convert_crystal` with defaults and it should work
for most cases.

By default this will write the ground state determinant from a Crystal run and
output several of the QWalk input files needed to get VMC or DMC started. 
There are also routines here to make excitations or multideterminant inputs out of the Crystal run.
'''

from __future__ import division,print_function
import numpy as np
import sys

# List of all atoms. Useful for translating elements from strings to integers and back. 
# Yes, I know it's neither periodic nor a table.
periodic_table = [
  "h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p","s","cl","ar",
  "k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br",
  "kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te",
  "i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm",
  "yb","lu","hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn",
  "fr","ra","ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr",
  "rf","db","sg","bh","hs","mt","ds","rg","cp","uut","uuq","uup","uuh","uus","uuo"
]

###############################################################################
# Reads in the geometry, basis, and pseudopotential from GRED.DAT.
def read_gred(gred="GRED.DAT"):
  ''' Read GRED.DAT file to read in Hamiltonian specs.

  Args:
    gred (str): Name of GRED.DAT file produced by properties using CRYAPI_OUT
  Returns:
    info (dict): Misc. information used by Crystal to organize the data.
    lat_parm (dict): Lattice parameters.
    ions (dict): Ion charges and positions, etc.
    basis (dict): Basis set information.
    pseudo (dict): Pseudopotential information.
  '''
    
    
  lat_parm = {}
  ions = {}
  basis = {}
  pseudo = {}

  gred = open(gred,'r').read()

  # Fix numbers with no space between them.
  gred = gred.replace("-"," -")
  gred = gred.replace("E -","E-") 

  gred_words = gred.split()
  nparms = [int(w) for w in gred_words[1:4]]
  cursor = 4

  # These follow naming of cryapi_inp (but "inf" -> "info").
  info = [int(w) for w in gred_words[cursor          :cursor+nparms[0]]]
  itol = [int(w) for w in gred_words[cursor+nparms[0]:cursor+nparms[1]]]
  par  = [int(w) for w in gred_words[cursor+nparms[1]:cursor+nparms[2]]]
  cursor += sum(nparms)

  lat_parm['struct_dim'] = int(info[9])

  # Lattice parameters.
  lat_parm['latvecs'] = \
      np.array(gred_words[cursor:cursor+9],dtype=float).reshape(3,3).T.round(15)
  if (lat_parm['latvecs'] > 100).any():
    print("Lattice parameter larger than 100 A! Reducing to 100.")
    print("If this is a dimension < 3 system, there is no cause for alarm.")
    print("Otherwise if this is a problem for you, please generalize crystal2qmc.")
    lat_parm['latvecs'][lat_parm['latvecs']>100] = 100.
  cursor += 9
  prim_trans= np.array(gred_words[cursor:cursor+9],dtype=float).reshape(3,3)
  cursor += 9
  lat_parm['conv_cell'] = prim_trans.dot(lat_parm['latvecs'])
  cursor += info[1] + 48*48 + 9*info[1] + 3*info[1] # Skip symmetry part.

  # Lattice "stars" (?) skipped.
  cursor += info[4]+1 + info[78]*3 + info[4]+1 + info[4]+1 + info[78] + info[78]*3

  # Some of ion information.
  natoms = info[23]
  ions['charges'] = [float(w) for w in gred_words[cursor:cursor+natoms]]
  cursor += natoms
  # Atom positions.
  atom_poss = np.array(gred_words[cursor:cursor+3*natoms],dtype=float)
  ions['positions'] = atom_poss.reshape(natoms,3)
  cursor += 3*natoms

  # Basis information (some ion information mixed in).
  nshells = info[19]
  nprim   = info[74]
  # Formal charge of shell.
  basis['charges'] = np.array(gred_words[cursor:cursor+nshells],dtype=float)
  cursor += nshells
  # "Adjoined gaussian" of shells.
  basis['adj_gaus'] = np.array(gred_words[cursor:cursor+nshells],dtype=float)
  cursor += nshells
  # Position of shell.
  shell_poss = np.array(gred_words[cursor:cursor+3*nshells],dtype=float)
  basis['positions'] = shell_poss.reshape(nshells,3)
  cursor += 3*nshells
  # Primitive gaussian exponents.
  basis['prim_gaus'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  # Coefficients of s, p, d, and (?).
  basis['coef_s'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  basis['coef_p'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  basis['coef_dfg'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  basis['coef_max'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  # Skip "old normalization"
  cursor += 2*nprim
  # Atomic numbers.
  ions['atom_nums'] = np.array(gred_words[cursor:cursor+natoms],dtype=int)
  cursor += natoms
  # First shell of each atom (skip extra number after).
  basis['first_shell'] = np.array(gred_words[cursor:cursor+natoms],dtype=int)
  cursor += natoms + 1
  # First primitive of each shell (skips an extra number after).
  basis['first_prim'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells + 1
  # Number of prims per shell.
  basis['prim_shell'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells
  # Type of shell: 0=s,1=sp,2=p,3=d,4=f.
  basis['shell_type'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells
  # Number of atomic orbtials per shell.
  basis['nao_shell'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells
  # First atomic orbtial per shell (skip extra number after).
  basis['first_ao'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells + 1
  # Atom to which each shell belongs.
  basis['atom_shell'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells

  # Pseudopotential information.
  # Pseudopotential for each element.
  pseudo_atom = np.array(gred_words[cursor:cursor+natoms],dtype=int)
  cursor += natoms
  cursor += 1 # skip INFPOT
  ngauss = int(gred_words[cursor])
  cursor += 1
  headlen = int(gred_words[cursor])
  cursor += 1
  # Number of pseudopotentials.
  numpseudo = int(gred_words[cursor])
  cursor += 1
  # Exponents of r^l prefactor.
  r_exps = -1*np.array(gred_words[cursor:cursor+ngauss],dtype=int)
  cursor += ngauss
  # Number of Gaussians for angular momenutum j
  n_per_j = np.array(gred_words[cursor:cursor+headlen],dtype=int)
  cursor += headlen
  # index of first n_per_j for each pseudo.
  pseudo_start = np.array(gred_words[cursor:cursor+numpseudo],dtype=int)
  cursor += numpseudo + 1
  # Actual floats of pseudopotential.
  exponents = np.array(gred_words[cursor:cursor+ngauss],dtype=float)
  cursor += ngauss
  prefactors = np.array(gred_words[cursor:cursor+ngauss],dtype=float)
  cursor += ngauss
  # Store information nicely.
  npjlen = int(headlen / len(pseudo_start))
  for aidx,atom in enumerate(ions['atom_nums']):
    psidx = pseudo_atom[aidx]-1
    start = pseudo_start[psidx]
    if psidx+1 >= len(pseudo_start): end = ngauss
    else                           : end = pseudo_start[psidx+1]
    if atom not in pseudo.keys():
      pseudo[atom] = {}
      pseudo[atom]['prefactors'] = prefactors[start:end]
      pseudo[atom]['r_exps'] = r_exps[start:end]
      pseudo[atom]['n_per_j'] = n_per_j[npjlen*psidx:npjlen*(psidx+1)]
      pseudo[atom]['exponents'] = exponents[start:end]

  ## Density matrix information.
  ## This is impossible to figure out.  See `cryapi_inp.f`.
  ## I think the order is overlap, Fock matrix, then DM in real space. 
  ## I'm not yet sure where the real space grid for the DM is defined.
  #atomic_charges = np.array(gred_words[cursor:cursor+natoms],dtype=float)
  #cursor += natoms
  ## A hacky way to get past all the symmetry stuff:
  ## Keep reading until you hit the floats.
  #nspin=info[63]+1
  #npgt=info[11]
  #nfgt=info[10]
  #print(npgt,nfgt)
  #for c in range(cursor,len(gred_words)):
  #  if 'E' in gred_words[c]:
  #    cursor=c
  #    break
  #print(gred_words[cursor-1:cursor+2])


  ## Skip symmetry information.
  #cursor += mvlaf*4 + info[19]*info[1] + 
  #print("atomic_charges",atomic_charges)

  return info, lat_parm, ions, basis, pseudo

###############################################################################
# Reads in kpoints and eigen{values,vectors} from KRED.DAT.
def read_kred(info,basis,kred="KRED.DAT"):
  ''' Read KRED.DAT file to read in the eigenstates and eigenvalues of a Crystal run.

  Args:
    info (dict): info dictionary output from read_gred()
    basis (dict): basis dictionary output from read_gred()
  Returns: 
    eigsys (dict): eigenvalues and vectors. 
  '''

  eigsys = {}

  kred = open(kred)
#  print(kred.readline())
#  kred=kred.read()
  kred_words = [] #kred.split()
  for lin in kred:
    kred_words += lin.split()
  cursor = 0

  # Number of k-points in each direction.
  eigsys['nkpts_dir'] = np.array([int(w) for w in kred_words[cursor:cursor+3]])
  cursor += 3
  # Total number of inequivilent k-points.
  nikpts = int(kred_words[cursor])
  cursor += 1
  # Reciprocal basis.
  recip_vecs = np.array(kred_words[cursor:cursor+9],dtype=float)
  eigsys['recip_vecs'] = recip_vecs.reshape(3,3)
  cursor += 9
  # Inequivilent k-point coord in reciprocal basis.
  ikpt_coords = np.array(kred_words[cursor:cursor+3*nikpts],int)
  ikpt_coords = list(map(tuple,ikpt_coords.reshape(nikpts,3)))
  # Useful to compare to old output format.
  eigsys['kpt_index'] = dict(zip(ikpt_coords,range(len(ikpt_coords))))
  cursor += 3*nikpts
  # is complex (0) or not (1), converted to True (if complex) or False
  ikpt_iscmpx = \
    np.array([int(w) for w in kred_words[cursor:cursor+nikpts]]) == 0
  eigsys['ikpt_iscmpx'] = dict(zip(ikpt_coords,ikpt_iscmpx))
  cursor += nikpts
  # Skip symmetry information.
  cursor += 9*48
  # Geometric weight of kpoints.
  eigsys['kpt_weights'] = np.array(kred_words[cursor:cursor+nikpts],dtype=float)
  cursor += nikpts
  # Eigenvalues: (how many) = (spin) * (number of basis) * (number of kpoints)
  eigsys['nspin'] = info[63]+1
  nevals = eigsys['nspin']*info[6]*nikpts
  eigsys['eigvals'] = np.array(kred_words[cursor:cursor+nevals],dtype=float)
  cursor += nevals
  # Weights of eigenvales--incorperating Fermi energy cutoff.
  nbands = int(round(nevals / nikpts / eigsys['nspin']))
  eigsys['eig_weights'] = np.array(kred_words[cursor:cursor+nevals],dtype=float)\
      .reshape(nikpts,eigsys['nspin'],nbands)
  cursor += nevals

  # Read in eigenvectors at inequivilent kpoints. Can't do all kpoints because we 
  # don't know if non-inequivilent kpoints are real or complex (without symmetry
  # info)
  nbands = int(round(nevals / nikpts / eigsys['nspin']))
  nkpts  = np.prod(eigsys['nkpts_dir'])
  nao = sum(basis['nao_shell'])
  ncpnts = int(nbands * nao)
  kpt_coords   = []
  # Format: eigvecs[kpoint][<real/imag>][<spin up/spin down>]
  eigvecs = {}
  for kpt in range(nkpts*eigsys['nspin']):
    try:
      new_kpt_coord = tuple([int(w) for w in kred_words[cursor:cursor+3]])
    except IndexError: # End of file.
      raise AssertionError("ERROR: KRED.DAT seems to have ended prematurely.\n" + \
            "Didn't find all {0} kpoints.".format(nikpts))
    cursor += 3

    # If new_kpt_coord is an inequivilent point...
    if new_kpt_coord in ikpt_coords:
      # If complex...
      if eigsys['ikpt_iscmpx'][new_kpt_coord]:
        eig_k = np.array(kred_words[cursor:cursor+2*ncpnts],dtype=float)
        cursor += 2*ncpnts
        eig_k = eig_k.reshape(ncpnts,2)
        kpt_coords.append(new_kpt_coord)
        if new_kpt_coord in eigvecs.keys():
          eigvecs[new_kpt_coord]['real'].append(
              eig_k[:,0].reshape(int(round(ncpnts/nao)),nao)
            )
          eigvecs[new_kpt_coord]['imag'].append(
              eig_k[:,1].reshape(int(round(ncpnts/nao)),nao)
            )
        else:
          eigvecs[new_kpt_coord] = {}
          eigvecs[new_kpt_coord]['real'] = \
            [eig_k[:,0].reshape(int(round(ncpnts/nao)),nao)]
          eigvecs[new_kpt_coord]['imag'] = \
            [eig_k[:,1].reshape(int(round(ncpnts/nao)),nao)]
      else: # ...else real.
        eig_k = np.array(kred_words[cursor:cursor+ncpnts],dtype=float)
        cursor += ncpnts
        kpt_coords.append(new_kpt_coord)
        if new_kpt_coord in eigvecs.keys():
          eigvecs[new_kpt_coord]['real'].append(
              eig_k.reshape(int(round(ncpnts/nao)),nao)
            )
          eigvecs[new_kpt_coord]['imag'].append(
              np.zeros((int(round(ncpnts/nao)),nao))
            ) # Not efficient, but safe.
        else:
          eigvecs[new_kpt_coord] = {}
          eigvecs[new_kpt_coord]['real'] = \
            [eig_k.reshape(int(round(ncpnts/nao)),nao)]
          eigvecs[new_kpt_coord]['imag'] = \
            [np.zeros((int(round(ncpnts/nao)),nao))]
    else: # ...else, skip.
      skip = True
      while skip:
        try: # If there's an int, we're at next kpoint.
          int(kred_words[cursor])
          skip = False
        except ValueError: # Keep skipping.
          cursor += ncpnts
        except IndexError: # End of file.
          skip = False
          break

  # It's probably true that kpt_coords == ikpt_coords, with repitition for spin
  # up and spin down, because we only read in inequivilent kpoints. However,
  # ordering might be different, and the ordering is correct for kpt_coords.
  # If there are bugs, this might be a source.
  eigsys['kpt_coords'] = ikpt_coords # kpt_coords
  eigsys['eigvecs'] = eigvecs

  return eigsys

###############################################################################
# Reads total spin from output file. 
# TODO Is there a way around this? Yes.
# Alternatively, this can read the CRYSTAL output file and still works!
def read_outputfile(fname = "prop.in.o"):
  fin = open(fname,'r')
  for line in fin:
    if "SUMMED SPIN DENSITY" in line:
      spin = float(line.split()[-1])
  if abs(round(spin) - spin) > 1e-8:
    print("Warning: spin %f is not close to integer!"%spin)
    print("  I'm rounding this to %d."%int(round(spin)))
  spin = int(round(spin))
  return spin

###############################################################################
def convert_crystal(base='qwalk',propoutfn='crys.in.o',kfmt='coord',kset='complex'):
  ''' You can just call this function to create the object and produce the files.
  This means you won't get an object to use for other things, but maybe you don't want it.
  It does require GRED.DAT, KRED.DAT, and crys.in.o to in the directory.

  Args:
    base (str): files are written with names [base].[filetype] and [base]_[kpt].[filetype].
    kfmt (str): 'coord' or 'int'. 'coord' would name the gamma point 000, int names it 0.
    kset (str): 'complex' or 'real'. Do you want only real or all the complex kpoints?
  Returns:
    dict: contains lists of each file produced organized by type.
  '''

  qw_writer=QWalkWriter()
  qw_writer.read_gredkred()
  outfiles=qw_writer.write_files(base,kfmt,kset)

  return outfiles

###############################################################################
class QWalkWriter:
  # This format should allow this writer to convert any program's results into QWalk input files.
  def __init__(self,nvirtual=50,gred=None,kred=None,xml=None):
    ''' This object can handle generation of input files for QWalk.

    Args:
      gred (str): GRED.DAT path, if you want to populate the object with that.
      kred (str): KRED.DAT path, if you want to populate the object with that.
      xml (str): Crystal CRYAPI_OUT XML if you want to populate the object with that.

    '''
    self.info, self.lat_parm, self.ions, self.basis, self.pseudo = {},{},{},{},{}
    self.eigsys = {}

    if xml is not None:
      raise NotImplementedError("Crystal17 XML reader not implemented.")
    elif gred is not None and kred is not None:
      self.read_gredkred(gred,kred,nvirtual=nvirtual)

  ###
  def read_gredkred(self,gred='GRED.DAT',kred='KRED.DAT',spinfn='crys.in.o',nvirtual=50):
    '''
    Populate the QWalkWriter with data from `gred` and `kred`, which are paths to the output of CRYAPI_OUT.

    Args:
      gred (str): GRED.DAT path, if you want to populate the object with that.
      kred (str): KRED.DAT path, if you want to populate the object with that.
      spinfn (str): Crystal output file (either crystal or properties) from which to find the total spin.
    '''
    # Currently also needs `spinfn` to find the total spin, but this shouldn't be required if you can find
    # the density matrix somewhere else in the output.

    self.info, self.lat_parm, self.ions, self.basis, self.pseudo = read_gred()
    self.eigsys = read_kred(self.info,self.basis)

    if self.eigsys['nspin'] > 1:
      self.eigsys['totspin'] = read_outputfile(spinfn)
    else:
      self.eigsys['totspin'] = 0

    # Useful quantities based on the read data.
    self.basis['ntot'] = int(round(sum(self.basis['charges'])))
    self.basis['nmo']  = sum(self.basis['nao_shell'])
    self.eigsys['nup'] = int(round(0.5 * (self.basis['ntot'] + self.eigsys['totspin'])))
    self.eigsys['ndn'] = int(round(0.5 * (self.basis['ntot'] - self.eigsys['totspin'])))
    self.maxmo_spin=min(max(self.eigsys['nup'],self.eigsys['ndn'])+nvirtual,self.basis['nmo'])

    # Change Crystal normalization to QWalk normalization.
    self.normalize_eigvec()

  ###
  def find_basis_cutoff(self):
    if self.lat_parm['struct_dim'] > 0:
      latvecs = self.lat_parm['latvecs']
      cutoff_divider = 2.000001
      cross01 = np.cross(latvecs[0], latvecs[1])
      cross12 = np.cross(latvecs[1], latvecs[2])
      cross02 = np.cross(latvecs[0], latvecs[2])

      heights = [0,0,0]
      heights[0]=abs(np.dot(latvecs[0], cross12)/np.dot(cross12,cross12)**.5)
      heights[1]=abs(np.dot(latvecs[1], cross02)/np.dot(cross02,cross02)**.5)
      heights[2]=abs(np.dot(latvecs[2], cross01)/np.dot(cross01,cross01)**.5)
      return min(heights)/cutoff_divider
    else:
      return 7.5

  ###
  def generate_slater(self,kpt,orbfile,basisfile,swaporbs=None):
    ''' Generate the string for a slater file.
    Args: 
      kpt (int): kpoint index.
      swaporbs (dict): Map orb numbers to each other to do excitations. 
        Example: {10:11,102:105} will excite 10 to 11 and 102 to 105.
    Returns: 
      str: contents of slater file.
    '''

    # Figure out orbital numbers.
    nmo  = self.basis['nmo']
    nup  = self.eigsys['nup']
    ndn  = self.eigsys['ndn']
    if self.maxmo_spin < 0:
      self.maxmo_spin=nmo
    uporbs = np.arange(nup)+1
    dnorbs = np.arange(ndn)+1
    if self.eigsys['nspin'] > 1:
      dnorbs += self.maxmo_spin

    # Modify for excitations.
    if swaporbs is not None:
      for idx,up in enumerate(uporbs):
        if up in swaporbs.keys(): uporbs[idx]=swaporbs[up]
      for idx,dn in enumerate(dnorbs):
        if dn in swaporbs.keys(): dnorbs[idx]=swaporbs[dn]

    if self.eigsys['ikpt_iscmpx'][kpt]: orbstr = "corbitals"
    else:                               orbstr = "orbitals"

    uporblines = ["{:5d}".format(orb) for orb in uporbs]
    width = 10
    for i in reversed(range(width,len(uporblines),width)):
      uporblines.insert(i,"\n ")
    dnorblines = ["{:5d}".format(orb) for orb in dnorbs]
    for i in reversed(range(width,len(dnorblines),width)):
      dnorblines.insert(i,"\n ")

    outlines = [
        "slater",
        "{0} {{".format(orbstr),
        "blas2_mo",
        "  magnify 1",
        "  nmo {0}".format(dnorbs.max()),
        "  orbfile {0}".format(orbfile),
        "  include {0}".format(basisfile),
        "  centers { useglobal }",
        "}",
        "detwt { 1.0 }",
        "states {",
        "  # Spin up orbitals.", 
        "  " + " ".join(uporblines),
        "  # Spin down orbitals.",
        "  " + " ".join(dnorblines),
        "}"
      ]
    return '\n'.join(outlines)

  ###
  def normalize_eigvec(self):
    ''' Normalize the eigenstate at `kpt`. '''
    snorm = 1./(4.*np.pi)**0.5
    pnorm = snorm*(3.)**.5
    dnorms = [
        .5*(5./(4*np.pi))**.5,
        (15./(4*np.pi))**.5,
        (15./(4*np.pi))**.5,
        .5*(15./(4.*np.pi))**.5,
        (15./(4*np.pi))**.5
      ]
    # f orbital normalizations are from 
    # <http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html>
    fnorms = [
        ( 7./(16.*np.pi))**.5,
        (21./(32.*np.pi))**.5,
        (21./(32.*np.pi))**.5,
        (105./(16.*np.pi))**.5, # xyz
        (105./(4.*np.pi))**.5,
        (35./(32.*np.pi))**.5,
        (35./(32.*np.pi))**.5
      ]

    # Duplicate coefficients for complex, and if multiple basis elements are d.
    # This is to align properly with the d-components of eigvecs.
    tmp = [[f for f in dnorms] for i in range(sum(self.basis['shell_type']==3))]
    dnorms = []
    for l in tmp: dnorms += l
    dnorms = np.array(dnorms)
    # Likewise for f.
    tmp = [[f for f in fnorms] for i in range(sum(self.basis['shell_type']==4))]
    fnorms = []
    for l in tmp: fnorms += l
    fnorms = np.array(fnorms)

    ao_type = []
    for sidx in range(len(self.basis['shell_type'])):
      ao_type += \
        [self.basis['shell_type'][sidx] for ao in range(self.basis['nao_shell'][sidx])]
    ao_type = np.array(ao_type)

    if any(ao_type==1):
      raise NotImplementedError("sp orbtials not implemented in normalize_eigvec(...)")

    for part in ['real','imag']:
      for spin in range(self.eigsys['nspin']):
        for kpt in self.eigsys['kpt_coords']:
          self.eigsys['eigvecs'][kpt][part][spin][:,ao_type==0] *= snorm
          self.eigsys['eigvecs'][kpt][part][spin][:,ao_type==2] *= pnorm
          self.eigsys['eigvecs'][kpt][part][spin][:,ao_type==3] *= dnorms
          self.eigsys['eigvecs'][kpt][part][spin][:,ao_type==4] *= fnorms
      
  ###
  def write_orb(self,kpt,orbfilename):
    ''' Write the orb file. The orbitals are quite large and not worth storing in RAM.
    Args: 
      kpt (int): kpoint index.
      orbfile (str): output file for orb.
    '''
    outf=open(orbfilename,'w')

    if self.maxmo_spin < 0:
      self.maxmo_spin=self.basis['nmo']
    eigvecs_real = self.eigsys['eigvecs'][kpt]['real']
    eigvecs_imag = self.eigsys['eigvecs'][kpt]['imag']
    atidxs = np.unique(self.basis['atom_shell'])-1
    nao_atom = np.zeros(atidxs.size,dtype=int)
    for shidx in range(len(self.basis['nao_shell'])):
      nao_atom[self.basis['atom_shell'][shidx]-1] += self.basis['nao_shell'][shidx]
    #nao_atom = int(round(sum(self.basis['nao_shell']) / len(ions['positions'])))
    coef_cnt = 1
    totnmo = self.maxmo_spin*self.eigsys['nspin'] #self.basis['nmo'] * self.eigsys['nspin']
    for moidx in np.arange(totnmo)+1:
      for atidx in atidxs+1:
        for aoidx in np.arange(nao_atom[atidx-1])+1:
          outf.write(" {:5d} {:5d} {:5d} {:5d}\n"\
              .format(moidx,aoidx,atidx,coef_cnt))
          coef_cnt += 1

    outf.write("COEFFICIENTS\n")
    eigreal_flat = [e[0:self.maxmo_spin,:].flatten() for e in eigvecs_real]
    eigimag_flat = [e[0:self.maxmo_spin,:].flatten() for e in eigvecs_imag]
    print_cnt = 0
    if self.eigsys['ikpt_iscmpx'][kpt]: #complex coefficients
      for eigr,eigi in zip(eigreal_flat,eigimag_flat):
        for r,i in zip(eigr,eigi):
          outf.write("({:<.12e},{:<.12e}) "\
              .format(r,i))
          print_cnt+=1
          if print_cnt%5==0: outf.write("\n")
    else: #Real coefficients
      for eigr in eigreal_flat:
        for r in eigr:
          outf.write("{:< 15.12e} ".format(r))
          print_cnt+=1
          if print_cnt%5==0: outf.write("\n")

    outf.close()
    return None

  ###
  def generate_sys(self,kpt):
    ''' Write the system file for QWalk. 

    Args:
      kpt (int): Which kpoint to write it for. 
      base (str): Base name to append the kpoint name to. 
    Returns: 
      str: contents of sys file.
    '''
    outlines = []
    min_exp = min(self.basis['prim_gaus'])
    cutoff_length = (-np.log(1e-8)/min_exp)**.5
    basis_cutoff = self.find_basis_cutoff()
    cutoff_divider = basis_cutoff*2.0 / cutoff_length
    if self.lat_parm['struct_dim'] != 0:
      outlines += [
          "system { periodic",
          "  nspin {{ {} {} }}".format(self.eigsys['nup'],self.eigsys['ndn']),
          "  latticevec {",
        ]
      for i in range(3):
        outlines.append("    {:< 15} {:< 15} {:< 15}".format(*self.lat_parm['latvecs'][i]))
      outlines += [
          "  }",
          "  origin { 0 0 0 }",
          "  cutoff_divider {0}".format(cutoff_divider),
          "  kpoint {{ {:4}   {:4}   {:4} }}".format(
              *(np.array(kpt)/self.eigsys['nkpts_dir']*2.)
            )
        ]
    else: # is molecule.
      outlines += [
          "system { molecule",
          "  nspin {{ {} {} }}".format(self.eigsys['nup'],self.eigsys['ndn']),
        ]
    for aidx in range(len(self.ions['positions'])):
      if self.ions['atom_nums'][aidx]-200-1 < 0:
        raise NotImplementedError("All-electron calculations not implemented yet.")
      outlines.append(
        "  atom {{ {0} {1} coor {2} }}".format(
          periodic_table[self.ions['atom_nums'][aidx]-200-1], # Assumes ECP.
          self.ions['charges'][aidx],
          "{:< 15} {:< 15} {:< 15}".format(*self.ions['positions'][aidx])
        )
      )
    outlines.append("}")
    done = []
    # TODO Generalize to no pseudopotential.
    for elem in self.pseudo.keys():
      atom_name = periodic_table[elem-200-1]
      n_per_j = self.pseudo[elem]['n_per_j']
      numL = sum(n_per_j>0)

      for i in range(1,len(n_per_j)):
        if (n_per_j[i-1]==0)and(n_per_j[i]!=0):
          raise NotImplementedError("ERROR: Weird pseudopotential, please generalize generate_sys(...).")

      n_per_j = n_per_j[n_per_j>0]
      order = list(np.arange(n_per_j[0],sum(n_per_j))) + \
              list(np.arange(n_per_j[0])) 
      exponents   = self.pseudo[elem]['exponents'][order]
      prefactors  = self.pseudo[elem]['prefactors'][order]
      r_exps      = self.pseudo[elem]['r_exps'][order]
      if numL > 2: aip = 12
      else:        aip =  6
      npjline = n_per_j[1:].tolist()+[n_per_j[0]]
      outlines += [
          "pseudo {",
          "  {}".format(atom_name),
          "  aip {:d}".format(aip),
          "  basis {{ {}".format(atom_name),
          "    rgaussian",
          "    oldqmc {",
          "      0.0 {:d}".format(numL),
          "      "+' '.join(["{}" for i in range(numL)]).format(*npjline)
        ]
      cnt = 0
      for eidx in range(len(exponents)):
        outlines.append("      {:d}   {:<12} {:< 12}".format(
          r_exps[cnt]+2,
          float(exponents[cnt]),
          float(prefactors[cnt])
        ))
        cnt += 1
      outlines += ["    }","  }","}"]
    return '\n'.join(outlines)

  ###
  def generate_jast2(self,optimizebasis=True):
    ''' Write a starting 2-body jastrow file for QWalk. 

    Args:
      optimizebasis (bool): Whether you want QWalk to optimize any basis functions present in the expansion.
    Returns: 
      str: contents of jast2 file.
    '''
    basis_cutoff = self.find_basis_cutoff()
    atom_types = [periodic_table[eidx-200-1] for eidx in self.ions['atom_nums']]
    atom_types=set(atom_types)
    outlines = [
        "jastrow2",
        "group {"
      ]
    if optimizebasis:
      outlines+=[
          "  optimizebasis"
        ]
    outlines+=[
        "  eebasis {",
        "    ee",
        "    cutoff_cusp",
        "    gamma 24.0",
        "    cusp 1.0",
        "    cutoff {0}".format(basis_cutoff),
        "  }",
        "  eebasis {",
        "    ee",
        "    cutoff_cusp",
        "    gamma 24.0",
        "    cusp 1.0",
        "    cutoff {0}".format(basis_cutoff),
        "  }",
        "  twobody_spin {",
        "    freeze",
        "    like_coefficients { 0.25 0.0 }",
        "    unlike_coefficients { 0.0 0.5 }",
        "  }",
        "}",
        "group {",
      ]
    if optimizebasis:
      outlines+=[
          "  optimizebasis"
        ]
    for atom_type in atom_types:
      outlines += [
        "  eibasis {",
        "    {0}".format(atom_type),
        "    polypade",
        "    beta0 0.2",
        "    nfunc 3",
        "    rcut {0}".format(basis_cutoff),
        "  }"
      ]
    outlines += [
        "  onebody {",
      ]
    for atom_type in atom_types:
      outlines += [
        "    coefficients {{ {0} 0.0 0.0 0.0}}".format(atom_type),
      ]
    outlines += [
        "  }",
        "  eebasis {",
        "    ee",
        "    polypade",
        "    beta0 0.5",
        "    nfunc 3",
        "    rcut {0}".format(basis_cutoff),
        "  }",
        "  twobody {",
        "    coefficients { 0.0 0.0 0.0 }",
        "  }",
        "}"
      ]
    return '\n'.join(outlines)

  ###
  def generate_basis(self):
    ''' Write a starting 2-body jastrow file for QWalk. 

    Args:
      optimizebasis (bool): Whether you want QWalk to optimize any basis functions present in the expansion.
    Returns: 
      str: contents of jast2 file.
    '''
    hybridized_check = 0.0
    hybridized_check += sum(abs(self.basis['coef_s'] * self.basis['coef_p']))
    hybridized_check += sum(abs(self.basis['coef_p'] * self.basis['coef_dfg']))
    hybridized_check += sum(abs(self.basis['coef_s'] * self.basis['coef_dfg']))
    if hybridized_check > 1e-10:
      raise NotImplementedError("Hybridized AOs (like sp) not implmemented in generate_basis(...)")

    # If there's no hybridization, at most one of coef_s, coef_p, and coef_dfg is
    # nonzero. Just add them, so we have one array.
    done_atoms = []
    coefs = self.basis['coef_s'] + self.basis['coef_p'] + self.basis['coef_dfg']

    shell_type = np.tile("Unknown...",self.basis['shell_type'].shape)
    typemap = ["S","SP","P","5D","7F_crystal","G","H"]
    for i in range(5): shell_type[self.basis['shell_type']==i] = typemap[i]

    cnt = 0
    aidx = 0
    atom_type = self.ions['atom_nums'][aidx]
    done_atoms.append(atom_type)
    outlines = [
        "basis {",
        "  {0}".format(periodic_table[atom_type-200-1]),
        "  aospline",
        "  normtype CRYSTAL",
        "  gamess {"
      ]
    for sidx in range(len(shell_type)):
      new_aidx = self.basis['atom_shell'][sidx]-1

      new_atom_type = self.ions['atom_nums'][new_aidx]
      if aidx != new_aidx:
        if new_atom_type in done_atoms:
          cnt+=self.basis['prim_shell'][sidx]
          continue
        else:
          outlines += ["  }","}"]
          atom_type = new_atom_type
          done_atoms.append(atom_type)
          aidx = new_aidx
          outlines += [
              "basis {",
              "  {0}".format(periodic_table[atom_type-200-1]),
              "  aospline",
              "  normtype CRYSTAL",
              "  gamess {"
            ]

      nprim = self.basis['prim_shell'][sidx]
      outlines.append("    {0} {1}".format(shell_type[sidx],nprim))
      for pidx in range(nprim):
        outlines.append("      {0} {1} {2}".format(
          pidx+1,
          self.basis['prim_gaus'][cnt],
          coefs[cnt]
        ))
        cnt += 1
    outlines += ["  }","}"]

    return '\n'.join(outlines)

  ###
  def write_files(self,base='qwalk',kfmt='coord',kset='complex'):
    ''' Write out all the QWalk system and wave function definition files. 

    Args:
      base (str): files are written with names [base].[filetype] and [base]_[kpt].[filetype].
      kfmt (str): 'coord' or 'int'. 'coord' would name the gamma point 000, int names it 0.
      kset (str): 'complex' or 'real'. Do you want only real or all the complex kpoints?
    Returns:
      dict: contains lists of each file produced organized by type.
    '''

    if (np.array(self.eigsys['kpt_coords']) >= 10).any():
      print("Cannot use coord kpoint format when SHRINK > 10.")
      print("Falling back on int format (old style).")
      kfmt = 'int'

   
    outfiles={
        'orb':[],
        'sys':[],
        'basis':[],
        'jast2':[],
        'slater':[]
      }
    
    with open("%s.jast2"%base,'w') as outf:
      outf.write(self.generate_jast2())
      outfiles['jast2'].append(outf.name)

    with open("%s.basis"%base,'w') as outf:
      outf.write(self.generate_basis())
      outfiles['basis'].append(outf.name)


    for kpt in self.eigsys['kpt_coords']:
      if self.eigsys['ikpt_iscmpx'][kpt] and kset=='real': continue
      if kfmt == 'int': 
        kbase = base + '_' + "{}".format(self.eigsys['kpt_index'][kpt])
      else:
        kbase = base + '_' + "{}{}{}".format(*kpt)

      with open("%s.sys"%kbase,'w') as outf:
        outf.write(self.generate_sys(kpt))
        outfiles['sys'].append(outf.name)

      orbfilename="%s.orb"%kbase
      self.write_orb(kpt,orbfilename)
      outfiles['orb'].append(orbfilename)

      with open("%s.slater"%kbase,'w') as outf:
        outf.write(
            self.generate_slater(
              kpt,orbfile=orbfilename,basisfile=outfiles['basis'][-1] ))
        outfiles['slater'].append(outf.name)

    return outfiles

###############################################################################
# Testing functions.
def test_basis(
    base="qwalk",
    propoutfn="prop.in.o",
    kfmt='coord',
    kset='complex'):
  print("Writing basis only.")

  info, lat_parm, ions, basis, pseudo = read_gred()
  write_basis(basis,ions,base=base)

def test_orb(
    base="qwalk",
    propoutfn="prop.in.o",
    kfmt='coord',
    kset='complex',
    kpt=0):
  print("Writing orb only.")

  info, lat_parm, ions, basis, pseudo = read_gred()
  eigsys = read_kred(info,basis)
  kpt=eigsys['kpt_coords'][kpt]
  basis['nmo']  = sum(basis['nao_shell']) # = nao
  write_orb(eigsys,basis,ions,kpt,base,kfmt)

###############################################################################
if __name__ == "__main__":
  # TODO clean up with proper arguement parsing.
  if len(sys.argv) > 1:
    base = sys.argv[1]
  else: 
    base = "qwalk"
  if len(sys.argv) > 2:
    propoutfn = sys.argv[2]
  else:
    propoutfn = "prop.in.o"
  if len(sys.argv) > 3:
    kfmt = sys.argv[3]
  else:
    kfmt="coord"
  if len(sys.argv) > 4:
    kset = sys.argv[4]
  else:
    kset="complex"
  print("Converting crystal with base {},".format(base))
  print("system spin drawn from {},".format(propoutfn))
  print("using {} kpoint naming convention,".format(kfmt))
  print("and using {} kpoint set.".format(kset))

  qw_writer=QWalkWriter()
  qw_writer.read_gredkred()
  qw_writer.write_files(base,kfmt,kset)
