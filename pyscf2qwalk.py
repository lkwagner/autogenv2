import sys
import numpy as np
from pyscf import gto 

##TODO:
## Find out what order the d and f atomic orbitals are in and check vs QWalk's ordering
## Fix the orbital coefficients (see gamess2qmc for how to do it) for normalization of the basis functions.
## CI coefficients and occupations.
## Periodic boundary conditions.
###########################################################
def print_orb(mol,m,f):
  aos_atom=mol.offset_nr_by_atom()
  coeff=m.mo_coeff
  nmo=len(coeff)
  count=0
  for ni in range(nmo):
    for ai,a in enumerate(aos_atom):
      for bi,b in enumerate(range(a[2],a[3])):
        f.write("%i %i %i %i\n"%(ni+1,bi+1,ai+1,count+1))
        count += 1

  count=0
  f.write("COEFFICIENTS\n")
  print(coeff.shape)
  for a in coeff.T:
    for b in a:
      f.write(str(b)+" ")
      count+=1
      if count%10==0:
        f.write("\n")

  f.write("\n")
  f.close() 
  print("count",count)
  return 

###########################################################

def print_basis(mol, f):
  aos_atom= mol.offset_nr_by_atom()
  mom_to_letter ={0:'s', 1:'p', 2:'d', 3:'f'}

  for at, li  in enumerate(aos_atom):
    f.write('''Basis { \n%s  \nAOSPLINE \nGAMESS {\n ''' %(mol.atom_pure_symbol(at) )) 
 
    for bas in range(li[0], li[1]):
      letter=mom_to_letter[mol.bas_angular(bas)]
      ex=mol.bas_exp(bas)
      c=mol.bas_ctr_coeff(bas)
      if len(c[0]) >1:
         print("skipping some coefficients")

      nbas=len(ex)
      f.write("%s  %i\n" %(letter,nbas))

      for i in range(nbas):
          f.write("  %d %f %f \n" %(i+1,ex[i],c[i][0]))
    f.write("}\n")
    f.write("}\n")
  f.close() 
  return 


###########################################################

def print_sys(mol, f):
  coords = mol.atom_coords()
  symbols=[]
  charges=[]
  coords_string=[]
    
  # write coordinates
  for i in range(len(coords)):
    symbols.append(mol.atom_pure_symbol(i))
    charges.append(mol.atom_charge(i))
    c = list(map(str, coords[i]))
    coords_string.append(' '.join(c))

  T_charge = sum(charges)
  T_spin = mol.spin
  spin_up =(T_charge + T_spin)/2
  spin_down = (T_charge - T_spin)/2

  f.write('SYSTEM { MOLECULE \nNSPIN { %d  %d }\n' %(spin_up, spin_down))
  for i in range(len(coords)):
    f.write('ATOM { %s  %d  COOR  %s }\n'   %(symbols[i], charges[i], coords_string[i]))
  f.write('}\n')

  # write pseudopotential 
  for i in range(len(coords)):
    ecp =  gto.basis.load_ecp(mol.ecp,symbols[i])
    f.write('''PSEUDO { \n %s \nAIP %d \nBASIS { \n%s \nRGAUSSIAN \nOLDQMC {\n  ''' %(symbols[i], charges[i], symbols[i]))
    coeff =ecp[1]
    numofpot =  len(coeff)  
    data=[]     
    N=[]        
    for i in range(1, len(coeff)):
      num= coeff[i][0]  
      eff =coeff[i][1]  
      count=0           
      for j, val in enumerate(eff):
        if(len(eff[j])>0) :     
          count+=1                      
          data.append([j]+eff[j])       
      N.append(str(count))

    num= coeff[0][0] 
    eff =coeff[0][1]
    count=0     
    for j, val in enumerate(eff):
      if(len(eff[j])>0) :
        count+=1                
        data.append([j]+eff[j]) 
    N.append(str(count))

    C=0.0       
    f.write("%f  %d\n" %(C, numofpot))
    f.write(' '.join(N)+ '\n') 
    for val in data:
      f.write(str(val[0]) +' ')
      S= list(map(str, val[1]))
      f.write(' '.join(S) + '\n')

    f.write("} \n")
    f.write("} \n")
    f.write("} \n")
  f.close() 
  return 
###########################################################

def print_slater(mol, mf, orbfile, basisfile, f):
  coeffs=(mf.mo_coeff)
  occ=mf.get_occ()

  if (isinstance(coeffs[0][0], np.float64)):
    orb_type = 'ORBITALS'
  else:
    orb_type = 'CORBITALS'

  te=sum(occ)
  ts=mol.spin
  us= (ts+te)/2
  ds= (te-ts)/2
  us_orb=[]
  ds_orb=[]
  max_orb=0
  for i in range(len(occ)):
    if(len(us_orb) < us and occ[i]>0):
      us_orb.append(i+1)
      max_orb= max(max_orb, i+1)
      occ[i]-=1 
      
    if(len(ds_orb)< ds and  occ[i]>0):
      ds_orb.append(i+1)
      max_orb= max(max_orb, i+1)

  up_orb=' '.join(list(map(str, us_orb)))
  down_orb= ' '.join(list(map(str, ds_orb)))

  f.write('''SLATER
  %s  { 
	CUTOFF_MO
  MAGNIFY 1.0 
  NMO %d 
  ORBFILE %s
  INCLUDE %s
  CENTERS { USEATOMS }
  }

  DETWT { 1.0 }
  STATES {
  #Spin up orbitals 
  %s 
  #Spin down orbitals 
  %s 
  }''' %(orb_type, max_orb,orbfile, basisfile, up_orb, down_orb ))
  f.close()
  return 
###########################################################


def print_qwalk(mol,mf, basename='qw'):
  orbfile=basename+".orb"
  basisfile=basename+".basis"

  print_orb(mol,mf,open(orbfile,'w'))
  print_basis(mol,open(basisfile,'w'))
  print_sys(mol,open(basename+".sys",'w'))
  print_slater(mol,mf,orbfile,basisfile,open(basename+".slater",'w'))
