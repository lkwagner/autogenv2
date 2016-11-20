from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

from pymatgen.io.cif import CifParser
from pymatgen.core.periodic_table import Element
from xml.etree.ElementTree import ElementTree


class CrystalWriter:
  def __init__(self,cif):
    self.kmesh=[8,8,8]
    self.spin_polarized=True
    self.initial_charges={}
    self.basis_params=[0.2,2,3]
    self.cutoff=0.2
    self.xml_name="BFD_Library.xml"
    self.struct=CifParser.from_string(cif).get_structures(primitive=False)[0]
    

########################################################
  def crystal_input(self):
    totinput=""
    geomlines = self.cif2geom()
    totinput+="\n".join(geomlines)
    basislines=self.basis_section()
    totinput+="\n".join(basislines)
    
    return totinput

    

########################################################

  def cif2geom(self):
    lat=self.struct.as_dict()['lattice']
    sites=self.struct.as_dict()['sites']
    geomlines=["CRYSTAL","0 0 0","1","%g %g %g %g %g %g"%\
          (lat['a'],lat['b'],lat['c'],lat['alpha'],lat['beta'],lat['gamma'])]

    geomlines+=["%i"%len(sites)]
    for v in sites:
      nm=v['species'][0]['element']
      nm=str(Element(nm).Z+200)
      geomlines+=[nm+" %g %g %g"%(v['abc'][0],v['abc'][1],v['abc'][2])]

    return geomlines
########################################################

  def basis_section(self):
    sites=self.struct.as_dict()['sites']
    elements=set()
    for s in sites:
      nm=s['species'][0]['element']
      elements.add(nm)
    basislines=[]
    elements = sorted(list(elements)) # Standardize ordering.
    for e in elements:
      basislines+=self.generate_basis(e)
    return basislines


########################################################
  def generate_basis(self,symbol):
    """
    Author: "Kittithat (Mick) Krongchon" <pikkienvd3@gmail.com> and Lucas K. Wagner
    Returns a string containing the basis section.  It is modified according to a simple recipe:
    1) The occupied atomic orbitals are kept, with exponents less than 'cutoff' removed.
    2) These atomic orbitals are augmented with uncontracted orbitals according to the formula 
        e_i = params[0]*params[2]**i, where i goes from 0 to params[1]
        These uncontracted orbitals are added for every occupied atomic orbital (s,p for most elements and s,p,d for transition metals)

    Args:
        symbol (str): The symbol of the element to be specified in the
            D12 file.

    Returns:
        str: The pseudopotential and basis section.


    Uses the following member variables:
        xml_name (str): The name of the XML pseudopotential and basis
            set database.
        cutoff: smallest allowed exponent
        params: parameters for generating the augmenting uncontracted orbitals
        initial_charges    
    """
    
    maxorb=3
    basis_name="vtz"
    nangular={"s":1,"p":1,"d":1,"f":1,"g":0}
    maxcharge={"s":2,"p":6,"d":10,"f":15}
    basis_index={"s":0,"p":2,"d":3,"f":4}
    transition_metals=["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"]
    if symbol in transition_metals:
      maxorb=4
      nangular['s']=2
    
    tree = ElementTree()
    tree.parse(self.xml_name)
    element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
    atom_charge = int(element.find('./Effective_core_charge').text)
    if symbol in self.initial_charges.keys():
      atom_charge-=self.initial_charges[symbol]
    basis_path = './Basis-set[@name="{}"]/Contraction'.format(basis_name)
    found_orbitals = []
    totcharge=0
    ret=[]
    ncontract=0
    for contraction in element.findall(basis_path):
        angular = contraction.get('Angular_momentum')
        if found_orbitals.count(angular) >= nangular[angular]:
            continue

        #Figure out which coefficients to print out based on the minimal exponent
        nterms = 0
        basis_part=[]
        for basis_term in contraction.findall('./Basis-term'):
            exp = basis_term.get('Exp')
            coeff = basis_term.get('Coeff')
            if float(exp) > self.cutoff:
                basis_part += ['  {} {}'.format(exp, coeff)]
                nterms+=1
        #now write the header 
        if nterms > 0:
          found_orbitals.append(angular)          
          charge=min(atom_charge-totcharge,maxcharge[angular])
          #put in a special case for transition metals:
          #depopulate the 4s if the atom is charged
          if symbol in transition_metals and symbol in self.initial_charges.keys() \
                  and self.initial_charges[symbol] > 0 and found_orbitals.count(angular) > 1 \
                  and angular=="s":
            charge=0
          totcharge+=charge
          ret+=["0 %i %i %g 1"%(basis_index[angular],nterms,charge)] + basis_part
          ncontract+=1

    #Add in the uncontracted basis elements
    angular_uncontracted=['s','p']
    if symbol in transition_metals:
        angular_uncontracted.append('d')

    for angular in angular_uncontracted:
        for i in range(0,self.basis_params[1]):
            exp=self.basis_params[0]*self.basis_params[2]**i
            line='{} {}'.format(exp,1.0)
            ret+=["0 %i %i %g 1"%(basis_index[angular],1,0.0),line]
            ncontract+=1

    return ["%i %i"%(Element(symbol).number+200,ncontract)] +\
            self.pseudopotential_section(symbol) +\
            ret
########################################################

  def pseudopotential_section(self,symbol):
    """
    Author: "Kittithat (Mick) Krongchon" <pikkienvd3@gmail.com>
    Returns a string of the pseudopotential section which is to be written
    as part of the basis set section.

    Args:
        symbol (str): The symbol of the element to be specified in the
            D12 file.
        xml_name (str): The name of the XML pseudopotential and basis
            set database.

    Returns:
        list of lines of pseudopotential section (edit by Brian Busemeyer).
    """
    tree = ElementTree()
    tree.parse(self.xml_name)
    element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
    eff_core_charge = element.find('./Effective_core_charge').text
    local_path = './Gaussian_expansion/Local_component'
    non_local_path = './Gaussian_expansion/Non-local_component'
    local_list = element.findall(local_path)
    non_local_list = element.findall(non_local_path)
    nlocal = len(local_list)
    m = [0, 0, 0, 0, 0]
    proj_path = './Gaussian_expansion/Non-local_component/Proj'
    proj_list = element.findall(proj_path)
    for projector in proj_list:
        m[int(projector.text)] += 1
    strlist = []
    strlist.append('INPUT')
    strlist.append(' '.join(map(str,[eff_core_charge,nlocal,
                                     m[0],m[1],m[2],m[3],m[4]])))
    for local_component in local_list:
        exp_gaus = local_component.find('./Exp').text
        coeff_gaus = local_component.find('./Coeff').text
        r_to_n = local_component.find('./r_to_n').text
        strlist.append(' '.join([exp_gaus, coeff_gaus,r_to_n]))
    for non_local_component in non_local_list:
        exp_gaus = non_local_component.find('./Exp').text
        coeff_gaus = non_local_component.find('./Coeff').text
        r_to_n = non_local_component.find('./r_to_n').text
        strlist.append(' '.join([exp_gaus, coeff_gaus,r_to_n]))
    return strlist


if __name__=="__main__":
  cwriter=CrystalWriter(open("si.cif").read())
  print(cwriter.crystal_input())
