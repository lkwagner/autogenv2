import sys

autogen_path="../../../"
sys.path.append(autogen_path)

from Crystal import CrystalWriter,CrystalReader


########  Set options. These are stored within the CrystalWriter object
options={'spin_polarized':False,
         'xml_name':autogen_path+"/BFD_PBC.xml",
         'functional':{'exchange':'PBE','correlation':'PBE','hybrid':0},
         'diis':True,
         'fmixing':95,
         'cutoff':0.1,
         'basis_params':[0.2,3,2],
         'restart':False,
         'kmesh':[2,2,2]
         }

######### Write input file

cwriter=CrystalWriter()
cwriter.set_options(options)
cwriter.set_struct_fromcif(open("si.cif").read())
f=open("crystal.inp",'w')
f.write(cwriter.crystal_input())
f.close()


###########  Check on job

import os
outname='crystal.inp.o'
if os.path.isfile(outname):
  creader=CrystalReader()
  status=creader.check_outputfile(outname)
  print("Status", status)
  if status=='not_finished':
    print("not finished")
  creader.collect(outname)
  print("Last few energies:",creader.out['etots'])
    
  
