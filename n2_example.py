import submitter,Manager
import CrystalWriter,CrystalRunner,CrystalReader
import PropertiesRunner,PropertiesReader
import os,json


setup={'id':'n2'}
setup['crystal']=CrystalWriter.CrystalWriter(xyz=open("n2.xyz").read())
setup['crystal'].set_options({
    'xml_name':"../BFD_Library.xml",
    'basis_params':[0.2,0,3],
    'cutoff':0.0,
    'dftgrid':'LGRID',
    'spin_polarized':False
  })

testjob = Manager.CrystalManager(
    writer=setup['crystal'],
    crys_reader=CrystalReader.CrystalReader(),
    crys_runner=CrystalRunner.LocalCrystalRunner(),
    prop_reader=PropertiesReader.PropertiesReader(),
    prop_runner=PropertiesRunner.LocalPropertiesRunner()
  )

currwd=os.getcwd()
d=setup['id']
try:
  os.mkdir(d)
except:
  pass
os.chdir(d)

testjob.nextstep()


