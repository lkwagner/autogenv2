import submission_tools,CrystalWriter,CrystalRun,PropertiesWriter,PropertiesRun
import local
import os,json


setup={'id':'h2',
        'job_record':submission_tools.JobRecord()
        }
setup['crystal']=CrystalWriter.CrystalWriter(xyz=open("h2.xyz").read())
setup['properties']=PropertiesWriter.PropertiesWriter(setup['crystal'])

setup['crystal'].xml_name="../BFD_Library.xml"
setup['crystal'].basis_params=[0.2,0,3]
setup['crystal'].cutoff=0.0    
setup['crystal'].dftgrid='LGRID'
setup['crystal'].spin_polarized=False

runcrys=CrystalRun.CrystalRun(local.LocalCrystal(),setup['crystal'])
runprop=PropertiesRun.PropertiesRun(setup['properties'])

currwd=os.getcwd()
d=setup['id']
try:
  os.mkdir(d)
except:
  pass
os.chdir(d)

runcrys.run(setup['job_record'])
runprop.run(setup['job_record'])

print(runcrys.check_status(setup['job_record']))

setup['crystal_output']=runcrys.output()

