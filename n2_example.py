import submission_tools,CrystalWriter,CrystalRunner,PropertiesRun
import local
import os,json


setup={'id':'n2',
        'job_record':submission_tools.JobRecord()
        }
setup['crystal']=CrystalWriter.CrystalWriter(xyz=open("n2.xyz").read())
setup['crystal'].set_options({
    'xml_name':"../BFD_Library.xml",
    'basis_params':[0.2,0,3],
    'cutoff':0.0,
    'dftgrid':'LGRID',
    'spin_polarized':False
  })

runcrys=CrystalRunner.CrystalRunner(local.LocalCrystal(),setup['crystal'])
runprop=PropertiesRun.PropertiesRun(setup['crystal'])

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

