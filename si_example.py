import submission_tools,CrystalWriter,CrystalRun,PropertiesWriter,PropertiesRun
import local
import os,json


setup={'id':'si',
        'job_record':submission_tools.JobRecord(),
       'crystal':CrystalWriter.CrystalWriter(open("si.cif").read() ),
       'properties':PropertiesWriter.PropertiesWriter() 
         }

setup['crystal'].set_options( {'xml_name':"../BFD_Library.xml",
                               'kmesh':[2,2,2],
                               'basis_params':[0.3,1,3],
                               'tolinteg':[8,8,8,8,16],
                               'spin_polarized':False,
                               'dftgrid':'LGRID'
                               }
                               ) 

runcrys=CrystalRun.CrystalRun(local.LocalCrystal(),setup['crystal'])
runprop=PropertiesRun.PropertiesRun()

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

