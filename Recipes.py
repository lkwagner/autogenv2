import Manager as mgmt
from copy import deepcopy
from Crystal import CrystalWriter,CrystalReader
from CrystalRunner import LocalCrystalRunner,CrystalRunnerPBS
from PropertiesReader import PropertiesReader
from PropertiesRunner import LocalPropertiesRunner
import average_tools as avg
import copy
import numpy as np

def separate_jastrow(f):
  tokens=f.readlines()
  in_jastrow=False
  nopen=0
  nclose=0
  ret=""
  for line in tokens:
    if line.find("JASTROW2") != -1:
      in_jastrow=True
    if in_jastrow:
      nopen+=line.count("{")
      nclose+=line.count("}")
    if in_jastrow and nopen >= nclose:
      ret+=line
  return ret
  

#########################################################
class Recipe:
  """ Contains DFT and QMC steps in the order they need to be performed.

  * One job per directory. 
  * For most users, a child class of this should be used."""
  
  def __init__(self,jobid,jobplans,managers):
    self.jobid=jobid
    self.managers=managers
    self.picklefn="%s.pickle"%jobid
#---------------------------------------
  def is_consistent(self,other):
    result=True
    if len(other.managers)!=len(self.managers):
      print('You have added or removed tasks for this job.')
      result=False


    for rec_manager,plan_manager in zip(other.managers,self.managers):
      plancheck=plan_manager.is_consistent(rec_manager)
      if plancheck==False:
        print('You have modified a job.')
        result=False
    return result
#---------------------------------------

  def nextstep(self):
    for manager in self.managers:
      manager.nextstep()
      if manager.status()!='ok':
        break
#---------------------------------------
  def write_summary(self):
    for manager in self.managers:
      if manager.status()=='ok':
        manager.write_summary()
    
#---------------------------------------
  def generate_report(self):
    print("generate_report not implemented for this Recipe")
    return {'id':self.jobid}
    

##########################################################
class LocalCrystalDFT(Recipe):
  """ An example of a Recipe that perfoms a crystal DFT calculation """
  
  def __init__(self,jobid,struct,crystal_opts,structtype='cif'):
    # May have it automatically detect file type? Probably wouldn't be too hard.
    inpcopy=deepcopy(crystal_opts)
    self.jobid=jobid

    #TODO primitive option.
    cwriter=CrystalWriter()
    if structtype=='cif':
      cwriter.set_struct_fromcif(struct)
    elif structtype=='xyz':
      cwriter.set_struct_fromxyz(struct)
    else:
      raise ValueError("structtype not recognized.")
    cwriter.set_options(crystal_opts)


    # For this simple case, only one Manager is needed.
    self.managers=[mgmt.CrystalManager(
        cwriter,
        CrystalReader(),
        LocalCrystalRunner(),
        PropertiesReader(),
        LocalPropertiesRunner()
      )]
    self.picklefn="%s.pickle"%jobid

##########################################################

from Crystal2QMCRunner import LocalCrystal2QMCRunner
from Crystal2QMCReader import Crystal2QMCReader
from Variance import VarianceWriter,VarianceReader
from Linear import LinearWriter,LinearReader
from PostProcess import PostProcessWriter,PostProcessReader
from QWalkRunner import LocalQWalkRunner,QWalkRunnerPBS
from DMC import DMCWriter,DMCReader
class LocalCrystalQWalk(Recipe):
  """ In this we will perform the following recipe:
    1) A Crystal calculation. 
    2) Convert the Crystal calculation to QWalk, form a Slater determinant trial function.
    3) Run variance optimization on a Jastrow factor for the gamma point.
    4) Remove OPTIMIZEBASIS from Jastrow, run energy optimization using LINEAR. 
    5) Run DMC on all k-points, saving configurations to a .trace file.
    6) Run properties on the .trace file.
    
    """
  
  def __init__(self,jobid,struct,
               crystal_opts={},
               variance_opts={},
               energy_opts={},
               dmc_opts={},
               structtype='cif',
               crystalrunner=CrystalRunnerPBS(),
               qwalkrunner=QWalkRunnerPBS(np=6)):
    # May have it automatically detect file type? Probably wouldn't be too hard.
    inpcopy=deepcopy(crystal_opts)
    self.jobid=jobid

    cwriter=CrystalWriter()
    if structtype=='cif':
      cwriter.set_struct_fromcif(struct)
    elif structtype=='xyz':
      cwriter.set_struct_fromxyz(struct)
    else:
      raise ValueError("structtype not recognized.")
    cwriter.set_options(crystal_opts)

    self.managers=[mgmt.CrystalManager(
        cwriter,
        crystalrunner,
        CrystalReader(),
        LocalPropertiesRunner(),
        PropertiesReader()
      ),
      mgmt.QWalkfromCrystalManager(
        LocalCrystal2QMCRunner(),
        Crystal2QMCReader()        
        ),
      mgmt.QWalkRunManager(
        VarianceWriter(variance_opts),
        copy.deepcopy(qwalkrunner),
        VarianceReader()
        ),
      mgmt.QWalkRunManager(
        LinearWriter(energy_opts),
        copy.deepcopy(qwalkrunner),
        LinearReader()
        ),
      mgmt.QWalkRunManager(
        DMCWriter(dmc_opts),
        copy.deepcopy(qwalkrunner),
        DMCReader()
        )
      ]
    self.picklefn="%s.pickle"%jobid

  #--------------------------------------------
  def nextstep(self):
    cry=0 #crystal index
    con=1 #converter index
    var=2 #variance index
    en=3 #energy index
    dmc=4 

    self.managers[cry].nextstep()
    if self.managers[cry].status()!='ok':
      return 

    self.managers[con].nextstep()
    if self.managers[con].status()!='ok':
      print("I think crystal is still running")
      return

    bases=self.managers[con].reader.out['basenames']
    ind=bases.index('qw_000')
    files={}
    for key in ['sysfiles','slaterfiles','jastfiles','basenames']:
      files[key]=[self.managers[con].reader.out[key][ind]]
    self.managers[var].writer.set_options(files)
    self.managers[var].nextstep()
    if self.managers[var].status()!='ok':
      return
   
    files={'basenames':[],
           'sysfiles':[],
           'wffiles':[] } 
    
    for base in [bases[ind]]:#self.managers[con].reader.out['basenames']:
      files['basenames'].append(base)
      files['wffiles'].append(base+".energy.wfin")
      files['sysfiles'].append(base+".sys")
      with open(base+".variance.wfout") as fin:
        fout=open(base+".energy.wfin",'w')
        for line in fin:
          fout.write(line.replace("OPTIMIZEBASIS",''))
        fout.close()
      

    self.managers[en].writer.set_options(files)
    self.managers[en].nextstep()
    if self.managers[en].status()!='ok':
      return


    jast=separate_jastrow(open("qw_000.energy.wfout"))
    files={'basenames':[],
           'sysfiles':[],
           'wffiles':[] } 
    for i in bases:
      wfname=i+'.dmc.wf'
      with open(wfname,'w') as f:
        f.write("slater-jastrow  \n" +\
            "wf1 { include %s.slater }\n"%i +\
            "wf2 { " + jast + "} \n ")
        f.close()
      files['wffiles'].append(wfname)
      files['sysfiles'].append(i+".sys")
      files['basenames'].append(i)

    self.managers[dmc].writer.set_options(files)
    self.managers[dmc].nextstep()
    if self.managers[dmc].status()!='ok':
      return
    
  #---------------------------------------
  def generate_report(self):
    cry=0 #crystal index
    con=1 #converter index
    var=2 #variance index
    en=3 #energy index
    dmc=4 
    ret={'id':self.jobid}
    
    if self.managers[cry].status()=='ok':
      ret['crystal_energy']=self.managers[cry].creader.out['total_energy']

    if self.managers[var].status()=='ok':
      varopt={}
      for f,out in self.managers[var].reader.output.items():
        sigma=[]
        for run in out:
          sigma.extend(run['sigma'])
        varopt[f]=sigma
      ret['variance_optimization']=varopt

    if self.managers[en].status()=='ok':
      enopt={}
      for f,out in self.managers[en].reader.output.items():
        en=[]
        err=[]
        for run in out:
          en.extend(run['energy'])
          err.extend(run['energy_err'])
        enopt[f]={'energy':copy.deepcopy(en),
                 'energy_err':copy.deepcopy(err)}
      ret['energy_optimization']=enopt
      
    if self.managers[dmc].status()=='ok':
      #here we average over k-points
      dmcret={'timestep':[],'energy':[],'energy_err':[]}
      basenames=self.managers[con].reader.out['basenames']
      timesteps=self.managers[dmc].writer.timesteps
      for t in timesteps:
        ens=[]
        errs=[]
        for base in basenames:
          nm=base+'t'+str(t)+".dmc.log"
          ens.append(self.managers[dmc].reader.output[nm]['properties']['total_energy']['value'][0])
          err.append(self.managers[dmc].reader.output[nm]['properties']['total_energy']['error'][0])
        dmcret['timestep'].append(t)
        dmcret['energy'].append(np.mean(ens))
        dmcret['energy_err'].append(np.sqrt(np.mean(np.array(err)**2)))
      ret['dmc']=dmcret
    return ret
        
    
    
#######################################################
from PySCF import PySCFWriter,PySCFReader
from PySCFRunner import PySCFRunnerPBS

class PySCFQWalk(Recipe):
  """ Use PySCF to generate a QWalk run. """
  
  def __init__(self,jobid,
               pyscf_opts={},
               variance_opts={},
               energy_opts={},
               dmc_opts={},
               post_opts={},
               pyscfrunner=PySCFRunnerPBS(np=4),
               qwalkrunner=QWalkRunnerPBS(np=4)):
    self.jobid=jobid
    self.picklefn="%s.pickle"%jobid

    assert post_opts=={} or dmc_opts['savetrace'],"""
      You need to save the trace (dmc_opts['savetrace']=True) to use postprocess options."""

    self.managers=[mgmt.PySCFManager(
                                     PySCFWriter(pyscf_opts),
                                     copy.deepcopy(pyscfrunner),
                                     PySCFReader()
                                    ),
      mgmt.QWalkRunManager(
                           VarianceWriter(variance_opts),
                           copy.deepcopy(qwalkrunner),
                           VarianceReader()
                          ),
      mgmt.QWalkRunManager(
                           LinearWriter(energy_opts),
                           copy.deepcopy(qwalkrunner),
                           LinearReader()
                           ),
      mgmt.QWalkRunManager(
                           DMCWriter(dmc_opts),
                           copy.deepcopy(qwalkrunner),
                           DMCReader()
                          ),
      mgmt.QWalkRunManager(
                           PostProcessWriter(post_opts),
                           copy.deepcopy(qwalkrunner),
                           PostProcessReader()
                          )
      ]
  #-----------------------------
  def nextstep(self):
    pyscf=0 #crystal index
    var=1 #variance index
    en=2 #energy index
    dmc=3 
    post=4 

    # PySCF.
    self.managers[pyscf].nextstep()
    if self.managers[pyscf].status()!='ok':
      return 

    # Variance minimization.
    base='qw'
    files={}
    files['sysfiles']=[base+'.sys']
    files['slaterfiles']=[base+'.slater']
    files['basenames']=[base]
    files['jastfiles']=[base+'.jast2']
    self.managers[var].writer.set_options(files)
    self.managers[var].nextstep()
    if self.managers[var].status()!='ok':
      return
   
    # Energy minimization.
    files={'basenames':[],
           'sysfiles':[],
           'wffiles':[] } 
    
    files['basenames'].append(base)
    files['wffiles'].append(base+".energy.wfin")
    files['sysfiles'].append(base+".sys")
    with open(base+".variance.wfout") as fin:
      fout=open(base+".energy.wfin",'w')
      for line in fin:
        fout.write(line.replace("OPTIMIZEBASIS",'').replace(" SLATER\n","SLATER OPTIMIZE_DET\n"))
      fout.close()
      

    self.managers[en].writer.set_options(files)
    self.managers[en].nextstep()
    if self.managers[en].status()!='ok':
      return

    # DMC.

    # Why does it need to seperate the jastrow?
    #jast=separate_jastrow(open("qw.energy.wfout"))
    files={'basenames':[],
           'sysfiles':[],
           'wffiles':[],
           'tracefiles':[]} 
    for i in [base]:
      wfname=i+'.dmc.wf'
      # Just copy over the results from the VMC energy minimization.
      with open(wfname,'w') as outf:
        outf.write(open(base+'.energy.wfout').read())
      files['wffiles'].append(wfname)
      files['sysfiles'].append(i+".sys")
      files['basenames'].append(i)
      files['tracefiles'].append(i+".trace")

    self.managers[dmc].writer.set_options(files)
    self.managers[dmc].nextstep()
    if self.managers[dmc].status()!='ok':
      return

    # Post process.

    self.managers[post].writer.set_options(files)
    self.managers[post].nextstep()
    if self.managers[post].status()!='ok':
      return

    
  #-----------------------------
  def generate_report(self):
    pyscf=0 #pyscf index
    var=1 #variance index
    en=2 #energy index
    dmc=3 
    post=4
    ret={'id':self.jobid}
    
    # Collect from PySCF.
    if self.managers[pyscf].status()=='ok':
      pyout={} 
      for f, out in self.managers[pyscf].reader.output.items(): 
        pyout[f]=out['energy']  
      ret['pyscf_energy']=pyout

    # Collect from VMC variance optimization.
    if self.managers[var].status()=='ok':
      varopt={}
      for f,out in self.managers[var].reader.output.items():
        sigma=[]
        for run in out:
          sigma.extend(run['sigma'])
        varopt[f]=sigma
      ret['variance_optimization']=varopt

    # Collect from VMC energy optimization.
    if self.managers[en].status()=='ok':
      enopt={}
      for f,out in self.managers[en].reader.output.items():
        en=[]
        err=[]
        for run in out:
          en.extend(run['energy'])
          err.extend(run['energy_err'])
        enopt[f]={'energy':en,'energy_err':err}
      ret['energy_optimization']=enopt
      
    # Collect from DMC. 
    if self.managers[dmc].status()=='ok':
      extra_obs=self.managers[dmc].writer.extra_observables
      basenames=self.managers[dmc].writer.basenames
      timesteps=self.managers[dmc].writer.timesteps

      dmcret={'timestep':[],'energy':[],'energy_err':[]}
      for obs in extra_obs:
        dmcret[obs['name']]=[]

      for t in timesteps:
        dmcret['timestep'].append(t)

        # Energy results.  ens=[] errs=[]
        for base in basenames:
          nm=base+'t'+str(t)+".dmc.log"
          en=self.managers[dmc].reader.output[nm]['properties']['total_energy']
          ens.append(en['value'][0])
          errs.append(en['error'][0])
        # k-average.
        dmcret['energy'].append(np.mean(ens))
        dmcret['energy_err'].append(np.sqrt(np.mean(np.array(errs)**2)))

        # Property results (if any).
        for obs in extra_obs:
          dmcret[obs['name']]=deepcopy(obs)
          fnames=[base+'t'+str(t)+".dmc.log" for base in basenames]
          allk=[self.managers[dmc].reader.output[nm]['properties'][avg.gosling_key(obs['name'])]
              for nm in fnames]
          dmcret[obs['name']].update(avg.kaverage(obs['name'],allk))

      ret['dmc']=dmcret
    return ret
