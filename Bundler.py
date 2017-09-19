import numpy as np
import subprocess as sub

class Bundler:
  ''' Class for handling the bundling of several jobs of approximately the same 
  length, but possibly in different locations. ''' 
  def __init__(self,queue='normal',
                    walltime='48:00:00',
                    jobname='AGBundler',
                    npb=16,ppn=32,
                    prefix=None,
                    postfix=None
                    ):
    ''' npb is the number of nodes desired per bundle. '''
    self.npb=npb
    self.ppn=ppn
    self.jobname=jobname
    self.jobs=[]
    self.queue=queue
    self.walltime=walltime
    if prefix is None: self.prefix=[]
    else:              self.prefix=prefix
    if postfix is None: self.postfix=[]
    else:               self.postfix=postfix
    self.queueid=[]

  def add_job(self,mgr):
    ''' mgr is a Manager. Add Managers that have a script ready 
    to run in their current directory.'''
    if mgr._runready: self.jobs.append(mgr)

  def _submit_bundle(self,mgrs,jobname=None,nn=None):
    if nn is None:      nn=sum([mgr.runner.nn for mgr in mgrs])
    if jobname is None: jobname=self.jobname

    qsublines=[
        "#PBS -q %s"%self.queue,
        "#PBS -l nodes=%i:ppn=%i:xe"%(self.nn,self.ppn),
        "#PBS -l walltime=%s"%self.walltime,
        "#PBS -j oe ",
        "#PBS -A bahu",
        "#PBS -N %s "%jobname,
        "#PBS -o %s.out "%jobname,
      ]
    for mgr in mgrs:
      # This might be better without an error-out.
      assert mgr._runready, "One of the Managers is not prepped for run."
      qsublines+=[
          "cd %s"%mgr.location,
          "bash %s &"%mgr.scriptfile
        ]
    qsublines+=["wait"]

    qsubfile=jobname+".qsub"
    with open(qsubfile,'w') as f:
      f.write('\n'.join(qsublines))
    try:
      result=sub.check_output("qsub %s"%(qsubfile),shell=True)
      queueid=result.decode().split()[0].split('.')[0]
      print("Submitted as %s"%queueid)
    except sub.CalledProcessError:
      print("Error submitting job. Check queue settings.")

    for mgr in mgrs:
      mgr.update_queueid(queueid)


  def submit(self,jobname=None):
    ''' Submit all the jobs in the Managers that were added.'''
    if jobname is None: jobname=self.jobname

    assign=np.cumsum([mgr.runner.nn for mgr in self.jobs])
    assign=((assign-0.1)//self.npb).astype(int)

    print(assign)

    for bidx in range(assign[-1]+1):
      self._submit_bundle(np.array(self.jobs)[assign==bidx],"%s_%d"%(jobname,bidx))
