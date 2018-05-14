
# coding: utf-8

# In[ ]:


# Access the autogen libraries.
import sys
sys.path.append('..')


# # Running a PySCF file.
# Autogen can run any file and keep track of its progress in a queue. 
# ### The `PySCFRunnerPBS` object.
# The runner contains the information about the resources you want to allocate to the program each time you attempt a run.
# This object also allows you to interact with a queueing program.
# This particular one works with PBS Torque.
# 
# *Note: you need to run this demo on a workstation with PBS Torque set up.*

# In[ ]:


from autorunner import PySCFRunnerPBS
import sys
pyscf_runner = PySCFRunnerPBS(
    queue='batch',        # This is the default queue for most PBS queues.
    walltime='1:00:00',
    jobname='03-queuing', # what the queue will report.
    ppath=sys.path,       # To run, we need to know the path to PySCF. I'm assuming its in your python path.
    nn=1,np=4
)


# ### Running a job.
# To run a job let's first make a simple input file.
# I'll use nitrogen, because hydrogen runs nearly instantly, and so won't be in the queue long enough to check.
# I also set the self-consistency tolerance extreamly small to give you some time to check the queue while its running.

# In[ ]:


# Simple input file generation. See 01-intro for example. 
from autopyscf import PySCFWriter
h2='\n'.join([
    'N 0.0  0.0 0.0 ',
    'N 1.1 0.0 0.0 '
  ])  
pyscf_writer=PySCFWriter({'xyz':h2,'direct_scf_tol':1e-12,'max_cycle':200})
pyscf_writer.pyscf_input('03-scratch.py','03-scratch.chkfile')


# The first step is to tell the runner that you want to run this command.
# This won't submit the job to the queue. 
# The reason it doesn't is that you may want to bundle the job---run multiple PySCF executables with one submission, for example.

# In[ ]:


# Run this file, please.
pyscf_runner.add_task('python3 03-scratch.py &> 03-scratch.py.out')
pyscf_runner.exelines


# Now call `submit` to send this to the queue.
# The queue id is stored in the `queueid` attribute.
# Call `check_status` to check on the job (right away, before it finishes!)

# In[ ]:


pyscf_runner.submit()
print("Queue id:",pyscf_runner.queueid)
pyscf_runner.check_status()


# Note the `queueid` is a list. 
# Each time something is run, the queue ID is saved in this list.
# If something funny happens with the submission, you can reference this information for debugging.
# The last submission is saved at the end.

# ### Checking the status of a job.
# 
# `check_status` will return 'unknown', 'running', or 'finished' depending on what information the queue has. 

# Try waiting a while and executing this command again.
# 
# If we wait long enough, the queue id is no longer in the system, so the status becomes 'unknown'. 
# Don't worry, this is normal!

# In[ ]:


pyscf_runner.check_status()


# Does it say "unknown" or "complete"? This usually means the job completed successfully. So the job is finished? Let's check its results.

# In[ ]:


from autopyscf import PySCFReader
pyscf_reader = PySCFReader()
print("Status:",pyscf_reader.collect('03-scratch.py.out','03-scratch.chkfile'))


# *Note: if you call this routine before the run is done, it will think the job is killed.
# You should wait until the runner reports 'complete' or 'unknown' and then call the reader.
# The next section automates this process.*

# # Next: automating the write-run-read process.
