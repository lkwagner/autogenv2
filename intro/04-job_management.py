
# coding: utf-8

# In[ ]:


# Access the autogen libraries.
import sys
sys.path.append('..')


# # Run management.
# Up to this point, the tasks have been elementary: writing, running, and reading a single task.
# The next step is to combine some of these tasks, and manage the files for the run. 
# This is what the manager is for.
# 
# Explicitly, the manager:
# * Stores the writer, reader, and runner on disk for after the script finishes.
# * Keeps track of file names and locations.
# * Communicates between these objects so you don't have to keep track of what is named what.
# 
# Let's stick with the hydrogen example.

# ### Manager for hydrogen DFT run.
# The job of this manager is to give you DFT results for the hydrogen molecule. 
# It has three "subordinates" at its disposal: the writer, the reader, and the runner from the last three notebooks.
# 
# First, lets set these elementary parts up.

# In[ ]:


from autopyscf import PySCFWriter,PySCFReader
from autorunner import PySCFRunnerLocal

h2='\n'.join([
    'H 0.0 0.0 0.0 ',
    'H 0.74 0.0 0.0 '
])

# The three components we've seen already:
pyscf_writer = PySCFWriter({
    'xyz':h2,
    'method':'RKS',
    'dft':'pbe,pbe'
})
pyscf_reader = PySCFReader()
pyscf_runner = PySCFRunnerLocal()


# The manger needs to know what space it has to work. 
# This is the `path` arguement.
# Data will be stored on disk at this location.
# 
# The manager also needs to have a `name`.
# All filenames are generated based on the `name` you give it. 
# You can have multiple managers with the same `path` so long as they have different `name`. 

# In[ ]:


from manager import PySCFManager
pyscf_manager = PySCFManager(
    path='04-scratch',
    name='h2_pbe',
    writer=pyscf_writer,
    reader=pyscf_reader,
    runner=pyscf_runner
)


# The manager will also talk you though all the events that are happening. 

# ### Perfoming the next step.
# If you've made it this far, the rest is easy. 
# Just keep calling `nextstep`. 
# This automatically will call the routines of the reader, the runner, and the writer according to the current status of each of their results.

# In[ ]:


pyscf_manager.nextstep()


# The manager had the input file written, and ran it using the runner.
# Let's check out the path.

# In[ ]:


import os
os.listdir(pyscf_manager.path)


# You can see a bunch of file were made in the `path` named according to `name`. 
# `h2_pbe.pkl` is where all the objects are stored for future script runs.
# Any other instance of a manager with the same `name` and `path` will reload `h2_pbe.pkl` to recover all the properties and results that this run produces.
# 
# Ok, now on to the next step, whatever that may be.

# In[ ]:


pyscf_manager.nextstep()


# If the status was "running" go ahead and run that last cell again (the job was still in progress).
# 
# Once the last cell reads "ready_for_analysis" that means the job has finished, and you can read in the results.

# ### Accessing the results.
# Well, this is just the same as when working with the reader object.

# In[ ]:


pyscf_manager.reader.output.keys()


# ### Multiple managers, same path.
# Managers are flexible in how you organize your work. 
# Lets do another run in the same folder, but using the PBE0 DFT functional.

# In[ ]:


from copy import deepcopy

# Same thing,  but change the functional.
pbe0_writer = deepcopy(pyscf_writer)
pbe0_writer.set_options({'dft':'pbe0'})

# deepcopy everything to prevent accidentally overwriting anything.
pbe0_manager = PySCFManager(
    path=pyscf_manager.path, # Will work in the same directory.
    name='h2_pbe0',         # New name so that we can keep these results seperate.
    writer=pbe0_writer,
    reader=deepcopy(pyscf_reader),
    runner=deepcopy(pyscf_runner)
)


# In[ ]:


os.listdir(pbe0_manager.path)


# There are two `.pkl` files, one for each manager.
# The `name` is different, so they can live side-by-side.

# # Next: Starting a QMC calculation (05-qwalk.ipynb).
