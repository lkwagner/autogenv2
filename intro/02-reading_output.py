
# coding: utf-8

# In[ ]:


# Access the autogen libraries.
import sys
sys.path.append('..')


# # Reading a PySCF output file.
# 
# After a file is succesfully run, we can read the results into an object and check the results are complete.
# 
# ### The `PySCFReader` object.
# `PySCFReader` reads in the results of a `PySCFWriter` produced-input file that has been run. 
# It can do quick diagonostics of the results and store information.

# In[ ]:


from autopyscf import PySCFReader
pyscf_reader = PySCFReader()


# ### Reading in some sample output.
# Let's generate some sample output for reading in. 

# In[ ]:


# Simple input file generation. See 01-intro for example. 
from autopyscf import PySCFWriter
import subprocess as sub
h2='\n'.join([
    'H 0.0  0.0 0.0 ',
    'H 0.7  0.0 0.0 '
  ])  
pyscf_writer=PySCFWriter({'xyz':h2})
pyscf_writer.pyscf_input('02-scratch-1.py','02-scratch-1.chkfile')
sub.run('python3 02-scratch-1.py > 02-scratch-1.py.out',shell=True)


# Now we can collect the results

# In[ ]:


status=pyscf_reader.collect('02-scratch-1.py.out','02-scratch-1.chkfile')
print("Status of job:",status)


# ### Accessing results of interest.
# After collection, the reader holds the information about the run in its member variable `output`.
# `output` will consist of a nested python dictionary that is easy to explore.

# In[ ]:


pyscf_reader.output.keys()


# In[ ]:


pyscf_reader.output['scf'].keys()


# In[ ]:


# Here are the energies of the MOs, for example.
pyscf_reader.output['scf']['mo_energy']


# In[ ]:


# For posterity, it stores the name of the file where this came from.
pyscf_reader.output['file']


# ### Error detection.
# Here's an example of how the reader can detect a problem has occured.
# To generate erroneous output, I'll give PySCF some unreasonable settings.

# In[ ]:


#PySCF is only allowed to do two steps.
pyscf_writer.set_options({'max_cycle':2})
pyscf_writer.pyscf_input('02-scratch-2.py','02-scratch-2.chkfile')
sub.run('python3 02-scratch-2.py > 02-scratch-2.py.out',shell=True)

# The output file says "SCF not converged" which the PySCFReader can detect.
status=pyscf_reader.collect('02-scratch-2.py.out','02-scratch-2.chkfile')
print("Status of job:",status)


# The results of these checks can be used to automatically decide the next course of action. 
# The next section utilizes this information to do some basic automated error correction.
# For example, if the job is killed it will attempt to restart it, and continue the calcultion.

# # Next: automated running of jobs in batch (03-queuing.ipynb).
