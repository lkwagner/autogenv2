
# coding: utf-8

# In[ ]:


# Access the autogen libraries.
import sys
sys.path.append('..')


# # Generating a PySCF input.
# 
# One of the simplest uses of autogen is to generate input files in succinct ways and have an object store the options for posterity. 

# ### The `PySCFWriter` object
# The PySCFWriter writes input files and stores the input parameters within itself.

# In[ ]:


from autopyscf import PySCFWriter

pyscf_writer = PySCFWriter()
pyscf_writer.__dict__


# ### Importing a structure.
# The only required input (without a default) is the structure definition. 
# For molecules, this means specifying the xyz coordinates of the atoms.

# In[ ]:


# A hydrogen molecule.
h2='\n'.join([
    'H 0.0 0.0 0.0 ',
    'H 0.74 0.0 0.0 '
])
pyscf_writer.set_options({'xyz':h2})


# ### Changing an option.
# Options can be generally changed using `set_options`. 
# Equivilently, you can instead have passed this dictionary of options to the constructor.

# In[ ]:


# Let's do a restricted PBE calculation instead of ROHF.
pyscf_writer.set_options({
    'method':'RKS',
    'dft':'pbe,pbe'
})


# Options get checked by `set_options` and rejected if they have obvious problems.

# In[ ]:


# This keyword is misspelled:
pyscf_writer.set_options({'converge_tol':1e-4})


# In[ ]:


# This spin setting incompatible with the charge setting.
pyscf_writer.set_options({'charge':2,'spin':1})


# ### Creating the input file.
# Once all the input is set up correctly, use `pyscf_input` to direct it to export the settings to a file.
# You should specify the chkfile, which is where PySCF will store much of the results.
# `h2_pyscf_input.py` is now a runable PySCF file.

# In[ ]:


pyscf_writer.pyscf_input('h2_pyscf_input.py','h2_pyscf_input.chkfile')


# # Next: how to read the results of the run (02-reading_output.ipynb).
