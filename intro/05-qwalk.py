
# coding: utf-8

# In[ ]:


# Access the autogen libraries.
import sys
sys.path.append('..')


# # Trial wave functions and QMC job management.
# This notebook walks through the process of building a QWalk job. 
# A few notes:
# * A QWalk writer, reader, and runner will work exactly the same as corresponding PySCF objects.
# * A QWalk manager will handle waiting for other runs to finish, and converting the files into a trial wave function.

# ### Setting up the trial wave function calculation.
# For the hydrogen run, let's optimize a Slater-Jastrow wave function with QWalk variance optimization.
# First step is setting up the mean field calculation.

# In[ ]:


from autopyscf import PySCFWriter,PySCFReader
from autorunner import PySCFRunnerLocal
from manager import PySCFManager

h2='\n'.join([
    'H 0.0 0.0 0.0 ',
    'H 0.74 0.0 0.0 '
])
path='05-scratch'

# This is covered in 04-job_management.
pyscf_manager = PySCFManager(
    path=path,
    name='h2_pbe',
    writer=PySCFWriter({
        'xyz':h2,
        'method':'RKS',
        'dft':'pbe,pbe'}),
    runner=PySCFRunnerLocal()
)


# The PySCF calculation is defined, but so far no calculations have been started.
# 
# QWalk writers can take explicit trial wave function input (as a string), but autogen can also handle it for you.
# Make a `TrialFunc` object with the manager you want to generate the trial wave function.

# In[ ]:


from trialfunc import SlaterJastrow
sj=SlaterJastrow(pyscf_manager)


# ### Setting up the QWalk run.
# Now we can set up the corresponding objects for QWalk.

# In[ ]:


from variance import VarianceWriter,VarianceReader
from autorunner import RunnerLocal
from manager import QWalkManager
from imp import reload
import time

vopt_writer = VarianceWriter({'iterations':6})

vopt_manager = QWalkManager(
    name='h2_vopt',
    path=path,
    writer=vopt_writer,
    runner=RunnerLocal(np=4),
    reader=VarianceReader(),
    trialfunc=sj # This is new!
)
while not vopt_manager.completed:
    vopt_manager.nextstep()

print(vopt_manager.reader.output['sigma'])
import matplotlib.pyplot as plt

plt.plot(vopt_manager.reader.output['sigma'])
plt.ylabel("Standard deviation of energy")
plt.xlabel("Optimization iteration")
plt.savefig("05-qwalk_variance.pdf",bbox_inches='tight')
# So what just happened here? 
# After setting up the calculation, we called `nextstep` repeatedly and the QWalk manager:
# * Noticed the PySCF run was not complete, and ran the PySCF run.
# * The PySCF run exported the qwalk files needed by QWalk.
# * The `QWalkManager` generated its imput files using the files from PySCF.
# * `QWalkManager` ran QWalk and collected the QWalk results.
# 
# The TrialFunc objects are meant to be modular and flexible. 
# You can build TrialFunc objects in a variety of ways, which the next notebook will demonstrate.

# # Next: building different trial functions (06-trialfunctions.ipynb).
