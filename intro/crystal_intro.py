
# This script is meant to serve as a brief introduction to using crystal and autogen.

# Access the autogen libraries.
import sys
sys.path.append('..')

from crystal import CrystalWriter,CrystalReader
from autorunner import PySCFRunnerLocal,RunnerLocal
from trialfunc import SlaterJastrow
from variance import VarianceWriter, VarianceReader
from crystalmanager import CrystalManager
from qwalkmanager import QWalkManager

# Where the work is stored.
path='crystal-scratch'

# A hydrogen molecule.
h2='\n'.join([
    '2',
    'h2 molecule',
    'H 0.0 0.0 0.0 ',
    'H 0.74 0.0 0.0 '
])

# Input parameters for SCF calculation.
crystal_writer=CrystalWriter({
    'functional':{'exchange':'PBE','correlation':'PBE','hybrid':25},
    'xml_name':'../../BFD_Library.xml'
  })
crystal_writer.set_struct_fromxyz(h2)

# Manage the PySCF job.
scfmanager=CrystalManager(
    name='h2_pbe',
    path=path,
    writer=crystal_writer,
    runner=RunnerLocal()
  )

scfmanager.nextstep()
