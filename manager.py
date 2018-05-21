''' Manager classes handle the file maintainance of various runs.

They:
  - Call write, run, and read in the appropriate order.
  - Keep a record of progress on disk.
  - Hold a list of all file paths together in one place.
  - Report progress to the user.

All choices revolving around file names should be handled here.

Tips on avoiding bugs with the write/run/read cycle:
  - Every time data is updated, write the updates to the disk.
  - Every time the data is used, update the Manager from the disk.
'''
import crystal
import propertiesreader
import os
import autopyscf 
import shutil as sh
import numpy as np
import pyscf2qwalk
import crystal2qmc
import pickle as pkl
import autorunner
from autopaths import paths

from copy import deepcopy

#TODO slow queue checking might be made more efficient by saving the qstat results somewhere.


