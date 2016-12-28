from __future__ import print_function
import os
import numpy as np
import subprocess as sub
import shutil
from submitter import LocalSubmitter
import crystal2qmc

####################################################

class LocalCrystal2QMCRunner(LocalSubmitter):
  """ Converts from a CRYSTAL run to a QWalk format"""
  _name_='LocalCrystal2QMCRunner'
  def __init__(self,qwbase='qw',propoutfn='prop.in.o',
               kfmt='coord',kset='complex'):
    self.qwbase=qwbase
    self.propoutfn=propoutfn
    self.kfmt=kfmt
    self.kset=kset
    self.kweights=""


  #-------------------------------------------------      
  def check_status(self):
    #We always just run it or not.
    return 'ok'

  #-------------------------------------------------      
  def run(self):
    self.kweights=crystal2qmc.convert_crystal(self.qwbase,self.propoutfn,
                                self.kfmt,self.kset)
    print(self.kweights)
