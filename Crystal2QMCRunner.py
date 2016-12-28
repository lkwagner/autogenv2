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
  def is_consistent(self,other):
    skipkeys = ['kweights']
    for otherkey in other.__dict__.keys():
      if otherkey not in self.__dict__.keys():
        print('other is missing a key.')
        return False
    for selfkey in self.__dict__.keys():
      if selfkey not in other.__dict__.keys():
        print('self is missing a key.')
        return False
    for key in self.__dict__.keys():
      if self.__dict__[key]!=other.__dict__[key] and key not in skipkeys:
        print("Different keys [{}] = \n{}\n or \n {}"\
            .format(key,self.__dict__[key],other.__dict__[key]))
        return False
    return True
    
  #-------------------------------------------------      
  def check_status(self):
    #We always just run it or not.
    return 'ok'

  #-------------------------------------------------      
  def run(self):
    self.kweights=crystal2qmc.convert_crystal(self.qwbase,self.propoutfn,
                                self.kfmt,self.kset)
    print(self.kweights)
