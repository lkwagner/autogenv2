import os,glob

class Crystal2QMCReader:
  """ Tries to extract properties of crystal run, or else diagnose what's wrong. """
  def __init__(self):
    self.completed=False
    self.out={}
#-------------------------------------------------      
  def collect(self,basename='qw'):
    """ Collect results from output."""
    self.out['sysfiles']=glob.glob(basename+"*.sys")
    self.out['orbfiles']=glob.glob(basename+"*.orb")
    self.out['slaterfiles']=glob.glob(basename+"*.slater")
    if len(self.out['sysfiles'])!=len(self.out['orbfiles']) or\
        len(self.out['orbfiles'])!=len(self.out['slaterfiles']):
      print("Inconsistency in conversion!")
      print(self.out['sysfiles'])
      print(self.out['slaterfiles'])
      print(self.out['orbfiles'])
      self.completed=False
      raise Error

    self.out['basenames']=[x.replace('.sys','') for x in self.out['sysfiles']]
    print(self.out['basenames'])
    self.completed=True

    
#-------------------------------------------------      
  def check_status(self,basename='qw'):
    """ Decide status of conversion."""
    status='ok'
    print("status",status)
    return status
    
