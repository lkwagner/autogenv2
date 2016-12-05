import PropertiesWriter
import os

class PropertiesRun:
  _name_="PropertiesRun"
  def __init__(self,pwriter,submitter=None):
    self._submitter=submitter
    self._pwriter=pwriter

  def run(self,job_record):
    with open("prop.in",'w') as f:
      f.write(self._pwriter.properties_input())
    if self._submitter==None:
      os.system("properties < prop.in > prop.in.o")
      return 'ok'
    
