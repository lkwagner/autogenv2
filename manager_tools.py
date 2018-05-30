import numpy as np
import os 

def resolve_status(runner,reader,outfile):
  #Check if the reader is done
  if reader.completed:
    return 'done'

  #Check if the job is in the queue or running. If so, we just return that.
  currstat=runner.check_status()
  if currstat=='running':
    return currstat
  
  #Now we are in a state where either there was an error,
  #the job hasn't been run, or we haven't collected the results
  if not os.path.exists(outfile):
    return 'not_started'

  #We are in an error state or we haven't collected the results. 
  return "ready_for_analysis"

######################################################################
def deep_compare(d1,d2):
  '''I have to redo dict comparison because numpy will return a bool array when comparing.'''
  if type(d1)!=type(d2):
    return False
  if type(d1)==dict:
    if d1.keys()!=d2.keys():
      return False
    allsame=True
    for key in d1.keys():
      allsame=allsame and deep_compare(d1[key],d2[key])
    return allsame
  else:
    try:
      return np.array_equal(d1,d2)
    except TypeError:
      return d1==d2

######################################################################
def update_attributes(copyto,copyfrom,skip_keys=[],take_keys=[]):
  ''' Save update of class attributes. If copyfrom has additional attributes, they are ignored.

  Args:
    copyto (obj): class who's attributes are being updated.
    copyfrom (obj): class who's attributes will be copied from.
    skip_keys (list): list of attributes (str) not to update. 
    take_keys (list): list of attributes (str) that are ok to update. Others will raise warning and be skipped.
  Returns:
    bool: Whether any changes were made.
  '''
  updated=False
  for key in copyfrom.__dict__.keys():
    if key in skip_keys: 
      #print("Skipping key (%s)"%key)
      pass
    elif key not in copyto.__dict__.keys():
      print("Warning: Object update. An attribute (%s) was skipped because it doesn't exist in both objects."%key)
    elif not deep_compare(copyto.__dict__[key],copyfrom.__dict__[key]):
      if key not in take_keys:
        print("Warning: update to attribute (%s) cancelled, because it requires job to be rerun."%key)
      else:
        #print("Copy",key)
        copyto.__dict__[key]=copyfrom.__dict__[key]
        updated=True
    else:
      #print("Keys match (%s)"%key)
      pass
  return updated

def separate_jastrow(wffile,optimizebasis=False):
  ''' Seperate the jastrow section of a QWalk wave function file.'''
  # Copied from utils/separate_jastrow TODO: no copy, bad
  wff=open(wffile,'r')
  tokens=wff.read().split('\n')
  in_jastrow=False
  nopen=0
  nclose=0
  jastlines=[]
  for line in tokens:
    if 'jastrow2' in line.lower():
      in_jastrow=True
    if in_jastrow:
      if not optimizebasis and 'optimizebasis' in line.lower():
        continue
      nopen+=line.count("{")
      nclose+=line.count("}")
    if in_jastrow and nopen >= nclose:
      jastlines.append(line)
  return '\n'.join(jastlines)
