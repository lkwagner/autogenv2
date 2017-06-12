
######################################################################
def diff_keys(old,new,skip_keys=[]):
  ''' Check if two objects have different keys, and what those keys are. '''
  issame=True
  diff={'old':[],'new':[]}

  for newkey in new.__dict__.keys():
    if newkey not in old.__dict__.keys():
      issame=False
      diff['new'].append(newkey)
  for oldkey in old.__dict__.keys():
    if oldkey not in new.__dict__.keys():
      issame=False
      diff['old'].append(oldkey)
  for key in old.__dict__.keys():
    if (key not in diff['new']) and (key not in diff['old']) and \
        (old.__dict__[key]!=new.__dict__[key]) and (key not in skip_keys):
      issame=False
      diff['old'].append(key)
      diff['new'].append(key)
  return issame,diff

######################################################################
def update_attributes(old,new,skip_keys=[],safe_keys=[]):
  ''' Replace attributes that do not affect accuracy. 
  Raise an AssertionError if there's a problematic discrepancy. 
  By default, all keys are not safe, so this mainly checks consistency.
  skip_keys are not checked or replaced.'''

  issame,diff=diff_keys(old,new,skip_keys)
  if not issame:
    print("Key update: {} from one doesn't match {} from new."\
        .format(diff['old'],diff['new']))
    for key in diff['new']:
      if key in safe_keys:
        print("Keeping {} from the latter.".format(diff['new']))
        old.__dict__[key]=new.__dict__[key]
      else:
        raise AssertionError("Unsafe update; new setting affects accuracy.")
