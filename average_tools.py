

################################################
def average_section(opts):
  ''' Generate QWalk input for average generator section. '''
  outlines=[]
  if opts['name'].lower()=='average_derivative_dm':
    outlines=[
        "  average { average_derivative_dm",
        "    average { tbdm_basis",
        "      orbitals { ",
        "        cutoff_mo",
        "        magnify 1",
        "        nmo %d"%opts['nmo'],
        "        orbfile %s"%opts['orbfile'],
        "        include %s"%opts['basis'],
        "        centers { useglobal }",
        "      }",
        "      states { %s }"%' '.join(map(str,opts['states'])),
        "    }",
        "  }"
      ]
  elif opts['name'].lower()=='region_fluctuation':
    outlines=[
        "  average { region_fluctuation maxn %i } "%opts['maxn']
        ]
  else:
    raise NotImplementedError("""
    '%s' is not implemented in autogen yet: 
    You should implement it, it should be easy!"""%opts['name'])
  return outlines

################################################
def check_opts(opts):
  ''' Make sure options are set completely.'''
  check={'average_derivative_dm':
            ['nmo','orbfile','basis','states'],
         'region_fluctuation':['maxn']
         }
  for key in check[opts['name']]:
    assert key in opts.keys(),"""
      '%s' missing from 'extra_observables' settings!
      Make sure all of %s are set."""%(key,', '.join(check))

################################################
def kaverage(name,data):
  ''' kaverage the data for each property.'''
  if name=='average_derivative_dm':
    return _kaverage_deriv(data)
  elif name=='region_fluctuation':
    return [] # TODO
  else:
    raise NotImplementedError("""
    '%s' is not implemented in autogen yet: 
    You should implement it, it should be easy!"""%name)

################################################
def _kaverage_deriv(data):
  res={}
  nparm=len(data[0]['dpenergy']['vals'])
  nkpt=len(data)

  # Parameters with nparm values.
  for prop in ['dpenergy','dpwf']:
    res[prop]=[
        sum([
          data[i][prop]['vals'][j]
          for i in range(nkpt)
        ])/nparm
        for j in range(nparm)
      ]
    res['%s_err'%prop]=[
        (sum([
          data[i][prop]['err'][j]**2
          for i in range(nkpt)
        ])/nparm)**0.5
        for j in range(nparm)
      ]
  return res

################################################
def gosling_key(input_keyword):
  ''' Hack because gosling doesn't name json keys the same as the input keywords for QWalk '''
  # Should we change this in QWalk?
  propmap={
      'average_derivative_dm':'derivative_dm',
      'tbdm_basis':'tbdm'
    }
  if input_keyword in propmap:
    return propmap[input_keyword]
  else:
    return input_keyword
