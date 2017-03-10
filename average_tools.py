

################################################
def average_section(opts):
  ''' Generate QWalk input for average generator section. '''
  outlines=[]
  if opts['name'].lower()=='average_derivative_dm':
    # opts['basis'] should be a fully defined orbitals section. 
    # See the QWalk docs on tbdm_basis for details.
    outlines=[
        "  average { average_derivative_dm",
        "    average { tbdm_basis",
        "      include %s"%opts['basis'],
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
  check={'average_derivative_dm': ['basis'],
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
def _kaverage_tbdm(data):
  nkpt=len(data)
  nstates=len(data[0]['states'])
  res={'obdm':{},'tbdm':{}}

  # Wait, do we need this? This is already in tbdm, no?
  for key in ['up','down']:
    res['obdm'][key]=[
        [
          sum([
            data[i]['obdm'][key][k][l]
            for i in range(nkpt)
          ])/nkpt 
          for l in range(nstates)
        ] 
        for k in range(nstates)
      ] 

  for key in ['upup','updown','downup','downdown']:
    res['tbdm'][key]=[
        [
          [
            [
              sum([
                data[i]['tbdm'][key][k][l][m][n]
                for i in range(nkpt)
              ])/nkpt 
              for n in range(nstates)
            ] 
            for m in range(nstates)
          ]
          for l in range(nstates)
        ] 
        for k in range(nstates)
      ] 
  return res

################################################
def _kaverage_deriv(data):
  res={}
  nparm=len(data[0]['dpenergy']['vals'])
  nkpt=len(data)
  nstates=len(data[0]['tbdm']['states'])

  # Parameters with one value per parameter values.
  for prop in ['dpenergy','dpwf']:
    res[prop]=[
        sum([
          data[i][prop]['vals'][j]
          for i in range(nkpt)
        ])/nkpt
        for j in range(nparm)
      ]
    res['%s_err'%prop]=[
        (sum([
          data[i][prop]['err'][j]**2
          for i in range(nkpt)
        ])/nkpt)**0.5
        for j in range(nparm)
      ]

  res['tbdm']=_kaverage_tbdm([data[k]['tbdm'] for k in range(nkpt)])

  res['dprdm']=[
      _kaverage_tbdm(
        [data[k]['dprdm'][i]['tbdm'] for k in range(nkpt)]
      ) for i in range(nparm)
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
