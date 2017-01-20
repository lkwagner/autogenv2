

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
