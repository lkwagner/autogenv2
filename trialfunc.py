# You can easily define new wave functions here. 
# The only requirement is to define the export() method, which defines how to generate the QWalk input. 
import os

#######################################################################
class TrialFunction:
  ''' Skeleton class to define API for Trial Functions.'''
  def export(self):
    raise NotImplementedError("All trial functions must have export function defined.")

#######################################################################
class Slater(TrialFunction):
  def __init__(self,slatman,kpoint=0):
    ''' Generate a Slater-Jastrow wave function from a manager that generates a Slater determinant and
      a manager that generates a Jastrow factor.

    Args: 
      slatman (Manager): Manager with a Slater-determinant-generating result.
      jastman (Manager): Manager with a Jastrow-generating result. 
      kpoint (int): kpoint number (as determined by the slatman converter). None implies its a finite system.
    Returns:
      str or None: None if managers are not ready, QWalk section (str) if they are.
    '''
    self.slatman=slatman
    self.kpoint=kpoint

  #------------------------------------------------
  def export(self,qmcpath):
    ''' Export the wavefunction section for this trial wave function.
    Args: 
      path (str): QWalkManager.path
      kpoint: the kpoint to choose for exporting. 
    Returns:
      str: system and wave fumction section for QWalk. Empty string if not ready.
    '''
    # This assumes you're using 2-body, should be easy to make a new object or maybe an arg for 3body.

    # Ensure files are correctly generated.
    if not (self.slatman.export_qwalk()):
      return ''

    if type(self.slatman.qwfiles['slater'])==str:
      slater=self.slatman.qwfiles['slater']
      sys=self.slatman.qwfiles['sys']
    else:
      slater=self.slatman.qwfiles['slater'][self.kpoint]
      sys=self.slatman.qwfiles['sys'][self.kpoint]

    outlines=[
        'include %s'%os.path.relpath(self.slatman.path+sys,qmcpath),
        'trialfunc { ',
        '  include %s '%os.path.relpath(self.slatman.path+slater,qmcpath),
        '}'
      ]
    return '\n'.join(outlines)

#######################################################################
class SlaterJastrow(TrialFunction):
  def __init__(self,slatman,jastman=None,kpoint=0):
    ''' Generate a Slater-Jastrow wave function from a manager that generates a Slater determinant and
      a manager that generates a Jastrow factor.

    Args: 
      slatman (Manager): Manager with a Slater-determinant-generating result.
      jastman (Manager): Manager with a Jastrow-generating result. 
      kpoint (int): kpoint number (as determined by the slatman converter). None implies its a finite system.
    Returns:
      str or None: None if managers are not ready, QWalk section (str) if they are.
    '''
    self.slatman=slatman
    if jastman is None:
      self.jastman=slatman
    else:
      self.jastman=jastman

    self.kpoint=kpoint

  #------------------------------------------------
  def export(self,qmcpath):
    ''' Export the wavefunction section for this trial wave function.
    Args: 
      path (str): QWalkManager.path
      kpoint: the kpoint to choose for exporting. 
    Returns:
      str: system and wave fumction section for QWalk. Empty string if not ready.
    '''
    # This assumes you're using 2-body, should be easy to make a new object or maybe an arg for 3body.

    # Ensure files are correctly generated.
    if not (self.slatman.export_qwalk() and self.jastman.export_qwalk()):
      return ''

    if type(self.slatman.qwfiles['slater'])==str:
      slater=self.slatman.qwfiles['slater']
      sys=self.slatman.qwfiles['sys']
    else:
      slater=self.slatman.qwfiles['slater'][self.kpoint]
      sys=self.slatman.qwfiles['sys'][self.kpoint]
    jastrow=self.jastman.qwfiles['jastrow2']

    # There may be a use case for these two to be different, but I want to check the first time this happens. 
    # You can have weird bugs if you use different system files for each wave function term, I think.
    # Should find a way to check for this bug. 
    # This doesn't work because we may have different jastrows in same directory, for instance.
    #assert (self.jastman.path+self.jastman.name)==\
    #    (self.slatman.path+self.slatman.name),\
    #    'System file probably should be the same between Jastrow and Slater files. '

    outlines=[
        'include %s'%os.path.relpath(self.slatman.path+sys,qmcpath),
        'trialfunc { slater-jastrow ',
        '  wf1 { include %s }'%os.path.relpath(self.slatman.path+slater,qmcpath),
        '  wf2 { include %s }'%os.path.relpath(self.slatman.path+jastrow,qmcpath),
        '}'
      ]
    return '\n'.join(outlines)

