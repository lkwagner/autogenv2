# You can easily define new wave functions here. 
# The only requirement is to define the export() method, which defines how to generate the QWalk input. 
import os

#######################################################################
class TrialFunction:
  ''' Skeleton class to define API for Trial Functions.'''
  def export(self):
    raise NotImplementedError("All trial functions must have export function defined.")

#######################################################################
class SlaterJastrow(TrialFunction):
  def __init__(self,slatman,jastman=None,kpoint=0):
    ''' Generate a Slater-Jastrow wave function from a manager that generates a Slater determinant and
      a manager that generates a Jastrow factor.

    Args: 
      slatman (Manager): Manager with a Slater-determinant-generating result.
      jastman (Manager): Manager with a Jastrow-generating result. 
      kpoint (int): kpoint number (as determined by the slatman converter).
    Returns:
      str or None: None if managers are not ready, QWalk section (str) if they are.
    '''
    # TODO more transparent kpoint selection.
    self.slatman=slatman
    if jastman is None:
      self.jastman=slatman
    else:
      self.jastman=jastman

  #------------------------------------------------
  def export(self,qmcpath):
    ''' Export the wavefunction section for this trial wave function.
    Args: 
      path (str): QWalkManager.path
      kpoint (int): the kpoint to choose for exporting. 
    Returns:
      str: system and wave fumction section for QWalk. Empty string if not ready.
    '''

    # Ensure files are correctly generated.
    if not (self.slatman.export_qwalk() and self.jastman.export_qwalk()):
      return ''

    if type(self.slatman.qwfiles['slater'])==list:
      slater=self.slatman.qwfiles['slater'][kpoint]
      sys=self.slatman.qwfiles['sys'][kpoint]
    else:
      slater=self.slatman.qwfiles['slater']
      sys=self.slatman.qwfiles['sys']
    jastrow=self.jastman.qwfiles['jastrow2']

    # There may be a use case for these two to be different, but I want to check the first time this happens. 
    # You can have weird bugs if you use different system files for each wave function term, I think.
    # Should find a way to check for this bug. 
    assert self.jastman.path+self.jastman.qwfiles['sys']==self.slatman.path+self.slatman.qwfiles['sys'],\
        'System file probably should be the same between Jastrow and Slater files. '

    outlines=[
        'include %s/%s'%(os.path.relpath(qmcpath,self.slatman.path),sys),
        'trialfunc { slater-jastrow ',
        '  wf1 { include %s/%s }'%(os.path.relpath(qmcpath,self.slatman.path),slater),
        '  wf2 { include %s/%s }'%(os.path.relpath(qmcpath,self.slatman.path),jastrow),
        '}'
      ]
    return '\n'.join(outlines)

