"""Functions to import phase and interface definitions from CTI files."""

import solution
import Interface
import XML

def importPhase(file = '', name = '', loglevel = 0):
    """Import a phase from a CTI file."""
    return importPhases(file, [name], loglevel)[0]

def importPhases(file = '', names = [], loglevel = 0):
    """Import multiple phases from one file. The phase names should be
    entered as a list of strings.  """
    s = []
    for nm in names:
        s.append(solution.Solution(src=file,id=nm,loglevel=loglevel))
    return s

def importInterface(file = '', name = '', phases = []):
    """Import an interface definition from a CTI file."""
    if name:
        src = file+'#'+name
    else:
        src = file
    return Interface.Interface(src = src, phases = phases)    
    
