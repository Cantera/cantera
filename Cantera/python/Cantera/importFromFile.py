import solution
import Interface
import XML

def importPhase(file = '', name = ''):
    return importPhases(file, [name])[0]

def importPhases(file = '', names = []):
    """Import multiple phase definitions.
    """
    s = []
    for nm in names:
        s.append(solution.Solution(src=file,id=nm))
    return s

def importInterface(file = '', name = '', phases = []):
    if name:
        src = file+'#'+name
    else:
        src = file
    return Interface.Interface(src = src, phases = phases)    
    
