import solution
import Interface
import XML

def preprocess(f):
    fn = f.split('.')
    prepr = 1
    base = f
    if len(fn) == 2:
        base = fn[0]
        if fn[1] == '.xml' or fn[1] == '.ctml':
            prepr = 0
    if prepr:
        from Cantera import pip
        pip.process(f)                
        src = base+'.xml'
    else:
        src = f
    return src

def importPhase(file = '', name = ''):
    return importPhases(file, [name])[0]

def importPhases(file = '', names = []):
    """Import multiple phase definitions.
    By importing all required phases in one file with one function call,
    the preprocessor and CTML parser only need to run once.
    """
    s = []
    #root = XML.XML_Node(name = 'doc', src = file, preprocess = 1)
    for nm in names:
        src = file+'#'+nm
        s.append(solution.Solution(src))
    return s

def importInterface(file = '', name = '', phases = []):
    #file = preprocess(file)
    #root = XML.XML_Node(name = 'doc', src = file, preprocess = 1)
    if name:
        src = file+'#'+name
    else:
        src = file
    return Interface.Interface(src = src, phases = phases)    
    
