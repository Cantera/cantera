import solution

def importPhase(file = '', name = ''):
    if name:
        src = file+'#'+name
    else:
        src = file
    return solution.Solution(src)

