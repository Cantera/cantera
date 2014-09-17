"""
Atomic elements.

"""

def elementMoles(s, element):
    """Number of moles of an element in one mole of a solution.

    s       -- an object representing a solution.
    element -- the symbol for an element in 's'.
    """
    # see if 'element' corresponds to a symbol for one of the elements
    # in s. If it does not, return zero moles.
    try:
        m = s.elementIndex(element)
        if m < 0.0: return 0.0
    except:
        return 0.0

    x = s.moleFractions()
    moles = 0.0
    for k in range(s.nSpecies()):
        moles += x[k]*s.nAtoms(k,m)
    return moles
