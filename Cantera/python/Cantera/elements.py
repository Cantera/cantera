
def elementMoles(mix, element):
    """Number of moles of an element in one mole of a mixture.

    mix     -- a mixture object.
    element -- the symbol for an element in 'mix'.
    """
    
    nsp = mix.nSpecies()
    try:
        m = mix.elementIndex(element)
	if m < 0.0: return 0.0
    except:
        return 0.0
    x = mix.moleFractions()
    moles = 0.0
    for k in range(nsp):
        moles += x[k]*mix.nAtoms(k,m)
    return moles


