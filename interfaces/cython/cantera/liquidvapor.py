from . import PureFluid

def Water():
    return PureFluid('liquidvapor.xml','water')

def Nitrogen():
    return PureFluid('liquidvapor.xml','nitrogen')

def Methane():
    return PureFluid('liquidvapor.xml','methane')

def Hydrogen():
    return PureFluid('liquidvapor.xml','hydrogen')

def Oxygen():
    return PureFluid('liquidvapor.xml','oxygen')

def Hfc134a():
    return PureFluid('liquidvapor.xml','hfc134a')

def CarbonDioxide():
    return PureFluid('liquidvapor.xml','carbondioxide')

def Heptane():
    return PureFluid('liquidvapor.xml','heptane')
