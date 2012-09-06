from . import PureFluid

def Water():
    """Create a `PureFluid` object using the equation of state for water."""
    return PureFluid('liquidvapor.xml','water')

def Nitrogen():
    """Create a `PureFluid` object using the equation of state for nitrogen."""
    return PureFluid('liquidvapor.xml','nitrogen')

def Methane():
    """Create a `PureFluid` object using the equation of state for methane."""
    return PureFluid('liquidvapor.xml','methane')

def Hydrogen():
    """Create a `PureFluid` object using the equation of state for hydrogen."""
    return PureFluid('liquidvapor.xml','hydrogen')

def Oxygen():
    """Create a `PureFluid` object using the equation of state for oxygen."""
    return PureFluid('liquidvapor.xml','oxygen')

def Hfc134a():
    """Create a `PureFluid` object using the equation of state for HFC-134a."""
    return PureFluid('liquidvapor.xml','hfc134a')

def CarbonDioxide():
    """Create a `PureFluid` object using the equation of state for carbon dioxide."""
    return PureFluid('liquidvapor.xml','carbondioxide')

def Heptane():
    """Create a `PureFluid` object using the equation of state for heptane."""
    return PureFluid('liquidvapor.xml','heptane')
