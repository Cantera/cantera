"""Fluids with complete liquid/vapor equations of state."""

from importFromFile import importPhase

def Water():
    return importPhase('liquidvapor.cti','water')

def Nitrogen():
    return importPhase('liquidvapor.cti','nitrogen')

def Methane():
    return importPhase('liquidvapor.cti','methane')

def Hydrogen():
    return importPhase('liquidvapor.cti','hydrogen')

def Oxygen():
    return importPhase('liquidvapor.cti','oxygen')
