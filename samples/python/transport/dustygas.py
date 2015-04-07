"""
 Dusty Gas transport model.

 The Dusty Gas model is a multicomponent transport model for gas
 transport through the pores of a stationary porous medium. This
 example shows how to create a transport manager that implements the
 Dusty Gas model and use it to compute the multicomponent diffusion
 coefficients.

 """

from Cantera import *
from Cantera.DustyGasTransport import *

# create a gas-phase object to represent the gas in the pores
g = importPhase('h2o2.cti')

# set the gas state
g.set(T = 500.0, P = OneAtm, X = "OH:1, H:2, O2:3, O:1.0E-8, H2:1.0E-8, H2O:1.0E-8, H2O2:1.0E-8, HO2:1.0E-8, AR:1.0E-8")

# create a Dusty Gas transport manager for this phase
d = DustyGasTransport(g)

# set its parameters
d.set(porosity = 0.2, tortuosity = 4.0,
      pore_radius = 1.5e-7, diameter = 1.5e-6)  # lengths in meters

# print the multicomponent diffusion coefficients
print d.multiDiffCoeffs()

# compute molar species fluxes
state1 = g.saveState()

g.set(P = 1.2*OneAtm)
state2 = g.saveState()
delta = 0.001

print d.molarFluxes(state1, state1, delta)
print d.molarFluxes(state1, state2, delta)
