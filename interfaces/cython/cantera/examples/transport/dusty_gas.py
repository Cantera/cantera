"""
Dusty Gas transport model.

The Dusty Gas model is a multicomponent transport model for gas transport
through the pores of a stationary porous medium. This example shows how to
create a transport manager that implements the Dusty Gas model and use it to
compute the multicomponent diffusion coefficients.

Requires:  cantera >= 2.5.0
"""

import cantera as ct

# create a gas-phase object to represent the gas in the pores, with a
# dusty gas transport manager
g = ct.DustyGas('h2o2.yaml')

# set the gas state
composition = {"OH": 1, "H": 2, "O2": 3, "O": 1.0E-8, "H2": 1.0E-8, "H2O": 1.0E-8,
               "H2O2": 1.0E-8, "HO2": 1.0E-8, "AR": 1.0E-8}
g.TPX = 500.0, ct.one_atm, composition

# set its parameters
g.porosity = 0.2
g.tortuosity = 4.0
g.mean_pore_radius = 1.5e-7
g.mean_particle_diameter = 1.5e-6  # lengths in meters

# print the multicomponent diffusion coefficients
print(g.multi_diff_coeffs)

# compute molar species fluxes
T1, rho1, Y1 = g.TDY

g.TP = g.T, 1.2 * ct.one_atm
T2, rho2, Y2 = g.TDY
delta = 0.001

print(g.molar_fluxes(T1, T1, rho1, rho1, Y1, Y1, delta))
print(g.molar_fluxes(T1, T2, rho1, rho2, Y1, Y2, delta))
