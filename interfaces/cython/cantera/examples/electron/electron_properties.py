import cantera as ct
import numpy as np
import scipy.constants as co

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
T = 2000.0  # temperature [K]
reactants = 'N2:3.76, O2:1'  # air composition

air = ct.Solution(infile='gri30_ion.cti', efile='lxcat.yaml')
air.TPX = T, p, reactants

# We can use species index to find electron properties
i = air.species_index('E')
print(air.mix_diff_coeffs[i])  # electron diffusivity
print(air.mobilities[i])    # electron mobility

# Here, we demonstrate another way to initiate the gas for electron properties. 
ecs = ct.ElectronCrossSection.listFromFile('lxcat.yaml')
gas = ct.Solution('gri30.yaml', electron='MaxwellBoltzmann', electron_cross_sections=ecs)
gas.TPX = T, p, reactants

# For the case without species index, we can obtain electron properties as following
N = p / (co.k * T)   # Gas number density
print(gas.electron_diffusivity(N))
print(gas.electron_mobility(N))
# Note: This is convenient for scaling properties with N. You can simply use N=1.

