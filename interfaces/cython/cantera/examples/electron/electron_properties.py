import cantera as ct
import numpy as np
import scipy.constants as co

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
T = 2000.0  # temperature [K]
reactants = 'N2:3.76, O2:1'  # air composition

air = ct.Solution(infile='gri30_ion.cti', efile='lxcat.yaml')
air.TPX = T, p, reactants
air.electric_field = 1e5

# set up proper grid for the electric field
grid = np.linspace(0.0, 30, num=500)
air.set_electron_energy_grid(grid)

print('electron diffusivity = {0:7f} [m^2/s]'.format(air.electron_diffusivity))
print('electron mobility = {0:7f} [m^2/V/s]'.format(air.electron_mobility))
print('electron total collision frequency = {0:7f} [Hz]'.format(air.electron_total_collision_frequency))
print('electron power gain = {0:7f} [eV/s]'.format(air.electron_power_gain))
print('electron elastic power loss = {0:7f} [eV/s]'.format(air.electron_elastic_power_loss))
print('electron inelastic power loss = {0:7f} [eV/s]'.format(air.electron_inelastic_power_loss))
print('mean electron energy = {0:7f} [eV]'.format(air.mean_electron_energy))

# Here, we demonstrate another way to initiate the gas for electron properties. 
# This method enable we to use multiple files
ecs = ct.ElectronCrossSection.listFromFile('lxcat.yaml')
gas = ct.Solution('gri30_plasma.cti', electron='WeaklyIonizedGas', electron_cross_sections=ecs)
gas.TPX = T, p, reactants
gas.electric_field = 1e5

# set up proper grid for the electric field
grid = np.linspace(0.0, 30, num=500)
gas.set_electron_energy_grid(grid)

# We can also use species index to find electron transport properties
gas.electron_transport_enabled = True
i = gas.species_index('E')
print('electron diffusivity = {0:7f} [m^2/s]'.format(gas.mix_diff_coeffs[i]))
print('electron mobility = {0:7f} [m^2/V/s]'.format(gas.mobilities[i]))

# and reaction rate of electron reaction which depends on electron temperature
gas.electron_temperature_reactions_enabled = True
print(gas.reaction_equation(331))
print('rate coefficient = {0:7f} [m^3/kmol/s]'.format(gas.forward_rate_constants[331]))

