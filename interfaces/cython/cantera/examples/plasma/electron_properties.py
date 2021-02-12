import cantera as ct
import numpy as np
import scipy.constants as co

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
T = 1000.0  # temperature [K]
reactants = 'O2:1'  # gas composition

gas = ct.Plasma(infile='oxygen_plasma.yaml')
gas.TPX = T, p, reactants
gas.electric_field = 1e5

# set up proper grid for the electric field
grid = np.linspace(0.0, 50, num=500)
gas.set_electron_energy_grid(grid)

print('electron diffusivity = {0:7f} [m^2/s]'.format(gas.electron_diffusivity))
print('electron mobility = {0:7f} [m^2/V/s]'.format(gas.electron_mobility))
print('electron total collision frequency = {0:7f} [Hz]'.format(gas.electron_total_collision_frequency))
print('electron power gain = {0:7f} [eV/s]'.format(gas.electron_power_gain))
print('electron elastic power loss = {0:7f} [eV/s]'.format(gas.electron_elastic_power_loss))
print('electron inelastic power loss = {0:7f} [eV/s]'.format(gas.electron_inelastic_power_loss))
print('mean electron energy = {0:7f} [eV]'.format(gas.mean_electron_energy))
