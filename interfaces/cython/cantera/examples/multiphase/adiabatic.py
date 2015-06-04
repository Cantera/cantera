"""
Adiabatic flame temperature and equilibrium composition for a fuel/air mixture
as a function of equivalence ratio, including formation of solid carbon.
"""

import cantera as ct
import numpy as np
import sys
import csv

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the mixture.

T = 300.0
P = 101325.0

# phases
gas = ct.Solution('gri30.xml')
carbon = ct.Solution('graphite.xml')

# the phases that will be included in the calculation, and their initial moles
mix_phases = [(gas, 1.0), (carbon, 0.0)]

# gaseous fuel species
fuel_species = 'CH4'

# equivalence ratio range
npoints = 50
phi = np.linspace(0.3, 3.5, npoints)

##############################################################################

mix = ct.Mixture(mix_phases)

# create some arrays to hold the data
tad = np.zeros(npoints)
xeq = np.zeros((mix.n_species,npoints))

if gas.n_atoms(fuel_species,'O') > 0 or gas.n_atoms(fuel_species,'N') > 0:
    raise "Error: only hydrocarbon fuels are supported."

stoich_O2 = gas.n_atoms(fuel_species,'C') + 0.25*gas.n_atoms(fuel_species,'H')

for i in range(npoints):
    X = {fuel_species: phi[i] / stoich_O2, 'O2': 1.0, 'N2': 3.76}
    # set the gas state
    gas.TPX = T, P, X

    # create a mixture of 1 mole of gas, and 0 moles of solid carbon.
    mix = ct.Mixture(mix_phases)
    mix.T = T
    mix.P = P

    # equilibrate the mixture adiabatically at constant P
    mix.equilibrate('HP', solver='gibbs', max_steps=1000)

    tad[i] = mix.T
    print('At phi = {0:12.4g}, Tad = {1:12.4g}'.format(phi[i], tad[i]))
    xeq[:,i] = mix.species_moles

# write output CSV file for importing into Excel
csv_file = 'adiabatic.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['phi','T (K)'] + mix.species_names)
    for i in range(npoints):
        writer.writerow([phi[i], tad[i]] + list(xeq[:,i]))
print('Output written to {0}'.format(csv_file))

if '--plot' in sys.argv:
    import matplotlib.pyplot as plt
    plt.plot(phi, tad)
    plt.xlabel('Equivalence ratio')
    plt.ylabel('Adiabatic flame temperature [K]')
    plt.show()
