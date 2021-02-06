"""
An equilibrium example with charged species in the gas phase
and multiple condensed phases.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""

import cantera as ct
import csv

# create objects representing the gas phase and the condensed phases. The gas
# is a mixture of multiple species, and the condensed phases are all modeled
# as incompressible stoichiometric substances. See file KOH.yaml for more
# information.
phases = ct.import_phases('KOH.yaml', ['K_solid', 'K_liquid', 'KOH_a', 'KOH_b',
                                       'KOH_liquid', 'K2O2_solid', 'K2O_solid',
                                       'KO2_solid', 'ice', 'liquid_water',
                                       'KOH_plasma'])

# create the Mixture object from the list of phases
mix = ct.Mixture(phases)
equil_data = []

# loop over temperature
for n in range(100):
    t = 350.0 + 50.0*n
    print('T = {0}'.format(t))
    mix.T = t
    mix.P = ct.one_atm
    mix.species_moles = "K:1.03, H2:2.12, O2:0.9"

    # set the mixture to a state of chemical equilibrium holding
    # temperature and pressure fixed
    # mix.equilibrate("TP",maxsteps=10000,loglevel=1)
    mix.equilibrate("TP", max_steps=10000, log_level=0)

    # store the moles of each species
    equil_data.append([t] + list(mix.species_moles))

with open("equil_koh.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['T'] + mix.species_names)

    writer.writerows(equil_data)
