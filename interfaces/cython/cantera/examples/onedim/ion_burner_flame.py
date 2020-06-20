"""
A burner-stabilized premixed methane-air flame with charged species.

Requires: cantera >= 2.5.0
"""

import cantera as ct

p = ct.one_atm
tburner = 600.0
reactants = 'CH4:1.0, O2:2.0, N2:7.52'  # premixed gas composition
width = 0.5  # m
loglevel = 1  # amount of diagnostic output (0 to 5)

gas = ct.Solution('gri30_ion.yaml')
gas.TPX = tburner, p, reactants
mdot = 0.15 * gas.density

f = ct.IonBurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show_solution()

f.transport_model = 'Ion'
f.solve(loglevel, auto=True)
f.solve(loglevel=loglevel, stage=2, enable_energy=True)
try:
    # save to HDF container file if h5py is installed
    f.write_hdf('ion_burner_flame.h5', group='ion', mode='w',
                description='solution with ionized gas transport')
except ImportError:
    f.save('ion_burner_flame.xml', 'mix', 'solution with mixture-averaged transport')

f.write_csv('ion_burner_flame.csv', quiet=False)
