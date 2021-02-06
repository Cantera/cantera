"""
A burner-stabilized lean premixed hydrogen-oxygen flame at low pressure.

Requires: cantera >= 2.5.0
"""

import cantera as ct

p = 0.05 * ct.one_atm
tburner = 373.0
mdot = 0.06
reactants = 'H2:1.5, O2:1, AR:7'  # premixed gas composition
width = 0.5  # m
loglevel = 1  # amount of diagnostic output (0 to 5)

gas = ct.Solution('h2o2.yaml')
gas.TPX = tburner, p, reactants

f = ct.BurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show_solution()

f.transport_model = 'Mix'
f.solve(loglevel, auto=True)
try:
    # save to HDF container file if h5py is installed
    f.write_hdf('burner_flame.h5', group='mix', mode='w',
                description='solution with mixture-averaged transport')
except ImportError:
    f.save('burner_flame.xml', 'mix', 'solution with mixture-averaged transport')

f.transport_model = 'Multi'
f.solve(loglevel)  # don't use 'auto' on subsequent solves
f.show_solution()
try:
    f.write_hdf('burner_flame.h5', group='multi',
                description='solution with multicomponent transport')
except ImportError:
    f.save('burner_flame.xml', 'multi', 'solution with multicomponent transport')

f.write_csv('burner_flame.csv', quiet=False)
