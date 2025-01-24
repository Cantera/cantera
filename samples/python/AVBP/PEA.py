import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Parameters of the solver
loglevel = 1
refine_grid = 'refine'

if ct.__version__ >= '2.5.0':
    mechanism_list = ['gri30.yaml', './inputs/BFER_methane.yaml', './inputs/BFER_methane.yaml']
else:
    mechanism_list = ['gri30.cti', './inputs/BFER_methane.cti', './inputs/BFER_methane.cti']
gas_list = ['', 'BFER_nopea', 'BFER_pea']

for mechanism_name, gas_name in zip(mechanism_list, gas_list):
    # Setting up the gas object (operating conditions)
    gas = ct.Solution(mechanism_name, gas_name)
    p = 1e5
    tin = 473.0
    phi = 1.4
    fuel = {'CH4': 1}
    oxidizer = {'O2': 1, 'N2': 3.76}
    gas.TP = tin, p
    gas.set_equivalence_ratio(phi, fuel, oxidizer)

    # Setting up the flame object
    f = ct.FreeFlame(gas, width=0.2)  # Careful about the domain size if you want to retrieve the same profile !
    f.inlet.X = gas.X
    f.inlet.T = gas.T

    f.energy_enabled = True
    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05, prune=0.01)
    f.set_max_jac_age(10, 10)
    f.set_time_step(1e-8, [10, 20, 40, 80, 100, 100, 150])
    f.max_time_step_count = 5000

    f.solve(loglevel, refine_grid)

    # Checking the results
    if ct.__version__ >= '2.5.0':
        plt.plot(f.flame.grid, f.velocity, label=mechanism_name + ' - ' + gas_name)
    else:
        plt.plot(f.grid, f.u, label=mechanism_name + ' - ' + gas_name)

plt.xlabel('Flame grid [m]')
plt.ylabel('Veloctity [m/s]')
plt.legend()
plt.show()
