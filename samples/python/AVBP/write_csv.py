import cantera as ct
import matplotlib.pyplot as plt
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Parameters of the solver
loglevel = 1
refine_grid = 'refine'

# Setting up the gas object (operating conditions)
if ct.__version__ >= '2.5.0':
    gas = ct.Solution('gri30.yaml')
else:
    gas = ct.Solution('gri30.cti')
p = 1e5
tin = 300.0
phi = 1.0
fuel = {'CH4':1}
oxidizer = {'O2': 1, 'N2': 3.76}
gas.TP = tin, p
gas.set_equivalence_ratio(phi, fuel, oxidizer)

# Setting up the flame object
f = ct.FreeFlame(gas, width=0.02)
f.inlet.X = gas.X
f.inlet.T = gas.T

f.energy_enabled = True
f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05, prune=0.01)
f.set_max_jac_age(10, 10)
f.set_time_step(1e-8, [10, 20, 40, 80, 100, 100, 150])
f.max_time_step_count = 500

f.solve(loglevel, refine_grid)

# Write AVBP solution for can2av
f.write_AVBP("./RESULTS/AVBP-Solution.csv")