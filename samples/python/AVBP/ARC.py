import cantera as ct
import matplotlib.pyplot as plt
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

p = 1e5
tin = 400.0
phi = 1.0
fuel = {'CH4':1}
oxidizer = {'O2': 1, 'N2': 3.76}

ct.compile_fortran('./inputs/Lu_ARC.f90')

if ct.__version__ >= '2.5.0':
    gas = ct.Solution('./inputs/Lu_ARC.yaml')
else:
    gas = ct.Solution('./inputs/Lu_ARC.cti')

gas.TP = tin, p
gas.set_equivalence_ratio(phi, fuel, oxidizer)
# gas.transport_model = 'UnityLewis'

loglevel = 1 
refine_grid = 'refine'

f = ct.FreeFlame(gas, width=0.02)
f.inlet.X = gas.X
f.inlet.T = gas.T

f.energy_enabled = True
f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05, prune=0.01)
f.set_max_jac_age(10, 10)
f.set_time_step(1e-8, [10, 20, 40, 80, 100, 100, 150,])
f.max_time_step_count = 500

f.solve(loglevel, refine_grid)

if ct.__version__ >= '2.5.0':
    plt.plot(f.flame.grid, f.velocity)
else:
    plt.plot(f.flame.grid, f.u)
    
plt.xlabel('Flame grid [m]')
plt.ylabel('Velocity [m/s]')
plt.show()