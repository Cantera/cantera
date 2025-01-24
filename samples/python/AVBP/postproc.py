import cantera as ct
import numpy as np
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

# Compute thermal flame thickness
z = f.flame.grid
T = f.T
grad = np.gradient(T,z)
thickness = (max(T) -min(T)) / max(grad)
print(' ')
print('laminar thermal flame thickness (delta_T) = ', thickness)
print(' ')

# Compute HRR flame thickness
# Remark : slight difference when using np.gradient (5% in relative)
# Therefore manual computation of gradient is used
HRR = f.heat_release_rate
gradHRR = (HRR[2:]-HRR[0:-2])/(z[2:]-z[0:-2])
HRR = HRR[1:-1]
x = z[1:-1]
gradHRR2 = (gradHRR[2:]-gradHRR[0:-2])/(x[2:]-x[0:-2])
x = x[1:-1]
gradHRR = gradHRR[1:-1]
HRR = HRR[1:-1]
delta_HRR = -2.0*(-2.0*HRR*gradHRR2)**0.5/gradHRR2
where = np.where(HRR > 1e3)
delta_HRR = np.nanmin(delta_HRR[where])

print(' ')
print('laminar HRR flame thickness (delta_HRR) = ', delta_HRR)
print(' ')
print('delta_T / delta_HRR = ', thickness / delta_HRR)
print(' ')

# Compute schmidt numbers
print('SCHMIDT NUMBERS (in burnt gases):')
f.set_gas_state(f.flame.n_points-1)
list_species = f.gas.species_names
for species in list_species:
  print (species,gas.viscosity/gas.mix_diff_coeffs[gas.species_index(species)]/gas.density)
print(' ')

print('PRANDTL NUMBER (in fresh gases):')
f.set_gas_state(f.flame.n_points-f.flame.n_points)
Pr = gas.cp_mass * gas.viscosity / gas.thermal_conductivity
print('Prandtl = '+str(Pr))