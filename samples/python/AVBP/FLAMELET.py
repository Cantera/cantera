import cantera as ct
import time
import numpy as np
from scipy.special import erfinv, erfcinv, erf
from scipy import integrate
import matplotlib.pyplot as plt

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

#gas init

T = 300.0                                 # temperature
P = 101325.0                              # pressure

rho_o = P / (8.314 / 0.032 * T)           # oxidizer density
tin_o = 300.0                             # oxidizer inlet temperature
mdot_o = 0.8                              # kg/m^2/s
comp_o = 'O2:1'                           # oxidizer composition
u_o = mdot_o / rho_o                      # oxidizer speed

rho_f = P / (8.314 / 0.016 * T)           # fuel density
tin_f = 300.0                             # fuel inlet temperature
mdot_f = 0.2                              # kg/m^2/s
comp_f = 'CH4:1'                          # fuel composition
u_f = mdot_f / rho_f                      # fuel speed

width = 1.0                               # width
SR = (u_o + u_f) / (width)                # strain rate

if ct.__version__ >= '2.5.0':
    ct.compile_fortran('./inputs/Lu_ARC.f90')
    gas_1 = ct.Solution('./inputs/Lu_ARC.yaml')
else:
    gas_1 = ct.Solution('./inputs/Lu.cti')

gas_1.TP = T, P

#flamelet creation
plt.rcParams['figure.figsize'] = (14, 10)

initial_grid = np.linspace(0,width,11)                      # initial grid for the calculation

f_1 = ct.Flamelet(gas_1)                      # Flamelet object creation

f_1.fuel_inlet.mdot = mdot_f                                # Fuel inlet conditions
f_1.fuel_inlet.Y = comp_f
f_1.fuel_inlet.T = tin_f
                                                            # Oxidizer inlet conditions
f_1.oxidizer_inlet.mdot = mdot_o
f_1.oxidizer_inlet.Y = comp_o
f_1.oxidizer_inlet.T = tin_o
chist = SR/np.pi*(np.exp(-2*((erfinv(1-2*0.2))**2)))
f_1.ChiSt = chist    # Stochiometric scalar dissipation rate

f_1.set_initial_guess()

# plt.plot(f_1.flame.grid, f_1.T)

# Calculation properties

loglevel = 1
f_1.set_max_jac_age(1, 1)
f_1.set_time_step(1.e-9, [50, 100, 150])
f_1.energy_enabled = True
f_1.set_refine_criteria(ratio=4, slope=0.1, curve=0.1)

# Solver
t1_A = time.time()
f_1.solve(loglevel, refine_grid = 'refine')
t1_B = time.time()

if ct.__version__ >= '2.5.0':
    f_1.save('./RESULTS/CH4-O2-converged-flamelet-Lu.yaml', 'Lu.yaml-zspace',overwrite=True)
else:
    f_1.save('./RESULTS/CH4-O2-converged-flamelet-Lu.xml')

# Calculation of the heat release
hr_1 = -np.sum(f_1.standard_enthalpies_RT * f_1.net_production_rates, 0) * ct.gas_constant * f_1.T
# print('Heat release rate = ', hr_1)

plt.rcParams['figure.figsize'] = (14, 10)

# Get interesting indices for computation of species
fuel_species = 'CH4'
ifuel = gas_1.species_index(fuel_species)
io2 = gas_1.species_index('O2')
ico = gas_1.species_index('CO')

# Initiate interesting vectors
ch4_1 = np.zeros(f_1.flame.n_points,'d')
o2_1 = np.zeros(f_1.flame.n_points,'d')
co_1 = np.zeros(f_1.flame.n_points,'d')

# Computes interesting quantities for analyzing a counter-flow flame
for n in range(f_1.flame.n_points):
    f_1.set_gas_state(n)
    ch4_1[n]= gas_1.Y[ifuel]
    o2_1[n]= gas_1.Y[io2]
    co_1[n]= gas_1.Y[ico]

T_1 = f_1.T
    
plt.plot(f_1.flame.grid,T_1/np.max(T_1),f_1.flame.grid,ch4_1/np.max(ch4_1),
         f_1.flame.grid,o2_1/np.max(o2_1),f_1.flame.grid,co_1/np.max(co_1), f_1.flame.grid,hr_1/np.max(hr_1))
plt.title(r'$T_{adiabatic}, Y_{CH_4},  Y_{O_2}, Y_{CO}$, HR vs. Mixture fraction',fontsize=25)
plt.xlabel(r'Mixture fraction', fontsize=15)
plt.ylabel('Normalized values of different quantities',fontsize=15)
plt.legend(['Temperature','$Y_{CH_4}$', '$Y_{O_2}$', '$Y_{CO}$', 'HR'],fontsize=15)

plt.show()