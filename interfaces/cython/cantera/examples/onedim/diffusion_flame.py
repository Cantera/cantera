"""
An opposed-flow ethane/air diffusion flame
"""

import cantera as ct
import numpy as np
import csv

p = ct.OneAtm  # pressure
tin_f = 300.0  # fuel inlet temperature
tin_o = 300.0  # oxidizer inlet temperature
mdot_o = 0.72  # kg/m^2/s
mdot_f = 0.24  # kg/m^2/s

comp_o = 'O2:0.21, N2:0.78, AR:0.01'  # air composition
comp_f = 'C2H6:1'  # fuel composition

initial_grid = np.linspace(0, 0.02, 6)

tol_ss = [1.0e-5, 1.0e-12]  # [rtol, atol] for steady-state problem
tol_ts = [5.0e-4, 1.0e-9]  # [rtol, atol] for time stepping

loglevel = 1  # amount of diagnostic output (0 to 5)
refine_grid = 1   # 1 to enable refinement, 0 to disable

gas = ct.Solution('gri30.xml', 'gri30_mix')
gas.TP = gas.T, p

f = ct.CounterflowDiffusionFlame(gas, initial_grid)
f.fuel_inlet.mdot = mdot_f
f.fuel_inlet.X = comp_f
f.fuel_inlet.T = tin_f

f.oxidizer_inlet.mdot = mdot_o
f.oxidizer_inlet.X = comp_o
f.oxidizer_inlet.T = tin_o

f.flame.setSteadyTolerances(default=tol_ss)
f.flame.setTransientTolerances(default=tol_ts)

f.setInitialGuess(fuel='C2H6')
f.energyEnabled = False
f.solve(loglevel, refine_grid=False)

f.energyEnabled = True
f.setRefineCriteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)
f.solve(loglevel, refine_grid=refine_grid)

f.showSolution()

f.save('c2h6_diffusion.xml')

z = f.flame.grid
T = f.T
u = f.u
V = f.V

with open('c2h6_diffusion.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)', 'rho (kg/m3)'] +
                    list(gas.speciesNames))
    for n in range(f.flame.nPoints):
        f.setGasState(n)
        writer.writerow([z[n], u[n], V[n], T[n], gas.density] + list(gas.X))

print('solution saved to c2h6_diffusion.csv')
