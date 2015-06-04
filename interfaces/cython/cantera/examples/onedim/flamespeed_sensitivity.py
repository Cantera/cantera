"""
Sensitivity analysis for a freely-propagating, premixed methane-air
flame. Computes the sensitivity of the laminar flame speed with respect
to each reaction rate constant.
"""

from __future__ import print_function

import cantera as ct
import numpy as np

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:0.45, O2:1.0, N2:3.76'

initial_grid = np.linspace(0, 0.03, 5)  # m
tol_ss = [1.0e-9, 1.0e-14]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-14]  # [rtol atol] for time stepping

# IdealGasMix object used to compute mixture properties
gas = ct.Solution('gri30.xml', 'gri30_mix')
gas.TPX = Tin, p, reactants

# Flame object
f = ct.FreeFlame(gas, initial_grid)
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

# Solve with the energy equation disabled
f.energy_enabled = False
f.set_max_jac_age(10, 10)
f.set_time_step(1e-5, [2, 5, 10, 20])
f.solve(loglevel=1, refine_grid=False)

# Solve with the energy equation enabled
f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
f.energy_enabled = True
f.solve(loglevel=1, refine_grid=True)

Su0 = f.u[0]
print('\nmixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))

print('Initial Solution:')
f.show_stats()

# Perturbation size. This must be large compared to the steady-state relative
# tolerance (tol_ss[0]. Sensitivities less than approximately tol_ss[0] / dk
# are not reliable.
dk = 1e-2

print()
print('Rxn #   k/S*dS/dk    Reaction Equation')
print('-----   ----------   ----------------------------------')
for m in range(gas.n_reactions):
    gas.set_multiplier(1.0) # reset all multipliers
    gas.set_multiplier(1+dk, m) # perturb reaction m
    f.solve(loglevel=0, refine_grid=False)
    Su = f.u[0]
    print('{: 5d}   {: 10.3e}   {}'.format(
          m, (Su-Su0)/(Su0*dk), gas.reaction_equation(m)))

# Sensitivity analysis requires additional function evaluations on the final
# grid, but no additional Jacobian evaluations.
print('\nInitial Solution + Sensitivity calculations:')
f.show_stats()
