# -*- coding: utf-8 -*-
###############################################################################
#
#   Copyright (c) 2014 Thomas Fiala (fiala@td.mw.tum.de), Lehrstuhl für
#   Thermodynamik, TU München. For conditions of distribution and use, see
#   copyright notice in License.txt.
#
###############################################################################
"""
This example computes the extinction point of a counterflow diffusion flame.
A hydrogen-oxygen diffusion flame at 1 bar is studied.

The tutorial makes use of the scaling rules derived by Fiala and Sattelmayer
(doi:10.1155/2014/484372). Please refer to this publication for a detailed
explanation. Also, please don't forget to cite it if you make use of it.
"""

import cantera as ct
import numpy as np
import os

# Create directory for output data files
data_directory = 'diffusion_flame_extinction_data/'
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

# PART 1: INITIALIZATION

# Set up an initial hydrogen-oxygen counterflow flame at 1 bar and low strain
# rate (maximum axial velocity gradient = 2414 1/s)

# Initial grid: 18mm wide, 21 points
reaction_mechanism = 'h2o2.xml'
gas = ct.Solution(reaction_mechanism)
initial_grid = np.linspace(0.0, 18.e-3, 21)
f = ct.CounterflowDiffusionFlame(gas, initial_grid)

# Define the operating pressure and boundary conditions
f.P = 1.e5  # 1 bar
f.fuel_inlet.mdot = 0.5  # kg/m^2/s
f.fuel_inlet.X = 'H2:1'
f.fuel_inlet.T = 300  # K
f.oxidizer_inlet.mdot = 3.0  # kg/m^2/s
f.oxidizer_inlet.X = 'O2:1'
f.oxidizer_inlet.T = 500  # K

# Define relative and absolute error tolerances
f.flame.set_steady_tolerances(default=[1.0e-5, 1.0e-12])
f.flame.set_transient_tolerances(default=[5.0e-4, 1.0e-11])

# Enable refinement
refine = True
# Set refinement parameters
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
f.set_grid_min(1e-20)

# Define a limit for the maximum temperature below which the flame is
# considered as extinguished and the computation is aborted
temperature_limit_extinction = 500  # K

# Initialize and solve
f.set_initial_guess(fuel='H2')
print('Creating the initial solution')
f.solve(loglevel=0, refine_grid=refine)

# Save to data directory
file_name = 'initial_solution.xml'
f.save(data_directory + file_name, name='solution',
       description='Cantera version ' + ct.__version__ +
       ', reaction mechanism ' + reaction_mechanism)


# PART 2: COMPUTE EXTINCTION STRAIN

# Exponents for the initial solution variation with changes in strain rate
# Taken from Fiala and Sattelmayer (2014)
exp_d_a = - 1. / 2.
exp_u_a = 1. / 2.
exp_V_a = 1.
exp_lam_a = 2.
exp_mdot_a = 1. / 2.

# Set normalized initial strain rate
alpha = [1.]
# Initial relative strain rate increase
delta_alpha = 1.
# Factor of refinement of the strain rate increase
delta_alpha_factor = 50.
# Limit of the refinement: Minimum normalized strain rate increase
delta_alpha_min = .001
# Limit of the Temperature decrease
delta_T_min = 1  # K

# Iteration indicator
n = 0
# Indicator of the latest flame still burning
n_last_burning = 0
# List of peak temperatures
T_max = [np.max(f.T)]
# List of maximum axial velocity gradients
a_max = [np.max(np.abs(np.gradient(f.u) / np.gradient(f.grid)))]

# Simulate counterflow flames at increasing strain rates until the flame is
# extinguished. To achieve a fast simulation, an initial coarse strain rate
# increase is set. This increase is reduced after an extinction event and
# the simulation is again started based on the last burning solution.
# The extinction point is considered to be reached if the abortion criteria
# on strain rate increase and peak temperature decrease are fulfilled.
while True:
    n += 1
    # Update relative strain rates
    alpha.append(alpha[n_last_burning] + delta_alpha)
    strain_factor = alpha[-1] / alpha[n_last_burning]
    # Create an initial guess based on the previous solution
    # Update grid
    f.flame.grid *= strain_factor ** exp_d_a
    normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
    # Update mass fluxes
    f.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
    f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
    # Update velocities
    f.set_profile('u', normalized_grid, f.u * strain_factor ** exp_u_a)
    f.set_profile('V', normalized_grid, f.V * strain_factor ** exp_V_a)
    # Update pressure curvature
    f.set_profile('lambda', normalized_grid, f.L * strain_factor ** exp_lam_a)
    try:
        f.solve(loglevel=0, refine_grid=refine)
    except Exception as e:
        # Throw Exception if solution fails
        print('Error: Did not converge at n =', n, e)
    if np.max(f.T) > temperature_limit_extinction:
        # Flame still burning, so go to next strain rate
        n_last_burning = n
        file_name = 'extinction_{0:04d}.xml'.format(n)
        f.save(data_directory + file_name, name='solution', loglevel=0,
               description='Cantera version ' + ct.__version__ +
               ', reaction mechanism ' + reaction_mechanism)
        T_max.append(np.max(f.T))
        a_max.append(np.max(np.abs(np.gradient(f.u) / np.gradient(f.grid))))
        # If the temperature difference is too small and the minimum relative
        # strain rate increase is reached, abort
        if ((T_max[-2] - T_max[-1] < delta_T_min) &
                (delta_alpha < delta_alpha_min)):
            print('Flame extinguished at n = {0}.'.format(n),
                  'Abortion criterion satisfied.')
            break
    else:
        # Procedure if flame extinguished but abortion criterion is not satisfied
        print('Flame extinguished at n = {0}. Restoring n = {1} with alpha = {2}'.format(
              n, n_last_burning, alpha[n_last_burning]))
        # Reduce relative strain rate increase
        delta_alpha = delta_alpha / delta_alpha_factor
        # Restore last burning solution
        file_name = 'extinction_{0:04d}.xml'.format(n_last_burning)
        f.restore(data_directory + file_name, name='solution', loglevel=0)


# Print some parameters at the extinction point
print('----------------------------------------------------------------------')
print('Parameters at the extinction point:')
print('Pressure p={0} bar'.format(f.P / 1e5))
print('Peak temperature T={0:4.0f} K'.format(np.max(f.T)))
print('Mean axial strain rate a_mean={0:.2e} 1/s'.format(f.strain_rate('mean')))
print('Maximum axial strain rate a_max={0:.2e} 1/s'.format(f.strain_rate('max')))
print('Fuel inlet potential flow axial strain rate a_fuel={0:.2e} 1/s'.format(
      f.strain_rate('potential_flow_fuel')))
print('Oxidizer inlet potential flow axial strain rate a_ox={0:.2e} 1/s'.format(
      f.strain_rate('potential_flow_oxidizer')))
print('Axial strain rate at stoichiometric surface a_stoich={0:.2e} 1/s'.format(
      f.strain_rate('stoichiometric', fuel='H2')))

# Plot the maximum temperature over the maximum axial velocity gradient
import matplotlib.pyplot as plt
plt.figure()
plt.semilogx(a_max, T_max)
plt.xlabel(r'$a_{max}$ [1/s]')
plt.ylabel(r'$T_{max}$ [K]')
plt.savefig(data_directory + 'figure_T_max_a_max.png')
