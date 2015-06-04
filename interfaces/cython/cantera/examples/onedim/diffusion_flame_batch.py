# -*- coding: utf-8 -*-
###############################################################################
#
#   Copyright (c) 2014 Thomas Fiala (fiala@td.mw.tum.de), Lehrstuhl für
#   Thermodynamik, TU München. For conditions of distribution and use, see
#   copyright notice in License.txt.
#
###############################################################################
"""
This example creates two batches of counterflow diffusion flame simulations.
The first batch computes counterflow flames at increasing pressure, the second
at increasing strain rates.

The tutorial makes use of the scaling rules derived by Fiala and Sattelmayer
(doi:10.1155/2014/484372). Please refer to this publication for a detailed
explanation. Also, please don't forget to cite it if you make use of it.

This example can e.g. be used to iterate to a counterflow diffusion flame to an
awkward  pressure and strain rate, or to create the basis for a flamelet table.
"""

import cantera as ct
import numpy as np
import os

# Create directory for output data files
data_directory = 'diffusion_flame_batch_data/'
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

# Set refinement: False for fast simulations, True for smoother curves
refine = True

# PART 1: INITIALIZATION

# Set up an initial hydrogen-oxygen counterflow flame at 1 bar and low strain
# rate (maximum axial velocity gradient = 2414 1/s)

# Initial grid: 18mm wide, 21 points
reaction_mechanism = 'h2o2.xml'
gas = ct.Solution(reaction_mechanism)
initial_grid = np.linspace(0.0, 18e-3, 21)
f = ct.CounterflowDiffusionFlame(gas, initial_grid)

# Define the operating pressure and boundary conditions
f.P = 1.e5  # 1 bar
f.fuel_inlet.mdot = 0.5  # kg/m^2/s
f.fuel_inlet.X = 'H2:1'
f.fuel_inlet.T = 300  # K
f.oxidizer_inlet.mdot = 3.0  # kg/m^2/s
f.oxidizer_inlet.X = 'O2:1'
f.oxidizer_inlet.T = 300  # K

# Define relative and absolute error tolerances
f.flame.set_steady_tolerances(default=[1.0e-5, 1.0e-12])
f.flame.set_transient_tolerances(default=[5.0e-4, 1.0e-11])

# Set refinement parameters, if used
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
f.set_grid_min(1e-20)

# Define a limit for the maximum temperature below which the flame is
# considered as extinguished and the computation is aborted
# This increases the speed of refinement is enabled
temperature_limit_extinction = 900  # K
def interrupt_extinction(t):
    if np.max(f.T) < temperature_limit_extinction:
        raise Exception('Flame extinguished')
    return 0.
f.set_interrupt(interrupt_extinction)

# Initialize and solve
f.set_initial_guess(fuel='H2')
print('Creating the initial solution')
f.solve(loglevel=0, refine_grid=refine)

# Save to data directory
file_name = 'initial_solution.xml'
f.save(data_directory + file_name, name='solution',
       description='Cantera version ' + ct.__version__ +
       ', reaction mechanism ' + reaction_mechanism)


# PART 2: BATCH PRESSURE LOOP

# Compute counterflow diffusion flames over a range of pressures
# Arbitrarily define a pressure range (in bar)
p_range = np.round(np.logspace(0, 2, 50), decimals=1)

# Exponents for the initial solution variation with changes in pressure Taken
# from Fiala and Sattelmayer (2014). The exponents are adjusted such that the
# strain rates increases proportional to p^(3/2), which results in flames
# similar with respect to the extinction strain rate.
exp_d_p = -5. / 4.
exp_u_p = 1. / 4.
exp_V_p = 3. / 2.
exp_lam_p = 4.
exp_mdot_p = 5. / 4.

# The variable p_previous (in bar) is used for the pressure scaling
p_previous = f.P / 1.e5
# Iterate over the pressure range
for p in p_range:
    print('pressure = {0} bar'.format(p))
    # set new pressure
    f.P = p * 1.e5
    # Create an initial guess based on the previous solution
    rel_pressure_increase = p / p_previous
    # Update grid
    f.flame.grid *= rel_pressure_increase ** exp_d_p
    normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
    # Update mass fluxes
    f.fuel_inlet.mdot *= rel_pressure_increase ** exp_mdot_p
    f.oxidizer_inlet.mdot *= rel_pressure_increase ** exp_mdot_p
    # Update velocities
    f.set_profile('u', normalized_grid,
                  f.u * rel_pressure_increase ** exp_u_p)
    f.set_profile('V', normalized_grid,
                  f.V * rel_pressure_increase ** exp_V_p)
    # Update pressure curvature
    f.set_profile('lambda', normalized_grid,
                  f.L * rel_pressure_increase ** exp_lam_p)

    try:
        # Try solving the flame
        f.solve(loglevel=0, refine_grid=refine)
        file_name = 'pressure_loop_' + format(p, '05.1f') + '.xml'
        f.save(data_directory + file_name, name='solution', loglevel=1,
               description='Cantera version ' + ct.__version__ +
               ', reaction mechanism ' + reaction_mechanism)
        p_previous = p
    except Exception as e:
        print('Error occurred while solving:', e, 'Try next pressure level')
        # If solution failed: Restore the last successful solution and continue
        f.restore(filename=data_directory + file_name, name='solution',
                  loglevel=0)


# PART 3: STRAIN RATE LOOP

# Compute counterflow diffusion flames at increasing strain rates at 1 bar
# The strain rate is assumed to increase by 25% in each step until the flame is
# extinguished
strain_factor = 1.25

# Exponents for the initial solution variation with changes in strain rate
# Taken from Fiala and Sattelmayer (2014)
exp_d_a = - 1. / 2.
exp_u_a = 1. / 2.
exp_V_a = 1.
exp_lam_a = 2.
exp_mdot_a = 1. / 2.

# Restore initial solution
file_name = 'initial_solution.xml'
f.restore(filename=data_directory + file_name, name='solution', loglevel=0)

# Counter to identify the loop
n = 0
# Do the strain rate loop
while np.max(f.T) > temperature_limit_extinction:
    n += 1
    print('strain rate iteration', n)
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
        # Try solving the flame
        f.solve(loglevel=0, refine_grid=refine)
        file_name = 'strain_loop_' + format(n, '02d') + '.xml'
        f.save(data_directory + file_name, name='solution', loglevel=1,
               description='Cantera version ' + ct.__version__ +
               ', reaction mechanism ' + reaction_mechanism)
    except Exception as e:
        if e.args[0] == 'Flame extinguished':
            print('Flame extinguished')
        else:
            print('Error occurred while solving:', e)
        break


# PART 4: PLOT SOME FIGURES

import matplotlib.pyplot as plt

fig1 = plt.figure()
fig2 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
ax2 = fig2.add_subplot(1,1,1)
p_selected = p_range[::7]

for p in p_selected:
    file_name = 'pressure_loop_{0:05.1f}.xml'.format(p)
    f.restore(filename=data_directory + file_name, name='solution', loglevel=0)

    # Plot the temperature profiles for selected pressures
    ax1.plot(f.grid / f.grid[-1], f.T, label='{0:05.1f} bar'.format(p))

    # Plot the axial velocity profiles (normalized by the fuel inlet velocity)
    # for selected pressures
    ax2.plot(f.grid / f.grid[-1], f.u / f.u[0],
             label='{0:05.1f} bar'.format(p))

ax1.legend(loc=0)
ax1.set_xlabel(r'$z/z_{max}$')
ax1.set_ylabel(r'$T$ [K]')
fig1.savefig(data_directory + 'figure_T_p.png')

ax2.legend(loc=0)
ax2.set_xlabel(r'$z/z_{max}$')
ax2.set_ylabel(r'$u/u_f$')
fig2.savefig(data_directory + 'figure_u_p.png')

fig3 = plt.figure()
fig4 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax4 = fig4.add_subplot(1,1,1)
n_selected = range(1, n, 5)
for n in n_selected:
    file_name = 'strain_loop_{0:02d}.xml'.format(n)
    f.restore(filename=data_directory + file_name, name='solution', loglevel=0)
    a_max = f.strain_rate('max') # the maximum axial strain rate

    # Plot the temperature profiles for the strain rate loop (selected)
    ax3.plot(f.grid / f.grid[-1], f.T, label='{0:.2e} 1/s'.format(a_max))

    # Plot the axial velocity profiles (normalized by the fuel inlet velocity)
    # for the strain rate loop (selected)
    ax4.plot(f.grid / f.grid[-1], f.u / f.u[0],
             label=format(a_max, '.2e') + ' 1/s')

ax3.legend(loc=0)
ax3.set_xlabel(r'$z/z_{max}$')
ax3.set_ylabel(r'$T$ [K]')
fig3.savefig(data_directory + 'figure_T_a.png')

ax4.legend(loc=0)
ax4.set_xlabel(r'$z/z_{max}$')
ax4.set_ylabel(r'$u/u_f$')
fig4.savefig(data_directory + 'figure_u_a.png')
