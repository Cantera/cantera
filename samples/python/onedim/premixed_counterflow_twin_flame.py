"""
Symmetric premixed twin flame
=============================

Simulate two counter-flow jets of premixed reactants shooting into each other. This
simulation differs from the similar :doc:`premixed_counterflow_flame.py
<premixed_counterflow_flame>` example as the latter simulates a jet of reactants
shooting into products.

Requires: cantera >= 3.0, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, premixed flame, strained flame, plotting

An illustration of this configuration is shown in the figure below:

.. image:: /_static/images/samples/twin-premixed-flame.svg
   :width: 90%
   :alt: Twin flame established by two opposed jets of fuel + oxidizer mixtures
   :align: center

"""

import sys
from pathlib import Path
import numpy as np
import cantera as ct


# Differentiation function for data that has variable grid spacing Used here to
# compute normal strain-rate
def derivative(x, y):
    dydx = np.zeros(y.shape, y.dtype.type)

    dx = np.diff(x)
    dy = np.diff(y)
    dydx[0:-1] = dy/dx

    dydx[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])

    return dydx


def compute_strain_rates(opposed_flame):
    # Compute the derivative of axial velocity to obtain normal strain rate
    strain_rates = derivative(opposed_flame.grid, opposed_flame.velocity)

    # Obtain the location of the max. strain rate upstream of the pre-heat zone.
    # This is the characteristic strain rate
    max_strain_location = abs(strain_rates).argmax()
    min_velocity_point = opposed_flame.velocity[:max_strain_location].argmin()

    # Characteristic Strain Rate = K
    strain_rate_point = abs(strain_rates[:min_velocity_point]).argmax()
    K = abs(strain_rates[strain_rate_point])

    return strain_rates, strain_rate_point, K


def compute_consumption_speed(opposed_flame):
    Tb = max(opposed_flame.T)
    Tu = min(opposed_flame.T)
    rho_u = max(opposed_flame.density)

    integrand = opposed_flame.heat_release_rate / opposed_flame.cp

    total_heat_release = np.trapz(integrand, opposed_flame.grid)
    Sc = total_heat_release / (Tb - Tu) / rho_u

    return Sc


# This function is called to run the solver
def solve_opposed_flame(opposed_flame, mass_flux=0.12, loglevel=1,
                        ratio=2, slope=0.3, curve=0.3, prune=0.05):
    """
    Execute this function to run the Opposed Flow Simulation. This function
    takes a CounterFlowTwinPremixedFlame object as the first argument
    """

    opposed_flame.reactants.mdot = mass_flux
    opposed_flame.set_refine_criteria(ratio=ratio, slope=slope,
                                      curve=curve, prune=prune)

    opposed_flame.show()
    opposed_flame.solve(loglevel, auto=True)

    # Compute the strain rate, just before the flame. This is not necessarily
    # the maximum We use the max. strain rate just upstream of the pre-heat zone
    # as this is the strain rate that computations compare against, like when
    # plotting Su vs. K
    strain_rates, strain_rate_point, K = compute_strain_rates(opposed_flame)

    return np.max(opposed_flame.T), K, strain_rate_point

# %%
# Define parameters
# -----------------

# Select the reaction mechanism
gas = ct.Solution('gri30.yaml')

# Create a CH4/Air premixed mixture with equivalence ratio=0.75, and at room
# temperature and pressure.
gas.set_equivalence_ratio(0.75, 'CH4', {'O2': 1.0, 'N2': 3.76})
gas.TP = 300, ct.one_atm

# Set the velocity
axial_velocity = 2.0  # in m/s

# Domain half-width of 2.5 cm, meaning the whole domain is 5 cm wide
width = 0.025

# %%
# Set up the simulation
# ---------------------

# Compute the mass flux, as this is what the Flame object requires
mass_flux = gas.density * axial_velocity  # units kg/m2/s

# Create the flame object
opposed_flame = ct.CounterflowTwinPremixedFlame(gas, width=width)

# Uncomment the following line to use a Multi-component formulation. Default is
# mixture-averaged
# opposed_flame.transport_model = 'multicomponent'

# %%
# Run the solver
# --------------
#
# The solver returns the peak temperature, strain rate and the point which we ascribe to
# the characteristic strain rate.

T, K, strain_rate_point = solve_opposed_flame(opposed_flame, mass_flux, loglevel=1)

# %%
# Compute flame properties
# ------------------------
#
# You can plot/see all state space variables by calling `opposed_flame.foo`, where
# `foo` is `T`, `Y[i]`, and so forth. The spatial variable (distance in meters) is in
# `opposed_flame.grid`. Thus to plot temperature vs distance, use `opposed_flame.grid``
# and `opposed_flame.T`.

Sc = compute_consumption_speed(opposed_flame)

if "native" in ct.hdf_support():
    output = Path() / "premixed_counterflow_twin_flame.h5"
else:
    output = Path() / "premixed_counterflow_twin_flame.yaml"
output.unlink(missing_ok=True)

opposed_flame.save(output, name="mix")

print(f"Peak temperature: {T:.1f} K")
print(f"Strain Rate: {K:.1f} 1/s")
print(f"Consumption Speed: {Sc * 100:.2f} cm/s")
opposed_flame.save("premixed_counterflow_twin_flame.csv", basis="mole", overwrite=True)

# %%
# Plot results
# ------------
#
# Generate plots to see results, if desired.
if '--plot' in sys.argv:

    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 6), facecolor='white')

    # Axial Velocity Plot
    plt.subplot(1, 2, 1)
    plt.plot(opposed_flame.grid, opposed_flame.velocity, 'r', lw=2)
    plt.xlim(opposed_flame.grid[0], opposed_flame.grid[-1])
    plt.xlabel('Distance (m)')
    plt.ylabel('Axial Velocity (m/s)')

    # Identify the point where the strain rate is calculated
    plt.plot(opposed_flame.grid[strain_rate_point],
             opposed_flame.velocity[strain_rate_point], 'gs')
    plt.annotate('Strain-Rate point',
                 xy=(opposed_flame.grid[strain_rate_point],
                     opposed_flame.velocity[strain_rate_point]),
                 xytext=(0.001, 0.1),
                 arrowprops={'arrowstyle': '->'})

    # Temperature Plot
    plt.subplot(1, 2, 2)
    plt.plot(opposed_flame.grid, opposed_flame.T, 'b', lw=2)
    plt.xlim(opposed_flame.grid[0], opposed_flame.grid[-1])
    plt.xlabel('Distance (m)')
    plt.ylabel('Temperature (K)')

    plt.tight_layout()
    plt.show()

else:
    print('************')
    print('Plotting option not enabled. Re-run script with --plot to see key plots.')
    print('************')
